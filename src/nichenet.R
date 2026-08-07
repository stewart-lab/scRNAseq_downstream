# NicheNet analysis: ligand-receptor/target prioritization between two
# conditions, for one or more receiver cell types.
# Configured via the nichenet section of config.json.
# Ported from src/nichenet.Rmd (kept as a reference notebook -- this script
# is the configurable, non-interactive version run via run_downstream_toolkit.sh).
#
# Two things from the original notebook are deliberately not ported as-is:
# - The trailing "Requires specific ligand-target pair" block is dropped --
#   it duplicated sig_network_maps() below (same get_ligand_signaling_path/
#   diagrammer_format_signaling_graph/export_graph calls), which the main
#   per-receiver loop already calls for every receiver. The trailing block
#   only worked by accident, reusing whichever receiver's order_ligands/
#   order_targets happened to be left over from the loop's last iteration.
# - The random forest diagnostic (assess_rf_class_probabilities etc.) is
#   opt-in (run_random_forest_diagnostics, default FALSE) and, when enabled,
#   runs properly inside the per-receiver loop instead of once afterward
#   using stale leftover variables.

suppressPackageStartupMessages(library(nichenetr)) # tested with v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(viridis)
library(digest)
library(DiagrammeR)

# Record R session info (packages + versions) for run reproducibility --
# appended to the run's shared provenance file when invoked via
# run_downstream_toolkit.sh (PROVENANCE_FILE env var), else written
# standalone to the current directory.
.provenance_file <- Sys.getenv("PROVENANCE_FILE", unset = paste0("./sessionInfo_nichenet.R_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
cat(paste0("\n--- R sessionInfo (nichenet.R) ---\n", paste(capture.output(sessionInfo()), collapse = "\n"), "\n"), file = .provenance_file, append = TRUE)

### load config ###
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(GIT_DIR, "config.json"))
cfg <- config$nichenet
docker <- config$docker

if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
  DATA_DIR <- "./data/input_data/"
} else {
  DATA_DIR <- cfg$DATA_DIR
}
# NicheNet's prior network files are a shared reference resource (not
# per-run data), kept under ./data/ in the repo itself rather than DATA_DIR
# -- same relative path in both docker and conda mode, since ./data/ is
# already covered by run_downstream_toolkit.sh's blanket Docker mount for
# every method (no dedicated mount needed), and conda mode runs with cwd at
# the repo root either way.
NETWORKS_DIR <- cfg$NETWORKS_DIR

name <- cfg$name
condition_oi <- cfg$condition_oi # the "case" condition
condition_reference <- cfg$condition_reference # the "control" condition
condition_colname <- cfg$condition_colname
celltype_colname <- cfg$celltype_colname
receivers <- cfg$receivers
expression_pct <- cfg$expression_pct
padj_thresh <- cfg$padj_thresh
lfc_thresh <- cfg$lfc_thresh
run_rf_diagnostics <- isTRUE(as.logical(cfg$run_random_forest_diagnostics))

### set output directory ###
timestamp <- Sys.getenv("RUN_TIMESTAMP", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
output <- paste0("./shared_volume/output_nichenet_", name, "_", timestamp)
print(output)
dir.create(output, recursive = TRUE, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(file.path(GIT_DIR, "config.json"), file.path(output, "config.json"))
if (file.exists(.provenance_file)) {
  file.copy(.provenance_file, file.path(output, basename(.provenance_file)), overwrite = TRUE)
  file.remove(.provenance_file)
}

### load and merge the two Seurat objects (one per condition) ###
seuratObj1 <- readRDS(file.path(DATA_DIR, cfg$`SEURAT.file1`))
seuratObj2 <- readRDS(file.path(DATA_DIR, cfg$`SEURAT.file2`))
print(seuratObj1@meta.data %>% head())
print(seuratObj2@meta.data %>% head())

seuratObj <- merge(x = seuratObj1, y = seuratObj2, merge.data = TRUE)
seuratObj[["RNA"]] <- JoinLayers(seuratObj[["RNA"]])
# For older Seurat objects, you may need to run this
seuratObj <- UpdateSeuratObject(seuratObj)
print(seuratObj)
print(seuratObj@meta.data %>% head())
print(table(seuratObj@meta.data[[condition_colname]]))
print(table(seuratObj@meta.data[[celltype_colname]]))

### load NicheNet prior networks ###
ligand_target_matrix <- readRDS(file.path(NETWORKS_DIR, cfg$ligand_target_matrix_path))
lr_network <- readRDS(file.path(NETWORKS_DIR, cfg$lr_network_path))
weighted_networks <- readRDS(file.path(NETWORKS_DIR, cfg$weighted_networks_path))
ligand_tf_matrix <- readRDS(file.path(NETWORKS_DIR, cfg$ligand_tf_matrix_path))

### dimensionality reduction + overview UMAP ###
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(object = seuratObj))
seuratObj <- RunUMAP(seuratObj, dims = 1:15)
gg1 <- DimPlot(seuratObj, reduction = "umap", group.by = celltype_colname)
gg2 <- DimPlot(seuratObj, reduction = "umap", group.by = condition_colname)
pdf(file = paste0(output, "seurat_umap.pdf"), height = 6, width = 10)
print(patchwork::wrap_plots(plots = list(gg1, gg2)))
dev.off()

# set cell type annotation
Idents(seuratObj) <- seuratObj@meta.data[[celltype_colname]]
print(table(Idents(seuratObj)))

############### visualization functions ###############

# heatmap of ligand activity
lig_activity_heatmap <- function(vis_ligand_aupr, output, name, receiver) {
  pdf(file = paste0(output, name, "_", receiver, "_ligand_activity_heatmap", condition_oi, ".pdf"), width = 3, height = 8)
  print(make_heatmap_ggplot(vis_ligand_aupr,
    y_name = "Prioritized ligands",
    x_name = "Ligand activity", legend_title = "AUPR", color = "darkorange"
  ) +
    theme(axis.text.x.top = element_blank()))
  dev.off()
}

# heat map of target potential
lig_target_heatmap <- function(vis_ligand_target, output, name, receiver) {
  pdf(file = paste0(output, name, "_", receiver, "_ligand-target_reg_potential_heatmap", condition_oi, ".pdf"), width = 8, height = 5)
  print(make_heatmap_ggplot(vis_ligand_target,
    y_name = "Prioritized ligands", x_name = "Predicted target genes",
    color = "purple", legend_title = "Regulatory potential"
  ) +
    scale_fill_gradient2(low = "whitesmoke", high = "purple"))
  dev.off()
}

# heat map of ligand-receptor
LR_heatmap <- function(vis_ligand_receptor_network, output, name, receiver) {
  pdf(file = paste0(output, name, "_", receiver, "_ligand-receptor_interactions_heatmap", condition_oi, ".pdf"))
  print(make_heatmap_ggplot(t(vis_ligand_receptor_network),
    y_name = "Prioritized ligands",
    x_name = "Receptors", color = "mediumvioletred", legend_title = "Prior interaction potential"
  ))
  dev.off()
}

ligand_lfc_heatmap <- function(vis_ligand_lfc, output, name, receiver) {
  pdf(file = paste0(output, name, "_", receiver, "_ligand_logfc_heatmap", condition_oi, ".pdf"))
  print(make_threecolor_heatmap_ggplot(vis_ligand_lfc,
    y_name = "Prioritized ligands",
    x_name = "LFC in Sender", low_color = "midnightblue", mid_color = "white",
    mid = median(vis_ligand_lfc), high_color = "red", legend_title = "LFC"
  ))
  dev.off()
}

ligand_dotplot <- function(seuratObj, best_upstream_ligands, sender_celltypes, output, name, receiver) {
  cells_to_keep <- colnames(seuratObj)[seuratObj@meta.data[[celltype_colname]] %in% sender_celltypes]
  pdf(file = paste0(output, name, "_", receiver, "_ligand_expr_dotplot", condition_oi, ".pdf"))
  print(DotPlot(subset(seuratObj, cells = cells_to_keep),
    features = rev(best_upstream_ligands), cols = "RdYlBu"
  ) +
    coord_flip() + scale_y_discrete(position = "right") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)))
  dev.off()
}

ligand_lineplot <- function(ligand_activities_all, potential_ligands_focused, output, name, receiver) {
  pdf(file = paste0(output, name, "_", receiver, "_ligand_activities_lineplot", condition_oi, ".pdf"))
  print(make_line_plot(
    ligand_activities = ligand_activities_all,
    potential_ligands = potential_ligands_focused
  ))
  dev.off()
}

chord_diagram <- function(sender_celltypes, circos_links, plot_name, output, name, receiver) {
  # Check for any cell types in circos_links that are missing from sender_celltypes
  extra_types <- setdiff(unique(circos_links$ligand_type), sender_celltypes)
  # Remove "General" from extra_types if it's there (we will add it specifically)
  extra_types <- extra_types[extra_types != "General"]

  # Construct list: senders + extra + General
  sender_celltypes_gen <- c(sender_celltypes, extra_types, "General")
  # Assign colors
  colors <- viridis(length(sender_celltypes_gen), option = "D")
  ligand_colors <- setNames(colors, sender_celltypes_gen)
  target_colors <- setNames("#999999", receiver)
  # prepare circos object
  vis_circos_obj <- prepare_circos_visualization(circos_links,
    ligand_colors = ligand_colors, target_colors = target_colors
  )
  # make legend
  par(bg = "transparent")
  # Default celltype order
  celltype_order <- unique(circos_links$ligand_type) %>%
    sort() %>%
    .[. != "General"] %>%
    c(., "General")
  # Create legend
  circos_legend <- ComplexHeatmap::Legend(
    labels = celltype_order,
    background = ligand_colors[celltype_order],
    type = "point",
    grid_height = unit(3, "mm"),
    grid_width = unit(3, "mm"),
    labels_gp = grid::gpar(fontsize = 8)
  )
  circos_legend_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(circos_legend))

  # combine plot with legend and make pdf
  pdf(file = paste0(output, name, "_", receiver, "_", plot_name, "_chord", condition_oi, ".pdf"))
  print(cowplot::plot_grid(
    make_circos_plot(vis_circos_obj,
      transparency = TRUE,
      link.visible = TRUE, args.circos.text = list(cex = 0.5)
    ),
    circos_legend_grob,
    rel_widths = c(1, 0.1)
  )) # combine chord plot and legend
  dev.off()
}

mushroom_plot <- function(prioritized_table, output, name, receiver) {
  mp <- make_mushroom_plot(prioritized_table, top_n = 30, show_all_datapoints = TRUE, show_rankings = TRUE, condition_oi = condition_oi)
  pdf(file = paste0(output, name, "_", receiver, "_ligand-prioritize_mushroom_plot", condition_oi, ".pdf"))
  print(mp + theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    legend.position = "bottom"
  ))
  dev.off()
}

# make signaling network maps (ligand -> receptor -> target TF -> target)
sig_network_maps <- function(order_ligands, order_targets, weighted_networks,
                             ligand_tf_matrix, output, name, receiver) {
  # Signaling graph depicting the key mediators between a ligand-target pair
  ligands_oi <- unique(as.vector(order_ligands))
  print("maps for ligands: ")
  print(ligands_oi)
  targets_oi <- unique(as.vector(order_targets))
  active_signaling_network <- get_ligand_signaling_path(
    ligands_all = ligands_oi,
    targets_all = targets_oi, weighted_networks = weighted_networks,
    ligand_tf_matrix = ligand_tf_matrix, top_n_regulators = 4, minmax_scaling = TRUE
  )
  # DiagrammeR signaling graph
  signaling_graph <- diagrammer_format_signaling_graph(
    signaling_graph_list = active_signaling_network,
    ligands_all = ligands_oi, targets_all = targets_oi, sig_color = "indianred", gr_color = "steelblue"
  )
  tryCatch(
    {
      signaling_graph %>%
        export_graph(
          file_name = paste0(output, name, "_", receiver, "_", "signaling_graph.pdf"),
          title = paste0("Ligand-Target Diagram for ", receiver)
        )
    },
    error = function(msg) {
      print("Error, cannot make signaling graph")
      return(NA)
    }
  )

  ## loop to graph all ligands with targets
  for (i in seq_along(ligands_oi)) {
    ligand <- ligands_oi[i]
    print(ligand)
    active_signaling_network_ligand <- get_ligand_signaling_path(
      ligands_all = ligand,
      targets_all = targets_oi, weighted_networks = weighted_networks,
      ligand_tf_matrix = ligand_tf_matrix, top_n_regulators = 4, minmax_scaling = TRUE
    )
    signaling_graph_lig <- diagrammer_format_signaling_graph(
      signaling_graph_list = active_signaling_network_ligand,
      ligands_all = ligand, targets_all = targets_oi, sig_color = "indianred", gr_color = "steelblue"
    )
    tryCatch(
      {
        signaling_graph_lig %>%
          export_graph(
            file_name = paste0(output, name, "_", receiver, "_", ligand, "_signaling_graph", condition_oi, ".pdf"),
            title = paste0("Ligand-Target Diagram for ", receiver, "_", ligand)
          )
      },
      error = function(msg) {
        print(paste0("Error, cannot make signaling graph for ligand ", ligand))
        return(NA)
      }
    )
  }
}

# Optional diagnostic (run_random_forest_diagnostics config flag): does the
# given gene more likely belong to the gene set of interest or to background
# genes, using regulatory potential scores of the top 30 predicted genes as
# predictors? Off by default -- mainly useful as a follow-up check when a
# receiver has no active_ligand_target_links.
run_rf_diagnostics_for_receiver <- function(geneset_oi, background_expressed_genes,
                                             best_upstream_ligands, ligand_target_matrix,
                                             output, name, receiver) {
  k <- 3
  n <- 10
  predictions_list <- lapply(1:n, assess_rf_class_probabilities,
    folds = k,
    geneset = geneset_oi, background_expressed_genes = background_expressed_genes,
    ligands_oi = best_upstream_ligands, ligand_target_matrix = ligand_target_matrix
  )
  performances_cv <- bind_rows(lapply(predictions_list, classification_evaluation_continuous_pred_wrapper))
  print(colMeans(performances_cv))
  fraction_cv <- bind_rows(lapply(predictions_list, calculate_fraction_top_predicted,
    quantile_cutoff = 0.95
  ), .id = "round")
  print(mean(filter(fraction_cv, true_target)$fraction_positive_predicted))
  print(mean(filter(fraction_cv, !true_target)$fraction_positive_predicted))
  print(lapply(predictions_list, calculate_fraction_top_predicted_fisher, quantile_cutoff = 0.95))
  top_predicted_genes <- lapply(1:n, get_top_predicted_genes, predictions_list)
  top_predicted_genes <- reduce(top_predicted_genes, full_join, by = c("gene", "true_target"))
  write.table(as.data.frame(top_predicted_genes), file = paste0(
    output, name,
    "_", receiver, "_top_predicted_target_genes.txt"
  ), sep = "\t", quote = FALSE)
}

############### main per-receiver loop ###############
for (i in seq_along(receivers)) {
  receiver <- receivers[i]
  print(receiver)

  ## Find potential ligands
  expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = expression_pct)
  all_receptors <- unique(lr_network$to)
  expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
  potential_ligands <- lr_network[lr_network$to %in% expressed_receptors, ]
  potential_ligands <- unique(potential_ligands$from)
  # filter by sender celltypes -- potential ligands are expressed here
  sender_celltypes <- as.vector(unique(seuratObj@meta.data[[celltype_colname]]))
  sender_celltypes <- sender_celltypes[sender_celltypes != receiver]
  list_expressed_genes_sender <- lapply(sender_celltypes, function(celltype) {
    get_expressed_genes(celltype, seuratObj, pct = expression_pct)
  })
  expressed_genes_sender <- unique(unlist(list_expressed_genes_sender))
  potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)
  print(potential_ligands_focused)

  ## DE for receiver group between conditions
  seurat_obj_receiver <- subset(seuratObj, idents = receiver)
  DE_table_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = condition_oi, ident.2 = condition_reference,
    group.by = condition_colname, min.pct = expression_pct
  )

  # find targets as genes that are DE and in the target matrix, widening the
  # threshold slightly if nothing passes at the configured cutoff (same
  # relative widening as the original notebook: padj +0.01, |lfc| -0.05)
  geneset_oi <- DE_table_receiver[DE_table_receiver$p_val_adj <= padj_thresh &
    abs(DE_table_receiver$avg_log2FC) >= lfc_thresh, ]
  geneset_oi <- rownames(geneset_oi)[rownames(geneset_oi) %in%
    rownames(ligand_target_matrix)]
  if (length(geneset_oi) == 0) {
    geneset_oi <- DE_table_receiver[DE_table_receiver$p_val_adj <= padj_thresh + 0.01 &
      abs(DE_table_receiver$avg_log2FC) >= max(lfc_thresh - 0.05, 0), ]
    geneset_oi <- rownames(geneset_oi)[rownames(geneset_oi) %in%
      rownames(ligand_target_matrix)]
  }
  print(paste0("DE genes: ", paste(geneset_oi, collapse = ", ")))
  # background genes: all genes expressed in the receiver that are also in
  # the ligand-target matrix
  background_expressed_genes <- expressed_genes_receiver[expressed_genes_receiver %in%
    rownames(ligand_target_matrix)]

  ## Ligand activity analysis
  suppressWarnings({
    ligand_activities <- predict_ligand_activities(
      geneset = geneset_oi,
      background_expressed_genes = background_expressed_genes,
      ligand_target_matrix = ligand_target_matrix,
      potential_ligands = potential_ligands
    )
  })
  ligand_activities <- ligand_activities[order(ligand_activities$aupr_corrected,
    decreasing = TRUE
  ), ]
  ligand_activities_all <- ligand_activities
  ligand_activities <- ligand_activities[ligand_activities$test_ligand %in% potential_ligands_focused, ]

  best_upstream_ligands <- top_n(ligand_activities, 50, aupr_corrected)$test_ligand
  print(best_upstream_ligands)
  ligand_activities_df <- as.data.frame(ligand_activities)
  write.table(ligand_activities_df,
    file = paste0(output, name, "_", receiver, "_ligand_activities", condition_oi, ".txt"),
    sep = "\t", quote = FALSE
  )

  ## Ligand-Target
  active_ligand_target_links_df <- lapply(best_upstream_ligands, get_weighted_ligand_target_links,
    geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200
  )
  active_ligand_target_links_df <- drop_na(bind_rows(active_ligand_target_links_df))
  write.table(as.data.frame(active_ligand_target_links_df),
    file = paste0(output, name, "_", receiver, "_ligand-target_links", condition_oi, ".txt"),
    sep = "\t", quote = FALSE
  )

  ## Ligand-Receptor
  ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
    best_upstream_ligands,
    expressed_receptors, lr_network, weighted_networks$lr_sig
  )
  write.table(as.data.frame(ligand_receptor_links_df),
    file = paste0(output, name, "_", receiver, "_ligand-receptor_links", condition_oi, ".txt"),
    sep = "\t", quote = FALSE
  )

  ###### Visualizations ######

  ## heatmap of ligand activity
  ligand_aupr_matrix <- column_to_rownames(ligand_activities, "test_ligand")
  ligand_aupr_matrix <- ligand_aupr_matrix[rev(best_upstream_ligands),
    "aupr_corrected",
    drop = FALSE
  ]
  vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1)
  lig_activity_heatmap(vis_ligand_aupr, output, name, receiver)

  ## Heatmap of ligand-target regulatory potential
  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix, cutoff = 0.05
  )
  if (all(is.na(active_ligand_target_links)) == FALSE) {
    order_ligands <- rev(intersect(best_upstream_ligands, colnames(active_ligand_target_links)))
    order_targets <- intersect(unique(active_ligand_target_links_df$target), rownames(active_ligand_target_links))
    vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])
    write.table(as.data.frame(vis_ligand_target), file = paste0(
      output, name,
      "_", receiver, "_ligand-target_reg_potential", condition_oi, ".txt"
    ), sep = "\t", quote = FALSE)
    lig_target_heatmap(vis_ligand_target, output, name, receiver)
    ## Chord diagram of ligand-target interactions
    ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj,
      best_upstream_ligands[1:min(20, length(best_upstream_ligands))],
      celltype_col = celltype_colname
    )
    active_ligand_target_links_df$target_type <- receiver
    circos_links <- get_ligand_target_links_oi(ligand_type_indication_df,
      active_ligand_target_links_df,
      cutoff = 0.40
    )
    chord_diagram(
      sender_celltypes, circos_links, "ligand-target_interactions",
      output, name, receiver
    )
    print("Making signaling network maps")
    sig_network_maps(
      order_ligands, order_targets, weighted_networks,
      ligand_tf_matrix, output, name, receiver
    )
    if (run_rf_diagnostics) {
      print("Running random forest diagnostics")
      run_rf_diagnostics_for_receiver(
        geneset_oi, background_expressed_genes, best_upstream_ligands,
        ligand_target_matrix, output, name, receiver
      )
    }
  } else {
    print(paste0(receiver, ": No active_ligand_target_links"))
    ligand_type_indication_df <- assign_ligands_to_celltype(seuratObj,
      best_upstream_ligands[1:min(20, length(best_upstream_ligands))],
      celltype_col = celltype_colname
    )
  }

  ## Heatmap of LR interaction potential
  vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(ligand_receptor_links_df,
    best_upstream_ligands,
    order_hclust = "receptors"
  )
  LR_heatmap(vis_ligand_receptor_network, output, name, receiver)

  ## Heatmap of ligand log-fold change across conditions
  celltype_order <- levels(Idents(seuratObj))
  DE_table_top_ligands <- lapply(celltype_order[celltype_order %in% sender_celltypes], function(celltype) {
    tryCatch(
      {
        get_lfc_celltype(celltype,
          seurat_obj = seuratObj, condition_colname = condition_colname,
          condition_oi = condition_oi, condition_reference = condition_reference,
          celltype_col = celltype_colname, min.pct = 0, logfc.threshold = 0,
          features = best_upstream_ligands
        )
      },
      error = function(e) {
        print(paste0("Error in get_lfc_celltype for ", celltype, ": ", e$message))
        return(NULL)
      }
    )
  })
  DE_table_top_ligands <- DE_table_top_ligands[!sapply(DE_table_top_ligands, is.null)]
  DE_table_top_ligands <- reduce(DE_table_top_ligands, full_join)
  DE_table_top_ligands <- column_to_rownames(DE_table_top_ligands, "gene")
  vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), ])
  celltype_colnames <- colnames(vis_ligand_lfc)
  vis_ligand_lfc <- as.data.frame(vis_ligand_lfc) %>%
    arrange(desc(celltype_colnames[1]), celltype_colnames[2])
  vis_ligand_lfc <- as.matrix(vis_ligand_lfc)
  ligand_lfc_heatmap(vis_ligand_lfc, output, name, receiver)

  ## Dot plot of ligand expression and fraction of expressed cells
  ligand_dotplot(seuratObj, best_upstream_ligands, sender_celltypes, output, name, receiver)

  ## Lineplot comparing sender-agnostic vs sender-focused ligand rankings
  ligand_lineplot(ligand_activities_all, potential_ligands_focused, output, name, receiver)

  ## Chord diagram of LR interactions
  lr_network_top_df <- rename(ligand_receptor_links_df, ligand = from, target = to)
  lr_network_top_df$target_type <- receiver
  lr_network_top_df <- inner_join(lr_network_top_df, ligand_type_indication_df)
  chord_diagram(
    sender_celltypes, lr_network_top_df, "ligand-receptor_interactions",
    output, name, receiver
  )

  ## Prioritization of LR pairs
  lr_network_filtered <- filter(lr_network, from %in%
    potential_ligands_focused & to %in% expressed_receptors)[
    ,
    c("from", "to")
  ]
  info_tables <- generate_info_tables(seuratObj,
    celltype_colname = celltype_colname, senders_oi = sender_celltypes, receivers_oi = receiver,
    lr_network = lr_network_filtered, condition_colname = condition_colname,
    condition_oi = condition_oi, condition_reference = condition_reference,
    scenario = "case_control"
  )
  processed_DE_table <- info_tables$sender_receiver_de
  processed_expr_table <- info_tables$sender_receiver_info
  processed_condition_markers <- info_tables$lr_condition_de

  prioritized_table <- generate_prioritization_tables(
    sender_receiver_info =
      processed_expr_table, sender_receiver_de = processed_DE_table,
    ligand_activities = ligand_activities, lr_condition_de =
      processed_condition_markers, scenario = "case_control"
  )
  prioritized_table$sender <- factor(prioritized_table$sender,
    levels = sender_celltypes
  )
  write.table(as.data.frame(prioritized_table), file = paste0(
    output, name,
    "_", receiver, "_ligand-receptor_prioritized_table", condition_oi, ".txt"
  ), sep = "\t", quote = FALSE)

  mushroom_plot(prioritized_table, output, name, receiver)
}

############### prioritize LR pairs across all receiver cell types ###############
nichenet_outputs <- lapply(receivers, function(receiver_ct) {
  tryCatch(
    {
      out <- nichenet_seuratobj_aggregate(
        receiver = receiver_ct,
        seurat_obj = seuratObj,
        condition_colname = condition_colname,
        condition_oi = condition_oi,
        condition_reference = condition_reference,
        sender = as.vector(unique(seuratObj@meta.data[[celltype_colname]]))[
          as.vector(unique(seuratObj@meta.data[[celltype_colname]])) != receiver_ct
        ],
        ligand_target_matrix = ligand_target_matrix,
        lr_network = lr_network,
        weighted_networks = weighted_networks,
        expression_pct = expression_pct
      )
      out$ligand_activities$receiver <- receiver_ct
      return(out)
    },
    error = function(e) {
      print(paste("Error processing", receiver_ct, ":", e$message))
      return(NULL)
    }
  )
})
nichenet_outputs <- nichenet_outputs[!sapply(nichenet_outputs, is.null)]

if (length(nichenet_outputs) > 1) {
  sender_celltypes_all <- as.vector(unique(seuratObj@meta.data[[celltype_colname]]))
  info_tables <- lapply(nichenet_outputs, function(out) {
    lr_network_filtered <- filter(
      lr_network[, c("from", "to")],
      from %in% out$ligand_activities$test_ligand & to %in% out$background_expressed_genes
    )
    generate_info_tables(seuratObj,
      celltype_colname = celltype_colname,
      senders_oi = sender_celltypes_all,
      receivers_oi = unique(out$ligand_activities$receiver),
      lr_network_filtered = lr_network_filtered,
      condition_colname = condition_colname,
      condition_oi = condition_oi,
      condition_reference = condition_reference,
      scenario = "case_control"
    )
  })
  info_tables_combined <- purrr::pmap(info_tables, bind_rows)
  ligand_activities_combined <- purrr::map_dfr(nichenet_outputs, "ligand_activities")
  prior_table_combined <- generate_prioritization_tables(
    sender_receiver_info = distinct(info_tables_combined$sender_receiver_info),
    sender_receiver_de = info_tables_combined$sender_receiver_de,
    ligand_activities = ligand_activities_combined,
    lr_condition_de = distinct(info_tables_combined$lr_condition_de),
    scenario = "case_control"
  )
  mp <- make_mushroom_plot(prior_table_combined, top_n = 30, show_all_datapoints = TRUE, show_rankings = TRUE)
  pdf(file = paste0(output, name, "_", "ligand-prioritize_combined_mushroom_plot.pdf"))
  print(mp + theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    legend.position = "bottom"
  ))
  dev.off()
} else {
  print("Fewer than 2 receivers produced output -- skipping combined multi-receiver prioritization.")
}

writeLines(capture.output(sessionInfo()), paste0(output, "sessionInfo.txt"))
system(paste("chmod -R 777", output))
