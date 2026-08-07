# Cell-cell communication interaction comparison between two conditions,
# using two previously-computed CellChat objects (one per condition, from
# cellchat.R).
# Configured via the cellchat_interactions section of config.json.
# Ported from src/cellchat_interactions.Rmd (kept as a reference notebook --
# this script is the configurable, non-interactive version run via
# run_downstream_toolkit.sh). Original notebook author: Beth Moore.

library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(circlize)

# Record R session info (packages + versions) for run reproducibility --
# appended to the run's shared provenance file when invoked via
# run_downstream_toolkit.sh (PROVENANCE_FILE env var), else written
# standalone to the current directory.
.provenance_file <- Sys.getenv("PROVENANCE_FILE", unset = paste0("./sessionInfo_cellchat_interactions.R_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
cat(paste0("\n--- R sessionInfo (cellchat_interactions.R) ---\n", paste(capture.output(sessionInfo()), collapse = "\n"), "\n"), file = .provenance_file, append = TRUE)

### load config ###
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(GIT_DIR, "config.json"))
cfg <- config$cellchat_interactions
docker <- config$docker

if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
  DATA_DIR <- "./data/input_data/"
} else {
  DATA_DIR <- cfg$DATA_DIR
}
# Resolve to an absolute path now, before setwd(output) below changes the
# working directory -- otherwise this relative path (docker mode) would
# resolve against the wrong directory once we've cd'd into the output dir.
# normalizePath() strips the trailing slash, so add it back -- downstream
# code joins DATA_DIR with object1_path/object2_path via paste0(), which
# assumes one is present.
DATA_DIR <- paste0(normalizePath(DATA_DIR), "/")

name <- cfg$name
object1_label <- cfg$object1_label
object2_label <- cfg$object2_label
pos.dataset <- cfg$positive_dataset
neg.dataset <- setdiff(c(object1_label, object2_label), pos.dataset)[1]
group1 <- cfg$group1
group1_label <- cfg$group1_label
group2 <- cfg$group2
group2_label <- cfg$group2_label

# reticulate is required by CellChat's functional/structural similarity
# analysis (computeNetSimilarityPairwise -> UMAP). Point it at this
# environment's own Python (installed alongside CellChat in the Dockerfile's
# `cellchat` conda env) rather than a personal virtualenv path.
reticulate::use_python(Sys.which("python"), required = TRUE)

### set output directory ###
timestamp <- Sys.getenv("RUN_TIMESTAMP", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
output <- paste0("./shared_volume/output_cellchat_interactions_", name, "_", timestamp)
print(output)
dir.create(output, recursive = TRUE, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(file.path(GIT_DIR, "config.json"), file.path(output, "config.json"))
if (file.exists(.provenance_file)) {
  file.copy(.provenance_file, file.path(output, basename(.provenance_file)), overwrite = TRUE)
  file.remove(.provenance_file)
}
# All plots below use relative filenames (matching the original notebook,
# which cd's into its own output dir) -- do the same here.
setwd(output)

### load the two previously-computed CellChat objects ###
ptm <- Sys.time()
cellchat.obj1 <- readRDS(paste0(DATA_DIR, cfg$object1_path))
cellchat.obj2 <- readRDS(paste0(DATA_DIR, cfg$object2_path))
object.list <- setNames(list(cellchat.obj1, cellchat.obj2), c(object1_label, object2_label))
print(table(object.list[[object1_label]]@idents))
print(table(object.list[[object2_label]]@idents))

# not_group1/not_group2: every other cluster in the negative-dataset object
# (mirrors cellchat.R's own not_source1/not_source2 computation)
celltype.list <- unique(object.list[[neg.dataset]]@idents)
not_group1 <- as.vector(celltype.list[celltype.list != group1])
not_group2 <- as.vector(celltype.list[celltype.list != group2])
groups <- list(
  list(id = group1, not_id = not_group1, label = group1_label),
  list(id = group2, not_id = not_group2, label = group2_label)
)

### merge ###
cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))
print(cellchat_merge)
execution.time <- Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

### save merged objects ###
save(object.list, file = paste0("cellchat_object.list_", name, ".RData"))
save(cellchat_merge, file = paste0("cellchat_merged_", name, ".RData"))

### Identify altered interactions and cell populations ###
# comparison = c(2, 1): second dataset (in names(object.list) order) vs. first
gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(2, 1))
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(2, 1), measure = "weight")
pdf(file = paste0(name, "_global_interactions_pos", pos.dataset, ".pdf"))
print(gg1 + gg2)
dev.off()

pdf(file = paste0(name, "_cell_diffinteractions_circ_pos", pos.dataset, ".pdf"))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat_merge, weight.scale = T)
netVisual_diffInteraction(cellchat_merge, weight.scale = T, measure = "weight")
dev.off()

# heatmap: red = increased signal in the second dataset vs. the first
gg1 <- netVisual_heatmap(cellchat_merge, comparison = c(1, 2))
gg2 <- netVisual_heatmap(cellchat_merge, measure = "weight", comparison = c(1, 2))
pdf(file = paste0(name, "_cell_diffinteractions_heatmap_pos", pos.dataset, ".pdf"))
print(gg1 + gg2)
dev.off()

### Compare the major sources and targets ###
num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in seq_along(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
    title = names(object.list)[i], weight.MinMax = weight.MinMax
  )
}
pdf(file = paste0(name, "_sig_changes_each_dataset.pdf"))
print(patchwork::wrap_plots(plots = gg))
dev.off()

### signaling changes for the two highlighted groups ###
# positive values = increase in the second dataset of comparison=c(2,1),
# negative = increase in the first
gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = group1, comparison = c(2, 1))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = group2, comparison = c(2, 1))
pdf(file = paste0(name, "_sig_changes_", group1_label, "-", group2_label, ".pdf"), height = 6, width = 10)
print(patchwork::wrap_plots(plots = list(gg1, gg2)))
dev.off()

### signaling changes for every other cluster ###
celltype.listA <- as.vector(celltype.list[celltype.list != group1 & celltype.list != group2])
for (i in seq_along(celltype.listA)) {
  print(celltype.listA[i])
  gg1 <- netAnalysis_signalingChanges_scatter(cellchat_merge, idents.use = celltype.listA[i], comparison = c(2, 1))
  pdf(file = paste0(name, "_sig_changes_", celltype.listA[i], ".pdf"), height = 6, width = 6)
  print(gg1)
  dev.off()
}

### Identify altered signaling with distinct network architecture ###
cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge, type = "functional")
cellchat_merge <- netEmbedding(cellchat_merge, type = "functional")
cellchat_merge <- netClustering(cellchat_merge, type = "functional")

pdf(file = paste0(name, "_sig_networks_functional_2D_plot.pdf"))
netVisual_embeddingPairwise(cellchat_merge, type = "functional", label.size = 3.5)
dev.off()

pdf(file = paste0(name, "_sig_networks_functional_2D_plot_pathslabel.pdf"))
netVisual_embeddingPairwise(cellchat_merge, type = "functional", label.size = 3.5, top.label = 0.2)
dev.off()

pdf(file = paste0(name, "_sig_networks_functional_2D_plot_zoomin.pdf"))
netVisual_embeddingPairwiseZoomIn(cellchat_merge, type = "functional", nCol = 2)
dev.off()

pdf(file = paste0(name, "_pathway_distance.pdf"))
rankSimilarity(cellchat_merge, type = "functional")
dev.off()

### Identify altered signaling with distinct interaction strength ###
gg1 <- rankNet(cellchat_merge,
  mode = "comparison", measure = "weight",
  sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE, comparison = c(2, 1)
)
gg2 <- rankNet(cellchat_merge,
  mode = "comparison", measure = "weight",
  sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE, comparison = c(2, 1)
)
pdf(file = paste0(name, "_overall_info_flow.pdf"))
print(gg1 + gg2)
dev.off()

### outgoing/incoming/all signaling role heatmaps, per dataset ###
i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i + 1]]@netP$pathways)
ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
  pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i],
  width = 5, height = 6, font.size = 3, font.size.title = 7
)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
  pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i + 1],
  width = 5, height = 6, font.size = 3, font.size.title = 7
)
pdf(file = paste0(name, "_outgoing_sig_paths.pdf"))
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

ht3 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
  pattern = "incoming", signaling = pathway.union, title = names(object.list)[i],
  width = 5, height = 6, color.heatmap = "GnBu", font.size = 3, font.size.title = 7
)
ht4 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
  pattern = "incoming", signaling = pathway.union, title = names(object.list)[i + 1],
  width = 5, height = 6, color.heatmap = "GnBu", font.size = 3, font.size.title = 7
)
pdf(file = paste0(name, "_incoming_sig_paths.pdf"))
draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
dev.off()

ht5 <- netAnalysis_signalingRole_heatmap(object.list[[i]],
  pattern = "all", signaling = pathway.union, title = names(object.list)[i],
  width = 5, height = 6, color.heatmap = "OrRd", font.size = 3, font.size.title = 7
)
ht6 <- netAnalysis_signalingRole_heatmap(object.list[[i + 1]],
  pattern = "all", signaling = pathway.union, title = names(object.list)[i + 1],
  width = 5, height = 6, color.heatmap = "OrRd", font.size = 3, font.size.title = 7
)
pdf(file = paste0(name, "_all_sig_paths.pdf"))
draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))
dev.off()

### Differential outgoing/incoming signaling heatmap (dataset2 vs dataset1) ###
# CellChat's netAnalysis_signalingRole_heatmap() has no built-in "differential"
# mode: it always plots a single object, and internally row-normalizes each
# pathway to its own max, so the two side-by-side heatmaps above are only
# comparable in relative pattern, not raw magnitude, between datasets.
# netVisual_heatmap() *does* support a differential view, but only for the
# overall interaction count/weight network, not the pathway-level centrality
# scores used here.
#
# This reimplements netAnalysis_signalingRole_heatmap's internal
# outgoing/incoming centrality matrix construction, but for two objects
# instead of one, aligns them to the same pathway and cluster sets
# (zero-filling anything missing from one side), and takes the raw
# difference mat2 - mat1. Diverging color scale via colorRamp3(): negative
# (comparison[1] relatively higher) -> blue; positive (comparison[2]
# relatively higher) -> red.
build_diverging_ramp <- function(mat, color.diff, power = 0.5, n.steps = 10) {
  neg.max <- abs(min(min(mat), 0))
  pos.max <- max(max(mat), 0)
  frac <- (0:n.steps) / n.steps

  breaks <- numeric(0)
  colors <- character(0)

  if (neg.max > 0) {
    neg.breaks <- -rev(neg.max * frac^(1 / power))
    neg.colors <- rev(colorRampPalette(c("white", color.diff[1]))(n.steps + 1))
    breaks <- c(breaks, neg.breaks)
    colors <- c(colors, neg.colors)
  }
  if (pos.max > 0) {
    pos.breaks <- pos.max * frac^(1 / power)
    pos.colors <- colorRampPalette(c("white", color.diff[2]))(n.steps + 1)
    if (length(breaks) > 0) {
      pos.breaks <- pos.breaks[-1]
      pos.colors <- pos.colors[-1]
    }
    breaks <- c(breaks, pos.breaks)
    colors <- c(colors, pos.colors)
  }
  colorRamp3(breaks, colors)
}

netAnalysis_signalingRole_heatmap_diff <- function(object.list, comparison = c(1, 2),
                                                     pattern = c("outgoing", "incoming", "all"),
                                                     signaling = NULL, slot.name = "netP",
                                                     color.use = NULL,
                                                     color.diff = c("#2166ac", "#b2182b"),
                                                     contrast.power = 0.5,
                                                     title = NULL, width = 10, height = 8,
                                                     font.size = 8, font.size.title = 10,
                                                     cluster.rows = FALSE, cluster.cols = FALSE) {
  pattern <- match.arg(pattern)

  get_mat <- function(object) {
    if (length(slot(object, slot.name)$centr) == 0) {
      stop("Please run `netAnalysis_computeCentrality` on all objects first!")
    }
    centr <- slot(object, slot.name)$centr
    outgoing <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    incoming <- matrix(0, nrow = nlevels(object@idents), ncol = length(centr))
    dimnames(outgoing) <- list(levels(object@idents), names(centr))
    dimnames(incoming) <- dimnames(outgoing)
    for (j in seq_along(centr)) {
      outgoing[, j] <- centr[[j]]$outdeg
      incoming[, j] <- centr[[j]]$indeg
    }
    if (pattern == "outgoing") {
      mat <- t(outgoing)
    } else if (pattern == "incoming") {
      mat <- t(incoming)
    } else {
      mat <- t(outgoing + incoming)
    }
    mat
  }

  obj1 <- object.list[[comparison[1]]]
  obj2 <- object.list[[comparison[2]]]
  mat1 <- get_mat(obj1)
  mat2 <- get_mat(obj2)

  if (is.null(signaling)) {
    signaling <- union(rownames(mat1), rownames(mat2))
  }
  clusters.union <- union(colnames(mat1), colnames(mat2))

  reindex <- function(mat) {
    out <- matrix(0,
      nrow = length(signaling), ncol = length(clusters.union),
      dimnames = list(signaling, clusters.union)
    )
    rows.keep <- rownames(mat)[rownames(mat) %in% signaling]
    cols.keep <- colnames(mat)[colnames(mat) %in% clusters.union]
    out[rows.keep, cols.keep] <- mat[rows.keep, cols.keep, drop = FALSE]
    out
  }
  mat1 <- reindex(mat1)
  mat2 <- reindex(mat2)

  mat.diff <- mat2 - mat1

  legend.name <- switch(pattern,
    outgoing = "Outgoing",
    incoming = "Incoming",
    all = "Overall"
  )
  if (is.null(title)) {
    title <- paste0(
      "Differential ", tolower(legend.name), " signaling\n(",
      names(object.list)[comparison[2]], " vs ", names(object.list)[comparison[1]], ")"
    )
  }

  if (is.null(color.use)) {
    color.use <- scPalette(ncol(mat.diff))
  }
  names(color.use) <- colnames(mat.diff)

  color.heatmap.use <- build_diverging_ramp(mat.diff, color.diff, power = contrast.power)

  df <- data.frame(group = colnames(mat.diff))
  rownames(df) <- colnames(mat.diff)
  col_annotation <- HeatmapAnnotation(
    df = df, col = list(group = color.use),
    which = "column", show_legend = FALSE, show_annotation_name = FALSE,
    simple_anno_size = grid::unit(0.2, "cm")
  )

  ha2 <- HeatmapAnnotation(Strength = anno_barplot(colSums(mat.diff),
    border = FALSE, gp = gpar(fill = color.use, col = color.use)
  ), show_annotation_name = FALSE)
  ha1 <- rowAnnotation(
    Strength = anno_barplot(rowSums(mat.diff), border = FALSE),
    show_annotation_name = FALSE
  )

  mat.plot <- mat.diff
  mat.plot[mat.plot == 0] <- NA

  Heatmap(mat.plot,
    col = color.heatmap.use, na_col = "white",
    name = paste0(names(object.list)[comparison[2]], " - ", names(object.list)[comparison[1]]),
    bottom_annotation = col_annotation, top_annotation = ha2,
    right_annotation = ha1, cluster_rows = cluster.rows, cluster_columns = cluster.cols,
    row_names_side = "left", row_names_rot = 0,
    row_names_gp = gpar(fontsize = font.size), column_names_gp = gpar(fontsize = font.size),
    width = unit(width, "cm"), height = unit(height, "cm"),
    column_title = title, column_title_gp = gpar(fontsize = font.size.title),
    column_names_rot = 90,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontface = "plain"),
      title_position = "leftcenter-rot", border = NA, legend_height = unit(20, "mm"),
      labels_gp = gpar(fontsize = 8), grid_width = unit(2, "mm")
    )
  )
}

# comparison = c(1, 2) matches names(object.list)[comparison] ==
# c(object1_label, object2_label), i.e. diff = object2 - object1
ht_diff_out <- netAnalysis_signalingRole_heatmap_diff(object.list,
  comparison = c(1, 2), pattern = "outgoing",
  signaling = pathway.union, width = 8, height = 10, font.size = 5, font.size.title = 8
)
pdf(file = paste0(name, "_diff_outgoing_sig_paths.pdf"), width = 6, height = 9)
draw(ht_diff_out)
dev.off()

ht_diff_in <- netAnalysis_signalingRole_heatmap_diff(object.list,
  comparison = c(1, 2), pattern = "incoming",
  signaling = pathway.union, width = 8, height = 10, font.size = 5, font.size.title = 8
)
pdf(file = paste0(name, "_diff_incoming_sig_paths.pdf"), width = 6, height = 9)
draw(ht_diff_in)
dev.off()

### Up-regulated and down-regulated signaling ligand-receptor pairs ###
# bubble plots of overall communication probability for the two highlighted groups
for (g in groups) {
  p1 <- netVisual_bubble(cellchat_merge, sources.use = g$id, targets.use = g$not_id, comparison = c(2, 1), angle.x = 45)
  p2 <- netVisual_bubble(cellchat_merge, sources.use = g$not_id, targets.use = g$id, comparison = c(2, 1), angle.x = 45)
  pdf(file = paste0(name, "_sig_cc_L-Rpairs_celltypes_bubble_", g$label, ".pdf"), height = 11, width = 8.5)
  print(p1 + p2)
  dev.off()
}

# increased/decreased signaling to/from each highlighted group, based on
# communication probability between datasets
group_bubbles <- list()
for (g in groups) {
  gg_inc_to <- netVisual_bubble(cellchat_merge,
    sources.use = g$not_id, targets.use = g$id, comparison = c(2, 1), max.dataset = 2,
    title.name = paste0("Increased signaling in ", pos.dataset, " to ", g$label), angle.x = 45, remove.isolate = T
  )
  gg_dec_to <- netVisual_bubble(cellchat_merge,
    sources.use = g$not_id, targets.use = g$id, comparison = c(2, 1), max.dataset = 1,
    title.name = paste0("Decreased signaling in ", pos.dataset, " to ", g$label), angle.x = 45, remove.isolate = T
  )
  gg_inc_from <- netVisual_bubble(cellchat_merge,
    sources.use = g$id, targets.use = g$not_id, comparison = c(2, 1), max.dataset = 2,
    title.name = paste0("Increased signaling in ", pos.dataset, " from ", g$label), angle.x = 45, remove.isolate = T
  )
  gg_dec_from <- netVisual_bubble(cellchat_merge,
    sources.use = g$id, targets.use = g$not_id, comparison = c(2, 1), max.dataset = 1,
    title.name = paste0("Decreased signaling in ", pos.dataset, " from ", g$label), angle.x = 45, remove.isolate = T
  )
  pdf(file = paste0(name, "_inc-dec_sig_cc_L-Rpairs_celltypes_bubble_", g$label, "s.pdf"), height = 11, width = 8.5)
  print(gg_inc_to + gg_dec_to + gg_inc_from + gg_dec_from)
  dev.off()
  group_bubbles[[g$label]] <- list(inc_to = gg_inc_to, dec_to = gg_dec_to, inc_from = gg_inc_from, dec_from = gg_dec_from)

  write.table(as.data.frame(gg_inc_to$data), file = paste0(name, "_inc", pos.dataset, "_sig_L-Rpairs_", g$label, "s_target.txt"), sep = "\t", quote = FALSE)
  write.table(as.data.frame(gg_dec_to$data), file = paste0(name, "_dec", pos.dataset, "_sig_L-Rpairs_", g$label, "s_target.txt"), sep = "\t", quote = FALSE)
  write.table(as.data.frame(gg_inc_from$data), file = paste0(name, "_inc", pos.dataset, "_sig_L-Rpairs_", g$label, "s_source.txt"), sep = "\t", quote = FALSE)
  write.table(as.data.frame(gg_dec_from$data), file = paste0(name, "_dec", pos.dataset, "_sig_L-Rpairs_", g$label, "s_source.txt"), sep = "\t", quote = FALSE)
}

### Identify dysfunctional signaling via differential expression ###
features.name <- paste0(pos.dataset, ".merged")
cellchat_merge <- identifyOverExpressedGenes(cellchat_merge,
  group.dataset = "datasets",
  pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE,
  thresh.pc = 0.1, thresh.fc = 0.05, thresh.p = 0.05, group.DE.combined = FALSE
)
net <- netMappingDEG(cellchat_merge, features.name = features.name, variable.all = TRUE)
net.up <- subsetCommunication(cellchat_merge, net = net, datasets = pos.dataset, ligand.logFC = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat_merge, net = net, datasets = neg.dataset, ligand.logFC = -0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat_merge)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat_merge)

### Visualize up/down-regulated L-R pairs for each highlighted group ###
pairLR.use.up <- net.up[, "interaction_name", drop = F]
pairLR.use.down <- net.down[, "interaction_name", drop = F]
for (g in groups) {
  # netVisual_bubble errors out (rather than returning an empty plot) when a
  # gene subset has zero matching interactions for this source/target pair --
  # a legitimate outcome for real data, not a bug, so skip that half of the
  # plot instead of halting the whole script (same defensive pattern as the
  # per-pathway circle plots below).
  gg1 <- tryCatch(
    netVisual_bubble(cellchat_merge,
      pairLR.use = pairLR.use.up, sources.use = g$id,
      targets.use = g$not_id, comparison = c(2, 1), angle.x = 90, remove.isolate = T,
      title.name = paste0("Up-regulated signaling in ", object2_label)
    ),
    error = function(e) {
      message("No up-regulated interactions for ", g$label, " -- skipping: ", conditionMessage(e))
      NULL
    }
  )
  gg2 <- tryCatch(
    netVisual_bubble(cellchat_merge,
      pairLR.use = pairLR.use.down, sources.use = g$id,
      targets.use = g$not_id, comparison = c(2, 1), angle.x = 90, remove.isolate = T,
      title.name = paste0("Down-regulated signaling in ", object2_label)
    ),
    error = function(e) {
      message("No down-regulated interactions for ", g$label, " -- skipping: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(gg1) && is.null(gg2)) next
  pdf(file = paste0(name, "_up-down_reg_cc_L-Rpairs_celltypes_bubble_", g$label, ".pdf"))
  print(if (!is.null(gg1) && !is.null(gg2)) gg1 + gg2 else if (!is.null(gg1)) gg1 else gg2)
  dev.off()
}

### Chord diagrams of up/down-regulated L-R pairs (top N by probability) ###
# netVisual_chord_gene can still throw circlize's "gap.degree is too large"
# error even with >0 interactions, when the filtered set has too few sectors
# for the configured gap sizes to fit -- a circlize plotting limitation, not
# a bug in the data. All 6 call sites below go through this wrapper so a
# single sparse slice doesn't halt every other chord diagram in the script.
safe_chord_gene <- function(...) {
  tryCatch(
    netVisual_chord_gene(...),
    error = function(e) {
      message("Skipping chord diagram -- ", conditionMessage(e))
      NULL
    }
  )
}

# Wrapper: pre-filter a net data frame to the top N interactions by
# probability for the given source/target combination, then draw the chord
# diagram. Avoids circlize's gap.degree error when the full net has too many
# L-R pairs.
chord_top_net <- function(obj, sources, targets, net_df, n = 50, ...) {
  sub <- net_df[net_df$source %in% sources & net_df$target %in% targets, ]
  if (nrow(sub) == 0) {
    message("No interactions for this source/target combination -- skipping")
    return(invisible(NULL))
  }
  net_top <- sub[order(sub$prob, decreasing = TRUE)[1:min(n, nrow(sub))], ]
  safe_chord_gene(obj,
    sources.use = sources, targets.use = targets,
    slot.name = "net", net = net_top, ...
  )
}

for (g in groups) {
  pdf(file = paste0(name, "_up-down_reg_cc_L-Rpairs_celltypes_chord_", g$label, ".pdf"))
  par(mfrow = c(2, 2), xpd = TRUE)
  chord_top_net(object.list[[2]], g$id, g$not_id, net.up, n = 50, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", object2_label))
  chord_top_net(object.list[[1]], g$id, g$not_id, net.down, n = 50, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", object2_label))
  chord_top_net(object.list[[2]], g$not_id, g$id, net.up, n = 50, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", object2_label))
  chord_top_net(object.list[[1]], g$not_id, g$id, net.down, n = 50, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", object2_label))
  dev.off()
}

### Visually compare cell-cell communication for each highlighted group ###
for (g in groups) {
  pdf(file = paste0(name, "_signal_cc_L-Rpairs_celltypes_chord_from", g$label, ".pdf"), height = 10, width = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    safe_chord_gene(object.list[[i]],
      sources.use = g$id, targets.use = g$not_id,
      lab.cex = 0.4, title.name = paste0("Signaling from ", g$label, " - ", names(object.list)[i])
    )
  }
  dev.off()

  pdf(file = paste0(name, "_signal_cc_L-Rpairs_celltypes_chord_to", g$label, ".pdf"), height = 10, width = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    safe_chord_gene(object.list[[i]],
      sources.use = g$not_id, targets.use = g$id,
      lab.cex = 0.4, title.name = paste0("Signaling to ", g$label, " - ", names(object.list)[i])
    )
  }
  dev.off()
}

### Top-25 L-R chord diagrams (probability-filtered, to avoid overcrowding) ###
# CellChat's pairLR.use validation has an operator precedence bug -- it
# rejects a data frame with only interaction_name. Workaround: zero out
# non-top interactions in a copy of the prob array so pairLR.use isn't
# needed at all.
filter_top_LR <- function(obj, sources, targets, n = 25) {
  prob <- obj@net$prob
  lr_strength <- apply(prob[sources, targets, , drop = FALSE], 3, sum)
  n_keep <- min(n, sum(lr_strength > 0))
  if (n_keep == 0) {
    return(obj)
  }
  keep_idx <- order(lr_strength, decreasing = TRUE)[1:n_keep]
  zero_idx <- setdiff(seq_len(dim(prob)[3]), keep_idx)
  obj@net$prob[, , zero_idx] <- 0
  obj
}

for (g in groups) {
  pdf(file = paste0(name, "_signal_cc_L-Rpairs_celltypes_chord_from", g$label, "_top25.pdf"), height = 10, width = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    obj_f <- filter_top_LR(object.list[[i]], g$id, g$not_id, n = 25)
    safe_chord_gene(obj_f,
      sources.use = g$id, targets.use = g$not_id,
      lab.cex = 0.2, small.gap = 0.5, big.gap = 5,
      title.name = paste0("Top 25 L-R from ", g$label, " - ", names(object.list)[i])
    )
  }
  dev.off()

  pdf(file = paste0(name, "_signal_cc_L-Rpairs_celltypes_chord_to", g$label, "_top25.pdf"), height = 10, width = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    obj_f <- filter_top_LR(object.list[[i]], g$not_id, g$id, n = 25)
    safe_chord_gene(obj_f,
      sources.use = g$not_id, targets.use = g$id,
      lab.cex = 0.2, small.gap = 0.5, big.gap = 5,
      title.name = paste0("Top 25 L-R to ", g$label, " - ", names(object.list)[i])
    )
  }
  dev.off()
}

### show all significant signaling pathways to/from each highlighted group ###
for (g in groups) {
  pdf(file = paste0(name, "_signal_cc_paths_celltypes_chord_from", g$label, ".pdf"), height = 10, width = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    safe_chord_gene(object.list[[i]],
      sources.use = g$id, targets.use = g$not_id,
      slot.name = "netP", title.name = paste0("Signaling pathways sending from ", g$label, " - ", names(object.list)[i]),
      legend.pos.x = 10, lab.cex = 0.2, small.gap = 0.5, big.gap = 5
    )
  }
  dev.off()

  pdf(file = paste0(name, "_signal_cc_paths_celltypes_chord_to", g$label, ".pdf"), height = 10, width = 10)
  par(mfrow = c(1, 2), xpd = TRUE)
  for (i in seq_along(object.list)) {
    safe_chord_gene(object.list[[i]],
      sources.use = g$not_id, targets.use = g$id,
      slot.name = "netP", title.name = paste0("Signaling pathways sending to ", g$label, " - ", names(object.list)[i]),
      legend.pos.x = 10, lab.cex = 0.2, small.gap = 0.5, big.gap = 5
    )
  }
  dev.off()
}

### detailed circle plots for every significant pathway ###
pathways.show.all <- unique(c(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways))
for (j in seq_along(pathways.show.all)) {
  print(pathways.show.all[j])
  tryCatch(
    {
      weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show.all[j])
      pdf(file = paste0(pathways.show.all[j], "_path_sig_networks_circle2.pdf"))
      par(mfrow = c(1, 2), xpd = TRUE)
      for (i in seq_along(object.list)) {
        groupSize <- as.numeric(table(object.list[[i]]@idents))
        netVisual_aggregate(object.list[[i]],
          signaling = pathways.show.all[j],
          layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, vertex.weight = groupSize,
          signaling.name = paste(pathways.show.all[j], names(object.list)[i])
        )
      }
      dev.off()
    },
    error = function(e) {
      message("weight max error: ", conditionMessage(e))
      tryCatch(
        {
          groupSize <- as.numeric(table(object.list[[1]]@idents))
          pdf(file = paste0(pathways.show.all[j], "_path_sig_networks_circle2.pdf"))
          netVisual_aggregate(object.list[[1]],
            signaling = pathways.show.all[j],
            layout = "circle", vertex.weight = groupSize,
            signaling.name = paste(pathways.show.all[j], names(object.list)[1])
          )
          dev.off()
        },
        error = function(e) {
          message("fallback (dataset 1) also failed: ", conditionMessage(e))
          groupSize <- as.numeric(table(object.list[[2]]@idents))
          pdf(file = paste0(pathways.show.all[j], "_path_sig_networks_circle2.pdf"))
          netVisual_aggregate(object.list[[2]],
            signaling = pathways.show.all[j],
            layout = "circle", vertex.weight = groupSize,
            signaling.name = paste(pathways.show.all[j], names(object.list)[2])
          )
          dev.off()
        }
      )
    }
  )
}

### aggregated communication weights per cell type ###
df.net1 <- as.data.frame(object.list[[1]]@net$weight)
df.net2 <- as.data.frame(object.list[[2]]@net$weight)
write.table(df.net1, file = paste0(name, "_netaggregate", names(object.list)[1], "_weights.txt"), quote = FALSE, row.names = FALSE)
write.table(df.net2, file = paste0(name, "_netaggregate", names(object.list)[2], "_weights.txt"), quote = FALSE, row.names = FALSE)

### pathway-level communication, both datasets combined ###
df.net.p1 <- subsetCommunication(object.list[[1]], slot.name = "netP")
df.net.p2 <- subsetCommunication(object.list[[2]], slot.name = "netP")
write.table(df.net.p1, file = paste0(name, "_", names(object.list)[1], "_cellchat_df_net_signal_paths.txt"), sep = "\t", quote = FALSE)
write.table(df.net.p2, file = paste0(name, "_", names(object.list)[2], "_cellchat_df_net_signal_paths.txt"), sep = "\t", quote = FALSE)

df.net.p.merged <- merge(df.net.p1, df.net.p2,
  by = c("source", "target", "pathway_name"), all = TRUE,
  suffixes = paste0("_", names(object.list))
)
df.net.p.merged[is.na(df.net.p.merged)] <- 0
prob_col1 <- paste0("prob_", names(object.list)[1])
prob_col2 <- paste0("prob_", names(object.list)[2])
df.net.p.merged$net_prob <- df.net.p.merged[[prob_col1]] - df.net.p.merged[[prob_col2]]
write.table(df.net.p.merged, file = paste0(name, "_cellchat_df_net_signal_paths.txt"), sep = "\t", quote = FALSE)

### signaling gene expression distributions between datasets, per pathway ###
cellchat_merge@meta$datasets <- factor(cellchat_merge@meta$datasets, levels = c(neg.dataset, pos.dataset))
for (j in seq_along(pathways.show.all)) {
  p1 <- plotGeneExpression(cellchat_merge,
    signaling = pathways.show.all[j], split.by = "datasets",
    colors.ggplot = T, type = "violin"
  )
  pdf(file = paste0(pathways.show.all[j], "_geneExpr_violin.pdf"))
  print(p1)
  dev.off()
}

### save ###
save(object.list, file = paste0("cellchat_object.list_", name, ".RData"))
save(cellchat_merge, file = paste0("cellchat_merge_", name, ".RData"))
writeLines(capture.output(sessionInfo()), paste0(name, "_sessionInfo.txt"))
# "." not `output`: `output` is the relative path computed before
# setwd(output) above, so it no longer resolves correctly now that we've
# cd'd into it -- we're already sitting in the directory to chmod.
system("chmod -R 777 .")
