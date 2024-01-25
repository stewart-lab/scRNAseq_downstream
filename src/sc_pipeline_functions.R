read_aligned_data <- function(base_directory, project_name, output) {
  # Existing directories
  filtered_data_directory <- paste0(base_directory, "/filtered/")
  raw_data_directory <- paste0(base_directory, "/raw/")

  # Alignment directory and subdirectories
  alignment_directory <- paste0(output, "/alignment/")
  alignment_raw_directory <- paste0(alignment_directory, "raw/")
  alignment_filtered_directory <- paste0(alignment_directory, "filtered/")

  # Create directories
  dir.create(alignment_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(alignment_raw_directory, recursive = TRUE, showWarnings = FALSE)
  dir.create(alignment_filtered_directory, recursive = TRUE, showWarnings = FALSE)

  # Try to copy Summary.csv to alignment directory
  summary_file_path <- paste0(base_directory, "/Summary.csv")
  alignment_summary_file_path <- paste0(alignment_directory, "Alignment_Summary.csv")

  tryCatch({
    file.copy(summary_file_path, alignment_summary_file_path, overwrite = TRUE)
    message("Summary.csv has been copied successfully to alignment directory.")
  }, warning = function(w) {
    message("Warning: ", conditionMessage(w))
  }, error = function(e) {
    message("Summary.csv could not be copied: ", conditionMessage(e))
  }, finally = {
    message("Continuing with the rest of the function.")
  })

  # Function to copy .gz files from a directory to a target directory
  copy_gz_files_to_directory <- function(source_directory, target_directory) {
    gz_files <- list.files(source_directory, pattern = "\\.gz$", full.names = TRUE)
    sapply(gz_files, function(file) {
      file.copy(file, file.path(target_directory, basename(file)), overwrite = TRUE)
    })
  }

  # Copy .gz files from raw and filtered data directories
  copy_gz_files_to_directory(raw_data_directory, alignment_raw_directory)
  copy_gz_files_to_directory(filtered_data_directory, alignment_filtered_directory)

  list(
    filtered = Read10X(filtered_data_directory),
    raw = Read10X(raw_data_directory),
    project = project_name
  )
}


prep_seurat_and_soupX <- function(data.raw, data, project) {
  dims_umap <- 1:config$prep_seurat_and_soupX$dims
  umap.method <- config$prep_seurat_and_soupX$umap.method
  tfidfMin <- config$prep_seurat_and_soupX$tfidfMin
  # Create SoupChannel object
  sc <- SoupChannel(data.raw, data)

  # Remove 'data.raw' as it's no longer needed.
  rm(data.raw)
  gc()

  # Create a Seurat object without filtering
  seurat_obj <- CreateSeuratObject(counts = data, project = project)

  # Remove 'data' as it's no longer needed.
  rm(data)
  gc(full = TRUE)

  # Perform transformations and find clusters
  seurat_obj <- SCTransform(seurat_obj, verbose = F)
  seurat_obj <- RunPCA(seurat_obj, verbose = F)

  # Use the provided 'dims' parameter in the following functions
  if (umap.method == "uwot") {
    # Run UMAP
    seurat_obj <- RunUMAP(seurat_obj, dims = dims_umap, umap.method = umap.method)
  } else if (umap.method == "umap-learn") {
    # Run UMAP
    seurat_obj <- RunUMAP(seurat_obj, dims = dims_umap, umap.method = umap.method, metric = "correlation")
  } else {
    stop("Invalid UMAP method")
  }
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_umap, verbose = F)
  seurat_obj <- FindClusters(seurat_obj, verbose = T)

  # Extract clusters from Seurat object and add to SoupChannel
  meta <- seurat_obj@meta.data
  umap <- seurat_obj@reductions$umap@cell.embeddings
  sc <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))

  # Estimate contamination fractitfidfMin
  sc <- autoEstCont(sc, tfidfMin = tfidfMin)
  out <- adjustCounts(sc, roundToInt = TRUE)

  list(seurat_obj = seurat_obj, meta = meta, umap = umap, out = out)
}

process_lane <- function(lane) {
  aligned_data <- read_aligned_data(lane$base_directory, lane$name, output)
  soupX_obj <- prep_seurat_and_soupX(data.raw = aligned_data$raw, data = aligned_data$filtered, project = aligned_data$project)

  # Now we no longer need 'aligned_data', we can remove it.
  rm(aligned_data)
  gc()

  # Call 'create_seurat_and_sce()' while 'soupX_obj' is still in scope.
  feature_set1 <- list(feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  sce_obj <- create_seurat_and_sce(out = soupX_obj$out, project = lane$name, feature_set = feature_set1)

  # Now we no longer need 'soupX_obj', we can remove it.
  rm(soupX_obj)
  gc(full = TRUE)

  # Here, you may want to save 'sce_obj' to a file if it's a large object.
  # Otherwise, return it.
  return(sce_obj)
}

create_seurat_and_sce <- function(out, project, feature_set) {
  # Create singular Seurat object
  seu <- CreateSeuratObject(counts = out, project = project)

  # Convert Seurat objects to SingleCellExperiment object
  sce <- as.SingleCellExperiment(seu)

  list(seu = seu, sce = sce)
}
run_scDblFinder_and_merge <- function(
    sce_list, save_plot = TRUE, file_name = "after_dbl_removal_and_merge",
    path = output) {
  set.seed(1234)
  result_list <- list()
  dbl_table_list <- list()
  for (i in seq_along(sce_list)) {
    sce_list[[i]] <- scDblFinder(sce_list[[i]])
    dbl_table <- as.data.frame(table(call = sce_list[[i]]$scDblFinder.class))
    dbl_table_list[[i]] <- dbl_table
    sce_list[[i]] <- sce_list[[i]][, sce_list[[i]]$scDblFinder.class ==
      "singlet"]
    result_list[[paste0("seurat_obj", i)]] <- as.Seurat(sce_list[[i]])
  }
  merged_dbl_table <- do.call(rbind, dbl_table_list)
  write.table(merged_dbl_table,
    file = paste0(path, "/merged_doublet_table.txt"),
    sep = "\t", row.names = TRUE, quote = FALSE
  )
  # Merge all seurat objects in the list
  merged_seurat_obj <- Reduce(function(x, y) merge(x, y), result_list)
  if (save_plot) {
    create_feature_scatter_plot(
      obj = merged_seurat_obj,
      feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
      file_name = file_name, path = path
    )
  }
  return(merged_seurat_obj)
}


filter_cells <- function(seurat_obj, path = output, save_plots = TRUE) {
  # Get parameters from config
  lower.nFeature <- config$filter_cells$lower.nFeature
  upper.nFeature <- config$filter_cells$upper.nFeature
  max.percent.mt <- config$filter_cells$max.percent.mt
  species <- config$species
  if (is.null(species)) {
    stop("Species is NULL. Please check your config file.")
  }

  # Define mt features based on species
  if (species == "pig") {
    mt.list <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")
    if (is.null(mt.list)) {
      warning("mt.list is not provided for pig species")
    } else {
      # Check if all features are in the matrix
      if (!all(mt.list %in% rownames(seurat_obj))) {
        for (mt in mt.list) {
          print(paste(mt, " in rownames: ", all(mt %in% rownames(seurat_obj))))
          if (all(mt %in% rownames(seurat_obj)) == FALSE) {
            print(paste(mt, " is not in rownames"))
            # remove mt from mt.list
            mt.list <- mt.list[!mt.list %in% mt]
          }
        }
      }
    }
    percent_mt <- PercentageFeatureSet(seurat_obj, features = mt.list, assay = "RNA")
  } else if (species == "human") {
    percent_mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", assay = "RNA")
  } else {
    stop("Species not recognized: please input 'pig' or 'human'")
  }
  seurat_obj[["percent.mt"]] <- percent_mt

  # Create pre-filter plot
  pre_filter_plot <- create_feature_scatter_plot(seurat_obj, "nCount_RNA", "percent.mt", file_name = "percent_mt_unfiltered", save = save_plots, path = path)

  # Filter cells based on nFeature_RNA and percent.mt
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > lower.nFeature & nFeature_RNA < upper.nFeature & percent.mt < max.percent.mt)

  # Create post-filter plots
  post_filter_plot <- create_feature_scatter_plot(seurat_obj, "nCount_RNA", "percent.mt", file_name = "percent_mt_filtered", save = save_plots, path = path)


  return(seurat_obj)
}

normalize_data <- function(seurat_obj, path = output) {
  # Get parameters from config
  min_size <- config$normalize_data$min_size
  min_mean <- config$normalize_data$min_mean
  feature <- config$normalize_data$feature

  sce <- as.SingleCellExperiment(seurat_obj)

  # Cluster the cells for scran
  clusters <- quickCluster(sce, use.ranks = FALSE, min.size = min_size)

  # Calculate size factors per cell
  sce <- computeSumFactors(sce, clusters = clusters, min.mean = min_mean)

  # Apply size factors to generate log normalized data
  sce <- logNormCounts(sce)

  # Replace Seurat normalized values with scran
  seurat_obj[["RNA"]] <- SetAssayData(seurat_obj[["RNA"]],
    slot = "data",
    new.data = logcounts(sce)
  )

  seurat_obj$sizeFactors <- sizeFactors(sce)
  seurat_obj <- UpdateSeuratObject(seurat_obj)

  # Create violin plots pre and post normalization
  vin_pre <- VlnPlot(seurat_obj, feature, slot = "counts")
  vin_post <- VlnPlot(seurat_obj, feature, slot = "data")

  # Save violin plots
  pdf(file = paste0(path, "violin_pre_norm.pdf"), width = 8, height = 6)
  print(vin_pre)
  dev.off()

  pdf(file = paste0(path, "violin_post_norm.pdf"), width = 8, height = 6)
  print(vin_post)
  dev.off()

  return(seurat_obj)
}

feature_selection <- function(seurat_obj) {
  n_features <- config$feature_selection$n_features
  analysis_type <- config$feature_selection$analysis_type

  if (analysis_type == "Seurat") {
    # Seurat method
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = n_features)
  } else if (analysis_type == "Scry") {
    # Scry method
    m <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
    devi <- scry::devianceFeatureSelection(m)
    dev_ranked_genes <- rownames(seurat_obj)[order(devi, decreasing = TRUE)]
    topdev <- head(dev_ranked_genes, n_features)
    # replace variable features with the deviance ranked genes
    VariableFeatures(seurat_obj) <- topdev
    seurat_obj <- UpdateSeuratObject(seurat_obj)
  } else {
    stop("Invalid analysis_type. Please choose 'Seurat' or 'Scry'.")
  }
  return(seurat_obj)
}

scale_data <- function(seurat_obj, path = output) {
  vars.2.regress <- config$scale_data$vars.2.regress
  marker.path.s <- config$scale_data$marker.path.s
  marker.path.g2m <- config$scale_data$marker.path.g2m
  species <- config$species # Get species information

  # Get all gene names
  all.genes <- rownames(seurat_obj)

  # Scale the data
  if (vars.2.regress == "cell.cycle") {
    # Read cell cycle markers
    cell.cycle.markers.s <- read.csv2(marker.path.s,
      sep = "\t", header = TRUE, row.names = 1
    )
    cell.cycle.markers.g2m <- read.csv2(marker.path.g2m,
      sep = "\t", header = TRUE, row.names = 1
    )
    varslist <- c(cell.cycle.markers.s, cell.cycle.markers.g2m)

    # Select species-specific cell cycle markers
    if (species == "human") {
      s.genes <- varslist[1]$human.gene.name
      g2m.genes <- varslist[1]$human.gene.name
    } else if (species == "pig") {
      s.genes <- varslist[4]$pig.gene.name
      g2m.genes <- varslist[8]$pig.gene.name
    } else {
      stop("Unsupported species")
    }
    # Check if the genes in the list are present in the dataset
    missing_genes <- setdiff(g2m.genes, rownames(seurat_obj))
    if (length(missing_genes) > 0) {
      print(paste("Missing genes: ", paste(missing_genes, collapse = ", ")))
    }

    # Perform cell cycle scoring
    seurat_obj <- CellCycleScoring(seurat_obj,
      s.features = s.genes,
      g2m.features = g2m.genes, set.ident = TRUE
    )

    # visualize before regressing out cell cycle
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    pdf(paste0(path, "pca_before_cc_regression.pdf"), width = 8, height = 6)
    print(DimPlot(seurat_obj))
    dev.off()
    # scale data and regress out cell cycle
    seurat_obj <- ScaleData(seurat_obj,
      vars.to.regress = c("S.Score", "G2M.Score"),
      features = all.genes
    )
    # visualize after regressing out cell cycle
    # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    pdf(paste0(path, "pca_after_cc_regression.pdf"), width = 8, height = 6)
    print(DimPlot(seurat_obj))
    dev.off()
  } else {
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  }
  return(seurat_obj)
}

sc_transform <- function(seurat_obj, path = output) {
  vars.2.regress <- config$sc_transform$vars.2.regress
  species <- config$species # Get species information
  marker.path.s <- config$sc_transform$marker.path.s
  marker.path.g2m <- config$sc_transform$marker.path.g2m

  # Get all gene names
  all.genes <- rownames(seurat_obj)

  # Scale the data
  if (vars.2.regress == "cell.cycle") {
    # Read cell cycle markers
    cell.cycle.markers.s <- read.csv2(marker.path.s,
      sep = "\t", header = TRUE, row.names = 1
    )
    cell.cycle.markers.g2m <- read.csv2(marker.path.g2m,
      sep = "\t", header = TRUE, row.names = 1
    )
    varslist <- c(cell.cycle.markers.s, cell.cycle.markers.g2m)

    # Select species-specific cell cycle markers
    if (species == "human") {
      s.genes <- varslist[1]$human.gene.name
      g2m.genes <- varslist[1]$human.gene.name
    } else if (species == "pig") {
      s.genes <- varslist[4]$pig.gene.name
      g2m.genes <- varslist[8]$pig.gene.name
    } else {
      stop("Unsupported species")
    }
    # Check if the genes in the list are present in the dataset
    missing_genes <- setdiff(g2m.genes, rownames(seurat_obj))
    if (length(missing_genes) > 0) {
      print(paste("Missing genes: ", paste(missing_genes, collapse = ", ")))
    }

    # Perform cell cycle scoring
    seurat_obj <- CellCycleScoring(seurat_obj,
      s.features = s.genes,
      g2m.features = g2m.genes, set.ident = TRUE
    )

    # visualize before regressing out cell cycle
    seurat_obj <- Seurat::SCTransform(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    pdf(paste0(path, "pca_before_cc_regression.pdf"), width = 8, height = 6)
    print(DimPlot(seurat_obj))
    dev.off()
    # scale data and regress out cell cycle
    seurat_obj <- Seurat::SCTransform(seurat_obj,
      vars.to.regress = c("S.Score", "G2M.Score")
    )
    # visualize after regressing out cell cycle
    # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
    seurat_obj <- RunPCA(seurat_obj, features = c(s.genes, g2m.genes))
    pdf(paste0(path, "pca_after_cc_regression.pdf"), width = 8, height = 6)
    print(DimPlot(seurat_obj))
    dev.off()
  } else {
    seurat_obj <- Seurat::SCTransform(seurat_obj)
  }
  return(seurat_obj)
}

run_and_visualize_pca <- function(seurat_obj, path = output) {
  top_n_dims <- config$run_and_visualize_pca$top_n_dims
  heatmap_dims <- 1:config$run_and_visualize_pca$heatmap_dims
  num_cells <- config$run_and_visualize_pca$num_cells
  dims <- 1:config$run_and_visualize_pca$dims
  num_replicate <- config$run_and_visualize_pca$num.replicate

  # Perform PCA
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

  # Visualize feature loadings
  pdf(paste0(path, "top_n_dims_with_genes.pdf"), width = 8, height = 6)
  print(VizDimLoadings(seurat_obj, dims = 1:top_n_dims, reduction = "pca"))
  dev.off()

  # Generate scatter plot of PCA results
  pdf(paste0(path, "pca_scatter_plot.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "pca"))
  dev.off()

  # Save PCA heat map
  pdf(paste0(path, "pca_heat_map.pdf"), width = 8, height = 6)
  print(DimHeatmap(seurat_obj, nfeatures = 5, dims = heatmap_dims, cells = num_cells, balanced = TRUE, fast = FALSE, combine = TRUE))
  dev.off()
  # Save Elbow plots
  pdf(paste0(path, "elbow_pca.pdf"), width = 8, height = 6)
  elbow_pca <- ElbowPlot(seurat_obj, reduction = "pca")
  print(elbow_pca)
  dev.off()

  # Perform JackStraw
  if (num_replicate == "NA") {
    print("jackstraw not run")
  } else {
    seurat_obj <- JackStraw(seurat_obj, num.replicate = num_replicate)
    seurat_obj <- ScoreJackStraw(seurat_obj, dims = dims)

    # Save JackStraw plot
    pdf(paste0(path, "jack_straw.pdf"), width = 8, height = 6)
    jack_straw <- JackStrawPlot(seurat_obj, dims = dims)
    print(jack_straw)
    dev.off()
  }

  # Return the updated Seurat object
  return(seurat_obj)
}

perform_batch_correction <- function(seurat_obj, path = output) {
  dims.use <- 1:config$perform_batch_correction$dims.use
  max_iter <- config$perform_batch_correction$max_iter

  # Generate pre-batch correction PCA plot
  pdf(paste0(path, "batch_uncorrected_pca.pdf"), width = 8, height = 6)
  p1_pre <- DimPlot(object = seurat_obj, reduction = "pca", pt.size = .1, group.by = "orig.ident")
  p2_pre <- VlnPlot(object = seurat_obj, features = "PC_1", group.by = "orig.ident", pt.size = .1)
  print(p1_pre + p2_pre)
  dev.off()

  # Run Harmony
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "orig.ident", dims.use = dims.use, max.iter.harmony = max_iter)

  # Access Harmony embeddings
  harmony_embeddings <- Embeddings(seurat_obj, "harmony")

  # Generate post-batch correction PCA plot
  pdf(paste0(path, "batch_corrected_pca.pdf"), width = 8, height = 6)
  p1_post <- DimPlot(object = seurat_obj, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
  p2_post <- VlnPlot(object = seurat_obj, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
  print(p1_post + p2_post)
  dev.off()

  # Return the updated Seurat object and harmony embeddings
  return(list(seurat_obj = seurat_obj, harmony_embeddings = harmony_embeddings))
}

run_umap <- function(seurat_obj, path = output) {
  dims_umap <- 1:config$run_umap$dims_umap
  umap.method <- config$run_umap$umap.method
  umap.red <- config$run_umap$umap.red
  # Run UMAP
  seurat_obj <- RunUMAP(seurat_obj,
    dims = dims_umap, umap.method = umap.method,
    reduction = umap.red, group.by = "orig.ident"
  )
  # Generate UMAP plot
  pdf(paste0(path, "umap_plot.pdf"), width = 8, height = 6)
  print(DimPlot(seurat_obj, reduction = "umap"))
  dev.off()

  # Return the updated Seurat object
  return(seurat_obj)
}

perform_clustering <- function(seurat_obj, path = output) {
  resolution <- config$perform_clustering$resolution
  algorithm <- config$perform_clustering$algorithm
  reduction <- config$perform_clustering$reduction
  dims_snn <- 1:config$perform_clustering$dims_snn

  # Check if Harmony embeddings exist in the Seurat object
  batch_corrected <- "harmony" %in% names(Embeddings(seurat_obj))

  # If batch correction was not performed and reduction is set to "harmony", update it to "pca"
  if (!batch_corrected && reduction == "harmony") {
    message("Batch correction was skipped. Updating reduction to 'pca'.")
    reduction <- "pca"
  }

  # Perform K-nearest neighbor (KNN) graph
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims_snn, reduction = reduction)

  # Save UMAP lanes plot
  pdf(paste0(path, "umap_lanes.pdf"), width = 8, height = 6)
  umap_lanes <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", pt.size = .5)
  print(umap_lanes)
  dev.off()

  # Cluster cells
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution, algorithm = algorithm)

  # Save UMAP clusters plot
  pdf(paste0(path, "umap_clusters.pdf"), width = 8, height = 6)
  umap_clusters <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = .5)
  print(umap_clusters)
  dev.off()

  # Return the updated Seurat object
  return(seurat_obj)
}


find_differentially_expressed_features <- function(seurat_obj, path = output) {
  # Get parameters from the config file
  min_pct <- config$find_differentially_expressed_features$min_pct
  logfc_threshold <- config$find_differentially_expressed_features$logfc_threshold
  top_n <- config$find_differentially_expressed_features$top_n

  # Find all markers
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = min_pct, logfc.threshold = logfc_threshold)

  # Extract top_n markers
  markers %>%
    group_by(cluster) %>%
    top_n(n = top_n, wt = avg_log2FC) -> topMarkers

  # Write out top_n markers
  write.table(topMarkers, file = paste0(path, "/seurat_obj.DE.markers.top", top_n, "_S2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

  # Define the top 4 features for FeaturePlot based on topMarkers
  top_features <- topMarkers$gene[1:4]

  # UMAP plots
  plot1 <- UMAPPlot(seurat_obj, group.by = "orig.ident")
  plot2 <- UMAPPlot(seurat_obj, label = T)
  plot3 <- FeaturePlot(seurat_obj, features = top_features, ncol = 2, pt.size = 0.1)

  # Combine plots
  combined_plot <- ((plot1 / plot2) | plot3) + plot_layout(width = c(1, 2))

  # Save plot
  pdf(file = paste0(path, "/combined_plot.pdf"), width = 8, height = 6)
  print(combined_plot)
  dev.off()

  # Heatmap of top_n markers
  heatmap <- DoHeatmap(seurat_obj, features = topMarkers$gene) + NoLegend()

  # Save heatmap
  pdf(file = paste0(path, "/heatmap_top_n_markers.pdf"))
  print(heatmap)
  dev.off()

  # Write out all markers
  write.table(markers, file = paste0(path, "/seurat_obj.DE.markers_S2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

  return(list(markers = markers, topMarkers = topMarkers))
}

analyze_known_markers <- function(seurat_obj, de_results, output_path = output) {
  # read in known GAMM retinoid markers
  known_markers_path <- config$score_and_plot_markers$known_markers_path
  known.markers <- read.csv2(known_markers_path, sep = "\t", header = TRUE)
  de.markers <- de_results[[1]]
  # match with any DE markers from data by merging dataframes
  marker_df <- merge(de.markers, known.markers, by = "gene")
  # write out marker df with known DE markers
  write.table(marker_df, file = paste0(output_path, "seurat_obj.knownDE.markers_S2.txt"), quote = FALSE, sep = "\t", row.names = FALSE)

  # first get unique cell type vector
  cell.types <- unique(marker_df$Cell.type)
  print(cell.types)

  # check gene number between known markers and DE known markers
  genesk <- unique(known.markers$gene)
  k <- length(genesk)
  genesDEk <- unique(marker_df$gene)
  DEk <- length(genesDEk)
  percent.markers.de <- (DEk / k) * 100
  percent.markers.de

  # check a specific marker
  df <- subset(marker_df, Cell.type == "Synaptic marker", select = c("gene"))
  length(df$gene)
  df2 <- subset(known.markers, Cell.type == "Synaptic marker", select = c("gene"))
  length(df2$gene)

  # subset all and plot using for loop
  for (i in 1:length(cell.types)) {
    new_df <- subset(marker_df, Cell.type == cell.types[i], select = c("gene", "Cell.type", "cluster"))
    new_vec <- unique(as.vector(new_df$gene))

    pdf(paste0(output_path, cell.types[i], "_featureplot.pdf"), width = 8, height = 6, bg = "white")
    # umap plot highlighting gene expression
    print(FeaturePlot(seurat_obj, features = new_vec))
    dev.off()

    pdf(paste0(output_path, cell.types[i], "_dotplot.pdf"), width = 8, height = 6, bg = "white")
    # expression dot plot
    dot.plot <- DotPlot(object = seurat_obj, features = new_vec)
    print(dot.plot + labs(title = cell.types[i]))
    dev.off()
  }
}

score_and_plot_markers <- function(seurat_obj, output_path = output) {
  known_markers_path <- config$score_and_plot_markers$known_markers_path
  known_markers <- config$score_and_plot_markers$known_markers
  top_n_markers <- config$score_and_plot_markers$top_n_markers
  cluster_type <- config$score_and_plot_markers$cluster_type
  pairwise <- config$score_and_plot_markers$pairwise
  logFC_thresh <- config$score_and_plot_markers$logFC_thresh
  auc_thresh <- config$score_and_plot_markers$auc_thresh


  sce_obj <- as.SingleCellExperiment(seurat_obj)

  # Score markers
  marker_info <- score_markers(sce_obj, cluster_type)

  # Read in known markers
  known.markers.df <- if (known_markers) {
    read.csv2(known_markers_path, sep = "\t", header = TRUE, row.names = 1)
  } else {
    NULL
  }

  annot_df <- data.frame(Cluster = integer(), Cell.type = character())

  # Determine clusters based on cluster_type
  clusters <- get_clusters(seurat_obj, cluster_type)

  # Process each cluster
  for (i in 1:length(clusters)) {
    # Process pairwise comparisons
    if (pairwise) {
      process_pairwise_comparisons(clusters, i, marker_info, output_path, logFC_thresh, auc_thresh, known.markers.df, seurat_obj)
    } else {
      print("No pairwise comparison")
    }

    # Process top genes for each cluster
    clust <- as.data.frame(marker_info[[clusters[i]]])
    annot_df <- process_top_genes(clust, clusters, i, known.markers.df, output_path, seurat_obj, annot_df)
  }

  return(annot_df)
}

score_markers <- function(sce_obj, cluster_type) {
  # Score markers based on cluster type
  if (cluster_type %in% c("seurat_clusters", "orig.ident")) {
    marker_field <- cluster_type
    marker_info <- scoreMarkers(sce_obj, sce_obj@colData@listData[[marker_field]], full.stats = TRUE)
  } else {
    stop("Invalid cluster_type. Please choose 'seurat_clusters' or 'orig.ident'.")
  }
  return(marker_info)
}

get_clusters <- function(seurat_obj, cluster_type) {
  if (cluster_type %in% c("seurat_clusters", "orig.ident")) {
    return(unique(seurat_obj@meta.data[[cluster_type]]))
  } else {
    stop("Invalid cluster_type. Please choose 'seurat_clusters' or 'orig.ident'.")
  }
}

process_top_genes <- function(clust, clusters, i, known.markers.df, output_path, seurat_obj, annot_df) {
  # Extract relevant configuration settings
  top_n_markers <- config$score_and_plot_markers$top_n_markers
  logFC_thresh <- config$score_and_plot_markers$logFC_thresh
  auc_thresh <- config$score_and_plot_markers$auc_thresh
  known_markers <- config$score_and_plot_markers$known_markers

  # Order by median Cohen's d and subset
  ordered <- subset(clust[order(clust$median.logFC.cohen, decreasing = TRUE), ], median.logFC.cohen > logFC_thresh & median.AUC > auc_thresh)
  top100 <- head(ordered, n = top_n_markers)

  # Create a subdirectory for Top100 DE genes
  top100_dir <- file.path(output_path, "Top100_DE_Genes")
  if (!dir.exists(top100_dir)) {
    dir.create(top100_dir, recursive = TRUE)
  }

  # Write out top100 genes
  write.table(top100,
    file = file.path(top100_dir, paste0("Top100DEgenes_clust_", clusters[i], ".txt")),
    quote = FALSE, sep = "\t", row.names = TRUE
  )

  # Process known markers
  return(process_known_markers(top100, known_markers, known.markers.df, clusters, i, output_path, top_n_markers, seurat_obj, annot_df))
}


process_pairwise_comparisons <- function(clusters, i, marker.info, output_path, logFC_thresh, auc_thresh, known_markers_df, seurat_obj) {
  # Extract cluster marker information
  clust <- as.data.frame(marker.info[[clusters[i]]])

  # Prepare the data for pairwise comparison
  clust_tbl <- rownames_to_column(clust, var = "gene") %>% as_tibble()
  df_clust0 <- clust_tbl %>% dplyr::select(starts_with("full.logFC.cohen"))
  df_clust <- clust_tbl %>% dplyr::select(c(
    "gene", starts_with("full.logFC.cohen"),
    starts_with("full.AUC")
  ))

  # Create a subdirectory for pairwise DE genes
  pairwise_de_path <- file.path(output_path, paste0("Pairwise_DE_Cluster_", clusters[i]))
  if (!dir.exists(pairwise_de_path)) {
    dir.create(pairwise_de_path, recursive = TRUE)
  }

  # Iterate through each pairwise comparison
  for (j in colnames(df_clust0)) {
    print(j)
    k <- str_split_1(j, "cohen.")
    l <- paste0("full.AUC.", k[2])

    # Subset based on logFC threshold and AUC threshold
    df_clust1 <- subset(df_clust, df_clust[j] > logFC_thresh & df_clust[l] > auc_thresh)
    df_clust1 <- df_clust1[order(df_clust1[j], decreasing = TRUE), ]

    # Write table of all DE genes for the comparison
    write.table(df_clust1,
      file = file.path(pairwise_de_path, paste0("DEgenes_", clusters[i], "_vs_", j, ".txt")),
      quote = FALSE, sep = "\t", row.names = TRUE
    )

    # Merge with known markers and check if empty
    df_clust2 <- merge(df_clust1, known_markers_df, by.x = "gene", by.y = "row.names")
    if (nrow(df_clust2) == 0) {
      print(paste0("This data frame is empty: ", clusters[i], ".vs_", j))
    } else {
      # Write table with known DE markers
      write.table(df_clust2,
        file = file.path(pairwise_de_path, paste0("KnownDEgenes_", clusters[i], "_vs_", j, ".txt")),
        quote = FALSE, sep = "\t", row.names = TRUE
      )

      # Make feature plot of genes
      new_vec0 <- unique(as.vector(df_clust2$gene))
      pdf(file.path(pairwise_de_path, paste0(clusters[i], "_vs_", j, "_featureplot.pdf")), bg = "white")
      print(FeaturePlot(seurat_obj, features = new_vec0), label = TRUE)
      dev.off()
    }
  }
}



process_known_markers <- function(top100, known_markers_flag, known_markers_df, clusters, i, output_path, top_n_markers, seurat_obj, annot_df) {
  if (known_markers_flag) {
    marker_df <- merge(top100, known_markers_df, by = "row.names")

    if (nrow(marker_df) == 0) {
      print(paste0("This data frame is empty: ", clusters[i]))
      new_row <- data.frame(Cluster = clusters[i], Cell.type = "unknown")
      annot_df <- rbind(annot_df, new_row)
    } else {
      subdirectory_path <- file.path(output_path, "Known_DE_Markers")
      if (!dir.exists(subdirectory_path)) {
        dir.create(subdirectory_path, recursive = TRUE)
      }
      # Write out marker dataframe with known DE markers
      write.table(marker_df,
        file = paste0(output_path, "KnownDE.markers_clust_", clusters[i], ".txt"),
        quote = FALSE, sep = "\t", row.names = FALSE
      )

      # Subset data
      new_df <- marker_df[, c("Row.names", "rank.logFC.cohen", "Cell.type")]
      new_vec <- unique(as.vector(new_df$Row.names))

      # Get top ranked
      rank <- top_n_markers + 1
      new_df.ordered <- new_df[order(new_df$rank.logFC.cohen), ]
      new_df.ordered <- subset(new_df.ordered, rank.logFC.cohen < rank)
      new_vec2 <- unique(as.vector(new_df.ordered$Row.names))

      if (identical(new_vec2, character(0))) {
        print(paste0("This vector does not have any ranks in top ", top_n_markers, ": ", clusters[i]))
        new_row <- data.frame(Cluster = clusters[i], Cell.type = "unknown")
        annot_df <- rbind(annot_df, new_row)
      } else {
        # UMAP plot highlighting gene expression
        pdf(paste0(output_path, clusters[i], "_featureplot_top", top_n_markers, "ranks.pdf"), bg = "white")
        print(FeaturePlot(seurat_obj, features = new_vec2), label = TRUE)
        dev.off()
        allcelltypes <- unique(as.vector(new_df.ordered$Cell.type))
        result_string <- paste(allcelltypes, collapse = "-")
        new_row <- data.frame(Cluster = clusters[i], Cell.type = result_string)
        annot_df <- rbind(annot_df, new_row)
      }
    }
  } else {
    print("No known marker set")
  }
  return(annot_df)
}



annotate_clusters_and_save <- function(seurat_obj, new_cluster_ids, output_path = output) {
  # Rename the clusters based on the new IDs
  names(new_cluster_ids) <- levels(seurat_obj)
  seurat_obj <- RenameIdents(seurat_obj, new_cluster_ids)

  # Generate and plot the UMAP plot

  pdf(paste0(output_path, "labeled-clusters.pdf"), bg = "white")
  print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5))
  dev.off()
  # Save the Seurat object
  saveRDS(seurat_obj, file = paste0(output_path, "seurat_obj_labeled.rds"))

  return(seurat_obj)
}

create_feature_scatter_plot <- function(obj, feature1, feature2, save = TRUE, file_name = NULL, path = output) {
  plot <- FeatureScatter(object = obj, feature1 = feature1, feature2 = feature2)

  # If save is TRUE, save the plot
  if (save) {
    # If file_name is not provided, use default naming scheme
    if (is.null(file_name)) {
      file_name <- paste0(feature1, "_vs_", feature2, "_scatter")
    }

    pdf(paste0(path, file_name, ".pdf"), width = 8, height = 6)
    print(plot)
    dev.off()
  }

  return(plot)
}

combine_feature_plots <- function(seurat_objs_list, feature_set1, feature_set2 = NULL, same_feature_set = TRUE, file_name_base = "post_SoupX_plot", path = output) {
  # If same_feature_set is TRUE, use feature_set1 for both plots
  if (same_feature_set) {
    feature_set2 <- feature_set1
  }

  # Create individual plots and save them as PDFs
  plots_list <- lapply(seq_along(seurat_objs_list), function(i) {
    if (i %% 2 == 1) { # if i is odd
      create_feature_scatter_plot(seurat_objs_list[[i]], feature_set1$feature1, feature_set1$feature2, save = TRUE, file_name = paste0(file_name_base, i), path = path)
    } else { # if i is even
      create_feature_scatter_plot(seurat_objs_list[[i]], feature_set2$feature1, feature_set2$feature2, save = TRUE, file_name = paste0(file_name_base, i), path = path)
    }
  })

  # Add plots together and save as a new PDF
  pdf(paste0(path, file_name_base, "_combined.pdf"), width = 8 * length(plots_list), height = 6)
  combined_plot <- cowplot::plot_grid(plotlist = plots_list)
  print(combined_plot)
  dev.off()

  # Delete the individual plot files
  for (i in seq_along(seurat_objs_list)) {
    file.remove(paste0(path, file_name_base, i, ".pdf"))
  }

  return(combined_plot)
}

annotate_with_clustifyR <- function(clustered_seurat_obj, output) {
  # Access the markers path
  markers_path <- config$score_and_plot_markers$known_markers_path
  markers <- read.csv2(markers_path, sep = "\t", header = TRUE)
  markers_df <- data.frame(markers$gene, markers$Cell.type)
  colnames(markers_df) <- c("gene", "cluster")

  # Clustify lists
  list_res <- clustify_lists(
    input = clustered_seurat_obj,
    cluster_col = "seurat_clusters",
    marker = markers_df,
    metric = "pct",
    marker_inmatrix = FALSE,
    obj_out = FALSE
  )

  if (length(unique(list_res)) > 1) {
    p1 <- plot_cor_heatmap(
      cor_mat = list_res,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      legend_title = "pct"
    )
    pdf(paste0(output, "correlation_heatmap.pdf"), width = 8, height = 6)
    print(p1)
    dev.off()
  } else {
    print(paste0("Insufficient distinct values in list_res for heatmap plotting: ", unique(list_res)))
    message("Insufficient distinct values in list_res for heatmap plotting.")
  }

  # Call cell types
  list_res2 <- cor_to_call(
    cor_mat = list_res,
    cluster_col = "seurat_clusters"
  )

  # Add clustifyr calls as metadata
  clust_call <- call_to_metadata(
    res = list_res2,
    metadata = clustered_seurat_obj@meta.data,
    cluster_col = "seurat_clusters",
    rename_prefix = "clustifyr_call"
  )
  clustered_seurat_obj <- AddMetaData(clustered_seurat_obj, metadata = clust_call)

  # Plot with clustifyr annotations
  pc <- DimPlot(clustered_seurat_obj,
    reduction = "umap", group.by = "clustifyr_call_type", label = TRUE,
    label.size = 3, repel = TRUE
  ) + ggtitle("Clustifyr annotated labels") +
    guides(fill = guide_legend(label.theme = element_text(size = 8)))
  pdf(paste0(output, "clustifyr_marker_annotation_umap.pdf"), width = 11, height = 6)
  print(pc)
  dev.off()

  # Save object with clustifyr annotation
  saveRDS(clustered_seurat_obj, file = paste0(output, "seurat_obj_clustifyr.rds"))
}
