# ENVIRONMENT SETUP
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(BiocParallel)
library(cowplot)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
library(scater)
library(SingleCellExperiment)
library(tidyverse)
# use_python("~/miniconda3/envs/scRNAseq_new2/bin/python")
# set variables
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "config.json"))
docker <- config$docker
if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$recluster$DATA_DIR
}
SEURAT.FILE <- config$recluster$SEURAT.FILE
# set up environment and output
# use_condaenv("/w5home/bmoore/miniconda3/envs/scRNAseq_new", required=TRUE)
# setwd(WD)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_recluster_", timestamp)
print(output)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(file.path(paste0(GIT_DIR, "/config.json")), file.path(paste0(
    "./",
    output, "config.json"
)), overwrite = TRUE)
config <- jsonlite::fromJSON(file.path(output, "config.json"))
source(paste0(GIT_DIR, "/src/sc_pipeline_functions.R"))
packageVersion("Seurat")

# Load Data
print("loading data")
seurat.obj <- readRDS(file = paste0(DATA_DIR, SEURAT.FILE))
print(seurat.obj)
# recluster
print("reclustering")
seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(object = seurat.obj))
resolutions <- config$recluster$perform_clustering$resolution
if (!is.list(resolutions) && length(resolutions) == 1) resolutions <- list(resolutions)

for (res in resolutions) {
    resolution <- as.character(res)
    print(paste("Processing resolution:", resolution))

    # cluster at this resolution (adds a new cluster column each iteration)
    seurat.obj <- perform_clustering(seurat.obj, output, "recluster", name = "", resolution = res)

    # Find markers with Differential expressed features (genes)
    # analyse known markers
    # note: results overwrite
    # get single cell object
    print("convert to sce")
    cluster_name <- paste0(config$recluster$perform_clustering$cluster_name, "_res", resolution)
    counts_matrix <- seurat.obj[["RNA"]]$counts
    data_matrix <- seurat.obj[["RNA"]]$data
    # get clusters
    clusters <- seurat.obj[[cluster_name]]
    # get embeddings
    dim_data_pca <- Embeddings(seurat.obj, reduction = "pca")
    dim_data_umap <- Embeddings(seurat.obj, reduction = "umap")
    # make sce object
    sce <- SingleCellExperiment(list(counts = counts_matrix, logcounts = data_matrix),
        colData = DataFrame(label = clusters)
    )
    reducedDims(sce) <- list(PCA = dim_data_pca, UMAP = dim_data_umap)
    print(sce)

    print("run DE analysis")
    annot_df <- score_and_plot_markers(seurat.obj, sce, output, "recluster", name = "", resolution = res)

    # manual annotation
    if (config$recluster$score_and_plot_markers$known_markers == "TRUE" || config$recluster$score_and_plot_markers$known_markers == "true" || config$recluster$score_and_plot_markers$known_markers == "T" || config$recluster$score_and_plot_markers$known_markers == "t") {
        print("annotate and save")
        # order annot_df
        annot_df.ordered <- annot_df[order(as.numeric(annot_df$Cluster)), ]
        new.cluster.ids <- annot_df.ordered$Celltype
        annotation_column <- paste0("CellType1_res", resolution)
        # annotate
        seurat.obj <- annotate_clusters_and_save(seurat.obj, new.cluster.ids, output, "recluster",
            cluster_type = cluster_name, annotation_column = annotation_column)
        # print table of celltypes
        print(table(seurat.obj@meta.data[[annotation_column]]))
        table1 <- table(seurat.obj@meta.data[[annotation_column]])
        table1 <- as.data.frame(table1)
        write.table(table1,
            file = paste0(output, "table_celltype_counts_res", resolution, ".txt"),
            sep = "\t", quote = F, row.names = F
        )
        # write annot df
        annot_df.ordered <- as.data.frame(annot_df.ordered)
        write.table(annot_df.ordered,
            file = paste0(output, "table_celltype_cluster_res", resolution, ".txt"),
            sep = "\t", quote = F, row.names = F
        )
        # plot umap without legend
        pdf(paste0(output, "labeled-clusters2_res", resolution, ".pdf"), bg = "white")
        print(DimPlot(seurat.obj,
            reduction = "umap", label = TRUE,
            pt.size = 0.5, group.by = annotation_column
        ) + NoLegend())
        dev.off()
    }

    # save metadata for this resolution
    print("write metadata")
    seurat.obj.metadata <- as.data.frame(seurat.obj@meta.data)
    write.table(seurat.obj.metadata,
        file = paste0(output, "manual_annot_metadata_res", resolution, ".txt"), sep = "\t",
        quote = F, row.names = T
    )
}

# save final object with all resolution cluster columns
print("save final object")
# saveRDS(seurat.obj, file = paste0(output, "seuratobj_recluster_all_res.rds"))

# save session info
writeLines(capture.output(sessionInfo()), paste0(output, "sessionInfo.txt"))
system(paste("chmod -R 777", output))

# read back in metadata

# # add in metadata from previous analysis
# metadata <- read.csv(paste0(output,"manual_annot_metadata_",resolution,".txt"),
#             row.names = 1, header = TRUE, sep = "\t")

# # add metadata to query data
# seurat.obj_recluster <- AddMetaData(seurat.obj_recluster, metadata)

# # check again
# colnames(seurat.obj_recluster@meta.data)
