# ENVIRONMENT SETUP
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(BiocParallel)
#library(DropletUtils)
library(cowplot)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
library(scater)
library(SingleCellExperiment)
library(tidyverse)
use_python("~/miniconda3/envs/scRNAseq_new2/bin/python")
# set variables
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "config.json"))
WD <- config$recluster$DATA_DIR
SEURAT.FILE <- config$recluster$SEURAT.FILE
# set up environment and output
#use_condaenv("/w5home/bmoore/miniconda3/envs/scRNAseq_new", required=TRUE)
setwd(WD)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output_recluster_", timestamp)
dir.create(output, showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(file.path(paste0(GIT_DIR,"/config.json")), file.path(paste0("./", 
          output,"config.json")), overwrite = TRUE)
config <- jsonlite::fromJSON(file.path(output, "config.json"))
source(paste0(GIT_DIR,"/src/sc_pipeline_functions.R"))
packageVersion("Seurat")

# Load Data
print("loading data")
seurat.obj <- readRDS(file = SEURAT.FILE)

# recluster
print("reclustering")
seurat.obj_recluster <- perform_clustering(seurat.obj, output, "recluster")

# save object
print("saving")
resolution <- as.character(config$recluster$perform_clustering$resolution)
#saveRDS(seurat.obj_recluster, file = paste0(output, "seuratobj_recluster_res", resolution, ".rds"))

# Find markers with Differential expressed features (genes)
# analyse known markers
# note: results overwrite
# get single cell object
print("convert to sce")
cluster_name <- config$recluster$perform_clustering$cluster_name
counts_matrix <- seurat.obj_recluster[["RNA"]]$counts
data_matrix <- seurat.obj_recluster[["RNA"]]$data
# get clusters
clusters <-  seurat.obj_recluster[[cluster_name]]
# get embeddings
dim_data_pca <- Embeddings(seurat.obj_recluster, reduction = "pca")
dim_data_umap <- Embeddings(seurat.obj_recluster, reduction = "umap")
# make sce object
sce <- SingleCellExperiment(list(counts=counts_matrix, logcounts=data_matrix),
    colData=DataFrame(label=clusters))
reducedDims(sce) <- list(PCA=dim_data_pca, UMAP=dim_data_umap)
print(sce)

print("run DE analysis")
annot_df <- score_and_plot_markers(seurat.obj_recluster,sce, output,"recluster")

# manual annotation
print("annotate and save")
# order annot_df
annot_df.ordered <- annot_df[order(as.numeric(annot_df$Cluster)), ]
new.cluster.ids <- annot_df.ordered$Cell.type
# annotate
seurat.obj_recluster<-annotate_clusters_and_save(seurat.obj_recluster, new.cluster.ids, output, "recluster")
# print table of celltypes
print(table(seurat.obj_recluster@meta.data$CellType1))

# save annotation as dataframe
print("write metadata and visualize")
seurat.obj_recluster.metadata<- as.data.frame(seurat.obj_recluster@meta.data)
write.table(seurat.obj_recluster.metadata, file = paste0(output,"manual_annot_metadata_",resolution,".txt"), sep="\t", 
            quote=F, row.names = T)
table1 <- table(seurat.obj_recluster@meta.data$CellType1)
table1 <- as.data.frame(table1)
write.table(table1, file= paste0(output,"table_celltype_counts.txt"), 
            sep="\t", quote=F, row.names = F)
# write annot df
annot_df.ordered <- as.data.frame(annot_df.ordered)
write.table(annot_df.ordered, file= paste0(output,"table_celltype_cluster.txt"), 
            sep="\t", quote=F, row.names = F)
# optional- plots and metadata
# plot umap without legend
pdf(paste0(output, "labeled-clusters2.pdf"), bg = "white")
print(DimPlot(seurat.obj_recluster, reduction = "umap", label = TRUE, 
      pt.size = 0.5, group.by = "CellType1")+ NoLegend())
dev.off()
# save session info
writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))


# read back in metadata

# # add in metadata from previous analysis
# metadata <- read.csv(paste0(output,"manual_annot_metadata_",resolution,".txt"), 
#             row.names = 1, header = TRUE, sep = "\t")

# # add metadata to query data
# seurat.obj_recluster <- AddMetaData(seurat.obj_recluster, metadata)

# # check again
# colnames(seurat.obj_recluster@meta.data)
