
# load_packages
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

# set_variables
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(GIT_DIR, "config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$de$DATA_DIR
}
SEURAT.FILE <- config$de$SEURAT.FILE
outname <- config$de$outname
setwd(DATA_DIR)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output_DE_", outname, "_",timestamp)
print(output)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(file.path(paste0(GIT_DIR,"/config.json")), file.path(paste0("./", 
          output,"config.json")), overwrite = TRUE)
config <- jsonlite::fromJSON(file.path(output, "config.json"))
source(paste0(GIT_DIR,"/src/sc_pipeline_functions.R"))
packageVersion("Seurat")

# readin
seurat.obj <- readRDS(file = paste0(SEURAT.FILE))
print(seurat.obj)


# convert_to_sce
print("convert to sce")
cluster_name <- config$de$score_and_plot_markers$cluster_type
counts_matrix <- seurat.obj[["RNA"]]$counts
data_matrix <- seurat.obj[["RNA"]]$data
# get clusters
clusters <-  seurat.obj[[cluster_name]]
# get embeddings
dim_data_pca <- Embeddings(seurat.obj, reduction = "pca")
dim_data_umap <- Embeddings(seurat.obj, reduction = "umap")
# make sce object
sce <- SingleCellExperiment(list(counts=counts_matrix, logcounts=data_matrix),
    colData=DataFrame(label=clusters))
reducedDims(sce) <- list(PCA=dim_data_pca, UMAP=dim_data_umap)
print(sce)

# DE_analysis
print("run DE analysis")
annot_df <- score_and_plot_markers(seurat.obj,sce, output,"de")
