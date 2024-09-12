# load packages
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
use_condaenv(condaenv = '/w5home/bmoore/miniconda3/envs/scRNAseq_best/', required = TRUE)
## to install SeuratDisk
# first install hdf5r:
# conda install -c conda-forge r-hdf5r
# Next install SeuratDisk:
# conda install -c conda-forge r-seuratdisk
# or
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
# note for me conda install -c conda-forge r-seuratdisk did not work, but the github install did

library(SeuratDisk)
# variables
GIT_DIR <- "/w5home/bmoore/scRNAseq_downstream/"
config <- fromJSON(file.path(GIT_DIR, "config.json"))
WD <- config$seurat2ann$WD
SEURAT_OBJ <- config$seurat2ann$SEURAT_OBJ
METADATA <- config$seurat2ann$METADATA
OUTPUT_name <- config$seurat2ann$OUTPUT_name
# set working dir
setwd(WD)
# read in Seurat object
seurat <- readRDS(file = SEURAT_OBJ)

if(METADATA=="NA"){
  print("No metadata file provided")
} else {
    print("Metadata file provided")
    print(METADATA)
    print("Reading in metadata")
    # add in metadata from previous analysis
    metadata <- read.csv(METADATA, row.names = 1,
                            header = TRUE, sep = "\t")
    # check metadata with query data
    print(colnames(seurat@assays$RNA@data)[1:10])
    print(rownames(metadata)[1:10])
    print(colnames(seurat@meta.data))

    # add metadata to query data
    seurat <- AddMetaData(seurat, metadata)
    print(colnames(seurat@meta.data))
}

# first save seurat as h5 seurat file
SaveH5Seurat(seurat, filename = OUTPUT_name)
# then convert to h5ad
Convert(OUTPUT_name, dest = "h5ad")