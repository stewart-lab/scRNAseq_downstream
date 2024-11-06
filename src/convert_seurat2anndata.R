# load packages
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
library(SeuratObject)
#use_condaenv(condaenv = 'scRNAseq_best', required = TRUE)
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
library(hdf5r)
# variables
GIT_DIR <- get_wd()
config <- fromJSON(file.path("./config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$seurat2ann$DATA_DIR
}
SEURAT_OBJ <- config$seurat2ann$SEURAT_OBJ
METADATA <- config$seurat2ann$METADATA
OUTPUT_name <- config$seurat2ann$OUTPUT_name
DIM.RED <- config$seurat2ann$DIM.RED
# set working dir
setwd(GIT_DIR)
# create output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_seurat2ann_", timestamp)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
GIT_DIR <- paste0(GIT_DIR, "/")
file.copy(paste0(GIT_DIR,"config.json"), file.path(output, "config.json"))
# read in Seurat object
seurat <- readRDS(file = paste0(DATA_DIR,SEURAT_OBJ))

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
# convert to seurat v3
v <- Version(seurat)
v <- as.numeric(unlist(v))
if(v[1]>=5){
  seurat[["RNA"]]$scale.data <- NULL
  seurat[["RNA3"]] <- as(object = seurat[["RNA"]], Class = "Assay")
  DefaultAssay(seurat) <- "RNA3"
}
# first save seurat as h5 seurat file
SaveH5Seurat(seurat, filename = paste0(output,OUTPUT_name))
# add user- specified dimension reduction (like in the case of integration)
library(hdf5r)
if(DIM.RED != "NA" && DIM.RED != ""){
    print("Adding dimension reduction")
    h5 <- H5File$new(paste0(output,OUTPUT_name), mode = "r+")
    print(h5[["reductions/"]])
    red.path <- paste0("reductions/", DIM.RED, "/cell.embeddings")
    new_dimred <- h5[[red.path]]
    h5[[paste0(DIM.RED,"1")]] <- new_dimred
    h5$close()
}
# then convert to h5ad
Convert(paste0(output,OUTPUT_name), dest = "h5ad")
system(paste("chmod -R 777", output))