# clustifyr analysis for automatic cell type annotation with marker list or reference
# To run clustifyr, you need clustered seurat object

### load libraries ###
print("load libraries")
library(reticulate)
#use_condaenv("/w5home/bmoore/miniconda3/envs/scRNAseq_best")
library(clustifyr)
library(ggplot2)
library(cowplot)
library(Seurat)
library(ComplexHeatmap)
library(rmarkdown)
library(purrr)
library(jsonlite)
library(dplyr)
### set variables ###
print("set variables")
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$clustifyr$DATA_DIR
}
REF.SEURAT <-  config$clustifyr$REF.SEURAT # if NA, only marker list is used
QUERY.SEURAT <- config$clustifyr$QUERY.SEURAT
cluster_name <- config$clustifyr$cluster_name
### set working directory and output ###
setwd(GIT_DIR)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_clustifyr_", timestamp)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
GIT_DIR <- paste0(GIT_DIR, "/")
## copy config to output
file.copy(paste0(GIT_DIR,"config.json"), file.path(output, "config.json"))
## source config and functions
source(paste0(GIT_DIR,"src/sc_pipeline_functions.R"))

# load data set
print("load data")
seurat.obj <- readRDS(file = paste0(DATA_DIR, QUERY.SEURAT))
# load reference
if(REF.SEURAT!="NA"){
ref.seurat<- readRDS(file = paste0(DATA_DIR, REF.SEURAT))
} else {
  print("No reference, using marker list only")
}

# rename idents
print("rename idents and remove NAs")
# check metadata
colnames(seurat.obj@meta.data)
# stash idents
seurat.obj[["idents"]] <- Idents(object = seurat.obj)
# make clusters idents
Idents(object = seurat.obj) <- cluster_name
table(Idents(object = seurat.obj))
# convert to character
seurat.obj@meta.data[cluster_name] <- mutate_if(seurat.obj@meta.data[cluster_name], is.factor, as.character)
unique(seurat.obj[[cluster_name]])
# remove NAs
#seurat.obj<- subset(seurat.obj[[cluster_name]] != "NA")
if(cluster_name=="seurat_clusters"){
  seurat.obj<- subset(seurat.obj, subset = seurat_clusters != "NA")
} else if (cluster_name=="seurat_clusters2"){
  seurat.obj<- subset(seurat.obj, subset = seurat_clusters2 != "NA")
} else if (cluster_name=="cca_clusters"){
  seurat.obj<- subset(seurat.obj, subset = cca_clusters != "NA")
} else {
  print(paste0("Not able to remove NAs from clusters: ",cluster_name,". 
        Clustifyr may throw an error if NAs are present in clusters."))
}

# annotate with marker list
# convert to single cell experiment
SCE <- as.SingleCellExperiment(seurat.obj)
# run clustifyr
seurat.obj<- annotate_with_clustifyR(seurat.obj, SCE, output, "clustifyr")
# with reference
if(REF.SEURAT!="NA"){
  #visualize and subset
  ref.seurat.sub <- visualize_and_subset_ref(ref.seurat, output, "clustifyr")
  # annotate with clustifyr
  clustifyr.obj <- annotate_with_clustifyR_ref(ref.seurat.sub, seurat.obj, SCE, output, "clustifyr")

} else {
  print("no references, quitting")
  q()
}

writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))
system(paste("chmod -R 777", output))