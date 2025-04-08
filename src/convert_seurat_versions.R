# get environment
library(dplyr)
library(Seurat)
library(SeuratData)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)


setwd("/w5home/bmoore/scRNAseq/GAMM/DATA/GEO_submission/processed_data_files/")
SEURAT.FILE <- "seurat_obj_labeled_S2.rds"
seurat.obj <- readRDS(file = SEURAT.FILE)
print(seurat.obj)
all.genes <- rownames(seurat.obj)
seurat.obj <- ScaleData(seurat.obj, features = all.genes)
# convert a v5 assay to a v3 assay
seurat.obj[["RNA3"]] <- as(object = seurat.obj[["RNA"]], Class = "Assay")
print(seurat.obj)
saveRDS(seurat.obj, file = "seurat_obj_S2_v3.rds")