# phate
print("load packages")
# install.packages("phateR")
# pip install --user phate
#if (!require(viridis)) install.packages("viridis")
#if (!require(ggplot2)) install.packages("ggplot2")
#if (!require(readr)) install.packages("readr")
#if (!require(Rmagic)) install.packages("Rmagic")
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(tidyverse)
#use_condaenv(condaenv = 'phate_env', required = TRUE)
print("set variables")
# set variables
#GIT_DIR <- "/w5home/bmoore/scRNAseq_downstream/"
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "/config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$phate$DATA_DIR
}
#file.copy(paste0("config.json"), file.path(DATA_DIR, "config_phate.json"))
SEURAT.FILE <- config$phate$SEURAT_OBJ
gene <- config$phate$gene
knn  <- config$phate$knn
decay  <- config$phate$decay
tset  <- config$phate$t
gamma  <- config$phate$gamma
embed_key  <- config$phate$embed_key
ANNOT  <- config$phate$ANNOT
# set up environment and output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_phate_", timestamp)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(paste0("config.json"), file.path(output, "config.json"))


print("load seurat")
# read in seurat object
seurat.obj <- readRDS(file = paste0(DATA_DIR, SEURAT.FILE))
seurat.obj <- UpdateSeuratObject(seurat.obj)

print("get matrix")
# get matrix from seurat object
Xmatrix <- LayerData(seurat.obj, assay="RNA", layer ="data")

# transpose - phate requires cells as rownames and genes as colnames
Xmatrix_t <- t(Xmatrix)
print(Xmatrix[1:5,1:10])
print(Xmatrix_t[1:5,1:10])

# convert to regular matrix
Xmatrix_t <- as.matrix(Xmatrix_t)
storage.mode(Xmatrix_t) <- "double"

print("run phate")
# run PHATE
# options: knn : Number of nearest neighbors (default: 5)
# decay : Alpha decay (default: 40)
# t : Number of times to power the operator (default: ‘auto’)
# gamma : Informational distance constant between -1 and 1 (default: 1)
Xmatrix_PHATE <- phate(Xmatrix_t, knn=knn, decay=decay, t=tset, gamma=gamma)
Xmatrix_PHATE

# visualize phate with a gene
print("visualize1")
p1 <- ggplot(Xmatrix_PHATE) +
        geom_point(aes(PHATE1, PHATE2, color=Xmatrix_t[,gene])) +
        labs(color=gene) +
        scale_color_viridis(option="B")
pdf(paste0(output,"phate_plot_",gene,".pdf"),bg = "white")
print(p1)
dev.off()


# add embedding to seurat object
print("add to seurat obj")
new_reduction <- CreateDimReducObject(embeddings = Xmatrix_PHATE$embedding, 
                                key = embed_key)
seurat.obj[[embed_key]] <- new_reduction

# visualize samples
print("visualize2")
plot1 <- DimPlot(seurat.obj, reduction = embed_key, label = TRUE, 
        pt.size = 0.5, group.by = "orig.ident")
pdf(paste0(output,"orig_ident_",embed_key,".pdf"), bg = "white", width = 8, height = 6)
print(plot1)
dev.off()

# visualize clusters
p1 <- DimPlot(seurat.obj, reduction = embed_key, label = TRUE, 
        pt.size = 0.5, group.by = ANNOT)
p2 <- DimPlot(seurat.obj, reduction = embed_key, label = TRUE, 
        pt.size = 0.5, group.by = "cca_clusters")
pdf(paste0(output,"labels-clusters_",embed_key,".pdf"), bg = "white", width = 8, height = 6)
print(p1+p2)
dev.off()

p3 <- DimPlot(seurat.obj, reduction = embed_key, label = TRUE, 
        pt.size = 0.5, group.by = ANNOT) + NoLegend()
p4 <- DimPlot(seurat.obj, reduction = embed_key, label = TRUE, 
        pt.size = 0.5, group.by = "seurat_clusters") + NoLegend()
pdf(paste0(output,"labels-clusters_",embed_key,"_2.pdf"), bg = "white", width = 8, height = 6)
print(p3+p4)
dev.off()

print("save object")
saveRDS(seurat.obj, file=paste0(output,"seurat_obj_labeled_",embed_key,".rds"))
system(paste("chmod -R 777", output))