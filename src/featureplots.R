# get environment
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
library(viridis)
use_condaenv(condaenv = '/w5home/bmoore/miniconda3/envs/scRNAseq_best/', required = TRUE)
# set variables
# set variables
GIT_DIR <- "/w5home/bmoore/scRNAseq_downstream/"
config <- fromJSON(file.path(GIT_DIR, "config.json"))
WD <- config$featureplots$WD
SEURAT_OBJ <- config$featureplots$SEURAT_OBJ
GENE_LIST <- config$featureplots$GENE_LIST
ANNOT <- config$featureplots$ANNOT
INPUT_NAME <- config$featureplots$INPUT_NAME
reduction <- config$featureplots$reduction
# set working directory
setwd(WD)
# load seurat object
seurat.obj <- readRDS(file = SEURAT_OBJ)
# load list of marker genes to plot
features <- read.csv(GENE_LIST, header = FALSE, sep = "\t")

# make cluster plots
plot1 <- DimPlot(seurat.obj, reduction = reduction, label = FALSE, 
        pt.size = 0.5, group.by = ANNOT)
plot2 <- DimPlot(seurat.obj, reduction = reduction, label = TRUE, 
        pt.size = 0.5, group.by = "seurat_clusters")

# make for loop a function
plot_function <- function(features, input_name, plot1, plot2) {
  cell_types <- unique(features$V2)
  for (c in seq(1, length(cell_types))) {
    print(c)
    print(cell_types[c])
    # subset features
    features1 <- features[features$V2 == cell_types[c],]
    marker.genes <- as.vector(features1$V1)
    # count to distinguish each plot
    count=1
    # loop to subset and plot
    for (i in seq(1, length(marker.genes), by=12)) {
      j= i+11
      markers1 <- marker.genes[i:j]
      plot3 <- FeaturePlot(seurat.obj, features = markers1, ncol = 3,
                           pt.size = 0.1, reduction = reduction) &
                           scale_color_viridis(option="B") &
                           xlim(c(-0.03,0.04)) & ylim(c(-0.03,0.04))
      combined_plot <- ((plot1 | plot2) / plot3) + plot_layout(width = c(2, 3),
                                                               heights = c(1, 4))
      # make pdf
      pdf(file = paste0("feature_plot_",as.character(input_name),"_", as.character(cell_types[c]),
                        as.character(count),".pdf"), width = 8, height = 11)
      print(combined_plot)
      dev.off()
      count=count+1
    }
  }
} 
# run function: input gene list, name, and two previous plots
plot_function(features,INPUT_NAME,plot1, plot2)
