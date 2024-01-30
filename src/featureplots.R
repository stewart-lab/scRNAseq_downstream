# get environment
library(dplyr)
library(Seurat)
library(patchwork)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
use_condaenv(condaenv = 'scRNAseq_best', required = TRUE)
# set variables
WD <- "/w5home/bmoore/scRNAseq/GAMM/GAMM_S1/output_20230921_142919/"
SEURAT_OBJ <- "GAMM_S1_clabeled-clusters_0.5.rds"
GENE_LIST <- "../../known_markers/genes_for_featureplots.txt"
ANNOT <- "CellType_manual" # annotation to use for plotting, "orig.ident", "CellType_manual", "CellType"
INPUT_NAME <- "Gamm_S1"
# set working directory
setwd(WD)
# load seurat object
seurat.obj <- readRDS(file = SEURAT_OBJ)
# load list of marker genes to plot
features <- read.csv(GENE_LIST, 
                          header = FALSE, sep = "\t")

# make cluster plots
plot1 <- UMAPPlot(seurat.obj, group.by = ANNOT, label = F)
plot2 <- UMAPPlot(seurat.obj, group.by = "seurat_clusters", label = T)

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
    for (i in seq(1, length(marker.genes), by=15)) {
      #print(i)
      j= i+14
      #print(j)
      markers1 <- marker.genes[i:j]
      plot3 <- FeaturePlot(seurat.obj, features = markers1, ncol = 3, 
                           pt.size = 0.1)
      combined_plot <- ((plot1 | plot2) / plot3) + plot_layout(width = c(2, 3),
                                                               heights = c(1, 2))
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
