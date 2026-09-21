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

# Record R session info (packages + versions) for run reproducibility --
# appended to the run's shared provenance file when invoked via
# run_downstream_toolkit.sh (PROVENANCE_FILE env var), else written
# standalone to the current directory.
.provenance_file <- Sys.getenv("PROVENANCE_FILE", unset = paste0("./sessionInfo_featureplots.R_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
cat(paste0("\n--- R sessionInfo (featureplots.R) ---\n", paste(capture.output(sessionInfo()), collapse = "\n"), "\n"), file = .provenance_file, append = TRUE)
# use_condaenv(condaenv = '/w5home/bmoore/miniconda3/envs/scRNAseq_best/', required = TRUE)
# set variables
# set variables
GIT_DIR <- getwd()
config <- fromJSON(file.path("./config.json"))
docker <- config$docker
if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
  DATA_DIR <- "./data/input_data/"
} else {
  DATA_DIR <- config$featureplots$DATA_DIR
}
SEURAT_OBJ <- config$featureplots$SEURAT_OBJ
GENE_LIST <- config$featureplots$GENE_LIST
ANNOT <- config$featureplots$ANNOT
INPUT_NAME <- config$featureplots$INPUT_NAME
reduction <- config$featureplots$reduction
# set working dir
# setwd(GIT_DIR)
# create output
timestamp <- Sys.getenv("RUN_TIMESTAMP", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
output <- paste0("./shared_volume/output_featureplots_", timestamp)
print(output)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
GIT_DIR <- paste0(GIT_DIR, "/")
file.copy(paste0(GIT_DIR, "config.json"), file.path(output, "config.json"))
if (file.exists(.provenance_file)) {
  file.copy(.provenance_file, file.path(output, basename(.provenance_file)), overwrite = TRUE)
  file.remove(.provenance_file)
}
# load seurat object
seurat.obj <- readRDS(file = paste0(DATA_DIR, SEURAT_OBJ))
# load list of marker genes to plot
features <- read.csv(paste0(GENE_LIST), header = TRUE, sep = "\t")

# make cluster plots
plot1 <- DimPlot(seurat.obj,
  reduction = reduction, label = FALSE,
  pt.size = 0.5, group.by = ANNOT
)
plot2 <- DimPlot(seurat.obj,
  reduction = reduction, label = TRUE,
  pt.size = 0.5, group.by = "seurat_clusters"
)

# make for loop a function
plot_function <- function(features, input_name, plot1, plot2) {
  cell_types <- unique(features$Celltype)
  for (c in seq(1, length(cell_types))) {
    print(c)
    print(cell_types[c])
    # subset features
    features1 <- features[features$Celltype == cell_types[c], ]
    marker.genes <- as.vector(features1$gene)
    # count to distinguish each plot
    count <- 1
    # loop to subset and plot
    for (i in seq(1, length(marker.genes), by = 12)) {
      j <- min(i + 11, length(marker.genes))
      markers1 <- marker.genes[i:j]
      plot3 <- FeaturePlot(seurat.obj,
        features = markers1, ncol = 3,
        pt.size = 0.1, reduction = reduction
      ) # &
      # scale_color_viridis() # &
      # xlim(c(-0.03,0.04)) & ylim(c(-0.03,0.04))
      # pad out to a full 4-row x 3-col grid (12 panels) with blank
      # placeholders so a short final batch doesn't get stretched to
      # fill the space reserved for a full page of feature plots
      n_missing <- 12 - length(markers1)
      if (n_missing > 0) {
        for (p in seq_len(n_missing)) {
          plot3 <- plot3 + patchwork::plot_spacer()
        }
        plot3 <- plot3 + plot_layout(ncol = 3)
      }
      combined_plot <- ((plot1 | plot2) / plot3) + plot_layout(
        width = c(2, 3),
        heights = c(1, 4)
      )
      # make pdf
      pdf(file = paste0(
        output, "feature_plot_", as.character(input_name), "_", as.character(cell_types[c]),
        as.character(count), ".pdf"
      ), width = 8, height = 11)
      print(combined_plot)
      dev.off()
      count <- count + 1
    }
  }
}
# run function: input gene list, name, and two previous plots
plot_function(features, INPUT_NAME, plot1, plot2)
# save session info
print("save session info")
writeLines(capture.output(sessionInfo()), paste0(output, "sessionInfo.txt"))
system(paste("chmod -R 777", output))
