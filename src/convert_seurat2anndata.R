# load packages
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(Seurat)
library(SeuratObject)
# use_condaenv(condaenv = 'scRNAseq_new', required = TRUE)
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
GIT_DIR <- getwd()
config <- fromJSON(file.path("./config.json"))
docker <- config$docker
if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
  DATA_DIR <- "./data/input_data/"
} else {
  DATA_DIR <- config$seurat2ann$DATA_DIR
}
METADATA_DIR <- config$seurat2ann$METADATA_DIR
DIM.RED <- config$seurat2ann$DIM.RED
# set working dir
setwd(GIT_DIR)
# create output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_seurat2ann_", timestamp)
print(output)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
GIT_DIR <- paste0(GIT_DIR, "/")
file.copy(paste0(GIT_DIR, "config.json"), file.path(output, "config.json"))

# find all .rds files in DATA_DIR
rds_files <- list.files(DATA_DIR, pattern = "\\.rds$", full.names = FALSE)
if (length(rds_files) == 0) {
  stop(paste("No .rds files found in", DATA_DIR))
}
print(paste("Found", length(rds_files), "RDS file(s) in", DATA_DIR))
print(rds_files)

# loop over each .rds file
for (rds_file in rds_files) {
  obj_name <- tools::file_path_sans_ext(rds_file)
  OUTPUT_name <- paste0(obj_name, ".h5Seurat")
  print(paste("=== Processing:", rds_file, "==="))

  # read in Seurat object
  seurat <- readRDS(file = paste0(DATA_DIR, rds_file))

  # extract cluster number from filename (e.g. "seurat_subset_cluster_0_NR.rds" -> "0")
  cluster_num <- sub(".*cluster_(\\d+).*", "\\1", rds_file)

  if (METADATA_DIR == "NA") {
    print("No metadata directory provided")
  } else {
    # build per-cluster metadata path: METADATA_DIR/<cluster_num>/subset.seu.meta.csv
    meta_file <- file.path(METADATA_DIR, cluster_num, "subset.seu.meta.csv")
    if (!file.exists(meta_file)) {
      print(paste("WARNING: Metadata file not found for cluster", cluster_num, ":", meta_file))
      print("Skipping metadata addition for this file")
    } else {
      print(paste("Metadata file found for cluster", cluster_num))
      print(meta_file)
      print("Reading in metadata")
      # add in metadata from previous analysis
      metadata <- read.csv(meta_file,
        row.names = 1,
        header = TRUE, sep = ","
      )
      # check metadata with query data
      print(colnames(seurat@assays$RNA@data)[1:10])
      print(rownames(metadata)[1:10])
      print(colnames(seurat@meta.data))

      # add metadata to query data
      seurat <- AddMetaData(seurat, metadata)
      print(colnames(seurat@meta.data))
    }
  }
  # convert to seurat v3
  v <- Version(seurat)
  v <- as.numeric(unlist(v))
  if (v[1] >= 5) {
    seurat[["RNA"]]$scale.data <- NULL
    seurat[["RNA3"]] <- as(object = seurat[["RNA"]], Class = "Assay")
    DefaultAssay(seurat) <- "RNA3"
  }
  # first save seurat as h5 seurat file
  SaveH5Seurat(seurat, filename = paste0(output, OUTPUT_name))
  # add user- specified dimension reduction (like in the case of integration)
  library(hdf5r)
  if (DIM.RED != "NA" && DIM.RED != "") {
    print("Adding dimension reduction")
    h5 <- H5File$new(paste0(output, OUTPUT_name), mode = "r+")
    print(h5[["reductions/"]])
    red.path <- paste0("reductions/", DIM.RED, "/cell.embeddings")
    new_dimred <- h5[[red.path]]
    h5[[paste0(DIM.RED, "1")]] <- new_dimred
    h5$close()
  }
  # then convert to h5ad
  Convert(paste0(output, OUTPUT_name), dest = "h5ad")
  print(paste("=== Done:", obj_name, ".h5ad ==="))
}
system(paste("chmod -R 777", output))
