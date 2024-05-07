
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
WD <- "/w5home/bmoore/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_seurat_mapping_20240125_112552/"
SEURAT_OBJ <- "human_D205_subset_annot.rds"
METADATA_GAMM <- "NA" #"/w5home/bmoore/Gamm_scRNAseq/data/gamm_metadata/gammS2_manual_annot_metadata_c0.5.txt"
OUTPUT_name <- "Reh_human_ortho_labeled-clusters.h5Seurat"

setwd(WD)


# read in Seurat object
seurat <- readRDS(file = SEURAT_OBJ)

if(METADATA_GAMM=="NA"){
  print("No metadata file provided")
} else {
    print("Metadata file provided")
    print(METADATA_GAMM)
    print("Reading in metadata")
    # add in metadata from previous analysis
    metadata.gamm <- read.csv(METADATA_GAMM, row.names = 1,
                            header = TRUE, sep = "\t")
    # check metadata with query data
    print(colnames(seurat@assays$RNA@data)[1:10])
    print(rownames(metadata.gamm)[1:10])
    print(colnames(seurat@meta.data))

    # add metadata to query data
    seurat <- AddMetaData(seurat, metadata.gamm)
    print(colnames(seurat@meta.data))
}

# remove "_1" from metadata
#rownames(metadata.gamm) <- sub("_1", "", rownames(metadata.gamm))
#rownames(metadata.gamm)[1:10]

# first save seurat as h5 seurat file
SaveH5Seurat(seurat, filename = OUTPUT_name)
# then convert to h5ad
Convert(OUTPUT_name, dest = "h5ad")

# can now be read in by scanpy:
# import scanpy
# adata = scanpy.read_h5ad("pbmc3k.h5ad")