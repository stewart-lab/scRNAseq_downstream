library(CellChat)
library(patchwork)
library(Seurat)
options(stringsAsFactors = FALSE)

# load Seurat object
# for seurat version < 5.0.0 normalized data: seurat_object[["RNA"]]@data # normalized data matrix

seurat_obj <- readRDS(file="")