
# load packages
print("load packages")
library(Seurat)
library(patchwork)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(reticulate)
library(multtest)
library(metap)
library(harmony)
library(SeuratWrappers)
library(ggplot2)
library(scran)
library(clustifyr)
library(SeuratDisk)
library(viridis)

# set variables
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path("./config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$subset_seurat$DATA_DIR
}
SEURAT_OBJ <- config$subset_seurat$SEURAT_OBJ
SUBSET_MODE <- config$subset_seurat$SUBSET_MODE
GENE_LIST <- config$subset_seurat$GENE_LIST
METADATA_COLUMN <- config$subset_seurat$METADATA_COLUMN
IDENTS <- config$subset_seurat$IDENTS
DIM.RED <- config$subset_seurat$DIM.RED
ANNOT <- config$subset_seurat$ANNOT
if (is.null(SUBSET_MODE) || SUBSET_MODE == "") {
    SUBSET_MODE <- "gene"
}
# set up environment and output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_subset_", timestamp)
print(output)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(paste0("config.json"), file.path(output, "config.json"))
source(paste0("src/sc_pipeline_functions.R"))

# load data
print("load data")
seurat.obj <- readRDS(file = paste0(DATA_DIR, SEURAT_OBJ))
seurat.obj

# subset cells by gene expression or metadata
if (SUBSET_MODE == "gene") {
    print("subset for positively expressed genelist")
    if (length(GENE_LIST) == 0) {
        stop("SUBSET_MODE is 'gene' but GENE_LIST is empty")
    }
    missing_genes <- setdiff(GENE_LIST, rownames(seurat.obj))
    if (length(missing_genes) > 0) {
        stop(paste("GENE_LIST contains genes not found in the Seurat object:",
                   paste(missing_genes, collapse = ", ")))
    }
    new.obj <- GetAssayData(object = seurat.obj, layer = "data")[GENE_LIST, ] > 0
    df <- as.data.frame(new.obj)
    cellstokeep <- which(apply(df, 2, any))
    seurat.obj_subset <- subset(seurat.obj, cells = cellstokeep)
} else if (SUBSET_MODE == "metadata") {
    print(paste("subset by metadata column:", METADATA_COLUMN))
    if (is.null(METADATA_COLUMN) || METADATA_COLUMN == "") {
        stop("SUBSET_MODE is 'metadata' but METADATA_COLUMN is not set")
    }
    if (length(IDENTS) == 0) {
        stop("SUBSET_MODE is 'metadata' but IDENTS is empty")
    }
    if (!(METADATA_COLUMN %in% colnames(seurat.obj@meta.data))) {
        stop(paste("METADATA_COLUMN not found in Seurat object metadata:",
                   METADATA_COLUMN))
    }
    Idents(object = seurat.obj) <- METADATA_COLUMN
    print(table(Idents(seurat.obj)))
    available_idents <- unique(as.character(Idents(seurat.obj)))
    missing_idents <- setdiff(IDENTS, available_idents)
    if (length(missing_idents) > 0) {
        stop(paste("IDENTS contains values not found in", METADATA_COLUMN, ":",
                   paste(missing_idents, collapse = ", ")))
    }
    seurat.obj_subset <- subset(x = seurat.obj, idents = IDENTS)
} else {
    stop(paste("Unknown SUBSET_MODE:", SUBSET_MODE, "- use 'gene' or 'metadata'"))
}
seurat.obj_subset

# plot
print("plot")
p1 <- DimPlot(seurat.obj_subset, reduction = DIM.RED, label = TRUE, 
        pt.size = 0.5, group.by = ANNOT)
p2 <- DimPlot(seurat.obj, reduction = DIM.RED, label = TRUE, 
        pt.size = 0.5, group.by = ANNOT)
p3 <- DimPlot(seurat.obj_subset, reduction = DIM.RED, 
        pt.size = 0.5, group.by = "orig.ident")
p4 <- DimPlot(seurat.obj_subset, reduction = DIM.RED, label = TRUE, 
        pt.size = 0.5, group.by = "seurat_clusters")
if (SUBSET_MODE == "gene" && length(GENE_LIST) > 0) {
    p5 <- FeaturePlot(seurat.obj_subset, features = GENE_LIST, ncol = length(GENE_LIST),
                               pt.size = 0.5, reduction = DIM.RED) &
                               scale_color_viridis(option = "B")
    combined_plot <- ((p1 | p2) / (p3 | p4) / p5) + plot_layout(width = c(2, 2, 2),
                                                          heights = c(1, 1, 1))
} else {
    combined_plot <- (p1 | p2) / (p3 | p4)
}
pdf(file = paste0(output,"seurat_subset.pdf"), width = 8, height = 8)
    print(combined_plot)
dev.off()

# save
print("save")
saveRDS(seurat.obj_subset, file= paste0(output,"seurat.obj_subset.rds"))
# save session info
writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))
system(paste("chmod -R 777", output))