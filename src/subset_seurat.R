
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

# set variables
GIT_DIR <- get_wd()
config <- jsonlite::fromJSON(file.path("./config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$subset_seurat$DATA_DIR
}
SEURAT_OBJ <- config$subset_seurat$SEURAT_OBJ
GENE_LIST <- config$subset_seurat$GENE_LIST
DIM.RED <- config$subset_seurat$DIM.RED
ANNOT <- config$subset_seurat$ANNOT
# set up environment and output
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_subset_", timestamp)
dir.create(output, mode = "0777", showWarnings = FALSE)
output <- paste0(output, "/")
file.copy(paste0("config.json"), file.path(output, "config.json"))
source(paste0("src/sc_pipeline_functions.R"))

# load data
print("load data")
seurat.obj <- readRDS(file = paste0(DATA_DIR, SEURAT_OBJ))
seurat.obj

# subset for positively expressed genelist
print("subset for positively expressed genelist")
new.obj <- GetAssayData(object = seurat.obj, layer = "data")[GENE_LIST,]>0
#get the required cells
df <- as.data.frame(new.obj)
# get cells that are true
cellstokeep <- which(apply(df, 2, any))
#subset the required cells
seurat.obj_subset <- subset(seurat.obj, cells = cellstokeep)
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
p5 <- FeaturePlot(seurat.obj_subset, features = GENE_LIST, ncol = length(GENE_LIST),
                           pt.size = 0.5, reduction = DIM.RED) &
                           scale_color_viridis(option="B")
combined_plot <- ((p1 | p2) / (p3 | p4) / p5) + plot_layout(width = c(2, 2, 2),
                                                      heights = c(1, 1, 1))
pdf(file = paste0(output,"seurat_subset.pdf"), width = 8, height = 8)
    print(combined_plot)
dev.off()

# save
print("save")
saveRDS(seurat.obj_subset, file= paste0(output,"seurat.obj_subset.rds"))
# save session info
writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))
system(paste("chmod -R 777", output))