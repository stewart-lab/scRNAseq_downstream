
# load libraries
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
library(Azimuth)
library(ggplot2)
library(scran)
library(clustifyr)
library(SeuratDisk)
#use_python("~/miniconda3/envs/scRNAseq_new2/bin/python")
#use_condaenv("/w5home/bmoore/miniconda3/envs/scRNAseq_best")

# set variables
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "config.json"))
DATA_DIR <- config$seurat_integration$DATA_DIR
filename_list <- config$seurat_integration$filename_list
new_names <- config$seurat_integration$orig_ident_rename
# set up environment and output
setwd(GIT_DIR)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output_seuratintegrate_", timestamp)
dir.create(output, showWarnings = FALSE)
output <- paste0(output, "/")
GIT_DIR <- paste0(GIT_DIR, "/")
file.copy(file.path(paste0(GIT_DIR,"config.json")), file.path(paste0("./", 
          output,"config.json")), overwrite = TRUE)
source(paste0(GIT_DIR,"src/sc_pipeline_functions.R"))

# load data
setwd(DATA_DIR)
# Create an empty list to store the imported objects
imported_objects <- list()
print("loading data")
# Loop through the filenames
for (i in seq_along(filename_list)) {
  # Extract the base filename without extension
  base_name <- tools::file_path_sans_ext(basename(filename_list[i]))
  
  # Import the file
  imported_object <- readRDS(filename_list[i])
  
  # Assign the imported object to the list with the base filename as the name
  imported_objects[[base_name]] <- imported_object
}
# reset dir
setwd(GIT_DIR)
# merge into one object
# merge
# Convert the list to a vector
seurat_vector <- unlist(imported_objects)
vector<- c()
for(i in 2:length(seurat_vector)){
    vector <- c(vector, seurat_vector[[i]])
}
# This will be used to give unique identifiers to cells from each object
object_names <- names(imported_objects)
if (is.null(object_names)) {
  object_names <- paste0("Object", seq_along(imported_objects))
}

# Now, we'll merge the objects
print("merging seurat objects")
merged_seurat <- merge(seurat_vector[[1]], y = vector, 
                        add.cell.ids = object_names, 
                        project = config$seurat_integration$project_name)
# make sure all objects are class v5
print("convert to seurat 5")
merged_seurat[["RNA5"]] <- as(object = merged_seurat[["RNA"]], Class = "Assay5")
DefaultAssay(merged_seurat) <- "RNA5"
merged_seurat

# stash idents and reset identity
print("rename idents")
merged_seurat[["CellType_combined"]] <- Idents(object = merged_seurat)
table(merged_seurat$CellType_combined)
# Set identity classes to an existing column in meta data
Idents(object = merged_seurat) <- "orig.ident"
# rename
#new_names <- new_names
names(new_names) <- levels(merged_seurat)
merged_seurat <- RenameIdents(object = merged_seurat, 
                                new_names)

merged_seurat[["orig.ident2"]] <- Idents(object = merged_seurat)

# split layers
print("split layers")
merged_seurat[["RNA5"]] <- split(merged_seurat[["RNA5"]], f = merged_seurat$orig.ident2)
merged_seurat
table(merged_seurat$orig.ident2)

# feature selection, scale, dimensionality reduction on combined data
print("feature selection, scale, dimensionality reduction on combined data")
merged_seurat <- feature_selection(merged_seurat, "integration")
# scale data
merged_seurat<- scale_data(merged_seurat, output, "integration")
# run PCA
merged_seurat<- run_and_visualize_pca(merged_seurat, output, "integration")
# run UMAP
merged_seurat<- run_umap(merged_seurat, output, "integration")

# integrate data
print("integrate data")
integration_type <- config$seurat_integration$integration_type
if(integration_type == "cca"){
  merged_seurat <- IntegrateLayers(
    object = merged_seurat, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca",
    k.weight = 96,
    verbose = FALSE
  )
} else if(integration_type == "rpca"){
  merged_seurat <- IntegrateLayers(
    object = merged_seurat, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    verbose = FALSE
  )
} else if(integration_type == "harmony"){
  merged_seurat <- IntegrateLayers(
    object = merged_seurat, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE)
} else if(integration_type == "mnn"){
  merged_seurat <- IntegrateLayers(
    object = merged_seurat, method = FastMNNIntegration,
    new.reduction = "integrated.mnn",
    verbose = FALSE)
} else if(integration_type == "scvi"){
  merged_seurat <- IntegrateLayers(
    object = merged_seurat, method = scVIIntegration,
    new.reduction = "integrated.scvi",
    conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE
  )
} else {
  print("Need to specify integration type (cca, rpca, harmony, mnn, or scvi)")
}

# cluster integrated data and run umap again

# run clustering
## Note: you must have the reduction in the config file set to new reduction from above in order for this to work- like integrated.cca
print("clustering")
merged_seurat <- perform_clustering(merged_seurat, output, "integration")

# run umap again to visualize clusters in integratied space
cluster_name <- config$seurat_integration$perform_clustering$cluster_name
if(integration_type == "cca"){
  merged_seurat <- RunUMAP(merged_seurat, reduction = "integrated.cca", 
                        dims = 1:30, reduction.name = "umap.cca")
  # visualize
  p1 <- DimPlot(merged_seurat,reduction = "umap.cca",
    group.by = c("orig.ident2", "CellType_combined", cluster_name),
    combine = FALSE, label.size = 2
  )
  pdf(paste0(output,"umap_merged_cca.pdf"))
  print(p1)
  dev.off()
} else if(integration_type == "rpca"){
  merged_seurat <- RunUMAP(merged_seurat, reduction = "integrated.rpca", 
                        dims = 1:30, reduction.name = "umap.rpca")
  # visualize
  p1 <- DimPlot(merged_seurat,reduction = "umap.rpca",
    group.by = c("orig.ident", "CellType", cluster_name),
    combine = FALSE, label.size = 2
  )
  pdf(paste0(output,"umap_merged_rpca.pdf"))
  print(p1)
  dev.off()
} else if(integration_type == "harmony"){
  merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", 
                        dims = 1:30, reduction.name = "umap.harmony")
  # visualize
  p1 <- DimPlot(merged_seurat,reduction = "umap.harmony",
    group.by = c("orig.ident", "CellType", cluster_name),
    combine = FALSE, label.size = 2
  )
  pdf(paste0(output,"umap_merged_harmony.pdf"))
  print(p1)
  dev.off()
} else if(integration_type == "mnn"){
  merged_seurat <- RunUMAP(merged_seurat, reduction = "integrated.mnn", 
                        dims = 1:30, reduction.name = "umap.mnn")
  # visualize
  p1 <- DimPlot(merged_seurat,reduction = "umap.mnn",
    group.by = c("orig.ident", "CellType", cluster_name),
    combine = FALSE, label.size = 2
  )
  pdf(paste0(output,"umap_merged_mnn.pdf"))
  print(p1)
  dev.off()
} else if(integration_type == "scvi"){
  merged_seurat <- RunUMAP(merged_seurat, reduction = "integrated.scvi", 
                        dims = 1:30, reduction.name = "umap.scvi")
  # visualize
  p1 <- DimPlot(merged_seurat,reduction = "umap.scvi",
    group.by = c("orig.ident", "CellType", cluster_name),
    combine = FALSE, label.size = 2
  )
  pdf(paste0(output,"umap_merged_scvi.pdf"))
  print(p1)
  dev.off()
}

# change metadata to characters before writing
i <- sapply(merged_seurat@meta.data, is.factor)
merged_seurat@meta.data[i] <- lapply(merged_seurat@meta.data[i], as.character)

# join data before saving
print("join data and save")
merged_seurat <- JoinLayers(merged_seurat)
merged_seurat
# save
#saveRDS(merged_seurat, file= paste0(output,"merged_seurat.rds"))
# write as h5ad obj
# convert assay to V3 first
print("convert to v3 to save as h5 file")
merged_seurat[["RNA3"]] <- as(object = merged_seurat[["RNA5"]], Class = "Assay")
DefaultAssay(merged_seurat) <- "RNA3"

# Convert to SingleCellExperiment
# first get counts matrix
counts_matrix <- merged_seurat[["RNA3"]]$counts
data_matrix <- merged_seurat[["RNA3"]]$data
# get clusters
clusters <-  merged_seurat[[cluster_name]]
# get embeddings
dim_data_cca <- Embeddings(merged_seurat, reduction = "integrated.cca")
dim_data_umap <- Embeddings(merged_seurat, reduction = "umap.cca")
# make sce object
sce <- SingleCellExperiment(list(counts=counts_matrix, logcounts=data_matrix),
    colData=DataFrame(label=clusters))
reducedDims(sce) <- list(PCA=dim_data_cca, UMAP=dim_data_umap)

# first save seurat as h5 seurat file
#SaveH5Seurat(merged_seurat, filename = paste0(output,"merged_seurat2.h5Seurat"))
# # then convert to h5ad
#Convert(paste0(output,"merged_seurat2.h5Seurat"), dest = "h5ad")

# now perform DE analysis
# run score and plot markers
print("run DE analysis")
annot_df <- score_and_plot_markers(merged_seurat, sce, output, "integration")

# annotate
print("annotate")
if (config$seurat_integration$score_and_plot_markers$known_markers) {
    new_df_ordered <- annot_df[order(as.numeric(annot_df$Cluster)), ]
    # Get new cluster names ordered by cluster number
    new_cluster_ids <- new_df_ordered$Cell.type
    merged_seurat <- annotate_clusters_and_save(merged_seurat, new_cluster_ids, output, "integration")
  } else {
   print("no manual annotation because no known marker set")
  }

# annotate with clustifyr
# convert to single cell experiment
# SCE <- as.SingleCellExperiment(merged_seurat)
# # run clustifyr
# print("annotate with clustifyr")
# annotate_with_clustifyR(merged_seurat, SCE, output, "integration")
# # save session info
writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))