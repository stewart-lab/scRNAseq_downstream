# celltype gpt
## load libraries
library(Seurat)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(patchwork)
library(scran)
library(dplyr)

## load config and get output file
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(getwd(), "config.json"))
docker <- config$docker
if(docker=="TRUE"||docker=="true"||docker=="T"||docker=="t"){
    DATA_DIR <- "./data/input_data/"
} else {
    DATA_DIR <- config$celltypeGPT$DATA_DIR
}
SEURAT_OBJ <- config$celltypeGPT$seurat.obj
setwd(GIT_DIR)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("./shared_volume/output_celltypeGPT_", timestamp)
dir.create(output, showWarnings = FALSE)
output <- paste0(output, "/")
GIT_DIR <- paste0(GIT_DIR, "/")
file.copy(paste0(GIT_DIR,"config.json"), file.path(output, "config.json"))
source(paste0(GIT_DIR,"src/sc_pipeline_functions.R"))

#install packages if needed
#install.packages("openai")
#remotes::install_github("Winnie09/GPTCelltype")

## set api key
Sys.setenv(OPENAI_API_KEY = config$celltypeGPT$openAI_key)

## load gpt packages
library(GPTCelltype)
library(openai)

# load seurat object
seurat.obj <- readRDS(file = paste0(DATA_DIR, SEURAT_OBJ))
seurat.obj

## get markers
# to get markers
cluster_name <- config$celltypeGPT$cluster_name
#sce_obj <- as.SingleCellExperiment(seurat.obj)
# convert to single cell
counts_matrix <- seurat.obj[["RNA"]]$counts
data_matrix <- seurat.obj[["RNA"]]$data
# get clusters
clusters <-  seurat.obj[[cluster_name]]
# make sce object
sce_obj <- SingleCellExperiment(list(counts=counts_matrix, logcounts=data_matrix),
    colData=DataFrame(label=clusters))
# update cluster name
cluster_name2 <- paste0("label.",cluster_name)
# score markers
marker_info <- scoreMarkers(sce_obj, sce_obj@colData@listData[[cluster_name2]], full.stats = TRUE)

## get marker dataframe
clusters <- unique(seurat.obj[[cluster_name]])
combined_df <- data.frame(cluster = character(), gene = character(), 
                       rank.logFC.cohen = integer(), median.logFC.cohen = numeric(), 
                       median.AUC = numeric())
for (i in 1:length(clusters[[1]])){
  clust <- as.data.frame(marker_info[[clusters[[1]][i]]])
# Order by median Cohen's d and subset
  ordered <- subset(clust[order(clust$median.logFC.cohen, decreasing = TRUE), ], 
                    median.logFC.cohen > 0.25 & median.AUC > 0.5)
  # get top 100
  top100 <- head(ordered, n = 100)
  new_df <- top100[, c("rank.logFC.cohen", "median.logFC.cohen", "median.AUC")]
  # add cluster to dataframe
  new_df['cluster'] <- clusters[[1]][i]
  # add gene name from rownames
  new_df['gene'] <-  row.names(top100)
  # add new_df to combined_df
  combined_df <- rbind(combined_df, new_df)
}
names(combined_df)[2]<- paste("avg_log2FC")

## annotate with GPT4
# get tissue
tissue_name <- config$celltypeGPT$tissue_name
# Cell type annotation by GPT-4
if(tissue_name=="NA"){
  res <- gptcelltype(combined_df, model = 'gpt-4')
} else {
  res <- gptcelltype(combined_df, model = 'gpt-4', tissuename = tissue_name)
}

## add to object
# make sure clusters are idents
Idents(seurat.obj) <- cluster_name
print(table(Idents(seurat.obj)))
# Assign cell type annotation back to Seurat object
seurat.obj@meta.data$celltype_gpt4 <- as.factor(res[as.character(Idents(seurat.obj))])

## visualize
# Visualize cell type annotation on UMAP
pdf(paste0(output, "celltypeGPT_umap.pdf"), 
  width = 10, height = 8)
print(DimPlot(seurat.obj,group.by='celltype_gpt4'))
dev.off()

## save
saveRDS(seurat.obj, file= paste0(output,"seurat_celltypeGPT_annot.rds"))
writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))