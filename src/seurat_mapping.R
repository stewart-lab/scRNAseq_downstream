
# Seurat tansfer mapping for cross species analysis
# do preprocessing on human and pig data (preprocess_crossspecies.Rmd)
# both reference and query should undergo normalization and find variable features,
# then only reference is scaled and undegoes dimentionality reduction

# Load libraries, set wd, and source functions
# set Rterm: /w5home/bmoore/miniconda3/envs/scRNAseq_best/bin/R
#conda upgrade -c conda-forge --all
#BiocManager::install(version = "3.18")
#BiocManager::install("ComplexHeatmap")
library(dplyr)
library(Seurat)
library(patchwork)
library(scran)
library(BiocParallel)
library(DropletUtils)
library(cowplot)
library(harmony)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
library(scPred)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
# set variables
GIT_DIR <- getwd()
jsonlite::fromJSON(file.path(getwd(), "config.json"))
source(paste0("src/sc_pipeline_functions.R"))
WD <- config$seurat_mapping$DATA_DIR
REF.SEURAT <- config$seurat_mapping$REF.SEURAT
QUERY.SEURAT <- config$seurat_mapping$QUERY.SEURAT
# set up environment and output
use_condaenv("scRNAseq_new", required=TRUE)
setwd(WD)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output_seurat_mapping_", timestamp)
dir.create(output, showWarnings = FALSE)
file.copy(paste0(GIT_DIR, "/config.json"), file.path(output, "config.json"))
config <- fromJSON(file.path(output, "config.json"))
output <- paste0(output, "/")

# read in seurat objects that were preprocessed

query.seurat <- readRDS(file = QUERY.SEURAT)
ref.seurat <- readRDS(file = REF.SEURAT)
print(query.seurat)
print(ref.seurat)

# get metadata- if needed
# get metadata if cell type annotation is needed
get_meta <- config$seurat_mapping$get_metadata$get_meta
if (get_meta=="TRUE"){
    # set type "ref" or "query" ref= metadata.file1, query= metadata.file2
    ref.seurat <- get_metadata(ref.seurat, "ref", output, "seurat_mapping")
    query.seurat <- get_metadata(query.seurat, "query", output, "seurat_mapping")
}
print(colnames(ref.seurat@meta.data))
print(colnames(query.seurat@meta.data))

# some visualizations

# feature plots
count_and_feature_plots(ref.seurat,query.seurat)

# ref umap and subset metadata

ref.seurat.sub <- visualize_and_subset_ref(ref.seurat,output,"seurat_mapping")


# mapping and annotating query datasets

obj.list <- transfer_anchors(ref.seurat.sub, query.seurat)
print(obj.list)
query.seurat <- obj.list[[1]]
anchors <- obj.list[[2]]



# projection of a query onto the reference UMAP structure.

obj.list <- project_query_on_ref(ref.seurat.sub, query.seurat, anchors)
ref.seurat <- obj.list[[1]]
query.seurat <- obj.list[[2]]

# visualize probabilities

#  probabilities of each cell in the @meta.data slot of the query Seurat object.
# get labels
labels<- as.vector(unique(colnames(query.seurat@meta.data)))
library(stringr)
labels <- labels[startsWith(labels, "prediction.score.")]
labels<-labels[! labels %in% c('prediction.score.max')]
# visualize the probabilities over the UMAP plot:
pdf(paste0(output, "query_celltype_prediction_featureplot.pdf"), width = 10, height = 10)
print(FeaturePlot(query.seurat, labels, keep.scale = "all"))
dev.off()

# get prediction score histograms

for(l in seq_along(1:length(labels))){
    # get celltype from prediction score label
    l.list <-  strsplit(labels[l], split = "score.", fixed = TRUE)
    celltype <- sapply(l.list, tail, 1)
    print(celltype)
    # fix cones
    if(celltype=="L.M.cone"){
        celltype <- "L/M cone"
    }else if(celltype=="S.cone"){
        celltype <- "S cone"
    }else if(celltype=="Prog.Glia"){
        celltype <- "Prog/Glia"
    }
    # check if selltype is in predicted.id
    if(celltype %in% query.seurat$predicted.id){
        # subset seurat object for cell type
        query.seurat.sub <- subset(query.seurat, subset = predicted.id == celltype)
        # get prediction scores for label
        data <- FetchData(object = query.seurat.sub, vars =labels[l])
        # get mean and median
        x.mean <- mean(data[[1]])
        x.med <- round(median(data[[1]]), digits = 3)
        # remove unwanted characters before writing
        celltype <- str_replace_all(celltype," ","")
        celltype <- str_replace_all(celltype,"/","")
        # plot histogram
        pdf(file = paste0(output, celltype, "_pred.score_histplot.pdf"), width = 3.5, height = 3.5, 
        bg = "white")
        print(hist(data[[1]], plot=TRUE, xlim = c(0,1), 
            main=paste0(celltype,"_pred.score median: ",x.med)))
        dev.off()
        print(paste0("mean: " ,x.mean))
        print(paste0("median: " ,x.med))
    }else{
        print(paste0(celltype, " not in predicted ids"))
    }
    
}


# get comparison between predicted and manual annotations

query.ref.data.table <- get_manual_comparison(query.seurat)
print(query.ref.data.table)

# make heatmap

hm <- heatmap_func(query.ref.data.table, output)

saveRDS(query.seurat, file= paste0(output,"query_seurat_pred.rds"))
saveRDS(ref.seurat, file= paste0(output,"human_seurat_map.rds"))

# save session info

writeLines(capture.output(sessionInfo()), paste0(output,"sessionInfo.txt"))

#query.seurat <- readRDS(file = "/w5home/bmoore/scRNAseq/GAMM/human_data/reh_cellrep_2020/output_seurat_mapping_20230804_143108/gamms2_rpca_pred.rds")
