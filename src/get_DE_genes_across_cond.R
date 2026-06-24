library(Seurat)
library(SeuratData)
library(patchwork)
library(scran)
library(BiocParallel)
library(cowplot)
library(reticulate)
library(purrr)
library(jsonlite)
library(rmarkdown)
library(ggplot2)
library(scater)
library(SingleCellExperiment)
library(tidyverse)
library(DESeq2)
library(dplyr)

#### load config ####
# script_dir <- dirname(normalizePath(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][[1]] |> sub("--file=", "", x = _)))
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path("./config.json"))
docker <- config$docker
if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
  DATA_DIR <- "./data/input_data/"
} else {
  DATA_DIR <- config$de_cond$DATA_DIR
}

seurat_file    <- config$de_cond$seurat_object
cond_col       <- config$de_cond$group_by$condition
replicate_col  <- config$de_cond$group_by$sample_replicates
cluster_col    <- config$de_cond$group_by$cluster_column

#### output directory ####
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0(DATA_DIR, "DEseq2_", timestamp)
print(output)
dir.create(output, mode = "0777", showWarnings = FALSE)
# copy config "./",
file.copy(file.path(paste0(GIT_DIR,"/config.json")), file.path(paste0( 
          output,"/config.json")), overwrite = TRUE)
# source functions
source(paste0(GIT_DIR,"/src/sc_pipeline_functions.R"))

#### read in merged seurat object ####
print("loading data")
merged_seurat <- readRDS(file = file.path(DATA_DIR, seurat_file))
print(merged_seurat)
print(colnames(merged_seurat@meta.data))

#### plot umap ####
pdf(paste0(output, "/seurat-clusters_umap.pdf"), bg = "white", width = 10, height = 6)
print(DimPlot(merged_seurat,
  reduction = "umap", label = TRUE,
  pt.size = 0.5, group.by = cluster_col
))
dev.off()

pdf(paste0(output, "/condition_umap.pdf"), bg = "white", width = 10, height = 6)
print(DimPlot(merged_seurat,
  reduction = "umap", label = TRUE,
  pt.size = 0.5, group.by = cond_col
))
dev.off()

#### pseudobulk aggregation ####
bulk_seurat <- AggregateExpression(merged_seurat,
  assays = "RNA",
  group.by = c(cond_col, replicate_col, cluster_col),
  return.seurat = TRUE
)
Cells(bulk_seurat)

bulk_seurat$celltype.comp <- paste(bulk_seurat[[cluster_col, drop = TRUE]], bulk_seurat[[cond_col, drop = TRUE]], sep = "_")
Idents(bulk_seurat) <- "celltype.comp"

#### output expression matrices ####
df <- as.data.frame(bulk_seurat[["RNA"]]$data)
write.csv(df, file = paste0(output, "/bulk_norm_expr_data.csv"))

df.counts <- as.data.frame(bulk_seurat[["RNA"]]$counts)
write.csv(df.counts, file = paste0(output, "/bulk_counts_expr_data.csv"))

pseudo_bulk_matrix_avg <- AverageExpression(merged_seurat,
  assays = "RNA",
  group.by = c(cond_col, replicate_col, cluster_col),
  return.seurat = FALSE
)
df.avg <- as.data.frame(pseudo_bulk_matrix_avg$RNA)
df.avg <- cbind(Gene = rownames(df.avg), df.avg)
write.csv(df.avg, file = paste0(output, "/bulk_avg_expr_data.csv"))

#### loop: R vs NR comparisons ####
merged_seurat$celltype.comp <- paste(merged_seurat[[cluster_col, drop = TRUE]], merged_seurat[[cond_col, drop = TRUE]], sep = "_")
comparisons <- c("2Wk", "4Wk", "6Wk")
celltypes <- as.vector(unique(merged_seurat[[cluster_col, drop = TRUE]]))

for (i in seq_along(celltypes)) {
  celltype <- celltypes[i]
  celltype1 <- if (grepl("_", celltype, fixed = TRUE)) gsub("_", "-", celltype) else celltype
  print(celltype1)

  for (j in seq_along(comparisons)) {
    week <- comparisons[j]
    print(week)

    tryCatch({
      bulk.de <- FindMarkers(
        object = bulk_seurat,
        ident.1 = paste0(celltype1, "_g", week, "-R"),
        ident.2 = paste0(celltype1, "_g", week, "-NR"),
        test.use = "DESeq2", min.cells.feature = 2,
        min.cells.group = 2,
        slot = "counts"
      )
      bulk.de <- as.data.frame(bulk.de)
      bulk.de <- subset(bulk.de, p_val_adj < 0.05)
      bulk.de <- subset(bulk.de, avg_log2FC >= 0.25 | avg_log2FC <= -0.25)
      print(head(bulk.de))
      bulk.de <- bulk.de[with(bulk.de, order(p_val_adj, -avg_log2FC)), ]
      write.table(bulk.de, file = paste0(output, "/", celltype1, "_", week, "_R-NR_DEseq2_genes_lfc0.25.txt"), sep = "\t")
    },
    error = function(cond) {
      message(conditionMessage(cond))
      bulk.de <<- data.frame()
    })

    if (dim(bulk.de)[1] != 0) {
      n <- min(10, dim(bulk.de)[1])
      top.10.genes <- unique(na.omit(rownames(bulk.de)[1:n]))

      seurat.obj_subset <- subset(merged_seurat, subset = .data[[cluster_col]] == celltype)

      pdf(paste0(output, "/violinplot_", celltype, "_", week, "_R-NR_DEseq2_lfc0.25.pdf"), bg = "white")
      print(VlnPlot(seurat.obj_subset,
        features = top.10.genes,
        idents = c(paste0(week, "-NR"), paste0(week, "-R")),
        group.by = cond_col
      ))
      dev.off()

      pdf(paste0(output, "/violinplot_", celltype1, "_", week, "_R-NR_DEseq2_lfc0.25_bulk.pdf"), bg = "white")
      print(VlnPlot(bulk_seurat,
        features = top.10.genes,
        idents = c(paste0(celltype1, "_g", week, "-NR"), paste0(celltype1, "_g", week, "-R")),
        group.by = cond_col
      ))
      dev.off()

      pdf(paste0(output, "/dotplot_", celltype1, "_", week, "_R-NR_DEseq2_lfc0.25.pdf"), bg = "white")
      print(DotPlot(bulk_seurat,
        features = top.10.genes,
        idents = c(paste0(celltype1, "_g", week, "-NR"), paste0(celltype1, "_g", week, "-R")),
        group.by = cond_col
      ) + RotatedAxis())
      dev.off()
    } else {
      print(paste0(celltype, "_", week, " has no DE genes between R and NR"))
    }
  }
}

#### loop: R/NR vs Sham comparisons ####
celltypes <- as.vector(unique(merged_seurat[[cluster_col, drop = TRUE]]))

for (i in seq_along(celltypes)) {
  celltype <- celltypes[i]
  celltype1 <- if (grepl("_", celltype, fixed = TRUE)) gsub("_", "-", celltype) else celltype
  print(celltype1)

  seurat.obj_subset <- subset(merged_seurat, subset = .data[[cluster_col]] == celltype)

  for (j in seq_along(comparisons)) {
    week <- comparisons[j]
    print(week)

    # R vs Sham
    tryCatch({
      bulk.de <- FindMarkers(
        object = bulk_seurat,
        ident.1 = paste0(celltype1, "_g", week, "-R"),
        ident.2 = paste0(celltype1, "_Sham"),
        test.use = "DESeq2", min.cells.feature = 2,
        min.cells.group = 2,
        slot = "counts"
      )
      bulk.de <- subset(bulk.de, p_val_adj < 0.05)
      bulk.de <- subset(bulk.de, avg_log2FC >= 0.25 | avg_log2FC <= -0.25)
      bulk.de <- as.data.frame(bulk.de)
      bulk.de <- bulk.de[with(bulk.de, order(p_val_adj, -avg_log2FC)), ]
      write.table(bulk.de, file = paste0(output, "/", celltype1, "_", week, "_R-Sham_DE_genes_lfc0.25.txt"), sep = "\t")
    },
    error = function(cond) {
      message(conditionMessage(cond))
      bulk.de <<- data.frame()
    })

    if (dim(bulk.de)[1] != 0) {
      n <- min(10, dim(bulk.de)[1])
      top.10.genes <- unique(na.omit(rownames(bulk.de)[1:n]))

      pdf(paste0(output, "/violinplot_", celltype1, "_", week, "_R-Sham_DE_lfc0.25.pdf"), bg = "white")
      print(VlnPlot(seurat.obj_subset,
        features = top.10.genes,
        idents = c("Sham", paste0(week, "-R")),
        group.by = cond_col
      ))
      dev.off()

      pdf(paste0(output, "/violinplot_", celltype1, "_", week, "_R-Sham_DEseq2_lfc0.25_bulk.pdf"), bg = "white")
      print(VlnPlot(bulk_seurat,
        features = top.10.genes,
        idents = c(paste0(celltype1, "_Sham"), paste0(celltype1, "_g", week, "-R")),
        group.by = cond_col
      ))
      dev.off()
    } else {
      print(paste0(celltype, " has no DE genes between Sham and R", week))
    }

    # NR vs Sham
    tryCatch({
      bulk.de <- FindMarkers(
        object = bulk_seurat,
        ident.1 = paste0(celltype1, "_g", week, "-NR"),
        ident.2 = paste0(celltype1, "_Sham"),
        test.use = "DESeq2", min.cells.feature = 2,
        min.cells.group = 2
      )
      bulk.de <- subset(bulk.de, p_val_adj < 0.05)
      bulk.de <- subset(bulk.de, avg_log2FC >= 0.25 | avg_log2FC <= -0.25)
      bulk.de <- as.data.frame(bulk.de)
      bulk.de <- bulk.de[with(bulk.de, order(p_val_adj, -avg_log2FC)), ]
      write.table(bulk.de, file = paste0(output, "/", celltype1, "_", week, "_NR-Sham_DE_genes_lfc0.25.txt"), sep = "\t")
    },
    error = function(cond) {
      message(conditionMessage(cond))
      bulk.de <<- data.frame()
    })

    if (dim(bulk.de)[1] != 0) {
      n <- min(10, dim(bulk.de)[1])
      top.10.genes <- unique(na.omit(rownames(bulk.de)[1:n]))

      pdf(paste0(output, "/violinplot_", celltype1, "_", week, "_NR-Sham_DE_lfc0.25.pdf"), bg = "white")
      print(VlnPlot(seurat.obj_subset,
        features = top.10.genes,
        idents = c("Sham", paste0(week, "-NR")),
        group.by = cond_col
      ))
      dev.off()

      pdf(paste0(output, "/violinplot_", celltype1, "_", week, "_NR-Sham_DEseq2_lfc0.25_bulk.pdf"), bg = "white")
      print(VlnPlot(bulk_seurat,
        features = top.10.genes,
        idents = c(paste0(celltype1, "_Sham"), paste0(celltype1, "_g", week, "-NR")),
        group.by = cond_col
      ))
      dev.off()
    } else {
      print(paste0(celltype, " has no DE genes between Sham and NR", week))
    }
  }
}

capture.output(sessionInfo(), file = paste0(output, "/sessionInfo.txt"))
