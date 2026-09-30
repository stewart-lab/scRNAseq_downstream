library(Seurat)
library(ggplot2)

base_dir <- "/w5home/bmoore/scRNAseq/LiFangChu/10x_somite_data/output_20260916_181955/pTLS_D4_Mgel_F20CH5"
sub_dir <- file.path(base_dir, "output_subset-recluster_20260928_144119")

full_meta_file <- file.path(base_dir, "manual_annot_metadata_pTLS_D4_Mgel_F20CH5.txt")
sub_meta_file <- file.path(sub_dir, "manual_annot_metadata_6_res0.9.txt")
obj_file <- file.path(base_dir, "final_processed_obj.rds")

parent_cluster <- "6"
sub_col <- "seurat_clusters_res0.9" # subcluster labels in the subset metadata
new_col <- "merged_clusters"

# ---- Load ----
# Header has one fewer field than the rows, so the first column (barcodes) becomes row names
full_meta <- read.delim(full_meta_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
sub_meta <- read.delim(sub_meta_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
seu <- readRDS(obj_file)

# ---- Sanity checks ----
cluster6_cells <- rownames(full_meta)[as.character(full_meta$seurat_clusters) == parent_cluster]
stopifnot(
    all(rownames(sub_meta) %in% rownames(full_meta)), # subset cells exist in full metadata
    setequal(rownames(sub_meta), cluster6_cells), # subset is exactly cluster 6
    setequal(rownames(full_meta), colnames(seu)) # metadata matches the Seurat object
)
print(table(sub_meta[[sub_col]]))

# ---- Build merged cluster labels ----
# Every cluster keeps its label; cluster 6 cells become "6_1", "6_2", ...
merged <- as.character(full_meta$seurat_clusters)
names(merged) <- rownames(full_meta)
merged[rownames(sub_meta)] <- paste0(parent_cluster, "_", sub_meta[[sub_col]])

# Order levels numerically so 6_1/6_2 sit between 5 and 7
lvls <- unique(merged)
lvls <- lvls[order(as.numeric(sub("_.*", "", lvls)), as.numeric(sub("^[^_]*_?", "", lvls)), na.last = FALSE)]
merged <- factor(merged, levels = lvls)

print(table(merged))
stopifnot(sum(table(merged)) == nrow(full_meta))

# ---- Add to Seurat object (only the new column; existing metadata is untouched) ----
seu <- AddMetaData(seu, metadata = merged[colnames(seu)], col.name = new_col)
Idents(seu) <- new_col

# Cross-check against the object's own clusters (should only differ for cluster 6)
print(table(original = seu$seurat_clusters, merged = seu[[new_col, drop = TRUE]]))

# ---- Visualize ----
p <- DimPlot(seu, reduction = "umap", group.by = new_col, label = TRUE, repel = TRUE) +
    ggtitle("Clusters with cluster 6 subclustered (res 0.9)")
ggsave(file.path(sub_dir, "umap_merged_clusters_6_res0.9.pdf"), p, width = 8, height = 6)

# ---- Save (new files; does not overwrite final_processed_obj.rds) ----
write.table(seu@meta.data, file.path(sub_dir, "manual_annot_metadata_merged_6_res0.9.txt"),
    sep = "\t", quote = F, row.names = T
)
saveRDS(seu, file.path(sub_dir, "final_processed_obj_merged_6_res0.9.rds"))
