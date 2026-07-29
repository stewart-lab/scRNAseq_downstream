# CASSIA cell-type annotation for a (re)clustered Seurat object, with an
# optional cross-species marker remap via a biomart orthology table.
# Configured via the cassia section of config.json.
# Ported from cassiaRMSTestGammS2_clus1.Rmd (GAMM S2 cluster-1 Docker test)
# so it can run as a regular method through run_downstream_toolkit.sh.

# Bind reticulate to the cassia_env conda Python BEFORE library(CASSIA).
# CASSIA's .onLoad otherwise initializes a uv-managed Python via py_require(),
# and setup_cassia_env(method = "conda") then fails with "another version
# of Python has already been initialized".
cassia_python <- Sys.getenv(
  "RETICULATE_PYTHON",
  unset = "/opt/conda/envs/cassia_env/bin/python"
)
if (!file.exists(cassia_python)) {
  stop("CASSIA Python not found at: ", cassia_python)
}
Sys.setenv(RETICULATE_PYTHON = cassia_python)
options(
  CASSIA.env_name = "cassia_env",
  CASSIA.env_method = "conda",
  CASSIA.conda_env = "cassia_env"
)

library(reticulate)
library(CASSIA)
library(Seurat)
library(dplyr)
library(stringr)
library(patchwork)
library(ggplot2)
library(jsonlite)

# Record R session info (packages + versions) for run reproducibility --
# appended to the run's shared provenance file when invoked via
# run_downstream_toolkit.sh (PROVENANCE_FILE env var), else written
# standalone to the current directory.
.provenance_file <- Sys.getenv("PROVENANCE_FILE", unset = paste0("./sessionInfo_cassia.R_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
cat(paste0("\n--- R sessionInfo (cassia.R) ---\n", paste(capture.output(sessionInfo()), collapse = "\n"), "\n"), file = .provenance_file, append = TRUE)

# Deadlock workaround: any print()/sys.stdout.write() call made from a
# background Python thread gets routed by reticulate through a queue that
# only the main R thread services. runCASSIA_pipeline() (and
# merge_annotations_all) blocks the main R thread synchronously in
# as_completed() while those same worker threads run -- so the callback
# queue never gets serviced, and whichever thread tries to print first
# (BatchProgressTracker's periodic render, or the plain print() in
# merge_annotations's exception handler, or anything else) hangs forever,
# often while holding a lock that then blocks every other worker behind it.
#
# First attempt (replacing sys.stdout directly) did NOT work: reticulate
# wraps every R->Python call boundary in RemapOutputStreams, which
# reassigns sys.stdout to a *fresh* rpytools.output.OutputRemap instance on
# each such call -- silently clobbering any sys.stdout object we set
# beforehand, including the whole runCASSIA_pipeline() call (confirmed via
# py-spy: the deadlock recurred at the exact same _render -> write call
# site after the sys.stdout patch was in place).
#
# Real fix: OutputRemap.write() is what actually calls back into R
# (self.handler(message)); patch that *class method* instead of the
# sys.stdout object. Since RemapOutputStreams creates new OutputRemap
# instances on every call, but they're all instances of this same class,
# patching the method survives being re-wrapped indefinitely -- the same
# principle as the BatchProgressTracker._render patch below.
reticulate::py_run_string("
import rpytools.output as _ro
import os

def _direct_write(self, message):
    try:
        os.write(1, message.encode('utf-8', errors='replace'))
    except Exception:
        pass
    return len(message)

_ro.OutputRemap.write = _direct_write
")

# Do not call setup_cassia_env(method = "conda") here -- it would try to
# re-init Python after reticulate is already attached to cassia_env.
print(py_config())
py_cassia <- get("py_cassia", envir = asNamespace("CASSIA"))
if (is.null(py_cassia)) {
  stop("CASSIA Python module failed to load (py_cassia is NULL)")
}
message("CASSIA Python ready: ", py_config()$python)

### load config ###
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(GIT_DIR, "config.json"))
cfg <- config$cassia

docker <- config$docker
if (docker == "TRUE" || docker == "true" || docker == "T" || docker == "t") {
  DATA_DIR <- "./data/input_data/"
} else {
  DATA_DIR <- cfg$DATA_DIR
}
SEURAT_OBJ <- cfg$SEURAT_OBJ
# Optional: a precomputed FindAllMarkers-style CSV (columns include at least
# gene/cluster/avg_log2FC/p_val_adj), path relative to DATA_DIR. When set,
# the Seurat object load, recluster, and FindAllMarkers steps below are all
# skipped -- there's then no Seurat object to recluster or attach CASSIA's
# annotations back onto, only the marker-based CASSIA run itself.
MARKER_FILE <- cfg$MARKER_FILE
USE_MARKER_FILE <- !is.null(MARKER_FILE) && nzchar(MARKER_FILE)
# Optional: a tab-separated cross-species orthology table (pig.gene.name/
# human.gene.name columns) to remap marker gene symbols before annotation.
# Leave blank in config.json to skip remapping and annotate with the marker
# gene symbols as-is (e.g. when markers are already in the target species).
BIOMART_FILE <- cfg$BIOMART_FILE
USE_BIOMART <- !is.null(BIOMART_FILE) && nzchar(BIOMART_FILE)
TISSUE <- cfg$tissue
SPECIES <- cfg$species
CLUSTER_COL <- cfg$CLUSTER_COL
CASSIA_OUT_NAME <- cfg$CASSIA_OUT_NAME
DIM.RED <- cfg$DIM.RED

DO_RECLUSTER <- isTRUE(cfg$recluster$DO_RECLUSTER)
RECLUSTER_RESOLUTION <- cfg$recluster$RESOLUTION
RECLUSTER_DIMS <- seq_len(cfg$recluster$DIMS)
RECLUSTER_ALGORITHM <- cfg$recluster$ALGORITHM

DO_FIND_MARKERS <- isTRUE(cfg$DO_FIND_MARKERS)
DO_RUN_CASSIA <- isTRUE(cfg$DO_RUN_CASSIA)

### set output directory ###
timestamp <- Sys.getenv("RUN_TIMESTAMP", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
OUTPUT_DIR <- paste0("./shared_volume/output_cassia_", timestamp, "/")
dir.create(OUTPUT_DIR, recursive = TRUE, mode = "0777", showWarnings = FALSE)
file.copy(file.path(GIT_DIR, "config.json"), file.path(OUTPUT_DIR, "config.json"))
if (file.exists(.provenance_file)) {
  file.copy(.provenance_file, file.path(OUTPUT_DIR, basename(.provenance_file)), overwrite = TRUE)
  file.remove(.provenance_file)
}
message("Output directory: ", OUTPUT_DIR)

### api key ###
# Prefer cfg$openAI_key (config.json), falling back to an OPENAI_API_KEY
# already set in the environment (e.g. via docker --env-file/-e) so either
# workflow keeps working. persist = FALSE avoids writing the key into the
# container filesystem.
api_key <- cfg$openAI_key
if (is.null(api_key) || !nzchar(api_key)) {
  api_key <- Sys.getenv("OPENAI_API_KEY")
}
if (!nzchar(api_key)) {
  stop("OpenAI API key not set. Set cassia.openAI_key in config.json, or set OPENAI_API_KEY (e.g. via docker --env-file/-e).")
}
setLLMApiKey(api_key, provider = "openai", persist = FALSE)

markers_merged_file <- file.path(OUTPUT_DIR, "markers_merge_unique.csv")

if (USE_MARKER_FILE) {
  ### load a precomputed marker file instead of a Seurat object ###
  message("cassia.MARKER_FILE set (", MARKER_FILE, "); skipping Seurat object load, recluster, and FindAllMarkers")
  seurat.obj <- NULL
  marker_file_path <- file.path(DATA_DIR, MARKER_FILE)
  if (!DO_FIND_MARKERS && file.exists(markers_merged_file)) {
    message("Skipping marker file load; will load ", markers_merged_file, " later")
    all.markers <- NULL
  } else {
    if (!file.exists(marker_file_path)) {
      stop(paste("MARKER_FILE not found at:", marker_file_path))
    }
    all.markers <- read.csv(marker_file_path, stringsAsFactors = FALSE)
    write.csv(all.markers, file = file.path(OUTPUT_DIR, "markers_celltypes_all.csv"), row.names = FALSE)
  }
} else {
  ### read in seurat object ###
  seurat.obj <- readRDS(file.path(DATA_DIR, SEURAT_OBJ))
  print(seurat.obj)
  print(colnames(seurat.obj@meta.data))
  if (CLUSTER_COL %in% colnames(seurat.obj@meta.data)) {
    print(table(seurat.obj[[CLUSTER_COL]], useNA = "ifany"))
  }

  ### recluster (creates multiple subclusters for CASSIA/FindAllMarkers) ###
  if (DO_RECLUSTER) {
    message("Reclustering at resolution ", RECLUSTER_RESOLUTION)

    if (!"pca" %in% names(seurat.obj@reductions)) {
      message("No PCA found; running RunPCA on variable features")
      if (length(VariableFeatures(seurat.obj)) == 0) {
        seurat.obj <- FindVariableFeatures(seurat.obj)
      }
      seurat.obj <- RunPCA(seurat.obj, features = VariableFeatures(seurat.obj), verbose = FALSE)
    }

    seurat.obj <- FindNeighbors(
      seurat.obj,
      dims = RECLUSTER_DIMS,
      reduction = "pca",
      verbose = FALSE
    )

    CLUSTER_COL <- paste0("seurat_clusters_res", RECLUSTER_RESOLUTION)
    seurat.obj <- FindClusters(
      seurat.obj,
      resolution = RECLUSTER_RESOLUTION,
      algorithm = RECLUSTER_ALGORITHM,
      cluster.name = CLUSTER_COL,
      verbose = FALSE
    )

    # Optional: refresh UMAP on this subset's PCA for clearer plots
    if (!DIM.RED %in% names(seurat.obj@reductions) || DO_RECLUSTER) {
      seurat.obj <- RunUMAP(seurat.obj, dims = RECLUSTER_DIMS, reduction = "pca", verbose = FALSE)
      DIM.RED <- "umap"
    }

    message("New cluster column: ", CLUSTER_COL)
    print(table(seurat.obj[[CLUSTER_COL]]))

    pdf(file.path(OUTPUT_DIR, paste0("umap_recluster_", CLUSTER_COL, ".pdf")), width = 8, height = 6)
    print(DimPlot(seurat.obj, reduction = DIM.RED, group.by = CLUSTER_COL, label = TRUE, pt.size = 0.5))
    dev.off()

    saveRDS(seurat.obj, file = file.path(OUTPUT_DIR, "seurat.obj_reclustered.rds"))
    message("Saved reclustered object to ", file.path(OUTPUT_DIR, "seurat.obj_reclustered.rds"))
  } else {
    message("Skipping recluster; using CLUSTER_COL = ", CLUSTER_COL)
  }

  ### set idents to (re)cluster column ###
  if (!(CLUSTER_COL %in% colnames(seurat.obj@meta.data))) {
    stop(paste("CLUSTER_COL not found in metadata:", CLUSTER_COL))
  }
  Idents(object = seurat.obj) <- CLUSTER_COL
  print(table(Idents(seurat.obj)))
  n_clusters <- length(unique(as.character(Idents(seurat.obj))))
  if (n_clusters < 2) {
    stop(paste(
      "Only", n_clusters, "cluster(s) after setup. FindAllMarkers/CASSIA need multiple clusters.",
      "Try a higher cassia.recluster.RESOLUTION (e.g. 0.8 or 1.0) or check cassia.recluster.DO_RECLUSTER."
    ))
  }

  ### find markers for each cluster ###
  if (!DO_FIND_MARKERS && file.exists(markers_merged_file)) {
    message("Skipping FindAllMarkers; will load ", markers_merged_file, " later")
    all.markers <- NULL
  } else {
    all.markers <- FindAllMarkers(object = seurat.obj)
    all.markers <- as.data.frame(all.markers)
    write.csv(all.markers, file = file.path(OUTPUT_DIR, "markers_celltypes_all.csv"), row.names = FALSE)
  }
}

### remap marker gene symbols to the target species via biomart orthology
### table (optional -- skipped if BIOMART_FILE is blank), then dedupe ###
if (!DO_FIND_MARKERS && file.exists(markers_merged_file)) {
  markers_merge_unique <- read.csv(markers_merged_file, stringsAsFactors = FALSE)
  message("Loaded markers_merge_unique: ", nrow(markers_merge_unique), " rows")
} else if (USE_BIOMART) {
  if (!file.exists(BIOMART_FILE)) {
    stop(paste(
      "Biomart file not found at:", BIOMART_FILE,
      "\nFor docker, place it under ./data/ (mounted to /data) and point",
      "cassia.BIOMART_FILE at it, e.g. \"./data/biomart/<file>.txt\"."
    ))
  }
  gene_symbols_df <- read.csv(file = BIOMART_FILE, header = TRUE, sep = "\t")
  print(head(gene_symbols_df))
  gene_symbols_df$pig.gene.name <- toupper(gene_symbols_df$pig.gene.name)

  markers_merge <- merge(all.markers, gene_symbols_df,
    by.x = "gene", by.y = "pig.gene.name",
    all = FALSE
  )
  print(head(markers_merge))

  markers_merge$pig_gene <- markers_merge$gene
  markers_merge$gene <- markers_merge$human.gene.name
  markers_merge <- markers_merge[, 1:7]
  markers_merge_unique <- distinct(markers_merge)
  write.csv(markers_merge_unique, file = markers_merged_file, row.names = FALSE)
} else {
  message("cassia.BIOMART_FILE not set; annotating with marker gene symbols as-is (no cross-species remap)")
  markers_merge_unique <- distinct(all.markers)
  write.csv(markers_merge_unique, file = markers_merged_file, row.names = FALSE)
}

### run CASSIA ###
if (DO_RUN_CASSIA) {
  runCASSIA_pipeline(
    output_file_name = CASSIA_OUT_NAME,
    output_dir = OUTPUT_DIR,
    tissue = TISSUE,
    species = SPECIES,
    marker = markers_merge_unique,
    max_workers = 4,
    # Re-enabled: the merge step was defaulting to overall_provider =
    # "openrouter" (we never overrode merge_provider/merge_model), and we
    # only supply an OpenAI key -- so every merge LLM call failed auth,
    # which triggered the print()-from-worker-thread deadlock in
    # merge_annotations's exception handler. Pointing merge at the same
    # provider/credentials as every other stage avoids that failure path;
    # the sys.stdout patch above covers any print() call that still happens
    # (e.g. on a genuine transient error) so it no longer deadlocks either
    # way.
    #
    # merge_model is gpt-4.1, NOT gpt-5.1 like the other stages: CASSIA's
    # merge_annotations call sends a `reasoning.effort` request parameter
    # for reasoning-capable models (gpt-5.1 included) whose value it
    # computes incorrectly for this specific model -- OpenAI's API rejects
    # it with "Invalid value: '<truncated text>'. Supported values are:
    # 'none', 'minimal', ..." (a CASSIA-side bug building that request, not
    # anything on our end). gpt-4.1 isn't a reasoning model, so it never
    # takes that code path. If CASSIA fixes this for gpt-5.1 upstream,
    # merge_model can go back to matching the other stages.
    do_merge_annotations = TRUE,
    merge_provider = "openai",
    merge_model = "gpt-4.1",
    annotation_model = "gpt-5.1",
    annotation_provider = "openai",
    score_model = "gpt-5.1",
    score_provider = "openai",
    annotationboost_model = "gpt-5.1",
    annotationboost_provider = "openai"
  )
} else {
  message("Skipping runCASSIA_pipeline; will reuse existing FINAL_RESULTS.csv from OUTPUT_DIR")
}

if (USE_MARKER_FILE) {
  message(
    "cassia.MARKER_FILE was set; no Seurat object to annotate or plot. See ",
    "CASSIA_Pipeline_*/03_csv_files/*_FINAL_RESULTS.csv under ", OUTPUT_DIR, " for results."
  )
} else {
  ### add cassia annotations back into seurat object ###
  # Find the most recent CASSIA_Pipeline_*/03_csv_files/*_FINAL_RESULTS.csv under
  # OUTPUT_DIR rather than hardcoding a timestamped path.
  final_results_candidates <- Sys.glob(file.path(
    OUTPUT_DIR, "CASSIA_Pipeline_*", "03_csv_files", "*_FINAL_RESULTS.csv"
  ))
  cassia_results_path <- if (length(final_results_candidates) > 0) {
    final_results_candidates[order(file.info(final_results_candidates)$mtime, decreasing = TRUE)][1]
  } else {
    file.path(OUTPUT_DIR, "NO_FINAL_RESULTS_FOUND.csv")
  }
  if (!file.exists(cassia_results_path)) {
    warning(paste(
      "No FINAL_RESULTS.csv found under", OUTPUT_DIR,
      "-- run runCASSIA_pipeline() first (set cassia.DO_RUN_CASSIA to true)."
    ))
  } else {
    message("Using CASSIA results: ", cassia_results_path)
    # columns_to_include = 2: always add both the raw per-cluster predictions
    # (general/sub/mixed_celltype, score) AND the merged-grouping columns when
    # present. Default (1) silently drops the raw prediction columns entirely
    # whenever merged groupings exist, which caught us off guard once
    # do_merge_annotations was enabled -- the plotting step below was written
    # against the raw columns and found nothing to plot.
    seurat.obj <- add_cassia_to_seurat(
      seurat_obj = seurat.obj,
      cassia_results_path = cassia_results_path,
      cluster_col = CLUSTER_COL,
      cassia_cluster_col = "Cluster ID",
      prefix = "Cassia_",
      columns_to_include = 2
    )
  }

  saveRDS(seurat.obj, file = file.path(OUTPUT_DIR, paste0(CASSIA_OUT_NAME, "_cassia.rds")))

  ### plot cassia annotations, if present ###
  # Curated columns worth visualizing. Prefer the merged-grouping columns
  # (Cassia_merged_grouping_1/2/3) when do_merge_annotations produced them --
  # they're CASSIA's own consolidated broad/detailed/very-detailed labels,
  # shorter and cleaner than the raw per-cluster predictions, which is the
  # whole point of running the merge step. Fall back to the raw columns
  # otherwise (e.g. if merging was disabled or failed). Deliberately excludes
  # from the raw-column fallback:
  #  - Cassia_score: a numeric confidence score, not a cell-type label.
  #    Grouping by it is actively misleading -- unrelated clusters that
  #    happen to share a score (e.g. two different clusters both scoring 88)
  #    get painted the same color.
  #  - mixed_celltype / sub_celltype_all / sub_celltype_1-3: CASSIA's raw
  #    free-text variants, either redundant with general/sub_celltype or too
  #    verbose for a legend.
  merged_cols <- intersect(
    c("Cassia_merged_grouping_1", "Cassia_merged_grouping_2", "Cassia_merged_grouping_3"),
    colnames(seurat.obj@meta.data)
  )
  raw_cols <- intersect(
    c("Cassia_general_celltype", "Cassia_sub_celltype", "Cassia_combined_celltype"),
    colnames(seurat.obj@meta.data)
  )
  plot_cols <- if (length(merged_cols) > 0) merged_cols else raw_cols

  if (length(plot_cols) >= 1) {
    # CASSIA's labels can run to a full sentence. Wrapping them into a
    # separate display-only column (rather than plotting the raw column)
    # keeps the legend from ballooning and pushing the actual UMAP panels
    # out of the page.
    for (col in plot_cols) {
      seurat.obj@meta.data[[paste0(col, "_wrapped")]] <-
        str_wrap(seurat.obj@meta.data[[col]], width = 40)
    }
    wrapped_cols <- paste0(plot_cols, "_wrapped")

    plots <- lapply(wrapped_cols, function(col) {
      # label = FALSE: CASSIA's cell-type names are too long to render as
      # in-plot cluster labels without overlapping each other. The legend
      # (already wrapped/sized above) conveys the full names instead.
      DimPlot(seurat.obj,
        reduction = DIM.RED, label = FALSE,
        pt.size = 0.5, group.by = col
      ) +
        theme(legend.text = element_text(size = 7)) +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
    })
    combined_plot <- wrap_plots(plots, ncol = length(plots))
    pdf(
      file = file.path(OUTPUT_DIR, paste0("Cassia_annot_", CASSIA_OUT_NAME, "_umap.pdf")),
      width = 7 * length(plots), height = 7
    )
    print(combined_plot)
    dev.off()
  } else {
    message("No curated Cassia_* cell-type columns found; skip plotting until add_cassia_to_seurat succeeds.")
  }
}
