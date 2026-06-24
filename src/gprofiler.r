# gets Go enrichment for list of De genes from single cell
# configured via the gprofiler section of config.json

suppressPackageStartupMessages(library("gprofiler2"))
library(jsonlite)
library(dplyr)
library("ggplot2")
library(reshape2)
library(devEMF)
library(ggrepel)
library(tidyr)
library(tibble)

### load config ###
GIT_DIR <- getwd()
config <- jsonlite::fromJSON(file.path(GIT_DIR, "config.json"))
cfg <- config$gprofiler

DATA_DIR    <- cfg$DATA_DIR
organism    <- cfg$organism
mthreshold  <- cfg$mthreshold
padj        <- cfg$padj
lfc         <- cfg$lfc
file_pattern <- cfg$file_pattern
output_name <- cfg$output_name
plot_title  <- cfg$title
lower       <- as.logical(cfg$lower)

### set output directory ###
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(GIT_DIR, "shared_volume", paste0("output_gprofiler_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
file.copy(file.path(GIT_DIR, "config.json"), file.path(output_dir, "config.json"))

cat("Output directory:", output_dir, "\n")

### find input files ###
gene_files <- list.files(DATA_DIR, pattern = file_pattern, full.names = TRUE)

if (length(gene_files) == 0) {
  stop(paste("No files matching", file_pattern, "found in:", DATA_DIR))
}

cat("Found", length(gene_files), "file(s) to process:\n")
cat(paste(" ", gene_files, collapse = "\n"), "\n")

for (gene_file in gene_files) {
  cat("\n--- Processing:", basename(gene_file), "---\n")

  output_file_name <- tools::file_path_sans_ext(basename(gene_file))
  output_file_name <- paste0(output_file_name, output_name)

  data <- read.table(gene_file, header = TRUE, stringsAsFactors = FALSE)
  data <- tibble::rownames_to_column(data, "gene")
  print(colnames(data))

  if (isTRUE(lower)) {
    data$gene <- tolower(data$gene)
  }

  # filter by lfc and padj
  if (is.null(lfc) && is.null(padj)) {
    print("using all genes")
  } else if (!is.null(lfc)) {
    if (sign(lfc) == 1) {
      data <- subset(data, avg_log2FC >= lfc & p_val_adj < padj)
    } else if (sign(lfc) == -1) {
      data <- subset(data, avg_log2FC <= lfc & p_val_adj < padj)
    }
  } else {
    data <- subset(data, p_val_adj < padj)
  }

  if (dim(data)[1] == 0) {
    cat("No DE genes to analyze for:", basename(gene_file), "\n")
    next
  }

  ############### enrichment analysis via gProfiler ################
  gostres <- gost(
    data$gene,
    organism = organism,
    measure_underrepresentation = FALSE,
    significant = TRUE
  )

  if (is.null(gostres) || is.null(gostres$result)) {
    cat("No significant GO terms found for:", basename(gene_file), "\n")
    next
  }

  # Flatten list columns in the results
  flattened_results <- gostres$result
  list_columns <- sapply(flattened_results, is.list)
  flattened_results[list_columns] <- lapply(flattened_results[list_columns], function(column) {
    sapply(column, function(entry) {
      if (is.null(entry)) return(NA)
      paste(entry, collapse = ", ")
    })
  })

  flattened_results <- flattened_results[order(flattened_results$p_value, decreasing = FALSE), ]
  flattened_results$p.adj <- p.adjust(flattened_results$p_value, method = "fdr",
                                      n = length(flattened_results$p_value))
  flattened_results <- flattened_results[flattened_results$term_size < mthreshold, ]

  print(flattened_results)

  if (nrow(flattened_results) == 0) {
    cat("No significant terms found for:", basename(gene_file), "\n")
    next
  }

  # if any term_name contains 'https://', swap term_id and term_name
  first_term_name <- flattened_results$term_name[1]
  if (grepl("https://", first_term_name)) {
    temp <- flattened_results$term_name
    flattened_results$term_name <- flattened_results$term_id
    flattened_results$term_id <- temp
  }

  write.csv(flattened_results, file.path(output_dir,
            paste0(output_file_name, ".csv")), row.names = FALSE)

  # make bar plot
  bar_data <- flattened_results[, c("term_name", "p.adj")]
  bar_data$neg.log.adj.pvalue <- -log(bar_data$p.adj)
  bar_data$neg.log.adj.pvalue %>% replace_na(300)
  bar_data <- subset(bar_data, neg.log.adj.pvalue >= 1.3)

  bar_data$term_name <- gsub("^GOBP_", "", as.character(bar_data$term_name))
  bar_data$term_name <- sapply(strsplit(as.character(bar_data$term_name), "_"), function(x) {
    last_parts <- tail(x, 3)
    paste(last_parts, collapse = "_")
  })

  bar_data <- bar_data %>%
    group_by(term_name) %>%
    slice(which.max(neg.log.adj.pvalue)) %>%
    ungroup()

  # use config title if not default, otherwise fall back to filename
  final_title <- if (!is.null(plot_title) && plot_title != "GO enrichment") {
    plot_title
  } else {
    tools::file_path_sans_ext(basename(gene_file))
  }

  p1 <- ggplot(bar_data, aes(x = reorder(term_name, neg.log.adj.pvalue), y = neg.log.adj.pvalue)) +
    geom_col(aes(fill = neg.log.adj.pvalue), position = "identity") +
    theme_minimal() +
    scale_fill_gradientn(colours = colorRampPalette(c("blue", "red"))(100)) +
    labs(title = final_title, x = "Term", y = "Neg Log adjusted P-value") +
    coord_flip()

  print("Checking for any duplicate terms or unusual values:")
  print(bar_data[, c("term_name", "neg.log.adj.pvalue")])

  nd <- file.path(output_dir, paste0(output_file_name, "_barplot.pdf"))
  pdf(file = nd, height = 11, width = 8.5)
  print(p1)
  dev.off()

  cat("Output written to:", output_file_name, "\n")
}
