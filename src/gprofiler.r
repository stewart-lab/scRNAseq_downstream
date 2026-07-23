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
library(rrvgo)
library(GO.db)

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

# Look up the Bioconductor OrgDb package matching `organism` (a gProfiler
# species code or custom GMT token) from a single shared mapping file,
# rather than a second config field -- config.json only ever needs
# `organism`, so organism and orgdb can't drift out of sync with each other.
# Add a row to this file whenever a new organism/custom GMT is introduced.
organism_map_file <- file.path(GIT_DIR, "data", "organism_orgdb_map.txt")
organism_map <- read.table(organism_map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
orgdb_match <- organism_map$orgdb[organism_map$organism == organism]
if (length(orgdb_match) == 0) {
  stop(paste0(
    "No orgdb mapping found for organism '", organism, "' in ", organism_map_file, ".\n",
    "Add a row mapping this organism/gmt-id to its Bioconductor OrgDb package name ",
    "(e.g. 'org.Hs.eg.db') before running rrvgo term reduction."
  ))
}
orgdb <- orgdb_match[1]
cat("Using orgdb:", orgdb, "for organism:", organism, "\n")

# Which GO ontologies to run rrvgo reduction over. gProfiler's GO results
# span Biological Process / Molecular Function / Cellular Component, but
# semantic similarity (and so rrvgo clustering) is only meaningful within a
# single ontology -- restricting to one (e.g. "BP") is a valid, deliberate
# choice, not a bug, but it does mean significant terms from the other
# ontologies won't be reduced (they still appear in the raw, unreduced
# output). Comma-separated in config.json, e.g. "BP" or "BP,MF,CC".
ontologies_raw <- cfg$ontologies
if (is.null(ontologies_raw) || !nzchar(ontologies_raw)) {
  ontologies <- "BP"
} else {
  ontologies <- toupper(trimws(strsplit(ontologies_raw, ",")[[1]]))
}
invalid_ontologies <- setdiff(ontologies, c("BP", "MF", "CC"))
if (length(invalid_ontologies) > 0) {
  stop(paste0(
    "Invalid value(s) in config.json's gprofiler.ontologies: ",
    paste(invalid_ontologies, collapse = ", "),
    ". Must be one or more of BP, MF, CC (comma-separated)."
  ))
}
cat("Reducing GO terms for ontologies:", paste(ontologies, collapse = ", "), "\n")

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
  # neg.log.adj.pvalue uses log10 (not natural log) so it's directly
  # comparable to the rrvgo-reduced plot's score, which is also -log10(p.adj)
  # (previously these used different log bases, making the two plots' bar
  # heights look inconsistent for the same term despite representing the
  # same p-value). The >= 1.3 threshold means p.adj <= ~0.05 either way,
  # since -log10(0.05) = 1.301.
  bar_data <- flattened_results[, c("term_name", "p.adj")]
  bar_data$neg.log.adj.pvalue <- -log10(bar_data$p.adj)
  bar_data$neg.log.adj.pvalue %>% replace_na(300)
  bar_data <- subset(bar_data, neg.log.adj.pvalue >= 1.3)

  bar_data$term_name <- gsub("^GOBP_", "", as.character(bar_data$term_name))
  bar_data$term_name <- sapply(strsplit(as.character(bar_data$term_name), "_"), function(x) {
    last_parts <- tail(x, 3)
    paste(last_parts, collapse = "_")
  })

  bar_data <- bar_data %>%
    group_by(term_name) %>%
    dplyr::slice(which.max(neg.log.adj.pvalue)) %>%
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
    labs(title = final_title, x = "Term", y = "-log10(adjusted P-value)") +
    coord_flip()

  print("Checking for any duplicate terms or unusual values:")
  print(bar_data[, c("term_name", "neg.log.adj.pvalue")])

  nd <- file.path(output_dir, paste0(output_file_name, "_barplot.pdf"))
  pdf(file = nd, height = 11, width = 8.5)
  print(p1)
  dev.off()

  ############### reduce similar GO terms via rrvgo ################
  # Same significance filter as the raw bar plot above (p.adj-based
  # neg.log.adj.pvalue >= 1.3, i.e. p.adj <= ~0.05), but keeps term_id
  # (needed by calculateSimMatrix) instead of the truncated/reformatted
  # term_name used for the raw plot.
  rrvgo_input <- flattened_results[, c("term_id", "term_name", "p.adj")]
  rrvgo_input$neg.log.adj.pvalue <- -log10(rrvgo_input$p.adj)
  rrvgo_input <- subset(rrvgo_input, neg.log.adj.pvalue >= 1.3)

  if (nrow(rrvgo_input) == 0) {
    cat("No terms pass the rrvgo significance threshold for:", basename(gene_file), "\n")
  } else {
    # Every term dropped from rrvgo reduction, and why, accumulates here and
    # gets written to _reduced_excluded.csv at the end of this block -- so
    # "why is term X in the raw CSV but missing from the reduced ones" is
    # answered by a file in the output directory, not just a console
    # message from whoever happened to run this.
    excluded_terms <- data.frame(
      term_id = character(0), term_name = character(0),
      p.adj = numeric(0), ontology = character(0), reason = character(0),
      stringsAsFactors = FALSE
    )

    # Semantic similarity (and so rrvgo clustering) is only meaningful
    # within a single GO ontology, so classify each term and restrict to
    # the ontologies requested in config.json (`ontologies`, default "BP").
    # Terms whose ontology can't be resolved at all (obsolete/merged/too-new
    # GO IDs not yet in the installed GO.db) can't be reduced regardless of
    # config -- they remain visible in the raw, unreduced output.
    rrvgo_input$ontology <- AnnotationDbi::Ontology(rrvgo_input$term_id)

    unclassified <- rrvgo_input[is.na(rrvgo_input$ontology), c("term_id", "term_name", "p.adj", "ontology")]
    if (nrow(unclassified) > 0) {
      unclassified$reason <- "no resolvable GO ontology (obsolete/merged/too-new GO ID not in installed GO.db)"
      excluded_terms <- bind_rows(excluded_terms, unclassified)
      message(
        nrow(unclassified), " term(s) have no resolvable GO ontology (obsolete/",
        "merged ID?) and are excluded from rrvgo reduction for ", basename(gene_file)
      )
    }

    not_requested <- rrvgo_input[
      !is.na(rrvgo_input$ontology) & !(rrvgo_input$ontology %in% ontologies),
      c("term_id", "term_name", "p.adj", "ontology")
    ]
    if (nrow(not_requested) > 0) {
      not_requested$reason <- paste0(
        "ontology '", not_requested$ontology, "' not in configured ontologies (",
        paste(ontologies, collapse = ", "), ")"
      )
      excluded_terms <- bind_rows(excluded_terms, not_requested)
    }

    rrvgo_input <- rrvgo_input[rrvgo_input$ontology %in% ontologies, ]

    if (nrow(rrvgo_input) == 0) {
      cat(
        "No terms in the configured ontologies (", paste(ontologies, collapse = ", "),
        ") for:", basename(gene_file), "\n"
      )
    } else {
      # run calculateSimMatrix/reduceSimMatrix separately per ontology,
      # then stack the results together with an `ontology` column.
      reduced_by_ontology <- lapply(intersect(ontologies, unique(rrvgo_input$ontology)), function(ont_i) {
        sub <- rrvgo_input[rrvgo_input$ontology == ont_i, ]
        candidate_ids <- unique(sub$term_id)

        sim_matrix <- calculateSimMatrix(
          sub$term_id,
          orgdb  = orgdb,
          ont    = ont_i,
          method = "Rel" # Relevance similarity -- generally recommended
        )
        scores <- setNames(-log10(sub$p.adj), sub$term_id)

        # reduceSimMatrix clusters terms via hclust, which needs >= 2 terms.
        # With too few significant terms in this ontology (or terms
        # calculateSimMatrix couldn't map in orgdb), there's nothing to
        # cluster -- pass what's mappable through unreduced instead.
        #
        # calculateSimMatrix's last line is `m[!out, !out]`: when exactly
        # one term survives its internal orgdb/ancestor filtering, that's a
        # 1x1 matrix subset, and R's default drop=TRUE indexing silently
        # collapses it to a bare scalar -- losing the dimnames that would
        # normally identify which term it was. A scalar NA means truly zero
        # terms mapped; a non-NA scalar means exactly one term mapped. With
        # a single candidate term there's no ambiguity about which one; with
        # several, re-check each individually (cheap: GOSemSim caches the
        # IC data) to find which one(s) actually mapped.
        has_matrix <- !is.null(dim(sim_matrix))
        if (has_matrix) {
          n_mapped_terms <- nrow(sim_matrix)
          mapped_ids <- rownames(sim_matrix)
        } else if (is.na(sim_matrix)) {
          n_mapped_terms <- 0
          mapped_ids <- character(0)
        } else {
          if (length(candidate_ids) == 1) {
            mapped_ids <- candidate_ids
          } else {
            mapped_ids <- Filter(function(id) {
              !is.na(suppressWarnings(calculateSimMatrix(
                id,
                orgdb = orgdb, ont = ont_i, method = "Rel"
              )))
            }, candidate_ids)
          }
          n_mapped_terms <- length(mapped_ids)
        }

        if (n_mapped_terms < 2) {
          message(
            "Only ", n_mapped_terms, " ", ont_i, " term(s) available for rrvgo ",
            "reduction -- skipping clustering for ", basename(gene_file)
          )
          mapped_ids <- as.character(mapped_ids)
          mapped_names <- sub$term_name[match(mapped_ids, sub$term_id)]
          reduced_i <- data.frame(
            go = mapped_ids,
            term = mapped_names,
            parentTerm = mapped_names,
            score = scores[mapped_ids],
            stringsAsFactors = FALSE
          )
        } else {
          reduced_i <- reduceSimMatrix(
            sim_matrix,
            scores,
            threshold = 0.7, # similarity cutoff; 0.7 = moderately aggressive collapsing
            orgdb     = orgdb
          )
        }
        # rep(), not a bare scalar: reduced_i can have 0 rows (e.g. the
        # ontology's only term(s) weren't found in orgdb at all), and
        # data.frame's `$<-` refuses to recycle a length-1 value onto a
        # 0-row data.frame.
        reduced_i$ontology <- rep(ont_i, nrow(reduced_i))

        # candidate_ids minus mapped_ids: terms with a valid, configured
        # ontology that calculateSimMatrix still couldn't map in orgdb
        # (rrvgo's own "Removed N terms not found in orgdb" warning, made
        # per-term and persisted instead of just a console warning).
        not_found_ids <- setdiff(candidate_ids, mapped_ids)
        excluded_i <- if (length(not_found_ids) > 0) {
          data.frame(
            term_id = not_found_ids,
            term_name = sub$term_name[match(not_found_ids, sub$term_id)],
            p.adj = sub$p.adj[match(not_found_ids, sub$term_id)],
            ontology = ont_i,
            reason = paste0("not found in orgdb (", orgdb, ") for ontology ", ont_i),
            stringsAsFactors = FALSE
          )
        } else {
          data.frame(
            term_id = character(0), term_name = character(0),
            p.adj = numeric(0), ontology = character(0), reason = character(0),
            stringsAsFactors = FALSE
          )
        }

        list(reduced = reduced_i, excluded = excluded_i)
      })
      # bind_rows, not rbind: the <2-mapped-terms passthrough branch above
      # produces fewer columns (go/term/parentTerm/score/ontology) than a
      # real reduceSimMatrix() result (which also has cluster/parent/size/
      # termUniqueness/etc), so different ontologies in this list can have
      # different column sets. bind_rows aligns by name and fills the rest
      # with NA; base rbind would error on the column mismatch.
      reduced <- bind_rows(lapply(reduced_by_ontology, `[[`, "reduced"))
      reduced <- as.data.frame(reduced)
      rownames(reduced) <- NULL
      excluded_terms <- bind_rows(excluded_terms, lapply(reduced_by_ontology, `[[`, "excluded"))

      print(head(reduced[, c("ontology", "go", "term", "parentTerm", "score")]))

      # full per-term detail (every significant term, its ontology, its
      # cluster's parent, and its own individual, untouched score) -- this
      # is also the file to consult if you want to see exactly which terms
      # were grouped under a given parent, row by row.
      write.csv(reduced, file.path(output_dir,
                paste0(output_file_name, "_reduced.csv")), row.names = FALSE)

      # one row per (ontology, parent cluster), listing its member (child)
      # terms -- a more direct answer to "what got grouped under this
      # parent" than scanning the full _reduced.csv for matching rows.
      grouping_summary <- reduced %>%
        group_by(ontology, parentTerm) %>%
        summarise(
          score = max(score),
          n_terms = n(),
          member_terms = paste(term, collapse = "; "),
          .groups = "drop"
        ) %>%
        arrange(ontology, desc(score))
      write.csv(grouping_summary, file.path(output_dir,
                paste0(output_file_name, "_reduced_grouping.csv")), row.names = FALSE)

      # hierarchical JSON view of the same data as _reduced.csv: ontology ->
      # parent cluster -> member terms, nested so a JSON viewer/editor's
      # fold/outline view shows parent-child structure directly, without
      # needing to filter or sort a spreadsheet to see what's grouped
      # under what (the treemap shows this visually but gets cramped with
      # many clusters; this is the same information as plain, browsable text).
      hierarchy <- lapply(split(reduced, reduced$ontology), function(ont_df) {
        clusters <- lapply(split(ont_df, ont_df$parentTerm), function(cluster_df) {
          cluster_df <- cluster_df[order(-cluster_df$score), ]
          is_parent_row <- cluster_df$term == cluster_df$parentTerm
          parent_go <- if (any(is_parent_row)) cluster_df$go[which(is_parent_row)[1]] else NA
          list(
            parentTerm = cluster_df$parentTerm[1],
            parent_go = parent_go,
            score = max(cluster_df$score),
            n_terms = nrow(cluster_df),
            terms = lapply(seq_len(nrow(cluster_df)), function(i) {
              list(
                go = cluster_df$go[i],
                term = cluster_df$term[i],
                score = cluster_df$score[i],
                is_parent = isTRUE(cluster_df$term[i] == cluster_df$parentTerm[i])
              )
            })
          )
        })
        # most significant cluster first, within each ontology
        clusters[order(-vapply(clusters, function(c) c$score, numeric(1)))]
      })
      write(
        jsonlite::toJSON(hierarchy, pretty = TRUE, auto_unbox = TRUE, na = "null"),
        file.path(output_dir, paste0(output_file_name, "_reduced_hierarchy.json"))
      )

      # one bar per (ontology, parent cluster) -- a term that got grouped
      # under a more significant parent no longer gets its own bar, since
      # which.max(score) within each parentTerm always selects the
      # parent's own row (parentTerm clusters are defined by their
      # highest-scoring member in the first place). Faceted by ontology
      # since scores are comparable across ontologies but clusters aren't.
      reduced_plot_data <- reduced %>%
        group_by(ontology, parentTerm) %>%
        dplyr::slice(which.max(score)) %>%
        ungroup()

      p2 <- ggplot(reduced_plot_data, aes(x = reorder(parentTerm, score), y = score)) +
        geom_col(aes(fill = score), position = "identity") +
        theme_minimal() +
        scale_fill_gradientn(colours = colorRampPalette(c("blue", "red"))(100)) +
        labs(title = final_title, x = "Term (parent)", y = "-log10(adjusted P-value)") +
        coord_flip() +
        facet_wrap(~ontology, scales = "free_y", ncol = 1)

      nd2 <- file.path(output_dir, paste0(output_file_name, "_reduced_barplot.pdf"))
      pdf(file = nd2, height = 11, width = 8.5)
      print(p2)
      dev.off()

      # treemap: rrvgo's own visualization for how terms were grouped --
      # each parent cluster is a labeled region, subdivided into its
      # member (child) terms, sized by score. A visual complement to
      # _reduced_grouping.csv above. One per ontology, since treemapPlot
      # doesn't facet.
      for (ont_i in unique(reduced$ontology)) {
        nd3 <- file.path(output_dir, paste0(output_file_name, "_reduced_treemap_", ont_i, ".pdf"))
        pdf(file = nd3, height = 8.5, width = 11)
        treemapPlot(reduced[reduced$ontology == ont_i, ], title = paste0(final_title, " (", ont_i, ")"))
        dev.off()
      }

      cat("Reduced-term output written to:", output_file_name, "_reduced\n")
    }

    # answers "why is term X in the raw CSV but missing from the reduced
    # ones" without needing to have captured this run's console output.
    if (nrow(excluded_terms) > 0) {
      write.csv(excluded_terms, file.path(output_dir,
                paste0(output_file_name, "_reduced_excluded.csv")), row.names = FALSE)
      cat(nrow(excluded_terms), "term(s) excluded from rrvgo reduction; see",
          paste0(output_file_name, "_reduced_excluded.csv\n"))
    }
  }

  cat("Output written to:", output_file_name, "\n")
}
