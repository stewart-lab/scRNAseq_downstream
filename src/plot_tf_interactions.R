# Load necessary libraries
suppressPackageStartupMessages({
    library(pheatmap)
    library(dplyr)
})

# ==========================================
# USER CONFIGURATION SECTION
# ==========================================

# 1. Path to your interaction matrices
# INPUT_TF_MATRIX: Rows = TFs, Columns = Ligands/Receptors
INPUT_TF_MATRIX_FILE <- "/w5home/bmoore/Pierre_sc_zebrafish/nichenet/Nichenet_ligand-tf_matrix.txt"

# INPUT_TARGET_MATRIX: Rows = Targets, Columns = Ligands/Receptors (Assumed)
INPUT_TARGET_MATRIX_FILE <- "/w5home/bmoore/Pierre_sc_zebrafish/nichenet/Nichenet_ligand-target_matrix.txt"

# 2. List of Transcription Factors (TFs) of interest
# These will be analyzed to find their upstream Ligands/Receptors
TARGET_TFS <- c(
    "Hoxb5",
    "Etv4"
)

# 3. List of Ligands/Receptors of interest
# These will be analyzed to find their downstream TFs and Targets
TARGET_LIGANDS <- c(
    "App", # REPLACE with your actual ligands
    "Jag2",
    "Sema4c",
    "Cntn1",
    "Agrn",
    "Sema3c",
    "Timp3",
    "Slit1",
    "Slit2"
)

# 4. Output directory for heatmaps
OUTPUT_DIR <- "/w5home/bmoore/Pierre_sc_zebrafish/nichenet/TF_Interaction_Heatmaps"

# 5. Top N interactions to show for Ligands (optional)
# Set to NULL to show all, or an integer (e.g., 50) to limit.
TOP_N_LIGAND_INTERACTIONS <- 50

# ==========================================
# HELPER FUNCTIONS
# ==========================================

read_matrix_file <- function(filepath) {
    if (!file.exists(filepath)) {
        stop(paste("File not found:", filepath))
    }
    message(paste("Reading matrix from:", filepath))

    if (grepl("\\.rds$", filepath, ignore.case = TRUE)) {
        mat <- readRDS(filepath)
    } else if (grepl("\\.csv$", filepath, ignore.case = TRUE)) {
        mat <- read.csv(filepath, row.names = 1, check.names = FALSE)
    } else {
        # Default to tab-separated for txt/tsv
        mat <- read.table(filepath, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
    }
    return(as.matrix(mat))
}

# Extract interactions for a specific feature (row or column)
# Returns a 1xN matrix where columns are the interactors
# Extract interactions for a specific feature by specifying axis
# axis = "row": Look for feature_name in rownames (e.g. TF in TF matrix), return Row (1xN).
# axis = "col": Look for feature_name in colnames (e.g. Ligand in any matrix), return Column (as 1xN transposed).
get_feature_interactions <- function(mat, feature_name, axis) {
    if (axis == "row") {
        if (feature_name %in% rownames(mat)) {
            return(mat[feature_name, , drop = FALSE])
        }
    } else if (axis == "col") {
        if (feature_name %in% colnames(mat)) {
            # Transpose so it returns as 1xN (Row = Feature, Cols = Interactors)
            return(t(mat[, feature_name, drop = FALSE]))
        }
    }
    return(NULL)
}

generate_sorted_heatmap <- function(data_matrix, feature_name, interaction_type, output_dir, top_n = NULL) {
    # data_matrix is 1 x N
    # Filter > 0
    valid_indices <- which(data_matrix[1, ] > 0)

    if (length(valid_indices) == 0) {
        warning(paste("No", interaction_type, "> 0 found for", feature_name))
        return(NULL)
    }

    # Subset
    data_subset <- data_matrix[, valid_indices, drop = FALSE]

    # Sort (Greatest to Smallest)
    sorted_indices <- order(data_subset[1, ], decreasing = TRUE)
    data_sorted <- data_subset[, sorted_indices, drop = FALSE]

    # Limit to top N if specified
    if (!is.null(top_n) && ncol(data_sorted) > top_n) {
        data_sorted <- data_sorted[, 1:top_n, drop = FALSE]
        message(paste("Limiting heatmap to top", top_n, "interactions for", feature_name))
    }

    # Define file path
    filename <- file.path(output_dir, paste0(feature_name, "_", interaction_type, "_heatmap.pdf"))

    # Dynamic width
    plot_width <- max(5, min(25, ncol(data_sorted) * 0.25))

    # Prepare arguments for pheatmap
    heatmap_args <- list(
        mat = data_sorted,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = FALSE,
        main = paste(feature_name, "->", interaction_type),
        filename = filename,
        width = plot_width,
        height = 3
    )

    # Handle low variance / constant values to prevent 'breaks are not unique' error
    unique_vals <- unique(as.vector(data_sorted))
    if (length(unique_vals) <= 1) {
        val <- unique_vals[1]
        # Create artificial breaks around the value to satisfy pheatmap
        if (abs(val) < 1e-9) {
            heatmap_args$breaks <- seq(-1, 1, length.out = 100)
        } else {
            # Create a small window around the value
            heatmap_args$breaks <- seq(val * 0.99, val * 1.01, length.out = 100)
        }
    }

    tryCatch(
        {
            do.call(pheatmap, heatmap_args)
            message(paste("Generated heatmap:", filename))
        },
        error = function(e) {
            message(paste("Error generating heatmap for", feature_name, ":", e$message))
        }
    )
}

# ==========================================
# MAIN SCRIPT
# ==========================================

# Create output directory
if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Global Data Cache (Load only if needed)
tf_matrix <- NULL
target_matrix <- NULL

# --- PART 1: Process Target TFs (Find Upstream Ligands) ---
if (length(TARGET_TFS) > 0) {
    message("\n=== Processing Target TFs ===")

    # Load TF Matrix if not loaded
    if (is.null(tf_matrix)) {
        tf_matrix <- read_matrix_file(INPUT_TF_MATRIX_FILE)
    }

    for (tf in TARGET_TFS) {
        message(paste("Analyzing TF:", tf))

        # Check if TF is in matrix
        # For TFs, we expect them to be Rows in the TF Matrix
        data <- get_feature_interactions(tf_matrix, tf, "row")

        if (!is.null(data)) {
            # If TF is Row, Interactors are Ligands
            generate_sorted_heatmap(data, tf, "Upstream_Ligands", OUTPUT_DIR)
        } else {
            warning(paste("TF", tf, "not found in TF Matrix."))
        }
    }
}

# --- PART 2: Process Target Ligands (Find Downstream TFs and Targets) ---
if (length(TARGET_LIGANDS) > 0 && !all(TARGET_LIGANDS == "Ligand1")) { # strict check to avoid running example
    message("\n=== Processing Target Ligands ===")

    # Load TF Matrix (for Ligand->TF) if not loaded
    if (is.null(tf_matrix)) {
        try({
            tf_matrix <- read_matrix_file(INPUT_TF_MATRIX_FILE)
        })
    }

    # Load Target Matrix (for Ligand->Target)
    if (is.null(target_matrix)) {
        try({
            target_matrix <- read_matrix_file(INPUT_TARGET_MATRIX_FILE)
        })
    }

    for (ligand in TARGET_LIGANDS) {
        message(paste("Analyzing Ligand:", ligand))

        # 1. Ligand -> TFs
        if (!is.null(tf_matrix)) {
            # Extract Ligand Interactions (Ligand is Column)
            data_tf <- get_feature_interactions(tf_matrix, ligand, "col")

            if (!is.null(data_tf)) {
                generate_sorted_heatmap(data_tf, ligand, "Downstream_TFs", OUTPUT_DIR, top_n = TOP_N_LIGAND_INTERACTIONS)
            } else {
                warning(paste("Ligand", ligand, "not found in TF Matrix."))
            }
        }

        # 2. Ligand -> Targets (Gene Expression)
        if (!is.null(target_matrix)) {
            # Extract Ligand Interactions (Ligand is Column)
            data_target <- get_feature_interactions(target_matrix, ligand, "col")

            if (!is.null(data_target)) {
                generate_sorted_heatmap(data_target, ligand, "Downstream_Targets", OUTPUT_DIR, top_n = TOP_N_LIGAND_INTERACTIONS)
            } else {
                warning(paste("Ligand", ligand, "not found in Target Matrix."))
            }
        }
    }
} else {
    message("\nSkipping Ligand analysis (TARGET_LIGANDS not set or is default).")
}

message("\nProcessing complete.")
