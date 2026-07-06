#!/usr/bin/env Rscript

# ==============================================================================
# CMageR: Hierarchical Filtering R Backend
# 
# Takes predicted labels (.obs) and multi-resolution probability matrices,
# maps them against structural cluster lineages, filters out low-confidence 
# assignments, and passes the consensus-filtered results back to Python.
# ==============================================================================

# --- PHASE 1: CLI ARGUMENT PARSING ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("❌ Not enough arguments provided to Rscript. Expected 5.")
}

meta_path        <- args[1] 
prob_coarse_path <- args[2] 
prob_mid_path    <- args[3] 
prob_fine_path   <- args[4] 
r_results_path   <- args[5] 

# ─── FIX 3: CHANGE THE R WORKING DIRECTORY TO YOUR CLEAN CHUNK STORAGE ──────
# If the sourced functions use relative paths like "R_temp", it will now generate
# safely inside your writable, isolated chunk folder.
setwd(dirname(r_results_path))


# --- PHASE 2: ENVIRONMENT & SCRIPT RESOLUTION ---
# Dynamically resolve where this script lives to source relative helper files robustly
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", initial_options, value = TRUE)

if (length(file_arg) > 0) {
  # Executed via Rscript command line
  script_dir <- dirname(normalizePath(substring(file_arg, 8)))
} else {
  # Fallback for manual RStudio/interactive execution
  script_dir <- file.path(getwd(), "cmager", "r_source")
}

# Source the helper validation and loading suites
source(file.path(script_dir, "functions", "CT_prediction_functions.r"))
source(file.path(script_dir, "functions", "CT_load_data.r")) 
# Note: 'cluster_grouping' variable is automatically initialized by the source scripts above

# --- PHASE 3: DATA INGESTION ---
# Load in raw probability matrices from CellTypist chunk exports
all_matrices <- list(
  "coarse_grain" = read.ct.res(prob_coarse_path),
  "mid_grain"    = read.ct.res(prob_mid_path),
  "fine_grain"   = read.ct.res(prob_fine_path)
)

# Load chunk metadata observations
meta_df <- read.ct.res(meta_path)

# Restructure current raw assignments into structured granular data frames
all_predictions <- list(
  "coarse_grain" = data.frame(rownames(meta_df), meta_df$cmager_celltype_coarse),
  "mid_grain"    = data.frame(rownames(meta_df), meta_df$cmager_celltype_mid),
  "fine_grain"   = data.frame(rownames(meta_df), meta_df$cmager_celltype_fine)
)

# --- PHASE 4: HIERARCHICAL CO-ASSIGNMENT MAPS ---
# Re-label predictions based on ancestral/parent identities across the lineage tree
all_predictions <- prediction.df(
  prediction_dataframe_list = all_predictions,
  probability_matrix_list   = all_matrices,
  cluster_grouping          = cluster_grouping
)

# --- PHASE 5: LOW-PROBABILITY CONSENSUS FILTERING ---
prob_thresh <- 10
use_odds    <- TRUE

# Flag low-confidence or conflicting nodes and transition them cleanly to 'unknown'
predictions_processed <- filter_low_probability(
  all_predictions,
  prob_thresh,
  odds_ratio = use_odds
)

# --- PHASE 6: RESULT EXPORT ---
target_col <- "predictions2" 

results_df <- data.frame(
  predictions_filtered_coarse = predictions_processed[["coarse_grain"]][[target_col]],
  predictions_filtered_mid    = predictions_processed[["mid_grain"]][[target_col]],
  predictions_filtered_fine   = predictions_processed[["fine_grain"]][[target_col]],
  row.names                   = rownames(meta_df)
)

write.csv(results_df, r_results_path, row.names = TRUE)
cat("   📤 [R-Filter] Consensus filtering success. Handing process back to Python.\n")
