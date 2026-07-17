#!/usr/bin/env Rscript

# ==============================================================================
# CMageR: Age Prediction R Backend
# 
# Ingests a streaming chunk of single-cell data (.h5ad) and predicted labels (.obs),
# runs GAM-based age predictions across coarse, mid, and fine granularities, 
# and exports the predicted ages back to Python via a temporary CSV.
# ==============================================================================

# --- PHASE 1: CLI ARGUMENT PARSING ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("❌ Not enough arguments provided to Rscript. Expected 4.")
}

meta_path        <- args[1] # The metadata containing predicted labels (.obs)
chunk_h5ad_path  <- args[2] # The physical AnnData chunk containing raw matrix (.h5ad)
r_results_path   <- args[3] # Target output CSV path for age predictions
age_models_dir   <- args[4] # Directory containing the pre-trained GAM models

# --- PHASE 2: ENVIRONMENT & SCRIPT RESOLUTION ---
suppressPackageStartupMessages({
  library(anndata)
  library(Matrix)
  library(matrixStats)
  library(gamsel)
})

# Dynamically resolve where this script lives so we can find the "functions" folder reliably
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("--file=", initial_options, value = TRUE)

if (length(file_arg) > 0) {
  # Executed via Rscript command line
  script_dir <- dirname(normalizePath(substring(file_arg, 8)))
} else {
  # Fallback for manual RStudio/interactive execution
  script_dir <- file.path(getwd(), "cmager", "r_source")
}

# Source the helper functions
source(file.path(script_dir, "functions", "gam_modelling_functions.r"))

# --- PHASE 3: MODEL & RESOURCE LOADING ---
# Load GAM models
model_list_coarse <- readRDS(file.path(age_models_dir, "gamsel_list_coarse_light.rds"))
model_list_mid    <- readRDS(file.path(age_models_dir, "gamsel_list_mid_light.rds"))
model_list_fine   <- readRDS(file.path(age_models_dir, "foetal_age_gamsel_list_fine.rds"))

model_list_all <- list(
  "coarse_grain" = model_list_coarse, 
  "mid_grain"    = model_list_mid, 
  "fine_grain"   = model_list_fine
)

# Load prior data objects
cell_adjust <- readRDS(file.path(script_dir, "resources", "cell_adjust_total.rds"))
tech_genes  <- readRDS(file.path(script_dir, "resources", "CvN.rds"))

# --- PHASE 4: DATA INGESTION ---
# Load the AnnData chunk
adata <- anndata::read_h5ad(chunk_h5ad_path)
input_data <- adata$X
meta_df <- adata$obs
  rownames(meta_df) = adata$obs_names # some weird issue with rownames not being taken from the obs df
modality <- meta_df[, "modality"]

# Format prior predictions into the list structure required by the GAM functions
experiment_predictions <- list(
  "coarse_grain" = data.frame(
    "barcode"     = rownames(meta_df), 
    "predictions" = meta_df$cmager_celltype_coarse,
    row.names     = rownames(meta_df)
  ),
  "mid_grain"    = data.frame(
    "barcode"     = rownames(meta_df), 
    "predictions" = meta_df$cmager_celltype_mid,
    row.names     = rownames(meta_df)
  ),
  "fine_grain"   = data.frame(
    "barcode"     = rownames(meta_df), 
    "predictions" = meta_df$cmager_celltype_fine,
    row.names     = rownames(meta_df)
  )
)

# --- PHASE 5: AGE PREDICTION EXECUTION ---

# Note on fine_grain: Input data from Python is log(cp10k). 
# The model was trained on standard cp10k. 
# The internal R function handles this reversal (expm1) specifically for fine_grain.

granularities <- c("coarse_grain", "mid_grain", "fine_grain")

for (gran in granularities) {
  experiment_predictions <- predict_age_fromCT_results(
    input_data                  = input_data,
    experiment_predictions      = experiment_predictions,
    model_list                  = model_list_all[[gran]],
    granularity                 = gran,
    modality                    = modality,
    good_barcodes               = NULL,
    is_already_cp10k_normalised = TRUE,
    library_sizes               = NULL,
    log_transform               = FALSE,
    add_to_experimentRDS        = TRUE
  )
}

# --- PHASE 6: RESULT EXPORT ---
target_col <- "predicted_age" 

results_df <- data.frame(
  cmager_age_coarse = experiment_predictions[["coarse_grain"]][[target_col]],
  cmager_age_mid    = experiment_predictions[["mid_grain"]][[target_col]],
  cmager_age_fine   = experiment_predictions[["fine_grain"]][[target_col]],
  row.names         = rownames(meta_df)
)

write.csv(results_df, r_results_path, row.names = TRUE)
cat("   📤 [R-Age] Age predictions success. Handing process back to Python.\n")
