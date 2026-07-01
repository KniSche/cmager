"""
CMageR: Cell Type and Cell Age prediction for cardiac cells scRNA-seq between 4-15 PCW.
Cell Age Prediction Module.

Subprocess written in R:
1. Pulls in the CellTypist predictions  at coarse, mid, and fine_grain.
2. Applies coarse, mid, and fine_grain multivariate generalised additive models
3. sends results to python / metadata
"""

import os
import sys
import logging
import subprocess
import pandas as pd
import anndata as ad

logger = logging.getLogger("cmager")

def predict_age(
    adata: ad.AnnData,
    model_dir: str,
    sample_id: str, 
    output_dir: str, 
    chunk_id: str = "0",
    keep_temp_files: bool = False,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Executes developmental age prediction on single-cell data using R.
    
    Args:
        adata (ad.AnnData): The input single-cell dataset or chunk.
        model_dir (srt): The path to the models.
        sample_id (str): Unique identifier for the dataset (filename).
        output_dir (str): Base directory - needed to save temporary intermediate files.
        chunk_id (str): Identifier for the current streaming chunk.
        keep_temp_files (bool): If True, prevents deletion of the bridge CSVs.
        verbose (bool): Controls logging verbosity.
        
    Returns:
        pd.DataFrame: A dataframe containing the predicted cell ages.
    """
    logger.debug(f"   ⏳ Predicting cell ages for {sample_id} chunk {chunk_id}...")

    # --- PHASE 1: WORKSPACE & PATH SETUP ---
    CURRENT_FILE_DIR = os.path.dirname(os.path.abspath(__file__)) 
    PACKAGE_ROOT = os.path.dirname(CURRENT_FILE_DIR)             
    R_SCRIPT_PATH = os.path.join(PACKAGE_ROOT, "r_source", "run_age_predictions.R")
    
    temp_files_dir = os.path.join(output_dir, "temp_files")
    os.makedirs(temp_files_dir, exist_ok=True)
    
    # --- PHASE 2: DATA PREPARATION ---
    # Assign unique temp files using sample_id + chunk_id to avoid multi-process race conditions
    prefix = f"temp_{sample_id}_chk{chunk_id}"
    meta_path = os.path.join(temp_files_dir, f"{prefix}_obs.csv")
    chunk_h5ad_path = os.path.join(temp_files_dir, f"temp_chunk_{sample_id}_{chunk_id}.h5ad")
    r_results_path = os.path.join(temp_files_dir, f"{prefix}_age_outputs.csv")
    
    try:
        # Export metadata for R ingestion (aligning with hierarchical filtering pattern)
        adata.obs.to_csv(meta_path)
        
        # --- PHASE 3: EXTERNAL R EXECUTION ---
        # Fix an issue/warning on usage of "Poetry" virtual environments in Reticulate / R.
        custom_env = os.environ.copy()
        custom_env["RETICULATE_PYTHON"] = sys.executable 

        cmd = [
            "Rscript", R_SCRIPT_PATH,
            meta_path, chunk_h5ad_path, r_results_path, model_dir
        ]

        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            check=True, 
            env=custom_env 
        )
    
        # --- PHASE 4: RESULT RETRIEVAL ---
        age_prediction_outputs = pd.read_csv(r_results_path, index_col=0)
        return age_prediction_outputs

    except subprocess.CalledProcessError as e:
        logger.error(f"\n💥 R Age Predictor Failed on {sample_id} (Chunk {chunk_id})!")
        logger.error("=================== R STDOUT LOGS ===================")
        logger.error(e.stdout if e.stdout else "(No stdout output)")
        logger.error("=================== R STDERR LOGS ===================")
        logger.error(e.stderr if e.stderr else "(No stderr output)")
        logger.error("=====================================================")
        raise e
        
    finally:
        # --- PHASE 5: CLEANUP ---
        if not keep_temp_files:
            for temp_file in [meta_path, r_results_path]:
                if os.path.exists(temp_file): 
                    os.remove(temp_file)
        else:
            logger.debug(f"   📌 Retaining R interface CSV tables for chunk {chunk_id} inside 'temp_files'.")
