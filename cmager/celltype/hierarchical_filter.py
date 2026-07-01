"""
CMageR: Cell Type and Cell Age prediction for cardiac cells scRNA-seq between 4-15 PCW.
Hierarchical Filtering Module.

Subprocess written in R:
1. Pulls in the CellTypist predictions  at coarse, mid, and fine_grain.
2. Subsets probability matrices of mid_grain and fine_grain predictions based on the coarse_grain label
3. Additionally runs probability filters using odds ratio between first and second predicted class
4. spits out "unknown" class in _filtered prediction columns and sends back to python metadata.
"""

import os
import logging
import subprocess
import pandas as pd
import anndata as ad

logger = logging.getLogger("cmager")

def hierarchical_filtering(
    adata: ad.AnnData, 
    sample_id: str, 
    output_dir: str, 
    chunk_id: str = "0",
    keep_temp_files: bool = False,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Executes hierarchical cell-type re-labelling.
    
    Args:
        adata (ad.AnnData): The input single-cell dataset or chunk.
        sample_id (str): Unique identifier for the dataset (filename).
        output_dir (str): Base directory - needed to save temporary intermediate files.
        chunk_id (str): Identifier for the current streaming chunk.
        keep_temp_files (bool): If True, prevents deletion of the bridge CSVs.
        verbose (bool): Controls logging verbosity.
        
    Returns:
        pd.DataFrame: A dataframe containing the consensus-filtered predictions.
    """
    logger.debug(f"   📊 [R-Filter] Running consensus filtering for sample {sample_id}, chunk {chunk_id}...")
    
    # --- PHASE 1: WORKSPACE & PATH SETUP ---
    CURRENT_FILE_DIR = os.path.dirname(os.path.abspath(__file__)) 
    PACKAGE_ROOT = os.path.dirname(CURRENT_FILE_DIR)             
    R_SCRIPT_PATH = os.path.join(PACKAGE_ROOT, "r_source", "run_hierarchical_filter.R")
    
    temp_files_dir = os.path.join(output_dir, "temp_files")
    os.makedirs(temp_files_dir, exist_ok=True)
    
    # Assign unique temp files using sample_id + chunk_id to avoid multi-process race conditions
    prefix = f"temp_{sample_id}_chk{chunk_id}"
    meta_path = os.path.join(temp_files_dir, f"{prefix}_obs.csv")
    prob_coarse_path = os.path.join(temp_files_dir, f"{prefix}_probs_coarse.csv")
    prob_mid_path = os.path.join(temp_files_dir, f"{prefix}_probs_mid.csv")
    prob_fine_path = os.path.join(temp_files_dir, f"{prefix}_probs_fine.csv")
    r_results_path = os.path.join(temp_files_dir, f"{prefix}_filter_outputs.csv")
    
    try:
        # --- PHASE 2: DATA EXPORT FOR R INGESTION ---
        adata.obs.to_csv(meta_path)
        
        # Reconstruct assignment tables using class maps from .uns
        df_coarse = pd.DataFrame(adata.obsm["celltypist_probs_coarse"], index=adata.obs_names)
        if "celltypist_classes_coarse" in adata.uns:
            df_coarse.columns = adata.uns["celltypist_classes_coarse"]
        df_coarse.to_csv(prob_coarse_path)
        
        df_mid = pd.DataFrame(adata.obsm["celltypist_probs_mid"], index=adata.obs_names)
        if "celltypist_classes_mid" in adata.uns:
            df_mid.columns = adata.uns["celltypist_classes_mid"]
        df_mid.to_csv(prob_mid_path)
        
        df_fine = pd.DataFrame(adata.obsm["celltypist_probs_fine"], index=adata.obs_names)
        if "celltypist_classes_fine" in adata.uns:
            df_fine.columns = adata.uns["celltypist_classes_fine"]
        df_fine.to_csv(prob_fine_path)
        
        # --- PHASE 3: EXTERNAL R EXECUTION ---
        cmd = [
            "Rscript", R_SCRIPT_PATH,
            meta_path, prob_coarse_path, prob_mid_path, prob_fine_path, r_results_path
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # --- PHASE 4: RESULT RETRIEVAL ---
        filter_outputs = pd.read_csv(r_results_path, index_col=0)
        return filter_outputs

    except subprocess.CalledProcessError as e:
        # Changed to logger.error so crashes prominently break through terminal noise
        logger.error(f"\n💥 R Script Failed on {sample_id} (Chunk {chunk_id})!")
        logger.error("=================== R STDOUT LOGS ===================")
        logger.error(e.stdout if e.stdout else "(No stdout output)")
        logger.error("=================== R STDERR LOGS ===================")
        logger.error(e.stderr if e.stderr else "(No stderr output)")
        logger.error("=====================================================")
        raise e
        
    finally:
        # --- PHASE 5: CLEANUP ---
        if not keep_temp_files:
            for temp_file in [meta_path, prob_coarse_path, prob_mid_path, prob_fine_path, r_results_path]:
                if os.path.exists(temp_file): 
                    os.remove(temp_file)
        else:
            logger.debug(f"   📌 Retaining R interface CSV tables for chunk {chunk_id} inside 'temp_files'.")
