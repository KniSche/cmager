"""
CMageR: Cell Type and Cell Age prediction for cardiac cells scRNA-seq between 4-15 PCW.
CellTypist Prediction Module.

Handles data normalization and multi-granularity cell type classification 
using pre-trained CellTypist models from the fetal heart atlas (Teichmann lab, Sinha lab 2026)
"""

import os
import sys
import logging
from contextlib import contextmanager

import scanpy as sc
import pandas as pd
import celltypist
import anndata as ad

logger = logging.getLogger("cmager")

@contextmanager
def suppress_stdout(enable: bool):
    """
    Optionally redirects standard output to devnull.
    Used to muzzle hardcoded third-party prints (like CellTypist) during clean runs.
    """
    if enable:
        with open(os.devnull, 'w') as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:
                yield
            finally:
                sys.stdout = old_stdout
    else:
        yield


def predict_celltypes(
    adata: ad.AnnData, 
    model_dir: str, 
    low_ram: bool = False,
    verbose: bool = False
) -> ad.AnnData:
    """
    Normalizes data and runs multi-resolution CellTypist predictions 
    directly on the provided AnnData object or chunk slice.
    
    Args:
        adata (ad.AnnData): The input single-cell dataset or matrix slice.
        model_dir (str): Absolute or relative path to the pre-trained .pkl models.
        low_ram (bool): Inherited flag for downstream behavior mapping (unused here but preserved for API consistency).
        verbose (bool): If False, suppresses internal CellTypist terminal outputs.
                          
    Returns:
        ad.AnnData: The updated AnnData object with cell types stored in .obs 
                    and probability matrices stored in .obsm.
    """
    
    # --- PHASE 1: NORMALIZATION ---
    if 'log1p' not in adata.uns:
        logger.debug("   ⚙️ [CellTypist] Applying target normalization and log1p...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
    # --- PHASE 2: MODEL LOADING ---
    coarse_path = os.path.join(model_dir, "foetal_heart_coarse_v3.pkl")
    mid_path = os.path.join(model_dir, "foetal_heart_mid_v3.pkl")
    fine_path = os.path.join(model_dir, "foetal_heart_fine_v3.pkl")
    
    model_coarse = celltypist.models.Model.load(model=coarse_path)
    model_mid = celltypist.models.Model.load(model=mid_path)
    model_fine = celltypist.models.Model.load(model=fine_path)

    # --- PHASE 3: MULTI-RESOLUTION PREDICTION ---
    logger.debug("   🧠 Running multi-resolution CellTypist predictions...")
    
    # Safely silence the noisy internal print statements from celltypist
    with suppress_stdout(enable=not verbose):
        pred_coarse = celltypist.annotate(adata, model=model_coarse)
        pred_mid = celltypist.annotate(adata, model=model_mid)
        pred_fine = celltypist.annotate(adata, model=model_fine)
    
    # --- PHASE 4: METADATA INJECTION ---
    
    # 4a. Assign discrete string labels cleanly to .obs
    adata.obs['cmager_celltype_coarse'] = pred_coarse.predicted_labels['predicted_labels']
    adata.obs['cmager_celltype_mid'] = pred_mid.predicted_labels['predicted_labels']
    adata.obs['cmager_celltype_fine'] = pred_fine.predicted_labels['predicted_labels']

    # 4b. Assign numerical probability matrices safely to .obsm
    adata.obsm['celltypist_probs_coarse'] = pred_coarse.probability_matrix.values
    adata.obsm['celltypist_probs_mid'] = pred_mid.probability_matrix.values
    adata.obsm['celltypist_probs_fine'] = pred_fine.probability_matrix.values

    # 4c. Track classification targets in unstructured metadata (.uns)
    adata.uns['celltypist_classes_coarse'] = list(model_coarse.classifier.classes_)
    adata.uns['celltypist_classes_mid']    = list(model_mid.classifier.classes_)
    adata.uns['celltypist_classes_fine']   = list(model_fine.classifier.classes_)
    
    return adata
