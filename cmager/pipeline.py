"""
CMageR: Cell Type and Cell Age prediction for cardiac cells scRNA-seq between 4-15 PCW.

1. Performs cell type annotation through semi-hierarchical classification with CellTypist
2. Performs cell age prediction through multivariate generalised additive model

In both cases, data were carefully curated from fetal heart single cell and single nuclei RNA sequencing.

Functions here handles dataset discovery, batch splitting, and chunk-based processing (for controlling RAM usage),
parallel processing, and final matrix aggregation.
"""

import os
import gc
import glob
import math
import time
import logging
import multiprocessing
import shutil
from typing import Dict, List, Optional

import pandas as pd
import anndata as ad
import scanpy as sc
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

# Local Package Imports
from cmager.utils import HGNCManifoldMapper 
from cmager.celltype import predict_celltypes, hierarchical_filtering
from cmager.cellage import predict_age

# Initialise the central package logger
logger = logging.getLogger("cmager")


def setup_pipeline_logging(verbose: bool):
    """
    For --verbose == True , mutes MOST of the messages 
    (although, CellTypist and some of the hierarchical filter steps still get through!)
    """
    logger.propagate = False
    logger.handlers.clear()
    
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(message)s") 
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    
    if verbose:
        logger.setLevel(logging.DEBUG)
        logging.getLogger("celltypist").setLevel(logging.INFO)
    else:
        logger.setLevel(logging.INFO)
        # Forcefully silence noisy upstream annotation dependencies
        logging.getLogger("celltypist").setLevel(logging.WARNING)
        logging.getLogger("scanpy").setLevel(logging.WARNING)
        logging.getLogger("anndata").setLevel(logging.WARNING)


def discover_datasets(input_dir: str) -> Dict[str, dict]:
    """
    Scans the input directory and identifies individual datasets.
    
    Args:
        input_dir: Path to the target directory, can contain h5ads or 10x sub-directories.
        
    Returns:
        Python dictionary that maps unique sample_ids to their respective paths and data types.
    """
    datasets = {}
    
    # 1. Discover all standalone .h5ad files
    h5ad_files = glob.glob(os.path.join(input_dir, "*.h5ad"))
    for filepath in h5ad_files:
        sample_id = os.path.splitext(os.path.basename(filepath))[0]
        datasets[sample_id] = {"path": filepath, "type": "h5ad"}
        
    # 2. Discover 10X subdirectories
    for entry in os.scandir(input_dir):
        if entry.is_dir():
            path = entry.path
            is_10x = all(
                os.path.exists(os.path.join(path, f)) or os.path.exists(os.path.join(path, f + ".gz"))
                for f in ['matrix.mtx', 'features.tsv', 'barcodes.tsv']
            )
            if is_10x:
                datasets[entry.name] = {"path": path, "type": "10x"}
                
    return datasets


def pre_split_datasets_by_batch(discovered: Dict[str, dict], batch_key: str, output_dir: str) -> Dict[str, dict]:
    """
    Scans discovered datasets and splits them into temporary multi-batch .h5ad files.
    For example, if you wanted to split by a prior known factor such as "donor", or "orig_ident"...
    Batch has no effect in 10X file inputs as there is no readable metadata.
    """
    if not batch_key:
        return discovered

    optimised_discovered = {}
    temp_dir = os.path.join(output_dir, "temp_files")
    os.makedirs(temp_dir, exist_ok=True)
    
    for sample_id, info in discovered.items():
        # 10X folders are single-capture libraries; pass them straight through
        if info["type"] == "10x":
            optimised_discovered[sample_id] = info
            continue

        logger.debug(f"✂️  Pre-splitting: Inspecting batch keys for dataset [{sample_id}]...")
        adata = sc.read_h5ad(info["path"])

        if batch_key not in adata.obs.columns:
            logger.info(f"⚠️ Warning: Batch key '{batch_key}' not found in {sample_id}. Processing as a single block.")
            optimised_discovered[sample_id] = info
            del adata; gc.collect()
            continue

        unique_batches = adata.obs[batch_key].dropna().unique()
        
        for batch_val in unique_batches:
            # Create a safe, filesystem-friendly tracking name
            sanitised_val = str(batch_val).replace(" ", "_").replace("/", "-")
            sub_sample_id = f"{sample_id}_batch_{sanitised_val}"
            sub_h5ad_path = os.path.join(temp_dir, f"temp_split_{sub_sample_id}.h5ad")
            
            logger.debug(f"    └── Creating temporary batch file: {sub_h5ad_path}")
            batch_adata = adata[adata.obs[batch_key] == batch_val].copy()
            batch_adata.write_h5ad(sub_h5ad_path)

            optimised_discovered[sub_sample_id] = {
                "path": sub_h5ad_path,
                "type": "h5ad",
                "is_temporary": True
            }

        del adata; gc.collect()

    return optimised_discovered
    

def standardise_gene_symbols(adata: ad.AnnData, ref_features_df: pd.DataFrame, min_match_threshold: float = 0.10) -> ad.AnnData:
    """
    Scans the adata.var_names and all var columns against the model reference list.
    Independently locates the best columns for Gene Symbols and Ensembl IDs.
    """
    # Extract reference sets using column positions (0 = Ensembl, 1 = Symbol)
    ref_ensembl = set(ref_features_df.iloc[:, 0].dropna().astype(str).tolist())
    ref_symbols = set(ref_features_df.iloc[:, 1].dropna().astype(str).tolist())

    if not ref_symbols and not ref_ensembl:
        return adata

    # Helper function to calculate both match rates simultaneously
    def get_match_rates(gene_list):
        ens_matches = sum(1 for g in gene_list if g in ref_ensembl)
        sym_matches = sum(1 for g in gene_list if g in ref_symbols)
        ens_rate = ens_matches / len(ref_ensembl) if ref_ensembl else 0
        sym_rate = sym_matches / len(ref_symbols) if ref_symbols else 0
        return ens_rate, sym_rate

    # 1. Evaluate current var_names index
    current_genes = adata.var_names.astype(str).tolist()
    best_ens_rate, best_sym_rate = get_match_rates(current_genes)
    best_ens_col, best_sym_col = 'index', 'index'
    
    logger.debug(f"Checking index. Initial rates -> Ensembl: {best_ens_rate:.2%} | Symbol: {best_sym_rate:.2%}")

    # 2. Scan available columns for superior matches
    for col in adata.var.columns:
        sample_genes = adata.var[col].astype(str).tolist()
        ens_r, sym_r = get_match_rates(sample_genes)
        
        if ens_r > best_ens_rate:
            best_ens_rate, best_ens_col = ens_r, col
            
        if sym_r > best_sym_rate:
            best_sym_rate, best_sym_col = sym_r, col

    # 3. Apply Ensembl IDs (Force into 'gene_ids' for downstream mapping)
    if best_ens_rate >= min_match_threshold:
        logger.debug(f"Rescuing Ensembl IDs from '{best_ens_col}' (Match rate: {best_ens_rate:.2%})")
        if best_ens_col == 'index':
            adata.var['gene_ids'] = adata.var_names
        else:
            adata.var['gene_ids'] = adata.var[best_ens_col].astype(str)

    # 4. Apply Gene Symbols (Switch var_names if a better column exists)
    if best_sym_rate >= min_match_threshold:
        if best_sym_col != 'index':
            logger.debug(f"Switching var_names to superior symbol column: '{best_sym_col}' (Match rate: {best_sym_rate:.2%})")
            
            if adata.var_names.name is None:
                adata.var["original_features"] = adata.var_names
            else:
                adata.var[adata.var_names.name] = adata.var_names
                
            adata.var_names = adata.var[best_sym_col].astype(str)
    else:
        logger.info(f"⚠️ Warning: Best gene symbol match rate is low ({best_sym_rate:.2%}). Relying on Ensembl IDs for mapping.")

    adata.var_names_make_unique()
    return adata
    


def process_single_sample(
    sample_id: str, 
    info: dict, 
    internal_ref_features: str, 
    celltype_model_dir: str, 
    cellage_model_dir: str, 
    output_dir: str,
    chunk_size: int, 
    skip_reductions: bool, 
    modality: str = None, 
    keep_temp_files: bool = False, 
    verbose: bool = False, 
    progress_queue=None
) -> str:
    """
    The core function, stringing together the different models and results.
    """
    setup_pipeline_logging(verbose)
    logger.debug(f"\n🚀 Starting pipeline for: {sample_id}")
    
    sample_output_h5ad = os.path.join(output_dir, f"{sample_id}_annotated.h5ad")
    sample_output_csv = os.path.join(output_dir, f"{sample_id}_metadata_table.csv")
    worker_temp_dir = os.path.join(output_dir, "temp_files")
    os.makedirs(worker_temp_dir, exist_ok=True)
    
    temp_chunk_paths = []
    
    try:
        # --- PHASE 1: GLOBAL RESOURCE INITIALIZATION ---
        if not os.path.exists(internal_ref_features):
            raise FileNotFoundError(f"❌ Reference feature file not found at: {internal_ref_features}")
            
        ref_features_df = pd.read_csv(internal_ref_features, sep=None, engine='python', header=None)
        mapper = HGNCManifoldMapper(cache_dir="./metadata_cache")
        audit_csv_path = os.path.join(output_dir, f"{sample_id}_gene_mapping_audit.csv")

        # --- PHASE 2: MEMORY-EFFICIENT INGESTION SETUP ---
        if info["type"] == "h5ad":
            source_stream = ad.read_h5ad(info["path"], backed="r")
            total_cells = source_stream.n_obs
        elif info["type"] == "10x":
            logger.debug(f"📦 Loading 10X matrix into memory once for slicing...")
            source_stream = sc.read_10x_mtx(info["path"], var_names='gene_symbols', cache=False)
            total_cells = source_stream.n_obs
        else:
            raise ValueError(f"Unknown dataset type for {sample_id}")

        # --- PHASE 3: CHUNK STREAMING & ANNOTATION ---
        first_chunk = True 
        
        for start_idx in range(0, total_cells, chunk_size):
            end_idx = min(start_idx + chunk_size, total_cells)
            chunk_idx = str(start_idx // chunk_size)
            
            logger.info(f"⏳ Processing sample {sample_id}, chunk {chunk_idx} (cells {start_idx} to {end_idx})...")
            
            # 1. Take the slice (creates a view)
            chunk_view = source_stream[start_idx:end_idx]
            
            # 2. Safely resolve the view into a standalone, in-memory object
            if info["type"] == "h5ad":
                chunk = chunk_view.to_memory()
            else:
                # 10X matrices are already loaded in memory, so we just copy the view
                chunk = chunk_view.copy()
            
            # store the original barcode (as supplied in file) for easily pulling it out   
            if "original_barcode" not in chunk.obs.columns:
                chunk.obs["original_barcode"] = chunk.obs_names.astype(str)

            prefix = f"{sample_id}_"
            chunk.obs_names = [
                bc if bc.startswith(prefix) else f"{prefix}{bc}"
                for bc in chunk.obs_names.astype(str)
            ]  
                
            # 3a. Modality validation
            if modality:
                if modality in ["cells", "nuclei"]:
                    chunk.obs["modality"] = modality
                elif modality in chunk.obs.columns:
                    chunk.obs["modality"] = chunk.obs[modality].astype(str)
                else:
                    raise KeyError(f"❌ Modality string '{modality}' is neither 'cells'/'nuclei' nor an existing column in .obs.")
            elif "modality" not in chunk.obs.columns:
                raise KeyError(f"❌ No modality supplied via CLI, and no fallback 'modality' column found in dataset metadata.")
                    
            # 3b. Gene Alignment & Standardization
            chunk = standardise_gene_symbols(chunk, ref_features_df)
            
            current_audit_path = audit_csv_path if start_idx == 0 else None
            chunk = mapper.align_by_ensembl_and_manifold(
                adata=chunk, 
                model_features_df=ref_features_df,
                output_csv_path=current_audit_path
            )
            
            
            ################################
            # 2d. CellTypist Predictions
            chunk = predict_celltypes(
              chunk, 
              model_dir=celltype_model_dir, 
              skip_reductions=True,
              verbose=verbose
            ) 
              
            # 2e. R Hierarchical Filtering
            filter_results = hierarchical_filtering(
              chunk, 
              sample_id=sample_id, 
              output_dir=output_dir, 
              chunk_id=chunk_idx,
              keep_temp_files=keep_temp_files,
              verbose=verbose,
            )
            chunk.obs['cmager_celltype_coarse'] = filter_results['predictions_coarse'].values
            chunk.obs['cmager_celltype_mid']    = filter_results['predictions_mid'].values
            chunk.obs['cmager_celltype_fine']   = filter_results['predictions_fine'].values
            chunk.obs['cmager_celltype_coarse_filtered'] = filter_results['predictions_filtered_coarse'].values
            chunk.obs['cmager_celltype_mid_filtered']    = filter_results['predictions_filtered_mid'].values
            chunk.obs['cmager_celltype_fine_filtered']   = filter_results['predictions_filtered_fine'].values
            chunk.obs['cmager_coarse_probability']       = filter_results['predictions_coarse_probability'].values
            
            
            ################################
            # Resolve AnnData 'uns' bugs with scanpy log1p
            if 'log1p' in chunk.uns:
              del chunk.uns['log1p']
              
            # Write chunks for the age prediction model
            chunk_path = os.path.join(worker_temp_dir, f"temp_chunk_{sample_id}_{chunk_idx}.h5ad")  
            chunk.write_h5ad(chunk_path)
             
            ################################
            # 2f. Age Predictions
            age_results = predict_age(
              chunk, 
              model_dir=cellage_model_dir, 
              sample_id=sample_id, 
              output_dir=output_dir, 
              chunk_id=chunk_idx,
              keep_temp_files=keep_temp_files,
              verbose=verbose
            )
            chunk.obs['cmager_age_coarse'] = age_results['cmager_age_coarse'].values
            chunk.obs['cmager_age_mid']    = age_results['cmager_age_mid'].values
            chunk.obs['cmager_age_fine']   = age_results['cmager_age_fine'].values
            ################################
            
            # 2g. Metadata Extraction & Storage
            chunk.obs["sample_id"] = sample_id
            chunk_metadata = chunk.obs.copy()
            
            # Write chunks again
            chunk.write_h5ad(chunk_path)
            temp_chunk_paths.append(chunk_path)

            # Iteratively build the CSV tracking log
            if first_chunk:
                chunk_metadata.to_csv(sample_output_csv, mode='w', index=True, index_label="cell_barcode")
                first_chunk = False
            else:
                chunk_metadata.to_csv(sample_output_csv, mode='a', header=False, index=True)
            
            del chunk, chunk_metadata
            gc.collect()
            
            if progress_queue is not None:
              progress_queue.put("chunk_complete")
            
        # --- PHASE 4: CHUNK REASSEMBLY ---
        logger.debug(f"📦 Assembling processed chunks for {sample_id}...")
        loaded_chunks = [ad.read_h5ad(p) for p in temp_chunk_paths]
        
        sample_adata = ad.concat(loaded_chunks, axis=0, join="outer", merge="same")
        
        
        
        # calculate UMAPs on the per-sample_ID basis
        if not skip_reductions:
            logger.debug(f"🗺️ Calculating UMAP for sample {sample_id}...")
            bdata = sample_adata.copy() 

            sc.pp.filter_genes(bdata, min_cells=3)
            sc.pp.highly_variable_genes(bdata, min_mean=0.0125, max_mean=3, min_disp=0.5)
            bdata = bdata[:, bdata.var.highly_variable].copy()

            sc.pp.scale(bdata, max_value=10)
            sc.tl.pca(bdata, svd_solver='arpack')
            
            sc.pp.neighbors(bdata, n_neighbors=10, n_pcs=40)
            sc.tl.umap(bdata)

            # Map coordinates back to the original object
            sample_adata.obsm["X_umap"] = bdata.obsm["X_umap"]

            del bdata
            gc.collect()
            logger.debug(f"🎉 CellType and UMAP annotation complete!")
        else:
            logger.debug(f"⏭️ skipping PCA and UMAP calculations (--skip-reductions)")
            
        
        # write / save adata
        sample_adata.write_h5ad(sample_output_h5ad)
        
        if not keep_temp_files:
            for p in temp_chunk_paths:
                if os.path.exists(p): os.remove(p)

        return sample_output_h5ad
        
    except Exception as e:
        # Fallback cleanup on failure
        if not keep_temp_files:
            for p in temp_chunk_paths:
                if os.path.exists(p): os.remove(p)
        raise e
        
    finally:
        # Clean up global stream pointer to release file locks safely
        if 'source_stream' in locals():
            del source_stream
        gc.collect()


def concatenate_and_save(processed_file_paths: List[str], output_dir: str, original_barcode_order: Optional[List[str]] = None):
    """
    Reads the processed files from temp output dir, aggregates 
    a master AnnData object and a master metadata data table,
    re-indexes rows to perfectly match the original sequence layout,
    and then neatens the workspace by sorting per-sample files into folders.
    """
    if not processed_file_paths:
        logger.info("⚠️ No datasets found to aggregate.")
        return

    master_h5ad_path = os.path.join(output_dir, "cmager_results_anndata.h5ad")
    master_csv_path = os.path.join(output_dir, "cmager_results_metadata.csv")
    
    logger.debug(f"🧩 Merging tracking records for {len(processed_file_paths)} datasets...")

    # --- PART 1: CONCATENATING LIGHTWEIGHT METADATA TABLES ---
    first_csv = True
    for h5ad_path in processed_file_paths:
        csv_path = h5ad_path.replace("_annotated.h5ad", "_metadata_table.csv")
        
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            if first_csv:
                df.to_csv(master_csv_path, mode='w', index=False)
                first_csv = False
            else:
                df.to_csv(master_csv_path, mode='a', header=False, index=False)
    
    logger.debug(f"💾 Master metadata table successfully pooled at: {master_csv_path}")

    # --- PART 2: MATRIX CONCATENATION ---
    logger.debug("📦 Pooling expression matrices sequentially...")
    
    first_adata = ad.read_h5ad(processed_file_paths[0])
    
    if len(processed_file_paths) == 1:
        aggregated_adata = first_adata
    else:
        remaining_paths = processed_file_paths[1:]
        
        # Yield one dataset at a time into the concat function to prevent RAM spikes
        def stream_adatas():
            yield first_adata
            for path in remaining_paths:
                current_adata = ad.read_h5ad(path)
                yield current_adata
                del current_adata
                gc.collect()

        aggregated_adata = ad.concat(
            stream_adatas(),
            axis=0,
            join="outer",
            label="sample_batch",
            merge="same"
        )

    # Ensure that the final object is in the same order as the input object(s)
    if original_barcode_order:
        logger.debug("🎯 Restoring original barcode sequence sorting across final object...")
        
        # 1. Align the expression matrix reference matrix-side
        aggregated_adata = aggregated_adata[original_barcode_order].copy()
        
        # 2. Align the master tracking table CSV structure
        master_df = pd.read_csv(master_csv_path)
        master_df.set_index("cell_barcode", inplace=True, drop=False)
        master_df = master_df.reindex(original_barcode_order)
        master_df.to_csv(master_csv_path, index=False)
        
    # Write aggregated matrix tracking object
    aggregated_adata.write_h5ad(master_h5ad_path)
    logger.debug(f"💾 All objects aggregated and saved to: {master_h5ad_path}")

    # --- PART 3: POST-RUN WORKSPACE TIDY UP ---
    logger.debug("🧹 Sorting individual sample output files into dedicated subfolders...")
    
    for h5ad_path in processed_file_paths:
        filename = os.path.basename(h5ad_path)
        sample_id = filename.replace("_annotated.h5ad", "")
        
        sample_folder = os.path.join(output_dir, sample_id)
        os.makedirs(sample_folder, exist_ok=True)
        
        file_targets = [
            f"{sample_id}_annotated.h5ad",
            f"{sample_id}_metadata_table.csv",
            f"{sample_id}_gene_mapping_audit.csv"
        ]
        
        for target_file in file_targets:
            source_path = os.path.join(output_dir, target_file)
            destination_path = os.path.join(sample_folder, target_file)
            
            if os.path.exists(source_path):
                shutil.move(source_path, destination_path)
                logger.debug(f"  Moved: {target_file} -> {sample_id}/")
        
        
def run_batch_pipeline(input_dir: str, output_dir: str, chunk_size: int, skip_reductions: bool, workers: int, 
                       keep_temp_files: bool = False, modality: str = None, batch_key: str = None, verbose: bool = False):
    """
    Batch pipeline, running through all of the steps in parallel.
    """
    PACKAGE_ROOT = os.path.dirname(os.path.abspath(__file__))
    CELLTYPE_MODEL_DIR = os.path.join(PACKAGE_ROOT, "celltype", "models")
    CELLAGE_MODEL_DIR = os.path.join(PACKAGE_ROOT, "cellage", "models")
    INTERNAL_REF_FEATURES = os.path.join(PACKAGE_ROOT, "celltype", "models", "model_features.tsv")
   
    # 1. Discover Inputs
    raw_discovered = discover_datasets(input_dir)
    if not raw_discovered:
      raise FileNotFoundError(f"No valid .h5ad files or 10X matrices found in {input_dir}")

    logger.info("\n🔎 Input directory datasets:")
    logger.info("-" * 65)
    for sample_id, info in raw_discovered.items():
        short_path = os.path.basename(info["path"]) if info["path"] else "N/A"
        logger.info(f"  • ID: {sample_id:<25} | Type: {info['type']:<8} | Location: {short_path}")
    logger.info("-" * 65 + "\n")
    
    # 2. Stage Data for Processing
    discovered = pre_split_datasets_by_batch(raw_discovered, batch_key, output_dir)
    
    logger.info(f"CMageR: Utilizing up to {workers} workers...")
    total_expected_chunks = 0
    
    # 📋 TRACKER DEFINED: Collect indices in original layout structure
    original_barcode_order = []
    
    for sample_id, info in discovered.items():
        prefix = f"{sample_id}_"
        if info["type"] == "10x":
            total_expected_chunks += 1
            # Record 10x cell names layout sequences
            adata_header = sc.read_10x_mtx(info["path"], var_names='gene_symbols', cache=False)
            # Prefix sample_id to guarantee global uniqueness
            prefixed_bcs = [
                bc if bc.startswith(prefix) else f"{prefix}{bc}"
                for bc in adata_header.obs_names.astype(str)
            ]
            original_barcode_order.extend(prefixed_bcs)
        else:
            adata_header = sc.read_h5ad(info["path"], backed="r")
            total_cells = adata_header.shape[0]
            chunks_for_sample = math.ceil(total_cells / chunk_size)
            total_expected_chunks += chunks_for_sample
            # Prefix sample_id to guarantee global uniqueness
            prefixed_bcs = [
                bc if bc.startswith(prefix) else f"{prefix}{bc}"
                for bc in adata_header.obs_names.astype(str)
            ]
            original_barcode_order.extend(prefixed_bcs)

    logger.debug(f"📊 Total workload calculated: {total_expected_chunks} streaming chunks across {len(discovered)} datasets.")
    
    manager = multiprocessing.Manager()
    progress_queue = manager.Queue()

    processed_file_paths = []
    
    # 3. Parallel Execution Engine
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                process_single_sample,
                sample_id=sample_id,
                info=info,
                internal_ref_features=INTERNAL_REF_FEATURES,
                celltype_model_dir=CELLTYPE_MODEL_DIR,
                cellage_model_dir=CELLAGE_MODEL_DIR,
                output_dir=output_dir,
                chunk_size=chunk_size,
                modality=modality,
                skip_reductions=skip_reductions,
                keep_temp_files=keep_temp_files,
                verbose=verbose,
                progress_queue=progress_queue
            ): sample_id for sample_id, info in discovered.items()
        }

        completed_futures = set()

        with tqdm(total=total_expected_chunks, desc="🧬 Pipeline Progress", unit="chunk", leave=True) as pbar:
            while len(completed_futures) < len(futures):
                
                # Drain the queue for chunk progress ticks completely
                while not progress_queue.empty():
                    try:
                        msg = progress_queue.get_nowait()
                        if msg == "chunk_complete":
                            pbar.update(1)
                    except Exception:
                        break
                
                # Check the status of each running dataset job
                for future, sample_id in futures.items():
                    if future.done() and future not in completed_futures:
                        try:
                            saved_path = future.result()
                            processed_file_paths.append(saved_path)
                            completed_futures.add(future)
                        except Exception as e:
                            logger.info(f"🚨 Master Pipeline: Aborting aggregation phase due to failure in worker tracking [{sample_id}]")
                            raise e
                
                time.sleep(0.1)

    # 4. Master Concatenation Step (Passing barcode array along)
    logger.debug("\n📦 All workers finished. Initiating Master Pool aggregation...")
    concatenate_and_save(processed_file_paths, output_dir, original_barcode_order=original_barcode_order)
    
    # 5. Pipeline Cleanup
    if keep_temp_files:
      logger.debug("\n📌 [Troubleshooting Mode] Preserving all split batch inputs, chunks, and CSVs in 'temp_files'.")
    else:
      logger.debug("\n🧹 Cleaning up temporary batch files...")
      for sample_id, info in discovered.items():
          if info.get("is_temporary") and os.path.exists(info["path"]):
              os.remove(info["path"])
              
      temp_dir_path = os.path.join(output_dir, "temp_files")
      if os.path.exists(temp_dir_path) and not os.listdir(temp_dir_path):
          os.rmdir(temp_dir_path)
          logger.debug("🗑️  Successfully removed empty 'temp_files' workspace directory.")                      
                      
    logger.info("🎉 Pipeline successfully executed and fully aggregated!")

