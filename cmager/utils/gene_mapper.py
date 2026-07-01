"""
HGNC Gene Feature Mapping Module.

Handles cross-dataset and cross-model genomic feature name drifts over the years by mapping 
query features to stable Ensembl IDs and canonical symbols using the HGNC database.
Not perfect... but does the job without a full realignment!
"""

import os
import logging
from typing import Optional
import pandas as pd
import anndata as ad

logger = logging.getLogger("cmager")


class HGNCManifoldMapper:
    """
    Manages genomic feature name drifts across datasets and models by 
    resolving modern/historical symbols and mapping them to stable Ensembl IDs.
    """

    def __init__(self, cache_dir: str = "./cache"):
        """
        Initializes the mapper and ensures lookup indexes are ready.
        
        Args:
            cache_dir (str): Folder where the HGNC reference set is saved.
        """
        self.cache_dir = cache_dir
        os.makedirs(cache_dir, exist_ok=True)
        self.db_path = os.path.join(cache_dir, "hgnc_complete_set.csv")
        
        self._load_hgnc_manifold()
        self._build_lookup_indexes()

    def _load_hgnc_manifold(self) -> None:
        """Downloads the authoritative HGNC database if missing from the local cache."""
        if not os.path.exists(self.db_path):
            logger.info("🌐 Downloading authoritative HGNC mapping database...")
            url = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
            df = pd.read_csv(url, sep="\t", low_memory=False)
            df.to_csv(self.db_path, index=False)
            self.df_hgnc = df
        else:
            self.df_hgnc = pd.read_csv(self.db_path, low_memory=False)

    def _build_lookup_indexes(self) -> None:
        """Indexes symbols, previous names, and aliases to their official Ensembl IDs."""
        self.symbol_to_current = {}
        self.symbol_to_ensembl = {}  # Maps historical/current symbols to stable Ensembl IDs
        
        for _, row in self.df_hgnc.iterrows():
            current_symbol = str(row['symbol']).strip()
            ensembl_id = str(row['ensembl_gene_id']).strip() if not pd.isna(row['ensembl_gene_id']) else None
            
            if pd.isna(row['symbol']) or current_symbol.lower() == "nan":
                continue
            
            # Base mappings for modern symbols
            current_lower = current_symbol.lower()
            self.symbol_to_current[current_lower] = current_symbol
            if ensembl_id:
                self.symbol_to_ensembl[current_lower] = ensembl_id
            
            # Map historical symbols to current symbols AND Ensembl IDs
            if not pd.isna(row['prev_symbol']):
                for prev in str(row['prev_symbol']).split("|"):
                    prev_clean = prev.strip().lower()
                    if prev_clean:
                        self.symbol_to_current[prev_clean] = current_symbol
                        if ensembl_id:
                            self.symbol_to_ensembl[prev_clean] = ensembl_id
                        
            # Map alias symbols to current symbols AND Ensembl IDs
            if not pd.isna(row['alias_symbol']):
                for alias in str(row['alias_symbol']).split("|"):
                    alias_clean = alias.strip().lower()
                    if alias_clean:
                        if alias_clean not in self.symbol_to_current:
                            self.symbol_to_current[alias_clean] = current_symbol
                        if ensembl_id and alias_clean not in self.symbol_to_ensembl:
                            self.symbol_to_ensembl[alias_clean] = ensembl_id

    def to_modern(self, gene: str) -> str:
        """Resolves an outdated gene symbol to its modern approved version."""
        return self.symbol_to_current.get(str(gene).strip().lower(), gene)

    def get_ensembl_via_manifold(self, gene: str) -> Optional[str]:
        """Retrieves the Ensembl ID mapped to a specific gene symbol alias."""
        return self.symbol_to_ensembl.get(str(gene).strip().lower(), None)

    def align_by_ensembl_and_manifold(
        self, 
        adata: ad.AnnData, 
        model_features_df: pd.DataFrame, 
        output_csv_path: Optional[str] = None
    ) -> ad.AnnData:
        """
        Executes hierarchical feature alignment matching:
        Query Symbol -> HGNC Manifold -> Ensembl ID Lookup -> Reference Ensembl Match -> Reference Symbol.
        
        Args:
            adata (ad.AnnData): The input dataset containing features to align.
            model_features_df (pd.DataFrame): Model reference with Ensembl (col 0) and Symbols (col 1).
            output_csv_path (str, optional): Target path to write alignment verification logs.
            
        Returns:
            ad.AnnData: A cleanly aligned AnnData copy containing model-compatible var_names.
        """
        logger.debug(f"🧬 Aligning query symbols ({len(adata.var_names)} features) using HGNC-EnsemblID manifold...")
        
        # --- PHASE 1: INITIALIZATION & TRACKING COPIES ---
        adata_aligned = adata.copy()
        adata_aligned.var["old_names"] = adata.var_names.tolist()
        
        final_query_names = []
        audit_log = []
        
        # --- PHASE 2: MODEL REFERENCE INDEXING ---
        ensembl_col = model_features_df.columns[0]
        symbol_col = model_features_df.columns[1]
        
        model_by_ensembl = {}
        model_modern_symbols = {}
        
        for _, row in model_features_df.dropna(subset=[ensembl_col, symbol_col]).iterrows():
            ens_id = str(row[ensembl_col]).strip()
            sym_id = str(row[symbol_col]).strip()
            model_by_ensembl[ens_id] = sym_id
            model_modern_symbols[self.to_modern(sym_id)] = sym_id

        # --- PHASE 3: MULTI-TIERED ALIGNMENT LOOP ---
        # Safeguard tracking for native query Ensembl IDs if they are pre-populated in .var
        query_ensembl_ids = adata.var.get('gene_ids', [None] * adata.n_vars)
        
        for q_symbol, q_native_ens in zip(adata.var_names, query_ensembl_ids):
            q_native_ens_clean = str(q_native_ens).strip() if pd.notna(q_native_ens) else None
            
            # Resolution Step: Infer Ensembl ID via manifold if native tracking is absent
            discovered_query_ensembl = q_native_ens_clean if q_native_ens_clean else self.get_ensembl_via_manifold(q_symbol)
            q_modern = self.to_modern(q_symbol)
            
            target_symbol = q_symbol
            mapping_method = "No Match (Retained Query)"
            
            # Tier 1 & 2: Structural Match via Ensembl ID (Native or Manifold-Discovered)
            if discovered_query_ensembl and discovered_query_ensembl in model_by_ensembl:
                target_symbol = model_by_ensembl[discovered_query_ensembl]
                mapping_method = "Provided EnsemblID Match" if q_native_ens_clean else "HGNC_EnsemblID Manifold Match"
            
            # Tier 3: Fallback Match via Canonical Symbol Intersection
            elif q_modern in model_modern_symbols:
                target_symbol = model_modern_symbols[q_modern]
                mapping_method = "Manifold Symbol Alias Match"
            
            final_query_names.append(target_symbol)
            is_remapped = (q_symbol != target_symbol)
            
            audit_log.append({
                "query_feature": q_symbol,
                "assigned_model_feature": target_symbol,
                "is_remapped": is_remapped,               
                "alignment_strategy": mapping_method,
                "query_native_ensembl": q_native_ens_clean if q_native_ens_clean else "None Provided",
                "resolved_ensembl_used": discovered_query_ensembl if discovered_query_ensembl else "NA",
                "manifold_modern_feature": q_modern,
            })
                    
        # --- PHASE 4: AUDIT LOGGING & INJECTION ---
        if output_csv_path:
            pd.DataFrame(audit_log).to_csv(output_csv_path, index=False)
            logger.debug(f"📊 Mapping audit table saved to: {output_csv_path}")
        
        adata_aligned.var_names = final_query_names
        adata_aligned.var["is_remapped"] = [row["is_remapped"] for row in audit_log]
        
        # Suffix collision handler for duplicate genes assigned to identical references
        if not adata_aligned.var_names.is_unique:
            logger.debug("⚠️ Resolving duplicate feature collisions...")
            adata_aligned.var_names_make_unique()
            
        return adata_aligned
