import os
import pytest
from pathlib import Path
from cmager.pipeline import discover_datasets, run_batch_pipeline

TEST_DIR = Path(__file__).parent
DATA_DIR = TEST_DIR / "data"

def test_10x_dataset_discovery():
    """Test if the pipeline correctly identifies 10X folders."""
    tenx_path = DATA_DIR / "10x"
    
    # Assuming discover_datasets returns a list of found valid paths/files
    discovered = discover_datasets(input_dir=str(tenx_path))
    
    # Assert that it actually found the expected data
    assert len(discovered) > 0, "Pipeline failed to discover 10X data!"

def test_h5ad_processing_end_to_end(tmp_path):
    """Test a small run using a temporary output directory."""
    input_h5ad = DATA_DIR / "h5ad_combined"
    
    # tmp_path is a built-in pytest feature that creates a safe, temporary folder
    run_batch_pipeline(
        input_dir=str(input_h5ad),
        output_dir=str(tmp_path),
        chunk_size=100,  # Keep it tiny for fast testing
        low_ram=True,
        workers=6
    )
    
    # Check if the expected output files were successfully created
    assert (tmp_path / "cmager_results_anndata.h5ad").exists()
