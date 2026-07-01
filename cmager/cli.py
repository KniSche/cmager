"""
Command Line Interface (CLI) for CMageR.
Handles user inputs, validates arguments, and triggers the master orchestration pipeline.
"""

import os
import click
from cmager.pipeline import run_batch_pipeline

@click.command()
@click.version_option(version="0.1.0")
@click.option(
    "--input-dir", "-i",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    required=True,
    help="Directory containing .h5ad files or 10X CellRanger outputs."
)
@click.option(
    "--output-dir", "-o",
    type=click.Path(dir_okay=True, file_okay=False, writable=True),
    required=True,
    help="Directory to save individual sample results and the final aggregated dataset."
)
@click.option(
    "--chunk-size", "-c",
    type=int,
    default=50000,
    show_default=True,
    help="Number of cells to process per chunk stream."
)
@click.option(
    "--workers", "-w",
    type=int,
    default=1,
    show_default=True,
    help="Number of parallel worker processes to spawn."
)
@click.option(
    "--modality", "-m",
    type=str, 
    default=None,
    required=True,
    help="Dataset modality ('cells' or 'nuclei'), or the .obs column containing this metadata."
)
@click.option(
    "--batch", "-b",
    type=str, 
    default=None, 
    help="Column name in .obs to split large datasets by for optimized parallel processing."
)
@click.option(
    "--low-ram", 
    is_flag=True, 
    help="Skip computationally heavy UMAP and PCA dimensional reductions to save memory."
)
@click.option(
    "--keep-temp-files",
    is_flag=True,
    help="Preserve intermediate processing files in 'output_dir/temp_files' for debugging."
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Enable detailed execution logs and chunk-level tracking in the terminal."
)
def main(input_dir: str, output_dir: str, chunk_size: int, workers: int, 
         modality: str, batch: str, low_ram: bool, keep_temp_files: bool, verbose: bool):
    """
    CMageR: Cell Age prediction for 4-15 PCW in cardiac cells. 
    
    1. Performs cell type annotation through semi-hierarchical classification with CellTypist
    2. Performs cell age prediction through multivariate generalised additive model
    
    In both cases, data were carefully curated from fetal heart single cell and single nuclei RNA sequencing.
    """
    
    # 1. Prepare the execution environment
    os.makedirs(output_dir, exist_ok=True)
    
    # 2. Hand off control to the master orchestrator
    run_batch_pipeline(
        input_dir=input_dir, 
        output_dir=output_dir, 
        chunk_size=chunk_size,
        workers=workers,
        low_ram=low_ram,
        modality=modality,
        batch_key=batch,
        keep_temp_files=keep_temp_files,
        verbose=verbose
    )

if __name__ == "__main__":
    main()
