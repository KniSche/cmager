"""
Command Line Interface (CLI) for CMageR.
"""

import os
import click
from cmager.pipeline import run_batch_pipeline

SPLASH_SCREEN = r"""
 ___________________________________________________________________________
 
  ###  #   # |
 #     ##### |   CMageR:
 #     # # # |   Cell Type and Cell Age prediction
  ###  #   # |   for cardiac cells between 4-15 PCW.
 ----------- |   1. Cell type annotation: 
             |      semi-hierarchical classification with CellTypist
 ----------- |   
          0  |   2. Cell age prediction:
 ----------- |      multivariate generalised additive model
             |   
 -----0----- |   Training data were carefully curated from fetal heart data
             |   single cell and single nuclei RNA sequencing
 ---0------- |   Publications: 
 ___________________________________________________________________________
"""

# Custom Command class to safely force the splash screen to display on --help
class SplashCommand(click.Command):
    def format_help(self, ctx, formatter):
        # 1. Print the splash screen with cyan color
        click.echo(click.style(SPLASH_SCREEN, fg="cyan", bold=True))
        # 2. Print the normal options (--input-dir, etc.) underneath
        super().format_help(ctx, formatter)


# Attach the custom help class directly to your main command
@click.command(cls=SplashCommand, context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version="0.1.3")
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
    "--skip-reductions", 
    is_flag=True, 
    help="Skips UMAP and PCA dimensional reductions to save memory."
)
@click.option(
    "--keep-temp-files",
    is_flag=True,
    help="Preserve intermediate processing files in 'output_dir/temp_files' for debugging."
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Prints status updates (very detailed) across all workers in the terminal."
)


def main(input_dir: str, output_dir: str, chunk_size: int, workers: int, 
         modality: str, batch: str, skip_reductions: bool, keep_temp_files: bool, verbose: bool):
   
    # 1. Prepare the execution environment
    os.makedirs(output_dir, exist_ok=True)
    
    # 2. Hand off control to the master orchestrator
    run_batch_pipeline(
        input_dir=input_dir, 
        output_dir=output_dir, 
        chunk_size=chunk_size,
        workers=workers,
        skip_reductions=skip_reductions,
        modality=modality,
        batch_key=batch,
        keep_temp_files=keep_temp_files,
        verbose=verbose
    )

if __name__ == "__main__":
    main()
