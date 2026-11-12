

[![Python](https://img.shields.io/badge/python-3.13-0dc7ed.svg?labelColor=000000&logo=python)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.4-blue&logo=R)](https://www.r-project.org/)

![CMageR_logo slim](assets/logoLong_v2_light.png)

> **CMageR** | *Cardiomyocyte (CM)* *ager*
>
> What `cardiomyocytes` are you making? What `age` or `maturity` are they transcriptionally? If you have scRNA-seq data, 
> CMageR is the first step towards a computational reference to help you find out the answer to these questions!

* Trained on human fetal single-cell transcriptomics, CMageR is a predictive framework for cardiac cell types and cardiac cell ages of cells in scRNA-seq data.
* It is intended for discrete classification of cell types within culture (in vitro, or in vivo), and their ages within **4 to 15 post-conceptional weeks (PCW)**.


The pipeline automates:
1. **Cell Type Annotation:** Labels cells using a trained fetal cardiac cell type annotation step comprised of three granularities of CellTypist models followed by a semi-hierarchical classification re-annotation layer.
2. **Cell Age Prediction:** Estimates post-conception week age of individual cells using cell-type specific multivariate generalized additive models each trained using the R package gamsel.

---

## Installation & Setup

Because this tool uses both R and Python, it requires a custom environment. It was build and tested in conda in Linux (Ubuntu 24.04) using the supplied environment file. Firstly, rebuild the working environment using the provided conda `environment.yml` file. This handles all specific installations, including Python with the packages Scanpy, AnnData, CellTypist, and others utilities.

### Prerequisites
Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) installed on your machine.

### Step-by-Step Installation using conda and pip
Installation should only take a few minutes for all steps to be complete. 

1. **Create the Environment from the YAML File**
   This command reads the `environment.yml` configuration from this repository to resolve cross-dependencies and build the environment called `cmager`, containing Python and R. This step takes a minute or two while conda solves the environment:
   ```bash
   conda env create -f https://raw.githubusercontent.com/knische/cmager/main/environment.yml
   ```
   
2. **Activate the conda environment:**
   ```bash
   conda activate cmager
   ```
      
3. **Use pip to install from github:**
   If you have git:
   ```bash
   pip install git+https://github.com/knische/cmager.git
   ```
   or, in all cases:
   ```bash
   pip install https://github.com/knische/cmager/archive/refs/heads/main.zip
   ```

4. **Confirm that the package is installed:**
   Check out the usage and help page. This step may take a moment the first time:
   ```bash
   cmager --help
   ```



## Usage
1. **quick usage**
   Typical CLI usage where input directory contains either multiple 10X matrix directories (matrix.mtx, features.tsv, barcodes.tsv), or multiple anndata.h5ad objects.
   ```bash
   cmager -i "input_directory" -o "output_directory" --modality "nuclei"
   ```

2. **Advanced example**
   Example where a single anndata.h5ad file is input, and its metadata .obs contains "donor", and "modality" for batch, and modality arguments respectively.
   ```bash
   cmager -i "input_directory/all_data" -o "output_directory" --modality "modality" --batch "donor" -workers 12 --verbose --skip-reductions
   ```

## Inputs and Parameters

### Supported Input Formats

`cmager` automatically detects and ingests any of the following directory structures passed to `--input-dir`:

```text
# Option 1: Standard 10X Genomics CellRanger directories
input_directory/
├── sample_id_01/
│   ├── matrix.mtx.gz
│   ├── barcodes.tsv.gz
│   └── features.tsv.gz
└── sample_id_02/
    ├── matrix.mtx.gz
    ├── barcodes.tsv.gz
    └── features.tsv.gz

# Option 2: Individual AnnData (.h5ad) files
input_directory/
├── sample_id_01.h5ad
└── sample_id_02.h5ad

# Option 3: Single merged AnnData (.h5ad) file
input_directory/
└── all_data/
    └── all_data.h5ad
```

### CLI Options

| Option | Requirement / Default | Description |
| :--- | :--- | :--- |
| `-i, --input-dir DIRECTORY` | **Required** | Directory containing `.h5ad` files or 10X CellRanger output folders. |
| `-o, --output-dir DIRECTORY` | **Required** | Directory to save individual sample results and the final aggregated dataset. |
| `-m, --modality TEXT` | **Required** | Dataset modality (`cells` or `nuclei`), or the `.obs` column containing this metadata. |
| `-c, --chunk-size INTEGER` | Default: `50000` | Number of cells to process per streaming chunk. |
| `-w, --workers INTEGER` | Default: `1` | Number of parallel worker processes to spawn. |
| `-b, --batch TEXT` | Optional | Column name in `.obs` to split large datasets by for optimized parallel processing. |
| `--skip-reductions` | Flag (`False`) | Skips PCA and UMAP dimensional reductions to save memory. |
| `--keep-temp-files` | Flag (`False`) | Preserves intermediate processing files in `output_dir/temp_files` for debugging. |
| `-v, --verbose` | Flag (`False`) | Prints detailed status updates across all workers in the terminal. |
| `--version` | Flag | Show the version and exit. |
| `-h, --help` | Flag | Show this help message and exit. |
	

## Results
The output directory will contain a concatenated anndata.h5ad object as well as a separate results.csv file.
Additionally, per batch directories are created each with its respective results.

```text
output_directory/
├── cmager_results_anndata.h5ad         # Anndata object
├── cmager_results_metadata.csv         # Complete pooled annotation flat table 
├── sample_id_01/                       # Neatly organized sample-specific directory
│   ├── sample_id_01_annotated.h5ad     # Individual processed AnnData matrix
│   ├── sample_id_01_metadata_table.csv # Sample specific metrics
│   └── sample_id_01_gene_mapping_audit.csv # Mapping table of the genes mapped back to the model
└── sample_id_02/
    ├── sample_id_02_annotated.h5ad
    ├── sample_id_02_metadata_table.csv
    └── sample_id_02_gene_mapping_audit.csv
```

