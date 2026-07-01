

[![Python](https://img.shields.io/badge/python-3.13-0dc7ed.svg?labelColor=000000&logo=python)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.3db7ed.svg?labelColor=000000&logo=R)](https://www.r-project.org/)

![CMageR_logo slim](assets/logoLong_v2.png)

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

Because this tool uses both R and Python, it requires a custom environment. It was build and tested in conda using the supplied environment file. Firstly, rebuild the working environment using the provided conda `environment.yml` file. This handles all specific installations, including Scanpy, AnnData, CellTypist, and tracking utilities.

### Prerequisites
Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) installed on your machine.

### Step-by-Step Installation using conda and pip

1. **Create the Environment from the YAML File**
   This command reads the `environment.yml` configuration from this repository to resolves cross-dependencies and build the environment called `cmager`:
   ```bash
   conda env create -f https://raw.githubusercontent.com/knische/cmager/main/environment.yml
   ```
   
2. **Activate the conda environment**
   ```bash
   conda activate cmager
   ```
      
2. **using pip to install from github**
   This command installs the package cmager from this repository:
   ```bash
   pip install git+https://github.com/yourusername/cmager.git
   ```

4. **install the package as a command line interface (CLI)**
   This command reads your `environment.yml` configuration, resolves all cross-dependencies, and builds an isolated workspace named `cmager`:
   ```bash
   pip install -e .
   ```



## Usage

1. **quick usage**
   Typical CLI usage where input directory contains either 10X matrix directories (matrix.mtx, features.tsv, barcodes.tsv), or anndata.h5ad objects.
   ```bash
   cmager -i "input/directory" -o "output/directory" modality "nuclei"
   ```




## Results
When a user runs the batch pipeline, the framework processes assets across a temporary workspace and automatically houses variables cleanly upon execution:

```text
output_directory/
├── cmager_results_anndata.h5ad         # Anndata object
├── cmager_results_metadata.csv         # Complete pooled annotation flat table 
├── sample_id_01/                       # Neatly organized sample-specific directory
│   ├── sample_id_01_annotated.h5ad     # Individual processed AnnData matrix
│   ├── sample_id_01_metadata_table.csv # Sample specific metrics
│   └── sample_id_01_gene_mapping_audit.csv
└── sample_id_02/
    ├── sample_id_02_annotated.h5ad
    └── sample_id_02_metadata_table.csv
```

