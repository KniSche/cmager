

[![Python](https://img.shields.io/badge/python-3.13-0dc7ed.svg?labelColor=000000&logo=python)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4.4db7ed.svg?labelColor=000000&logo=R)](https://www.r-project.org/)

![CMageR_logo slim](assets/logoLong_v2.png)

# CMageR: Fetal Cardiac Single-Cell Annotation & Age Prediction Pipeline

> **CMageR** | *Cardiomyocytes (CM)*, *ager*
>
> What `cardiomyocytes` are you making? What `age` or `maturity` are they transcriptionally? 
> CMageR is the first step towards a computational reference to help you find out the answer to these questions!

* Trained on human fetal single-cell transcriptomics, CMageR is a predictive framework for cardiac cell types and cardiac cell ages of cells in scRNA-seq data.
* It is intended for discrete classification of cell types within culture (in vitro, or in vivo), and their ages within **4 to 15 post-conceptional weeks (PCW)**.
* It uses a custom-built pipeline in python combining logistic regression (CellTypist), and multivariate generalised additive models in R (gamsel)


The pipeline automates:
1. **Cell Type Annotation:** Employs a semi-hierarchical classification layer built on top of CellTypist using carefully curated developmental human heart references.
2. **Cell Age Prediction:** Estimates continuous developmental coordinates using an optimized multivariate generalized additive model (GAM).

---

## 🚀 Installation & Setup

New users can reproduce your exact working environment using the provided Anaconda `environment.yml` file. This handles all specific installations, including Scanpy, AnnData, CellTypist, and tracking utilities.

### Prerequisites
Ensure you have [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) installed on your machine.

### Step-by-Step Installation

1. **Clone the Repository**
   ```bash
   git clone [https://github.com/yourusername/cmager.git](https://github.com/yourusername/cmager.git)
   cd cmager
   ```

2. **Create the Environment from the YAML File**
   This command reads your `environment.yml` configuration, resolves all cross-dependencies, and builds an isolated workspace named `cmager`:
   ```bash
   conda env create -f environment.yml
   ```

3. **Activate the Environment**
   ```bash
   conda activate cmager
   ```

4. **Verify the Installation**
   ```bash
   python -c "import cmager; print('CMageR environment ready!')"
   ```

---

## 🛠️ Pipeline Architecture & Key Functions

The `pipeline.py` script acts as the master orchestration layer. It is explicitly engineered for **low-RAM footprints**, leveraging parallel multiprocessing worker pools alongside streaming chunk configurations.

### Core Functions

* **`discover_datasets(input_dir)`**
  Scans the targeted input directory. Dynamically pairs and resolves raw inputs, natively supporting standalone `.h5ad` matrices and standard 10X Genomics matrix folders containing `matrix.mtx`, `features.tsv`, and `barcodes.tsv`.
* **`pre_split_datasets_by_batch(discovered, batch_key, output_dir)`**
  An optional performance step. If a metadata categorical column (e.g., `donor` or `orig_ident`) is provided via the `batch_key`, this function pre-slices complex multi-sample files into temporary, isolated files to lower peak execution RAM.
* **`standardise_gene_symbols(adata, ref_features_df)`**
  Scans matrix var features against internal model genes. It automatically shifts the index to alternate metadata mapping columns if a superior feature intersection configuration is detected.
* **`process_single_sample(...)`**
  The primary execution worker thread. Spawns chunked streams through custom `HGNCManifoldMapper` layers, calculates underlying CellTypist classifications, and applies R-backed hierarchical consensus modeling alongside multivariate age estimations.
* **`concatenate_and_save(processed_file_paths, output_dir, original_barcode_order)`**
  The master aggregator loop. Gathers distinct worker outputs, concatenates the single-cell expressions sequentially, **restores exact original sequence row indexing layout sequences** via `original_barcode_order`, and cleans up individual runs into organized sample folders.

---

## 📊 Directory Structure Workflow

When a user runs the batch pipeline, the framework processes assets across a temporary workspace and automatically houses variables cleanly upon execution:

```text
output_directory/
├── cmager_results_anndata.h5ad         # Unified master matrix (Preserves original input ordering)
├── cmager_results_metadata.csv         # Complete pooled annotation flat table
├── sample_id_01/                       # Neatly organized sample-specific directory
│   ├── sample_id_01_annotated.h5ad     # Individual processed AnnData matrix
│   ├── sample_id_01_metadata_table.csv # Sample specific metrics
│   └── sample_id_01_gene_mapping_audit.csv
└── sample_id_02/
    ├── sample_id_02_annotated.h5ad
    └── sample_id_02_metadata_table.csv
```

---

## 💡 Quick Usage Example

The pipeline can be triggered directly via your custom execution CLI entry point script:

```python
from cmager.pipeline import run_batch_pipeline

run_batch_pipeline(
    input_dir="/path/to/raw_data",
    output_dir="/path/to/results",
    chunk_size=5000,       # Process matrix rows in safe memory frames
    low_ram=True,
    workers=4,             # Dictates ProcessPoolExecutor execution scale
    batch_key="donor"      # Optional pre-splitting key variable
)
```
