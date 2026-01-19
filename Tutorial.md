## Project Overview: Single-Cell Analysis of 3D Bioprinted Tissues

This project investigates cellular transitions from 2D cell lines (Day 0) to 3D bioprinted tissues. The dataset comprises ~1.6 million cells processed via **Parse Biosciences Evercode** combinatorial barcoding.

### 1. System Requirements & Setup

* **R Version:** 4.5+
* **Key Libraries:** `Seurat` (v5), `vroom` (fast loading), `anndata` (Python interoperability).
* **Hardware:** Minimum 32GB RAM recommended for individual sample processing; 64GB+ for merged analysis.

### 2. Data Acquisition

Raw data is transferred from the sequencing core using `rsync` to ensure file integrity and the ability to resume interrupted transfers:

```bash
rsync -avzP [source_path] ./data/raw/

```

### 3. Version Control (Git)

We use Git to track code logic, excluding large data and image files to keep the repository lightweight.

* **Ignored Folders:** `data/`, `figures/`, `.Rproj.user/`
* **Workflow:** `git add [script]`, `git commit -m "description"`, `git push origin main`.

### 4. Quality Control (QC) Pipeline

The data is processed through a custom R script (`scripts/01_qc_and_filtering.R`) that accounts for the unique properties of Parse Biosciences data:

#### A. Matrix Transposition

Parse output is oriented as `Cells x Genes`. Seurat requires `Genes x Cells`. The script automatically transposes the matrix upon loading.

#### B. Adaptive Filtering

We avoid "one-size-fits-all" filters to prevent the loss of smaller immune cell populations.

* **Immune Cell Preservation:** We maintain a low `nFeature_RNA` floor (200 genes) for 3D samples.
* **D0 vs 3D:** Day 0 cell lines have a significantly higher RNA ceiling (150k UMIs) compared to 3D tissue samples (30k UMIs).

#### C. Mitochondrial DNA (Death Marker)

We use **Ridge Plots** to visualize the `percent.mt` across sub-libraries. A peak shifting to the right indicates high cellular stress or death in a specific well-pool.

### 5. Interoperability

Filtered data is exported as **`.h5ad`** files. This allows the data to be opened seamlessly in Python using `Scanpy`:

```python
import scanpy as sc
adata = sc.read_h5ad("sample_filtered.h5ad")

```

---


### Your Next Step

Now that your documentation is ready and your data is filtered, the next logical step is to **Merge the samples** and perform **Dimensionality Reduction (PCA/UMAP)**.

**Would you like me to provide the "Sketching" code for Seurat v5 to handle the 1.6 million cells in the next step of your analysis?**