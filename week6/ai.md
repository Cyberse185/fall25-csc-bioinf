# AI Usage Documentation (Week 6)

Claude Sonnet 4.5 from Antrhopic was used throughout this assignment for code, debugging, installations. 

### 1. Environment and Setup

* **Conda Installation:** Provided the initial `conda` commands for `simpleaf` and helped troubleshoot an `UnsatisfiableError` on Apple Silicon.
* **Jupyter Scoping:** Solved the primary execution issue (`!simpleaf command not found`) by providing a Python helper function (`run_in_conda`) to correctly run shell commands from within the notebook's `conda` environment.

### 2. Pipeline Debugging

* **`alevin-fry` / `simpleaf`:**
    * Helped debug `FileNotFoundError` by providing `ls`/`find` commands to locate the correct paths for the reference genome and FASTQ files.
    * Corrected the `simpleaf quant` command to handle the non-gzipped (`.fastq`) input files.
* **`scanpy` Preprocessing:**
    * Resolved a `ValueError` (0 samples) by suggesting more lenient filtering parameters (`min_genes=5`) suitable for the small toy dataset.
    * Fixed an `InvalidParameterError` in `sc.tl.pca` by providing code to dynamically set `n_comps` based on the data's actual dimensions.
* **`celltypist` Annotation:**
    * **Gene Name Conversion (Critical Fix):** Diagnosed a `ValueError` from CellTypist ("No features overlap"). The AI identified that my data used Ensembl IDs while the model required gene symbols. It then **generated the Python code to parse the `.gtf` file**, create a mapping dictionary, and update `adata.var_names` to use gene symbols.
    * **Data Normalization:** Fixed a `ValueError` (invalid matrix) by providing the correct workflow: saving the log-normalized data to `adata.raw = adata.copy()` *before* scaling the data with `sc.pp.scale`.