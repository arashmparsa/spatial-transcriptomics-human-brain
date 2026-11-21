"""01_walkthrough — a compact, runnable mirror of the pipeline.

Written as a jupytext-style ``# %%`` script (open in Jupyter/VS Code as a
notebook, or run top-to-bottom with ``python notebooks/01_walkthrough.py``).
It exercises the same modules as ``run_pipeline.py`` on whatever dataset
``config.yaml`` resolves to (spatialLIBD DLPFC by default, or the 10x
fallback).
"""

# %% [markdown]
# # Spatial Transcriptomics of Human DLPFC (spatialLIBD) — walkthrough
#
# QC -> normalization -> clustering -> spatially variable genes, using Scanpy,
# Squidpy, and SpatialDE. Deconvolution (cell2location) is demonstrated at the
# end but left commented out because it needs a scRNA-seq reference and is heavy.

# %%
import sys
from pathlib import Path

import scanpy as sc
import yaml

# Make the repo root importable when run from notebooks/.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import cluster as cluster_mod
from src import io_utils
from src import normalize as normalize_mod
from src import qc as qc_mod
from src import spatial_genes as svg_mod

sc.settings.verbosity = 1

# %%
# Load config and resolve the dataset (prefers local data, falls back to 10x).
with open(ROOT / "config.yaml", "r", encoding="utf-8") as fh:
    cfg = yaml.safe_load(fh)

source = io_utils.resolve_data_source(cfg)
print("Loading:", source)
adata = io_utils.load_visium(source)
adata

# %% [markdown]
# ## 1. Quality control
# Mitochondrial fraction, counts, and genes per spot; then filter.

# %%
adata = qc_mod.run_qc(adata, cfg)

# %% [markdown]
# ## 2. Normalization + dimensionality reduction
# normalize_total + log1p, HVGs, scale, PCA, neighbors, UMAP.

# %%
adata = normalize_mod.normalize_and_reduce(adata, cfg)

# %% [markdown]
# ## 3. Leiden clustering + putative cortical-layer annotation

# %%
adata, layer_table = cluster_mod.run_clustering(adata, cfg)
layer_table

# %% [markdown]
# ## 4. Spatially variable genes (Squidpy Moran's I + SpatialDE)

# %%
svg_results = svg_mod.run_spatial_genes(adata, cfg)
svg_results["moran"].head(10)

# %% [markdown]
# ## 5. (Optional) cell2location deconvolution
# Requires a paired snRNA-seq DLPFC reference `.h5ad`; enable in `config.yaml`.
#
# ```python
# from src import deconvolution as deconv_mod
# adata = deconv_mod.run_deconvolution(adata, cfg)
# ```
