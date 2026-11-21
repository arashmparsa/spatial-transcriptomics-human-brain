"""IO utilities: load 10x Visium / spatialLIBD data into an AnnData object.

Supported inputs
----------------
1. A Space Ranger output directory (contains ``filtered_feature_bc_matrix.h5``
   and a ``spatial/`` subfolder) -> ``sc.read_visium``.
2. A filtered feature-barcode matrix ``*.h5`` -> ``sc.read_10x_h5``.
3. A spatialLIBD DLPFC sample id (e.g. ``151673``) resolved under ``data_dir``.
4. A 10x-hosted Visium section id -> ``sc.datasets.visium_sge`` (fallback).

All routes return an AnnData with unique var names and (where available)
``adata.obsm['spatial']`` populated so downstream Squidpy / Scanpy spatial
functions work.
"""

from __future__ import annotations

import os
from pathlib import Path

import scanpy as sc
from anndata import AnnData


def _looks_like_spaceranger_dir(path: Path) -> bool:
    """True if ``path`` looks like a Space Ranger output directory."""
    if not path.is_dir():
        return False
    has_spatial = (path / "spatial").is_dir()
    has_matrix = (
        (path / "filtered_feature_bc_matrix.h5").exists()
        or (path / "raw_feature_bc_matrix.h5").exists()
        or (path / "filtered_feature_bc_matrix").is_dir()
    )
    return has_spatial and has_matrix


def load_visium(sample_or_path: str, count_file: str = "filtered_feature_bc_matrix.h5") -> AnnData:
    """Load a Visium dataset from a directory, an ``.h5`` file, or a 10x SGE id.

    Parameters
    ----------
    sample_or_path
        One of:
          * a path to a Space Ranger output directory,
          * a path to a ``*.h5`` filtered feature-barcode matrix,
          * a 10x ``visium_sge`` sample id (e.g. ``"V1_Human_Brain_Section_1"``).
    count_file
        Name of the counts ``.h5`` inside a Space Ranger directory.

    Returns
    -------
    AnnData
        With unique var names; ``obsm['spatial']`` present for Visium inputs.
    """
    path = Path(sample_or_path)

    # 1. Space Ranger output directory.
    if _looks_like_spaceranger_dir(path):
        adata = sc.read_visium(path, count_file=count_file)

    # 2. A standalone .h5 filtered matrix (no spatial coords).
    elif path.suffix == ".h5" and path.exists():
        adata = sc.read_10x_h5(path)

    # 3. An .h5ad already prepared.
    elif path.suffix == ".h5ad" and path.exists():
        adata = sc.read_h5ad(path)

    # 4. Fall back to a 10x-hosted Visium section id.
    else:
        adata = sc.datasets.visium_sge(sample_id=str(sample_or_path))

    adata.var_names_make_unique()
    return adata


def resolve_data_source(cfg: dict) -> str:
    """Decide which dataset to load from the config, preferring local data.

    Resolution order:
      1. ``data.explicit_path`` if set.
      2. A spatialLIBD Space Ranger dir at ``data_dir/sample_id``.
      3. Any Space Ranger dir / ``.h5`` auto-detected under ``local_data_dir``.
      4. The 10x ``visium_sge`` fallback id (guaranteed runnable).

    Returns the string to hand to :func:`load_visium`.
    """
    data = cfg["data"]

    explicit = data.get("explicit_path")
    if explicit:
        return explicit

    # spatialLIBD sample staged under data_dir (the titular dataset).
    sample_dir = Path(data["data_dir"]) / str(data["sample_id"])
    if _looks_like_spaceranger_dir(sample_dir):
        return str(sample_dir)

    # Auto-detect pre-existing local 10x data under an optional local_data_dir.
    local = data.get("local_data_dir")
    if local:
        found = find_local_10x(local)
        if found is not None:
            return found

    # Guaranteed-runnable 10x fallback section.
    return data["visium_sge_fallback"]


def find_local_10x(root: str) -> str | None:
    """Search ``root`` for a Space Ranger directory or filtered ``.h5`` matrix.

    Returns the first match (Space Ranger dir preferred), else ``None``.
    """
    root_path = Path(root)
    if not root_path.exists():
        return None

    # Prefer full Space Ranger output directories.
    for dirpath, _dirnames, _filenames in os.walk(root_path):
        if _looks_like_spaceranger_dir(Path(dirpath)):
            return dirpath

    # Otherwise a filtered feature-barcode matrix .h5.
    for h5 in root_path.rglob("*filtered_feature_bc_matrix.h5"):
        return str(h5)
    for h5 in root_path.rglob("*.h5"):
        return str(h5)

    return None


def attach_layer_annotations(adata: AnnData, cfg: dict) -> AnnData:
    """Add manual spatialLIBD cortical-layer labels as ``adata.obs['layer']``.

    Reads ``data_dir/<sample_id>/layer_annotations.tsv`` (written by
    ``scripts/download_data.py``). Spots without an annotation get NaN.
    No-op if the file is absent.
    """
    import pandas as pd

    path = Path(cfg["data"]["data_dir"]) / str(cfg["data"]["sample_id"]) / "layer_annotations.tsv"
    if not path.exists():
        return adata
    ann = pd.read_csv(path, sep="	", index_col="barcode")
    adata.obs["layer"] = pd.Categorical(
        adata.obs_names.map(ann["layer"]),
        categories=["L1", "L2", "L3", "L4", "L5", "L6", "WM"],
    )
    n = int(adata.obs["layer"].notna().sum())
    print(f"[io] attached layer annotations for {n}/{adata.n_obs} spots")
    return adata
