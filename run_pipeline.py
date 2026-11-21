"""Orchestrate the spatialLIBD DLPFC Visium pipeline.

Steps (any subset, in order):
    load       -> read the Visium dataset into AnnData
    qc         -> QC metrics, plots, filtering        (writes adata_qc.h5ad)
    normalize  -> normalize, HVGs, PCA, neighbors, UMAP(writes adata_norm.h5ad)
    cluster    -> Leiden + layer annotation           (writes adata_clustered.h5ad)
    svg        -> Moran's I + SpatialDE SVGs           (writes SVG tables/plots)
    deconv     -> cell2location deconvolution          (writes proportions/maps)

Steps are resumable: each reads the AnnData written by the previous stage from
``results/`` if it isn't already in memory.

Examples
--------
    python run_pipeline.py --steps load qc normalize cluster svg
    python run_pipeline.py --steps cluster --resolution 1.0
    python run_pipeline.py --steps deconv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import scanpy as sc
import yaml
from anndata import AnnData

from src import cluster as cluster_mod
from src import deconvolution as deconv_mod
from src import io_utils
from src import normalize as normalize_mod
from src import qc as qc_mod
from src import spatial_genes as svg_mod

STEP_ORDER = ["load", "qc", "normalize", "cluster", "svg", "deconv"]

# Where each stage persists its AnnData.
ARTIFACTS = {
    "qc": "adata_qc.h5ad",
    "normalize": "adata_norm.h5ad",
    "cluster": "adata_clustered.h5ad",
    "svg": "adata_clustered.h5ad",       # svg annotates the clustered object
    "deconv": "adata_deconv.h5ad",
}


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def _results_path(cfg: dict, name: str) -> Path:
    return Path(cfg["paths"]["results_dir"]) / name


def _ensure_loaded(adata: AnnData | None, cfg: dict, prefer: str | None) -> AnnData:
    """Return an in-memory AnnData, loading a prior artifact if needed."""
    if adata is not None:
        return adata
    if prefer is not None:
        path = _results_path(cfg, prefer)
        if path.exists():
            print(f"[pipeline] loading prior artifact: {path}")
            return sc.read_h5ad(path)
    # Fall back to loading raw data from the resolved source.
    source = io_utils.resolve_data_source(cfg)
    print(f"[pipeline] loading dataset from: {source}")
    return io_utils.attach_layer_annotations(io_utils.load_visium(source), cfg)


def run(cfg: dict, steps: list[str]) -> None:
    results_dir = Path(cfg["paths"]["results_dir"])
    results_dir.mkdir(parents=True, exist_ok=True)
    Path(cfg["paths"]["figures_dir"]).mkdir(parents=True, exist_ok=True)

    adata: AnnData | None = None

    if "load" in steps:
        source = io_utils.resolve_data_source(cfg)
        print(f"[pipeline] loading dataset from: {source}")
        adata = io_utils.load_visium(source)
        adata = io_utils.attach_layer_annotations(adata, cfg)
        print(f"[pipeline] loaded {adata.n_obs} spots x {adata.n_vars} genes")

    if "qc" in steps:
        adata = _ensure_loaded(adata, cfg, prefer=None)
        adata = qc_mod.run_qc(adata, cfg)
        adata.write_h5ad(_results_path(cfg, ARTIFACTS["qc"]))

    if "normalize" in steps:
        adata = _ensure_loaded(adata, cfg, prefer=ARTIFACTS["qc"])
        adata = normalize_mod.normalize_and_reduce(adata, cfg)
        adata.write_h5ad(_results_path(cfg, ARTIFACTS["normalize"]))

    if "cluster" in steps:
        adata = _ensure_loaded(adata, cfg, prefer=ARTIFACTS["normalize"])
        adata, layer_table = cluster_mod.run_clustering(adata, cfg)
        layer_table.to_csv(results_dir / "cluster_layer_annotation.csv")
        adata.write_h5ad(_results_path(cfg, ARTIFACTS["cluster"]))

    if "svg" in steps:
        adata = _ensure_loaded(adata, cfg, prefer=ARTIFACTS["cluster"])
        svg_mod.run_spatial_genes(adata, cfg)

    if "deconv" in steps:
        adata = _ensure_loaded(adata, cfg, prefer=ARTIFACTS["cluster"])
        adata = deconv_mod.run_deconvolution(adata, cfg)
        adata.write_h5ad(_results_path(cfg, ARTIFACTS["deconv"]))

    print("[pipeline] done.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="spatialLIBD DLPFC Visium analysis pipeline."
    )
    parser.add_argument("--config", default="config.yaml", help="Path to config.yaml")
    parser.add_argument(
        "--steps",
        nargs="+",
        default=["load", "qc", "normalize", "cluster", "svg"],
        choices=STEP_ORDER,
        metavar="STEP",
        help=f"Subset of steps to run, in order. Choices: {STEP_ORDER}",
    )
    parser.add_argument("--resolution", type=float, default=None,
                        help="Override Leiden resolution.")
    parser.add_argument("--n-top-genes", type=int, default=None,
                        help="Override number of highly variable genes.")
    parser.add_argument("--deconv", action="store_true",
                        help="Force-enable the deconvolution step.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)

    # CLI overrides.
    if args.resolution is not None:
        cfg["cluster"]["resolution"] = args.resolution
    if args.n_top_genes is not None:
        cfg["normalize"]["n_top_genes"] = args.n_top_genes

    steps = list(args.steps)
    if args.deconv and "deconv" not in steps:
        steps.append("deconv")
    # Preserve canonical order regardless of how the user listed them.
    steps = [s for s in STEP_ORDER if s in steps]

    run(cfg, steps)


if __name__ == "__main__":
    main()
