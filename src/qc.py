"""Quality control: metrics, filtering, and QC plots for Visium spots.

Follows the standard Scanpy spatial QC recipe: annotate mitochondrial genes,
compute per-spot metrics with ``sc.pp.calculate_qc_metrics``, visualise the
distributions, then filter low-quality spots and rarely detected genes.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless-safe backend for pipeline runs
import matplotlib.pyplot as plt
import scanpy as sc
from anndata import AnnData


def compute_qc_metrics(adata: AnnData, mt_prefix: str = "MT-") -> AnnData:
    """Flag mitochondrial genes and compute per-spot / per-gene QC metrics.

    Adds ``pct_counts_mt``, ``n_genes_by_counts``, ``total_counts`` to
    ``adata.obs`` (in place) and returns ``adata``.
    """
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    return adata


def plot_qc(adata: AnnData, figures_dir: str) -> None:
    """Save violin and scatter QC plots to ``figures_dir``."""
    fig_dir = Path(figures_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Violin distributions of the three core QC metrics.
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    for ax, key in zip(
        axes, ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    ):
        sc.pl.violin(adata, key, jitter=0.4, ax=ax, show=False)
    fig.tight_layout()
    fig.savefig(fig_dir / "qc_violin.png", dpi=150)
    plt.close(fig)

    # Scatter: counts vs mito fraction and counts vs genes.
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", ax=axes[0], show=False)
    sc.pl.scatter(
        adata, x="total_counts", y="n_genes_by_counts", ax=axes[1], show=False
    )
    fig.tight_layout()
    fig.savefig(fig_dir / "qc_scatter.png", dpi=150)
    plt.close(fig)


def filter_cells_and_genes(
    adata: AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    max_pct_mt: float = 20.0,
) -> AnnData:
    """Filter low-quality spots and rarely detected genes.

    * spots with < ``min_genes`` detected genes are dropped,
    * spots with mito fraction > ``max_pct_mt`` are dropped,
    * genes seen in < ``min_cells`` spots are dropped.
    """
    n_start = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mt].copy()
    print(
        f"[qc] spots {n_start} -> {adata.n_obs} "
        f"(min_genes={min_genes}, max_pct_mt={max_pct_mt}); "
        f"genes -> {adata.n_vars} (min_cells={min_cells})"
    )
    return adata


def run_qc(adata: AnnData, cfg: dict) -> AnnData:
    """End-to-end QC step driven by the config dict."""
    qc_cfg = cfg["qc"]
    figures_dir = cfg["paths"]["figures_dir"]

    adata = compute_qc_metrics(adata, mt_prefix=qc_cfg["mt_prefix"])
    plot_qc(adata, figures_dir)
    adata = filter_cells_and_genes(
        adata,
        min_genes=qc_cfg["min_genes"],
        min_cells=qc_cfg["min_cells"],
        max_pct_mt=qc_cfg["max_pct_mt"],
    )
    return adata
