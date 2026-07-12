"""Leiden clustering and putative cortical-layer annotation.

Clusters the neighbor graph with Leiden, plots the clusters in UMAP and in
tissue space, and provides a light-touch annotation that maps clusters to
putative cortical layers using canonical laminar marker genes (Maynard et al.,
2021). The mapping is a heuristic starting point — real layer assignment should
be validated against the histology / manual annotations.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData

# Canonical DLPFC laminar markers (Maynard et al. 2021, Nat Neurosci).
# Used to *suggest* a layer label for each Leiden cluster by mean expression.
LAYER_MARKERS: dict[str, list[str]] = {
    "L1": ["AQP4", "RELN", "NDNF"],
    "L2": ["HPCAL1", "CARTPT"],
    "L3": ["CUX2", "RASGRF2"],
    "L4": ["RORB", "PCP4"],
    "L5": ["BCL11B", "FEZF2", "PCP4"],
    "L6": ["TLE4", "SEMA3E", "NR4A2"],
    "WM": ["MBP", "MOBP", "PLP1", "GFAP"],
}


def run_leiden(adata: AnnData, resolution: float = 0.8, key_added: str = "leiden",
               seed: int = 0) -> AnnData:
    """Run Leiden community detection on the existing neighbor graph."""
    sc.tl.leiden(
        adata, resolution=resolution, key_added=key_added, random_state=seed
    )
    n_clusters = adata.obs[key_added].nunique()
    print(f"[cluster] Leiden res={resolution}: {n_clusters} clusters")
    return adata


def annotate_layers(adata: AnnData, cluster_key: str = "leiden") -> pd.DataFrame:
    """Assign each cluster a putative cortical layer by top marker score.

    Uses ``sc.tl.score_genes`` per layer marker set, then labels each cluster
    with the layer whose mean score is highest across that cluster's spots.

    Returns a table (cluster -> putative_layer, per-layer mean scores) and
    writes ``adata.obs['putative_layer']``.
    """
    # Score each layer signature (only markers actually present in the data).
    layer_score_cols: dict[str, str] = {}
    for layer, markers in LAYER_MARKERS.items():
        present = [g for g in markers if g in adata.raw.var_names] if adata.raw \
            else [g for g in markers if g in adata.var_names]
        if not present:
            continue
        col = f"score_{layer}"
        sc.tl.score_genes(adata, gene_list=present, score_name=col, use_raw=True)
        layer_score_cols[layer] = col

    # Mean score per cluster, then argmax layer.
    scores = adata.obs.groupby(cluster_key)[list(layer_score_cols.values())].mean()
    scores.columns = list(layer_score_cols.keys())
    cluster_to_layer = scores.idxmax(axis=1)

    adata.obs["putative_layer"] = (
        adata.obs[cluster_key].map(cluster_to_layer).astype("category")
    )

    table = scores.copy()
    table["putative_layer"] = cluster_to_layer
    return table


def plot_clusters(adata: AnnData, figures_dir: str, cluster_key: str = "leiden") -> None:
    """Save UMAP and in-tissue spatial cluster plots."""
    fig_dir = Path(figures_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    # UMAP colored by cluster and putative layer.
    color = [cluster_key]
    if "putative_layer" in adata.obs:
        color.append("putative_layer")
    fig = sc.pl.umap(adata, color=color, wspace=0.4, show=False, return_fig=True)
    fig.savefig(fig_dir / "umap_clusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Spatial (in-tissue) cluster map — only if spatial coords + image present.
    if "spatial" in adata.obsm:
        try:
            sc.pl.spatial(adata, color=color, show=False)
            plt.savefig(fig_dir / "spatial_clusters.png", dpi=150, bbox_inches="tight")
            plt.close()
        except Exception as exc:  # noqa: BLE001 - plotting is best-effort
            print(f"[cluster] spatial plot skipped: {exc}")


def run_clustering(adata: AnnData, cfg: dict) -> tuple[AnnData, pd.DataFrame]:
    """End-to-end clustering + annotation step driven by the config dict."""
    ccfg = cfg["cluster"]
    seed = cfg["project"].get("random_seed", 0)
    figures_dir = cfg["paths"]["figures_dir"]

    adata = run_leiden(
        adata,
        resolution=ccfg["resolution"],
        key_added=ccfg["key_added"],
        seed=seed,
    )
    layer_table = annotate_layers(adata, cluster_key=ccfg["key_added"])
    plot_clusters(adata, figures_dir, cluster_key=ccfg["key_added"])
    return adata, layer_table
