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
    if "layer" in adata.obs:
        color.append("layer")  # manual spatialLIBD annotation
    fig = sc.pl.umap(adata, color=color, wspace=0.4, show=False, return_fig=True)
    fig.savefig(fig_dir / "umap_clusters.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Spatial (in-tissue) cluster map — only if spatial coords + image present.
    if "spatial" in adata.obsm:
        try:
            sc.pl.spatial(adata, color=color, show=False, wspace=0.3)
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

    if "layer" in adata.obs:
        results_dir = Path(cfg["paths"]["results_dir"])
        sweep = ari_sweep(adata, ccfg.get("ari_sweep", []), ccfg["resolution"], seed)
        sweep.to_csv(results_dir / "ari_vs_layers.csv", index=False)
        pd.crosstab(adata.obs[ccfg["key_added"]], adata.obs["layer"]).to_csv(
            results_dir / "cluster_vs_layer_counts.csv")
        print(sweep.to_string(index=False))
        markers = marker_check(adata, figures_dir)
        markers.to_csv(results_dir / "marker_by_layer.csv")
        print(markers.round(3).to_string())
    return adata, layer_table


# --------------------------------------------------------------------------
# Evaluation against manual cortical-layer annotation (spatialLIBD DLPFC)
# --------------------------------------------------------------------------

def score_against_layers(adata: AnnData, cluster_key: str) -> dict:
    """ARI / NMI between ``cluster_key`` and manual ``obs['layer']``.

    Only spots with a manual annotation are scored.
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    mask = adata.obs["layer"].notna().values
    truth = adata.obs["layer"].astype(str).values[mask]
    pred = adata.obs[cluster_key].astype(str).values[mask]
    return {
        "n_scored_spots": int(mask.sum()),
        "n_clusters": int(adata.obs[cluster_key].nunique()),
        "ARI": float(adjusted_rand_score(truth, pred)),
        "NMI": float(normalized_mutual_info_score(truth, pred)),
    }


def ari_sweep(adata: AnnData, resolutions: list[float], main_resolution: float,
              seed: int = 0) -> pd.DataFrame:
    """Score Leiden at several resolutions against the manual layers."""
    rows = []
    for res in sorted(set(list(resolutions) + [main_resolution])):
        key = f"leiden_r{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key, random_state=seed)
        row = {"method": "leiden", "resolution": res, **score_against_layers(adata, key)}
        rows.append(row)
    return pd.DataFrame(rows)


# Canonical markers with the compartment where each is expected to peak.
CHECK_MARKERS: dict[str, str] = {
    "MOBP": "WM",     # oligodendrocyte myelin gene -> white matter
    "PCP4": "L5",     # deep-layer (L5) neurons
    "SNAP25": "L1-L6 (neuronal grey matter, low in WM)",
}


def marker_check(adata: AnnData, figures_dir: str) -> pd.DataFrame:
    """Mean log-normalized expression of check markers per manual layer.

    Writes a spatial panel of the markers and a layer x marker dot plot.
    Returns a layer x gene table of mean log1p(CP10k) expression.
    """
    fig_dir = Path(figures_dir)
    genes = [g for g in CHECK_MARKERS if g in adata.raw.var_names]
    annotated = adata[adata.obs["layer"].notna()]
    expr = annotated.raw[:, genes].X
    expr = expr.toarray() if hasattr(expr, "toarray") else expr
    table = (
        pd.DataFrame(expr, index=annotated.obs_names, columns=genes)
        .groupby(annotated.obs["layer"].values, observed=True).mean()
    )
    table.index.name = "layer"

    if "spatial" in adata.obsm:
        sc.pl.spatial(adata, color=genes, use_raw=True, show=False, ncols=3)
        plt.savefig(fig_dir / "marker_spatial.png", dpi=150, bbox_inches="tight")
        plt.close()
    sc.pl.dotplot(annotated, var_names=list(dict.fromkeys(
        [g for ms in LAYER_MARKERS.values() for g in ms] + genes)),
        groupby="layer", use_raw=True, show=False)
    plt.savefig(fig_dir / "marker_dotplot_by_layer.png", dpi=150, bbox_inches="tight")
    plt.close()
    return table
