"""Spatially variable genes (SVGs) via two complementary methods.

1. **Squidpy Moran's I** — build a spatial neighbor graph on the Visium hex grid
   (``sq.gr.spatial_neighbors``) and compute spatial autocorrelation
   (``sq.gr.spatial_autocorr`` with ``mode='moran'``). Genes with high positive
   Moran's I are spatially structured.

2. **SpatialDE** — Svensson et al. (2018) Gaussian-process test on
   (normalized-count, coordinate) pairs, returning per-gene q-values.

Both write ranked tables to ``results/`` and plot the top SVGs in tissue space.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from anndata import AnnData


def moran_svgs(adata: AnnData, coord_type: str = "grid", n_neighs: int = 6,
               n_perms: int = 100, seed: int = 0) -> pd.DataFrame:
    """Compute Moran's I spatial autocorrelation with Squidpy.

    Returns the ``adata.uns['moranI']`` table sorted by descending Moran's I.
    Operates on highly variable genes when that annotation is present, to keep
    the permutation test tractable.
    """
    # Spatial neighbor graph on the Visium hexagonal grid.
    sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)

    genes = None
    if "highly_variable" in adata.var:
        genes = adata[:, adata.var["highly_variable"]].var_names.tolist()

    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=genes,
        n_perms=n_perms,
        seed=seed,
    )
    moran = adata.uns["moranI"].sort_values("I", ascending=False)
    return moran


def spatialde_svgs(adata: AnnData) -> pd.DataFrame:
    """Run SpatialDE on normalized counts + spatial coordinates.

    SpatialDE expects a spots x genes DataFrame of (normalized) expression and a
    coordinate array. We reconstruct normalized counts from the ``counts`` layer
    via ``NaiveDE.stabilize`` + ``regress_out`` per the SpatialDE tutorial.

    Returns the per-gene results table sorted by q-value.
    """
    import NaiveDE  # bundled with SpatialDE
    import SpatialDE

    # Raw counts and per-spot total counts for the regression covariate.
    counts_mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    counts_mat = counts_mat.toarray() if hasattr(counts_mat, "toarray") else np.asarray(counts_mat)
    counts = pd.DataFrame(counts_mat, index=adata.obs_names, columns=adata.var_names)

    sample_info = pd.DataFrame(index=adata.obs_names)
    coords = adata.obsm["spatial"]
    sample_info["x"] = coords[:, 0]
    sample_info["y"] = coords[:, 1]
    sample_info["total_counts"] = np.asarray(counts.sum(axis=1)).ravel()

    # Variance-stabilize (Anscombe) then regress out library size.
    norm_expr = NaiveDE.stabilize(counts.T).T
    resid_expr = NaiveDE.regress_out(
        sample_info, norm_expr.T, "np.log(total_counts)"
    ).T

    results = SpatialDE.run(sample_info[["x", "y"]].values, resid_expr)
    results = results.sort_values("qval")
    return results


def plot_top_svgs(adata: AnnData, genes: list[str], figures_dir: str,
                  prefix: str = "moran") -> None:
    """Plot the top SVGs in tissue space (best-effort)."""
    fig_dir = Path(figures_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    genes = [g for g in genes if g in (adata.raw.var_names if adata.raw else adata.var_names)]
    if not genes:
        return
    try:
        sc.pl.spatial(adata, color=genes, show=False, use_raw=True)
        plt.savefig(fig_dir / f"{prefix}_top_svgs.png", dpi=150, bbox_inches="tight")
        plt.close()
    except Exception as exc:  # noqa: BLE001 - plotting is best-effort
        print(f"[svg] spatial SVG plot skipped: {exc}")


def run_spatial_genes(adata: AnnData, cfg: dict) -> dict[str, pd.DataFrame]:
    """End-to-end SVG step: Moran's I (+ optional SpatialDE), save tables/plots."""
    scfg = cfg["svg"]
    seed = cfg["project"].get("random_seed", 0)
    results_dir = Path(cfg["paths"]["results_dir"])
    figures_dir = cfg["paths"]["figures_dir"]
    n_svgs = scfg["n_svgs"]
    results_dir.mkdir(parents=True, exist_ok=True)

    out: dict[str, pd.DataFrame] = {}

    # 1. Squidpy Moran's I.
    moran = moran_svgs(
        adata,
        coord_type=scfg["coord_type"],
        n_neighs=scfg["n_neighs"],
        seed=seed,
    )
    moran.to_csv(results_dir / "moran_svgs.csv")
    out["moran"] = moran
    plot_top_svgs(adata, moran.head(n_svgs).index.tolist(), figures_dir, prefix="moran")
    print(f"[svg] Moran's I top genes: {list(moran.head(n_svgs).index)}")

    # 2. SpatialDE (optional; slower).
    if scfg.get("run_spatialde", True):
        try:
            sde = spatialde_svgs(adata)
            sde.to_csv(results_dir / "spatialde_svgs.csv", index=False)
            out["spatialde"] = sde
            top = sde.head(n_svgs)["g"].tolist() if "g" in sde.columns else []
            plot_top_svgs(adata, top, figures_dir, prefix="spatialde")
            print(f"[svg] SpatialDE top genes: {top}")
        except Exception as exc:  # noqa: BLE001
            print(f"[svg] SpatialDE skipped ({type(exc).__name__}): {exc}")

    return out
