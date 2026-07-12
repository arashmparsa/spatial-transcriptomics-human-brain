"""cell2location deconvolution: estimate cell-type abundances per Visium spot.

Two-stage workflow (Kleshchevnikov et al., 2022):

1. **Reference signatures** — fit a negative-binomial regression model
   (``RegressionModel``) to an annotated scRNA-seq / snRNA-seq reference to
   estimate per-gene, per-cell-type expression signatures.
2. **Spatial mapping** — fit ``Cell2location`` to the Visium spots using those
   signatures to estimate the abundance of every cell type in every spot.

Outputs the per-spot 5%-quantile cell-type abundances (the conservative point
estimate cell2location recommends) as a table and as spatial maps.

Obtaining a paired DLPFC reference
----------------------------------
Provide an ``.h5ad`` of human DLPFC / prefrontal-cortex single-nucleus RNA-seq
with a ``cell_type`` column in ``.obs`` (and optionally a donor/batch key).
Suitable public sources include the Allen Brain Map human cortex snRNA-seq
atlas and the snRNA-seq released alongside spatialLIBD follow-up studies. The
reference must share gene identifiers (symbols) with the Visium data; the code
intersects genes and applies cell2location's standard feature filter.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from anndata import AnnData


def _get_raw_counts(adata: AnnData) -> AnnData:
    """Return an AnnData whose ``.X`` holds raw integer counts.

    cell2location must be trained on unnormalized counts. If a ``counts`` layer
    exists we restore it into ``.X``.
    """
    a = adata.copy()
    if "counts" in a.layers:
        a.X = a.layers["counts"].copy()
    return a


def train_reference_signatures(
    ref: AnnData,
    cell_type_key: str = "cell_type",
    batch_key: str | None = None,
    max_epochs: int = 250,
):
    """Estimate reference cell-type signatures with the regression model.

    Returns ``(ref, inf_aver)`` where ``inf_aver`` is the genes x cell-types
    signature matrix used for spatial mapping.
    """
    from cell2location.models import RegressionModel
    from cell2location.utils.filtering import filter_genes

    ref = _get_raw_counts(ref)

    # cell2location's recommended permissive gene filter for the reference.
    selected = filter_genes(
        ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12
    )
    ref = ref[:, selected].copy()

    RegressionModel.setup_anndata(
        ref, batch_key=batch_key, labels_key=cell_type_key
    )
    model = RegressionModel(ref)
    model.train(max_epochs=max_epochs)

    ref = model.export_posterior(
        ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500}
    )

    # Per-cell-type expression signature (genes x cell types).
    if "means_per_cluster_mu_fg" in ref.varm:
        inf_aver = ref.varm["means_per_cluster_mu_fg"][
            [f"means_per_cluster_mu_fg_{ct}" for ct in ref.uns["mod"]["factor_names"]]
        ].copy()
    else:
        inf_aver = ref.var[
            [f"means_per_cluster_mu_fg_{ct}" for ct in ref.uns["mod"]["factor_names"]]
        ].copy()
    inf_aver.columns = ref.uns["mod"]["factor_names"]
    return ref, inf_aver


def map_spatial(
    adata_vis: AnnData,
    inf_aver,
    n_cells_per_location: int = 8,
    detection_alpha: float = 20.0,
    max_epochs: int = 30000,
) -> AnnData:
    """Fit cell2location to the Visium spots and export per-spot abundances."""
    import cell2location

    adata_vis = _get_raw_counts(adata_vis)

    # Restrict both matrices to shared genes.
    shared = [g for g in adata_vis.var_names if g in inf_aver.index]
    adata_vis = adata_vis[:, shared].copy()
    inf_aver = inf_aver.loc[shared, :].copy()

    cell2location.models.Cell2location.setup_anndata(adata_vis)
    model = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        N_cells_per_location=n_cells_per_location,
        detection_alpha=detection_alpha,
    )
    model.train(max_epochs=max_epochs, batch_size=None, train_size=1)

    adata_vis = model.export_posterior(
        adata_vis, sample_kwargs={"num_samples": 1000, "batch_size": model.adata.n_obs}
    )

    # cell2location recommends the 5% quantile as the conservative abundance.
    abund = adata_vis.obsm["q05_cell_abundance_w_sf"]
    abund.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in abund.columns]
    adata_vis.obs[abund.columns] = abund.values
    adata_vis.uns["cell_type_names"] = list(abund.columns)
    return adata_vis


def save_outputs(adata_vis: AnnData, results_dir: str, figures_dir: str) -> None:
    """Write per-spot proportions table and spatial abundance maps."""
    results_path = Path(results_dir)
    fig_dir = Path(figures_dir)
    results_path.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    cell_types = adata_vis.uns.get("cell_type_names", [])
    if cell_types:
        # Convert abundances to per-spot proportions and save.
        ab = adata_vis.obs[cell_types].to_numpy()
        props = ab / np.clip(ab.sum(axis=1, keepdims=True), 1e-8, None)
        import pandas as pd

        pd.DataFrame(
            props, index=adata_vis.obs_names, columns=cell_types
        ).to_csv(results_path / "cell_type_proportions.csv")

        # Spatial abundance maps (best-effort).
        try:
            sc.pl.spatial(adata_vis, color=cell_types, show=False, ncols=4)
            plt.savefig(fig_dir / "spatial_cell_abundance.png", dpi=150, bbox_inches="tight")
            plt.close()
        except Exception as exc:  # noqa: BLE001
            print(f"[deconv] spatial abundance plot skipped: {exc}")


def run_deconvolution(adata_vis: AnnData, cfg: dict) -> AnnData:
    """End-to-end cell2location deconvolution driven by the config dict."""
    dcfg = cfg["deconv"]
    results_dir = cfg["paths"]["results_dir"]
    figures_dir = cfg["paths"]["figures_dir"]

    ref_path = dcfg["reference_h5ad"]
    if not Path(ref_path).exists():
        raise FileNotFoundError(
            f"Reference not found: {ref_path}. Provide a snRNA-seq DLPFC "
            f"reference .h5ad with a '{dcfg['cell_type_key']}' column in .obs. "
            "See the module docstring for suggested sources."
        )

    ref = sc.read_h5ad(ref_path)
    print(f"[deconv] reference: {ref.n_obs} cells x {ref.n_vars} genes")

    ref, inf_aver = train_reference_signatures(
        ref,
        cell_type_key=dcfg["cell_type_key"],
        batch_key=dcfg.get("batch_key"),
        max_epochs=dcfg["ref_max_epochs"],
    )
    adata_vis = map_spatial(
        adata_vis,
        inf_aver,
        n_cells_per_location=dcfg["n_cells_per_location"],
        detection_alpha=dcfg["detection_alpha"],
        max_epochs=dcfg["map_max_epochs"],
    )
    save_outputs(adata_vis, results_dir, figures_dir)
    return adata_vis
