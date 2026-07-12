"""Normalization and dimensionality reduction.

Standard Scanpy recipe: total-count normalize to a fixed library size, log1p,
select highly variable genes, scale, PCA, build the neighbor graph, and compute
a UMAP embedding. The un-normalized counts are stashed in ``adata.raw`` and in a
``counts`` layer so that count-based tools (SpatialDE, cell2location) can recover
them downstream.
"""

from __future__ import annotations

import scanpy as sc
from anndata import AnnData


def normalize_and_reduce(adata: AnnData, cfg: dict) -> AnnData:
    """Normalize, find HVGs, scale, and run PCA / neighbors / UMAP.

    Parameters are taken from ``cfg['normalize']`` and ``cfg['project']``.
    """
    ncfg = cfg["normalize"]
    seed = cfg["project"].get("random_seed", 0)

    # Preserve raw counts before any transformation.
    adata.layers["counts"] = adata.X.copy()

    # Library-size normalization + log transform.
    sc.pp.normalize_total(adata, target_sum=ncfg["target_sum"])
    sc.pp.log1p(adata)

    # Highly variable genes. 'seurat'/'cell_ranger' expect log data;
    # 'seurat_v3' expects raw counts (pass layer='counts' in that case).
    flavor = ncfg.get("hvg_flavor", "seurat")
    if flavor == "seurat_v3":
        sc.pp.highly_variable_genes(
            adata, n_top_genes=ncfg["n_top_genes"], flavor="seurat_v3",
            layer="counts",
        )
    else:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=ncfg["n_top_genes"], flavor=flavor
        )

    # Keep the full log-normalized matrix in .raw for plotting / DE.
    adata.raw = adata

    # Restrict to HVGs for the embedding.
    adata = adata[:, adata.var["highly_variable"]].copy()

    # Scale (clip extreme values) then PCA.
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=ncfg["n_pcs"], svd_solver="arpack", random_state=seed)

    # Neighbor graph + UMAP.
    sc.pp.neighbors(
        adata,
        n_neighbors=ncfg["n_neighbors"],
        n_pcs=ncfg["n_pcs"],
        random_state=seed,
    )
    sc.tl.umap(adata, random_state=seed)

    return adata
