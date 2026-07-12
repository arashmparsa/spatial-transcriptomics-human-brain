"""Spatial transcriptomics analysis of human DLPFC (spatialLIBD) Visium data.

Modules
-------
io_utils        : load 10x Visium / spatialLIBD data into AnnData
qc              : quality-control metrics, filtering, and QC plots
normalize       : normalization, HVGs, scaling, PCA, neighbors, UMAP
cluster         : Leiden clustering and cortical-layer annotation
spatial_genes   : spatially variable genes (Squidpy Moran's I + SpatialDE)
deconvolution   : cell2location cell-type deconvolution
"""

__version__ = "0.1.0"
__author__ = "Arash Mehrabi Parsa"

__all__ = [
    "io_utils",
    "qc",
    "normalize",
    "cluster",
    "spatial_genes",
    "deconvolution",
]
