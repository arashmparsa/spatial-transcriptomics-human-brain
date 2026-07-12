"""Fetch / locate a Visium dataset for the pipeline.

Resolution strategy (honest about what's what):

* **Titular dataset — spatialLIBD DLPFC (Maynard et al. 2021).** This is the
  dataset the project is named for. Its Space Ranger outputs are primarily
  distributed as a Bioconductor ``SpatialExperiment`` object, so a plain
  programmatic download of the raw matrices is not always available. If you have
  the Space Ranger output for a section (e.g. sample ``151673``) place it under
  ``data/<sample_id>/`` and this script will detect it.

* **Local pre-existing 10x data.** If ``local_data_dir`` (e.g.
  ``F:/Latch_Bio_Data``) contains a Space Ranger directory or a filtered
  feature-barcode ``.h5``, it is detected and used. This lets already-downloaded
  Visium HD colon / lung / mouse sections drive the QC / clustering / SVG steps.

* **Guaranteed-runnable fallback.** Otherwise download a 10x-hosted human brain
  Visium section via ``sc.datasets.visium_sge`` so the full pipeline runs
  end-to-end out of the box. This is a different section from the DLPFC data and
  is used only as a fallback.

Nothing large is committed to git (see ``.gitignore``).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

# Allow running as a plain script (python scripts/download_data.py).
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.io_utils import find_local_10x  # noqa: E402


SPATIALLIBD_NOTE = """\
spatialLIBD DLPFC (the titular dataset) is distributed via Bioconductor:
    - R: BiocManager::install('spatialLIBD'); spe <- spatialLIBD::fetch_data('spe')
    - or download per-sample Space Ranger outputs and place them at
      data/<sample_id>/ (containing spatial/ and filtered_feature_bc_matrix.h5).
See Maynard et al. 2021, Nat Neurosci; https://doi.org/10.1038/s41593-020-00787-0
"""


def load_config(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def main() -> None:
    parser = argparse.ArgumentParser(description="Fetch/locate Visium data.")
    parser.add_argument("--config", default="config.yaml")
    parser.add_argument(
        "--force-fallback",
        action="store_true",
        help="Skip local detection and download the 10x visium_sge fallback.",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    data = cfg["data"]
    data_dir = Path(data["data_dir"])
    data_dir.mkdir(parents=True, exist_ok=True)

    sample_id = str(data["sample_id"])
    sample_dir = data_dir / sample_id

    # 1. spatialLIBD DLPFC section already staged under data/<sample_id>/ ?
    if not args.force_fallback and (sample_dir / "spatial").is_dir():
        print(f"[download] found spatialLIBD DLPFC sample at: {sample_dir}")
        return

    # 2. Pre-existing local 10x data (e.g. F:/Latch_Bio_Data) ?
    local = data.get("local_data_dir")
    if not args.force_fallback and local:
        found = find_local_10x(local)
        if found:
            print(f"[download] detected local 10x data: {found}")
            print("[download] set data.explicit_path in config.yaml to pin it.")
            return

    # 3. Guaranteed-runnable 10x-hosted fallback section.
    print(SPATIALLIBD_NOTE)
    fallback = data["visium_sge_fallback"]
    print(f"[download] no DLPFC/local data found; downloading 10x fallback: {fallback}")
    import scanpy as sc

    adata = sc.datasets.visium_sge(sample_id=fallback)
    adata.var_names_make_unique()
    out = data_dir / f"{fallback}.h5ad"
    adata.write_h5ad(out)
    print(f"[download] wrote fallback dataset -> {out} "
          f"({adata.n_obs} spots x {adata.n_vars} genes)")


if __name__ == "__main__":
    main()
