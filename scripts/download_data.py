"""Download the spatialLIBD human DLPFC Visium section used by this pipeline.

Source: Maynard et al. (2021), Nat Neurosci, https://doi.org/10.1038/s41593-020-00787-0
(the ``spatialLIBD`` / LieberInstitute HumanPilot release).

Files fetched for one section (default ``151673``) into ``data/<sample_id>/``
(``data/`` is git-ignored):

* ``filtered_feature_bc_matrix.h5``      - Space Ranger counts (public S3 bucket)
* ``spatial/tissue_hires_image.png``     - H&E image, hires (public S3 bucket)
* ``spatial/tissue_lowres_image.png``    - H&E image, lowres (HumanPilot GitHub)
* ``spatial/tissue_positions_list.csv``  - spot barcodes -> array/pixel coords
* ``spatial/scalefactors_json.json``     - image scale factors
* ``layer_annotations.tsv``              - manual cortical layer labels (L1-L6, WM)
                                           for this section, from
                                           ``10X/barcode_level_layer_map.tsv``

Usage::

    python scripts/download_data.py                 # sample 151673
    python scripts/download_data.py --sample 151507
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

S3 = "https://spatial-dlpfc.s3.us-east-2.amazonaws.com"
GH = "https://raw.githubusercontent.com/LieberInstitute/HumanPilot/master/10X"


def _fetch(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[download] exists, skipping: {dest}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"[download] {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)


def download(sample: str, data_dir: Path) -> Path:
    out = data_dir / sample
    spatial = out / "spatial"
    _fetch(f"{S3}/h5/{sample}_filtered_feature_bc_matrix.h5",
           out / "filtered_feature_bc_matrix.h5")
    _fetch(f"{S3}/images/{sample}_tissue_hires_image.png",
           spatial / "tissue_hires_image.png")
    _fetch(f"{GH}/{sample}/tissue_lowres_image.png",
           spatial / "tissue_lowres_image.png")
    _fetch(f"{GH}/{sample}/tissue_positions_list.txt",
           spatial / "tissue_positions_list.csv")
    _fetch(f"{GH}/{sample}/scalefactors_json.json",
           spatial / "scalefactors_json.json")

    # Layer annotations: keep only this section's rows.
    layer_out = out / "layer_annotations.tsv"
    if not layer_out.exists():
        all_map = data_dir / "barcode_level_layer_map.tsv"
        _fetch(f"{GH}/barcode_level_layer_map.tsv", all_map)
        rows = [ln for ln in all_map.read_text().splitlines()
                if ln.split("\t")[1] == sample]
        layer_out.write_text("barcode\tsample\tlayer\n" + "\n".join(rows) + "\n")
        print(f"[download] wrote {len(rows)} annotated spots -> {layer_out}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sample", default="151673")
    parser.add_argument("--data-dir", default="data")
    args = parser.parse_args()
    out = download(args.sample, Path(args.data_dir))
    print(f"[download] done: {out}")


if __name__ == "__main__":
    main()
