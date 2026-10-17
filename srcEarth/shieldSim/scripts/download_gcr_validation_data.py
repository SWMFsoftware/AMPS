#!/usr/bin/env python3
"""Download and verify the pinned AMS-02 and ACE/CRIS validation sources.

This script is intentionally not run by CTest: normal tests are offline and use
the committed immutable fixtures.  The downloader is an auditable acquisition
recipe for a clean checkout.  Every response is checked before installation;
an upstream change is an error that requires scientific review and a new
fixture revision, never an automatic checksum update.
"""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import tempfile
import urllib.request


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_DESTINATION = ROOT / "tests" / "reference_spectra" / "gcr"


# (relative destination, immutable URL/request, reviewed SHA-256).  The HAPI
# data request lists fields explicitly so a server-side default cannot change
# the CSV layout consumed by ACECRISHapiAdapter.
ARTIFACTS = (
    ("ams02_2025/raw/table-S1-S147.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S1-S147.csv",
     "6be27d8115c47b2102eec21bee7629751677641c8c808727dc3041ee0e2f030f"),
    ("ams02_2025/raw/table-S148-S294.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S148-S294.csv",
     "c2e78af4346acfb13206c6f9cada18e2c9e2b8473795453488d3ba50379e84de"),
    ("ams02_2025/raw/table-S295-S441.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S295-S441.csv",
     "28c62dad061f0d39182c2be4d4af371df0b94759a0734706c02305988ccf0e1d"),
    ("ams02_2025/raw/table-S442-S588.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S442-S588.csv",
     "fc49fe5605bcc27d923c4773dd54c30a0834b94c0a6bb31c1b233e2749446fd7"),
    ("ams02_2025/raw/table-S589-S735.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S589-S735.csv",
     "4e511736fe2f8688f042e528e7f97ecce22d721c87b154692cc3516f7ac67a37"),
    ("ams02_2025/raw/table-S736-S882.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S736-S882.csv",
     "89358efc6ad9dcc01343f7390bb6a37d7b850cc49b5ce51f2f1e02d68be4334c"),
    ("ams02_2025/raw/table-S883-S1029.csv",
     "https://ams02.space/sites/default/files/publication/202403/table-S883-S1029.csv",
     "0ef5277ae6ede83e32882b1182ad58705c8012490d09931cbe306bd94ebdffaa"),
    ("ace_cris/ac_h3_cris_hapi_info.json",
     "https://cdaweb.gsfc.nasa.gov/hapi/info?id=AC_H3_CRIS",
     "319400e7c2d3298807340cf6d6a627103e14e547e7208c955159fcef60bf6ede"),
    ("ace_cris/ac_h3_cris_20191201_all.csv",
     "https://cdaweb.gsfc.nasa.gov/hapi/data?id=AC_H3_CRIS&time.min=2019-12-01T00:00:00Z&time.max=2019-12-02T00:00:00Z&parameters=Time,flux_B,flux_C,flux_N,flux_O,flux_F,flux_Ne,flux_Na,flux_Mg,flux_Al,flux_Si,flux_P,flux_S,flux_Cl,flux_Ar,flux_K,flux_Ca,flux_Sc,flux_Ti,flux_V,flux_Cr,flux_Mn,flux_Fe,flux_Co,flux_Ni,cnt_B,cnt_C,cnt_N,cnt_O,cnt_F,cnt_Ne,cnt_Na,cnt_Mg,cnt_Al,cnt_Si,cnt_P,cnt_S,cnt_Cl,cnt_Ar,cnt_K,cnt_Ca,cnt_Sc,cnt_Ti,cnt_V,cnt_Cr,cnt_Mn,cnt_Fe,cnt_Co,cnt_Ni&format=csv",
     "33e018b7934971670dab40e883f03c0f4a3542cdd5a96103d2cf729957820cfd"),
    ("ace_cris/raw/ac_h3_cris_20190101_v04.cdf",
     "https://cdaweb.gsfc.nasa.gov/pub/data/ace/cris/level_2_cdaweb/cris_h3/2019/ac_h3_cris_20190101_v04.cdf",
     "3b43c0fab98d1e4947810fde32c1dbccb556d652affac330f67b9e97f04fd5a1"),
    ("ace_cris/raw/cris_data_Bartels.hdf",
     "https://izw1.caltech.edu/ACE/ASC/DATA/level2/cris/cris_data_Bartels.hdf",
     "85c39309100957466f2d47251d65e22b0760ea81f88c03d18656f750b4ad33a5"),
    ("ace_cris/raw/cris_factors.txt",
     "https://izw1.caltech.edu/ACE/ASC/DATA/level2/cris/cris_factors.txt",
     "fe48e9572fef3fbe6ae70e5fbb87041f336fbfb333403ad1e64ea2961d050d8c"),
    ("ace_cris/raw/cris_release_notes.txt",
     "https://izw1.caltech.edu/ACE/ASC/DATA/level2/cris/cris_release_notes",
     "49bf5e06c64f9b7f55d77603f141c8e94a87f62e4bc8e622e4ba2782b7b8096d"),
    ("ace_cris/raw/cris_l2desc.html",
     "https://izw1.caltech.edu/ACE/ASC/level2/cris_l2desc.html",
     "286713205255575bab0e0b3bb827c62eb94a4ce31f511df2d77dc4063be54304"),
)


def digest(path: pathlib.Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=pathlib.Path, default=DEFAULT_DESTINATION)
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="re-fetch already valid files and prove that upstream bytes still match",
    )
    args = parser.parse_args()
    destination = args.destination.resolve()

    with tempfile.TemporaryDirectory(prefix="shieldsim-gcr-download-") as temporary:
        staging = pathlib.Path(temporary)
        for relative, url, expected in ARTIFACTS:
            target = destination / relative
            if target.is_file() and digest(target) == expected and not args.force_download:
                print(f"VERIFIED {relative}")
                continue
            staged = staging / relative
            staged.parent.mkdir(parents=True, exist_ok=True)
            request = urllib.request.Request(url, headers={"User-Agent": "shieldSim-validation/1"})
            with urllib.request.urlopen(request, timeout=120) as response, staged.open("wb") as output:
                shutil.copyfileobj(response, output)
            actual = digest(staged)
            if actual != expected:
                raise RuntimeError(
                    f"upstream bytes changed for {relative}: expected {expected}, received {actual}"
                )
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and digest(target) != expected:
                raise RuntimeError(
                    f"refusing to overwrite unreviewed local bytes at {target}; move them manually first"
                )
            staged.replace(target)
            print(f"INSTALLED {relative}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
