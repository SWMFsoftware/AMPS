#!/usr/bin/env python3
"""Download and verify the public inputs recovered for CRaTER validation.

This is an acquisition recipe, not a Stage-8 validation run.  Normal tests are
offline and consume the committed, checksum-pinned files.  The script refuses
to update a checksum when an upstream object changes: such a change requires a
new reviewed fixture revision.

Two BON2020 spectra are embedded in a 383 MB GLACE response-kernel ZIP.  The
ZIP is an immutable Zenodo object, so this script retrieves only the exact
DEFLATE byte ranges for those 56 KB JSON members.  Output SHA-256 checks make
the ranged extraction equivalent to downloading and unpacking the full ZIP,
without storing hundreds of megabytes unrelated to shieldSim validation.
"""

from __future__ import annotations

import argparse
import hashlib
import pathlib
import shutil
import tempfile
import urllib.request
import zlib


ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_DESTINATION = ROOT / "tests" / "reference_crater"
DEFAULT_OBSERVATION = (
    ROOT / "tests" / "reference_spectra" / "observations"
    / "crater_gcr_let_spectra_v1_3.csv"
)
PDS = "https://pds-ppi.igpp.ucla.edu/data/LRO-L-CRAT-3_4-DDR-PROCESSED-V1.0"


# (relative destination, exact source URL, reviewed SHA-256).  Several PDS
# catalog files are living archive endpoints rather than content-addressed
# URLs.  Their reviewed bytes are pinned here; a future quarterly PDS update is
# deliberately reported as an error instead of being accepted silently.
ARTIFACTS = (
    ("zenodo/record_v1_3.json",
     "https://zenodo.org/api/records/22883077",
     "0512a5f0685459c6045cc17b7cc35c91a552e7818da40d351490e29677b4e76e"),
    ("zenodo/000README_GCR_LET_v1_3.pdf",
     "https://zenodo.org/api/records/22883077/files/000README_GCR_LET_v1_3.pdf/content",
     "f96cfbe4acf5c0d7efc326ff903de5a9a0835048fd36379f25570279c1d1d37a"),
    ("glace/record_v1_0.json",
     "https://zenodo.org/api/records/8343472",
     "58be6667b34532ab935e31a7719f0f5356323d06b99eb0b6fd5aa66a44a63310"),
    ("glace/000README_GLACE_v1_0.pdf",
     "https://zenodo.org/api/records/8343472/files/000README_GLACE_v1_0.pdf/content",
     "c9163e0fc9266a5bc3e058ae79d3ad30aec9dbe44e9e03ee1efe3327a605df33"),
    ("publications/CRaTER_SSR2010.pdf",
     "https://crater.unh.edu/pdfs/CRaTER_SSR2010.pdf",
     "206d909ae0ce6c2b4d25e75745df5652016ee438bd018eb49b09896430313f9c"),
    ("publications/deepblue_item_2027.42_163885.json",
     "https://backend.production.deepblue-documents.lib.umich.edu/server/api/pid/find?id=2027.42%2F163885",
     "ffa8aab4e62737ff47b4b43adab66bad34b07c9883728d9b68865fc1513fb415"),
    ("publications/Looper2020_CRaTER_GCR_LET.pdf",
     "https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/85659e4d-b599-4c2f-a24f-3befde5e7756/content",
     "ab0ec871a90b24942aafdde939fd3d3061f0e6638bbe758c9fff83d5d409c00e"),
    ("publications/swe21079-sup-0001-2020SW002543-ts1.pdf",
     "https://backend.production.deepblue-documents.lib.umich.edu/server/api/core/bitstreams/7a9ca2dd-ab81-40e4-a0b9-506ee27a7576/content",
     "76696c3e98dd3f6ae2470e23e34ad67795a694b8c20d88da4a847cb0a4610c69"),
    ("pds/AAREADME.TXT", f"{PDS}/AAREADME.TXT",
     "b8aabb344e08caabbf048c57dba781e488a17bf733cff6b202e53e95e61aa7fe"),
    ("pds/CHECKSUMS.TXT", f"{PDS}/CHECKSUMS.TXT",
     "b2ea8217c6eca2aadf4a192aa0e2f7889e1a368763fb244d63cf176297e2a37c"),
    ("pds/ERRATA.TXT", f"{PDS}/ERRATA.TXT",
     "49a87c0cb426abcfac03021ab0c1713141c74e78f4e4d2081c0ccd7f9f1d2ceb"),
    ("pds/VOLDESC.CAT", f"{PDS}/VOLDESC.CAT",
     "1d80dcb1effd227fbedc7e794b0525efb1afe141fcce46046228d9165ccc6d86"),
    ("pds/CALINFO.TXT", f"{PDS}/CALIB/CALINFO.TXT",
     "aa9f39de362d65e358511526db0d470f7aa2ad2b63da75beef09e36d09fcc269"),
    ("pds/CRAT_CAL_PLAN.LBL", f"{PDS}/CALIB/CRAT_CAL_PLAN.LBL",
     "91c27fb763e1f78ae00ff0cb314e191cc25439188bb74cbc645aa095ad96a8a9"),
    ("pds/CRAT_CAL_PLAN.PDF", f"{PDS}/CALIB/CRAT_CAL_PLAN.PDF",
     "620dc29900c4b797476d0b826fa800db8387e7a439c333f308e83987e103668d"),
    ("pds/CATINFO.TXT", f"{PDS}/CATALOG/CATINFO.TXT",
     "7702a1dfedd1452750f8d182f80d2ad90f4e4c424ccaef7be4fdf7901c574c8e"),
    ("pds/LRO_CRAT_INST.CAT", f"{PDS}/CATALOG/LRO_CRAT_INST.CAT",
     "374ae94992089b29b6544d54823112a7b7cf812e7f24a371e740d040b3b5fe11"),
    ("pds/LRO_CRAT_L2_DS.CAT", f"{PDS}/CATALOG/LRO_CRAT_L2_DS.CAT",
     "d26b4815e51e9bc1c6bcd6e1c826f2a5bee28af964abe273855b945cb4db9799"),
    ("pds/LRO_CRAT_REF.CAT", f"{PDS}/CATALOG/LRO_CRAT_REF.CAT",
     "5eeddfc8359ca3ae448685137133fb3325e2f2f12d818402fdeb5ea3c23c1ce8"),
    ("pds/LRO_SOFTWARE.CAT", f"{PDS}/CATALOG/LRO_SOFTWARE.CAT",
     "268c02a9751589f635fd65079b7d2fc5d3f154881613d7935fb7d90cb90271df"),
    ("pds/CRATDD.FUL", f"{PDS}/DOCUMENT/CRATDD.FUL",
     "39630b8ee03da5fc688aaeedcc6e3e78fdfe568c7e4318f7decf7a9b405369e9"),
    ("pds/DOCINFO.TXT", f"{PDS}/DOCUMENT/DOCINFO.TXT",
     "539ae0f713454075e452fce479ff76a1f3bb1c91a3f78496510dcc592a091d8b"),
    ("pds/GAPS_HK.LBL", f"{PDS}/DOCUMENT/GAPS_HK.LBL",
     "28102a11eb64acb9bd7e3688675faeadb64040847305d2ade8cf41ee4b6c581c"),
    ("pds/GAPS_HK.TAB", f"{PDS}/DOCUMENT/GAPS_HK.TAB",
     "fcdb1918471c471a1ecef700ece563a388da73846d731202bc9c17b8e81f21e4"),
    ("pds/GAPS_PRI.LBL", f"{PDS}/DOCUMENT/GAPS_PRI.LBL",
     "db2cff4488b430d2afc46e73a4fb05c05fa5c28e9faa07fa695219c23b715525"),
    ("pds/GAPS_PRI.TAB", f"{PDS}/DOCUMENT/GAPS_PRI.TAB",
     "da8e8839484c018bd1e019e5619aff56c1ffa4b72d1e2b84ecc62defa4c29afd"),
    ("pds/GAPS_SEC.LBL", f"{PDS}/DOCUMENT/GAPS_SEC.LBL",
     "6eba213700c6aeb2b24c68ea0761b7660bc80794e9d93343493819811f45ef2a"),
    ("pds/GAPS_SEC.TAB", f"{PDS}/DOCUMENT/GAPS_SEC.TAB",
     "1483563a5ca98a87de4347f8781a27419f94822e0a669eac2a4cc11e78b125e8"),
    ("pds/MODE_CHANGES.LBL", f"{PDS}/DOCUMENT/MODE_CHANGES.LBL",
     "0518269e65e6b539531a5f369903524f41157b3d4d8421ebf11a5b779ce9b976"),
    ("pds/MODE_CHANGES.TAB", f"{PDS}/DOCUMENT/MODE_CHANGES.TAB",
     "a62abdd9b04e898e66b8f2ae7802a18f854537c38ea80605bce3a80cd19c7bb3"),
    ("pds/VOLSIS.PDF", f"{PDS}/DOCUMENT/VOLSIS.PDF",
     "c3e7b416313d5a329f1ccbd99a79ea4dcfbe8984ee21d7991e439e071d23c4a8"),
    ("pds/CRAT_L2_PRI.FMT", f"{PDS}/LABEL/CRAT_L2_PRI.FMT",
     "e1d8b438f599ca0e95e0d70205a8c7c947705897e32716650dc1cf3bbb192bf2"),
    ("pds/CRAT_L2_PRI_2009259_V01.LBL",
     f"{PDS}/DATA/2009/2009259/CRAT_L2_PRI_2009259_V01.LBL",
     "27822bbd9dce676e2736145f0ce30c59743f2625612857aaab55919d2bee65a7"),
    ("pds/PIPE2.LBL", f"{PDS}/SOFTWARE/PIPE2.LBL",
     "5ce14527050ea7049fe60d9b2d277e4335330e06de16dd1a74ead78f2f308692"),
    ("pds/PIPE2.ZIP", f"{PDS}/SOFTWARE/PIPE2.ZIP",
     "2e03ba6e04b41aa5efba91048b1ed550be60151c4e9fd3581a4c4290e13b6adf"),
    ("pds/SOFTINFO.TXT", f"{PDS}/SOFTWARE/SOFTINFO.TXT",
     "49bc5d395ef4b2bd2446baa64e9fe34fa90d44fa7ebd28c78def884f081f5e4e"),
    ("pds/ZIPINFO.TXT", f"{PDS}/SOFTWARE/ZIPINFO.TXT",
     "c2e976ac931f1f9784d7e7650bdca243d8366e50e7da9950453545324567d434"),
)


OBSERVATION = (
    "https://zenodo.org/api/records/22883077/files/CRaTER_GCR_LET_spectra.txt/content",
    "3e6be9ca1f8ac9b39653a5cec25bc5284b0e7d2ed8832b3dece77334c1e503a3",
)


# The starts are compressed-data offsets, not local-header offsets.  They were
# read from the immutable archive central directory and independently checked
# against each extracted member's SHA-256.  A Range response must be HTTP 206;
# accepting HTTP 200 here could accidentally download the entire archive.
GLACE_ARCHIVE = (
    "https://zenodo.org/api/records/8343472/files/dry_10cm_protons.zip/content"
)
RANGED_DEFLATE_ARTIFACTS = (
    ("glace/solmin_gcr.json", 371009884, 15725,
     "a2c2e0125d36c5e6a13313bad8ab2972c9c9cd6951b9180d660c35dc1ee63093"),
    ("glace/solmax_gcr.json", 370738852, 15718,
     "dff78c45e7beb58b1f291de1cfff9e4c6d661e9b3e76472a65cb3a95066b29d1"),
)


def digest(path: pathlib.Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def request(url: str, *, byte_range: tuple[int, int] | None = None):
    headers = {"User-Agent": "shieldSim-validation/1"}
    if byte_range is not None:
        start, length = byte_range
        headers["Range"] = f"bytes={start}-{start + length - 1}"
    return urllib.request.Request(url, headers=headers)


def install_bytes(target: pathlib.Path, payload: bytes, expected: str,
                  staging: pathlib.Path) -> None:
    staged = staging / target.name
    staged.write_bytes(payload)
    actual = digest(staged)
    if actual != expected:
        raise RuntimeError(
            f"downloaded bytes changed for {target}: expected {expected}, received {actual}"
        )
    if target.exists() and digest(target) != expected:
        raise RuntimeError(
            f"refusing to overwrite unreviewed local bytes at {target}; move them manually first"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    staged.replace(target)


def fetch_direct(target: pathlib.Path, url: str, expected: str,
                 staging: pathlib.Path) -> None:
    with urllib.request.urlopen(request(url), timeout=180) as response:
        payload = response.read()
    install_bytes(target, payload, expected, staging)


def fetch_embedded(target: pathlib.Path, start: int, compressed_size: int,
                   expected: str, staging: pathlib.Path) -> None:
    with urllib.request.urlopen(
        request(GLACE_ARCHIVE, byte_range=(start, compressed_size)), timeout=180
    ) as response:
        if response.status != 206:
            raise RuntimeError(
                "GLACE server did not honor the byte-range request; refusing a full archive download"
            )
        compressed = response.read(compressed_size + 1)
    if len(compressed) != compressed_size:
        raise RuntimeError(
            f"short/long GLACE range for {target.name}: expected {compressed_size} bytes, "
            f"received {len(compressed)}"
        )
    install_bytes(target, zlib.decompress(compressed, -zlib.MAX_WBITS), expected, staging)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=pathlib.Path, default=DEFAULT_DESTINATION)
    parser.add_argument(
        "--observation-destination", type=pathlib.Path, default=DEFAULT_OBSERVATION,
        help="destination for the published LET table shared with Stage-2 fixtures",
    )
    parser.add_argument(
        "--force-download", action="store_true",
        help="re-fetch valid files and prove that upstream bytes still match",
    )
    args = parser.parse_args()
    destination = args.destination.resolve()
    observation = args.observation_destination.resolve()

    # User policy requires generated/acquisition scratch state to remain below
    # test_output rather than hidden in the system /tmp tree.
    scratch_root = ROOT / "test_output"
    scratch_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="crater-stage8-download-", dir=scratch_root
    ) as temporary:
        staging = pathlib.Path(temporary)
        for relative, url, expected in ARTIFACTS:
            target = destination / relative
            if target.is_file() and digest(target) == expected and not args.force_download:
                print(f"VERIFIED {relative}")
                continue
            fetch_direct(target, url, expected, staging)
            print(f"INSTALLED {relative}")

        for relative, start, size, expected in RANGED_DEFLATE_ARTIFACTS:
            target = destination / relative
            if target.is_file() and digest(target) == expected and not args.force_download:
                print(f"VERIFIED {relative}")
                continue
            fetch_embedded(target, start, size, expected, staging)
            print(f"INSTALLED {relative}")

        observation_url, observation_sha = OBSERVATION
        if (observation.is_file() and digest(observation) == observation_sha
                and not args.force_download):
            print(f"VERIFIED {observation.relative_to(ROOT)}")
        else:
            fetch_direct(observation, observation_url, observation_sha, staging)
            print(f"INSTALLED {observation.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
