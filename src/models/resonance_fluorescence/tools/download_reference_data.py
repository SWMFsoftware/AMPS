#!/usr/bin/env python3
"""Explicitly refresh upstream resonance-fluorescence reference data.

This maintenance tool is never used at AMPS runtime.  It refuses to replace
the reviewed snapshot unless ``--replace-reviewed-data`` is given.  After a
refresh, review the scientific diff, update manifest hashes, regenerate the
compiled include, and run ``make test``; a successful download alone does not
authorize a reference-data change.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from decimal import Decimal
import io
import json
from pathlib import Path
import tarfile
import tempfile
import time
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen


MODULE = Path(__file__).resolve().parents[1]
REFERENCE = MODULE / "reference" / "g_factors"
H_URL = "https://users.cbk.waw.pl/~ikowalska/LyaModel/files/table1v4.txt"
LOWELL_URL = "https://asteroid.lowell.edu/api/comet/calc/gfactor"
H_PATH = REFERENCE / "h_lya_ikl2020" / "table1v4.txt"
OH_DIR = REFERENCE / "oh_ax00_schleicher_ahearn1988"
OH_RAW = OH_DIR / "lowell_oh00_raw.jsonl"
OH_CSV = OH_DIR / "lowell_oh00_velocity_gfactor_1au.csv"

# Atomic oxygen 130.4 nm: Killen et al. (2022) Figure-14 data, the
# SOHO/SUMER disk-centre atlas served by BASS2000, and NIST ASD line data.
O_KVB_URL = (
    "https://content.cld.iop.org/journals/0067-0049/263/2/37/revision1/"
    "apjsac9eabf14.tar.gz"
)
O_KVB_DIR = REFERENCE / "o_i_1304_kvb2022"
O_KVB_MEMBERS = (
    "oI_1302_gvalues.txt",
    "oI_1304_gvalues.txt",
    "oI_1306_gvalues.txt",
    "ReadMe_dbf14.txt",
)
O_SUMER_DIR = REFERENCE / "o_i_sumer2001_bass2000"
O_SUMER_URL = "https://bass2000.obspm.fr/php/getSolarSpectrumDB.php"
O_SUMER_WINDOWS = {
    "sumer_disk_center_1297_1311A.txt": {"WL": "1297", "DW": "14"},
    "sumer_disk_center_1019_1032A.txt": {"WL": "1019", "DW": "13"},
}
O_NIST_DIR = REFERENCE / "o_i_nist_asd"
O_NIST_URL = "https://physics.nist.gov/cgi-bin/ASD/lines1.pl"
O_NIST_QUERIES = {
    "o_i_lines_1020_1310A_vacuum.csv": ("1020", "1310"),
    "o_i_lines_11270_11310A_air.csv": ("11270", "11310"),
}

# Retain the historical convenience column byte-for-byte.  It is not an input
# to the core evaluator; applications must state their own effective photon
# wavelength when converting the authoritative OH energy rate.
PLANCK_J_S = 6.62607015e-34
LIGHT_SPEED_M_S = 299792458.0
OH_CONVENIENCE_WAVELENGTH_M = 3090.0e-10


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as stream:
        stream.write(text)
        temporary = Path(stream.name)
    temporary.replace(path)


def _download_h() -> None:
    request = Request(H_URL, headers={"User-Agent": "AMPS-reference-data/1.0"})
    with urlopen(request, timeout=30) as response:
        text = response.read().decode("utf-8")
    required = {"AK", "muK", "sigmaK", "kappa", "AR", "dmu", "sigmaR", "bbkg", "abkg"}
    found = {
        fields[0]
        for line in text.splitlines()
        if len(fields := line.split()) == 3 and fields[0] in required
    }
    if found != required:
        raise RuntimeError(f"IKL table is incomplete; found {sorted(found)}")
    _atomic_text(H_PATH, text)


def _fetch_lowell(velocity: Decimal, attempts: int = 4) -> dict[str, Any]:
    url = LOWELL_URL + "?" + urlencode({"r": "1", "v": str(velocity)})
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": "AMPS-reference-data/1.0"})
            with urlopen(request, timeout=30) as response:
                payload = json.load(response)
            one_au = payload["gfactor"]["one_AU"]
            if payload.get("extrapolated") is not False:
                raise RuntimeError(f"Lowell extrapolated velocity {velocity} km/s")
            if abs(float(one_au["v"]) - float(velocity)) > 5.0e-10:
                raise RuntimeError(f"Lowell returned a mismatched velocity for {velocity}")
            if float(one_au["OH"]) <= 0.0:
                raise RuntimeError(f"Lowell returned non-positive OH g at {velocity}")
            return payload
        except (OSError, TimeoutError, json.JSONDecodeError) as error:
            last_error = error
            if attempt + 1 < attempts:
                time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"Lowell request failed at {velocity}: {last_error}")


def _download_oh(workers: int) -> None:
    velocities = [Decimal("-60") + Decimal("0.1") * index for index in range(1201)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        payloads = list(pool.map(_fetch_lowell, velocities))

    raw = "".join(
        json.dumps(
            {"requested_velocity_km_s": str(velocity), "response": payload},
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
        for velocity, payload in zip(velocities, payloads)
    )
    _atomic_text(OH_RAW, raw)

    photon_energy_erg = (
        PLANCK_J_S * 1.0e7 * LIGHT_SPEED_M_S / OH_CONVENIENCE_WAVELENGTH_M
    )
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=OH_DIR, delete=False
    ) as stream:
        fields = (
            "heliocentric_radial_velocity_km_s",
            "oh00_energy_g_1au_erg_s-1_molecule-1",
            "oh00_photon_g_1au_s-1_at_3090A",
            "api_extrapolated",
            "source",
        )
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for velocity, payload in zip(velocities, payloads):
            energy = float(payload["gfactor"]["one_AU"]["OH"])
            writer.writerow(
                {
                    fields[0]: str(velocity),
                    fields[1]: f"{energy:.12e}",
                    fields[2]: f"{energy / photon_energy_erg:.12e}",
                    fields[3]: "false",
                    fields[4]: LOWELL_URL,
                }
            )
        temporary = Path(stream.name)
    temporary.replace(OH_CSV)


def _get_bytes(url: str, timeout: int = 60) -> bytes:
    request = Request(url, headers={"User-Agent": "AMPS-reference-data/1.0"})
    with urlopen(request, timeout=timeout) as response:
        return response.read()


def _atomic_bytes(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=path.parent, delete=False) as stream:
        stream.write(data)
        temporary = Path(stream.name)
    temporary.replace(path)


def _download_o() -> None:
    archive = _get_bytes(O_KVB_URL)
    members: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:gz") as tar:
        for member in tar.getmembers():
            # Untrusted archive: accept only the expected regular files.
            if not member.isfile() or member.name not in O_KVB_MEMBERS:
                raise RuntimeError(f"unexpected Killen et al. archive member {member.name!r}")
            members[member.name] = tar.extractfile(member).read()
    if set(members) != set(O_KVB_MEMBERS):
        raise RuntimeError(f"incomplete Killen et al. archive: {sorted(members)}")
    _atomic_bytes(O_KVB_DIR / "apjsac9eabf14.tar.gz", archive)
    for name, data in members.items():
        _atomic_bytes(O_KVB_DIR / name, data)

    for name, window in O_SUMER_WINDOWS.items():
        query = urlencode({**window, "resol": "0.001", "fmt": "txt"})
        data = _get_bytes(O_SUMER_URL + "?" + query, timeout=120)
        if not data.startswith(b"Wavelength (Angstroms), Normalized intensity"):
            raise RuntimeError(f"unexpected BASS2000 response for {name}")
        _atomic_bytes(O_SUMER_DIR / name, data)

    for name, (low, high) in O_NIST_QUERIES.items():
        query = urlencode(
            {
                "spectra": "O I", "limits_type": "0", "low_w": low, "upp_w": high,
                "unit": "0", "de": "0", "format": "2", "line_out": "0",
                "remove_js": "on", "en_unit": "0", "output": "0", "bibrefs": "1",
                "page_size": "15", "show_obs_wl": "1", "show_calc_wl": "1",
                "unc_out": "1", "order_out": "0", "show_av": "2", "tsb_value": "0",
                "A_out": "0", "f_out": "on", "S_out": "on", "loggf_out": "on",
                "intens_out": "on", "allowed_out": "1", "forbid_out": "1",
                "conf_out": "on", "term_out": "on", "enrg_out": "on", "J_out": "on",
                "g_out": "on",
            }
        )
        data = _get_bytes(O_NIST_URL + "?" + query)
        if b"Aki(s^-1)" not in data.splitlines()[0]:
            raise RuntimeError(f"unexpected NIST ASD response for {name}")
        _atomic_bytes(O_NIST_DIR / name, data)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h-only", action="store_true")
    parser.add_argument("--oh-only", action="store_true")
    parser.add_argument("--o-only", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--replace-reviewed-data", action="store_true")
    args = parser.parse_args()
    if sum((args.h_only, args.oh_only, args.o_only)) > 1:
        parser.error("--h-only, --oh-only and --o-only are mutually exclusive")
    if not 1 <= args.workers <= 8:
        parser.error("--workers must be in 1..8")
    if not args.replace_reviewed_data:
        parser.error("refusing to replace reviewed data without --replace-reviewed-data")
    everything = not (args.h_only or args.oh_only or args.o_only)
    if everything or args.h_only:
        _download_h()
    if everything or args.oh_only:
        _download_oh(args.workers)
    if everything or args.o_only:
        _download_o()
    print("Reference files refreshed. Review diffs and hashes before regeneration.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
