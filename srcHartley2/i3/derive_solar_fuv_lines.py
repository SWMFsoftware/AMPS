#!/usr/bin/env python3
"""Derive the per-case 1-AU solar O I 130.4-nm and H Lyman-beta photon fluxes.

The O I 130.4-nm brightness models of ``src/models/resonance_fluorescence``
need two case inputs (see that module's README):

* the complete O I 130.4-nm triplet photon flux at 1 AU on the SORCE SOLSTICE
  absolute scale of Killen et al. (2022), for the two direct models;
* the 1-AU photon flux in 102.4--102.8 nm (H Lyman-beta), for the Bowen model.

Source: LASP/LISIRD FISM2 daily spectra (dataset ``fism_daily_hr``, 0.1-nm
bins).  Each case uses the same Earth-equivalent exposure interval as the
F10.7P and Lyman-alpha inputs (``reference/i3_cases.csv``); daily values are
weighted by their overlap with that interval.  The triplet flux is FISM2's
130.1--130.7 nm flux divided by its 2009-03-01 value -- the SOLSTICE day of
Killen et al. -- times the reference triplet flux that their tables imply
(module manifest).  The ratio cancels the FISM2/SOLSTICE calibration offset.
The Lyman-beta flux is the absolute FISM2 flux of the four bins that the
Bowen model's SUMER profile is normalized to.

By default the frozen snapshot in ``reference/fism2`` is used (no network);
``--download`` refreshes it from LISIRD first.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen


I3 = Path(__file__).resolve().parent
REPOSITORY = I3.parents[1]
CASES = I3 / "reference" / "i3_cases.csv"
SNAPSHOT = I3 / "reference" / "fism2"
OUTPUT = I3 / "reference" / "solar_fuv_lines.csv"
MANIFEST = (
    REPOSITORY / "src/models/resonance_fluorescence/reference/g_factors/manifest.json"
)
LISIRD = "https://lasp.colorado.edu/lisird/latis/dap/fism_daily_hr.csv"

PLANCK_J_S = 6.62607015e-34
LIGHT_SPEED_M_S = 299792458.0
REFERENCE_DAY = "2009060"  # 2009-03-01, the SORCE SOLSTICE day of Killen et al.
# (name, first, last bin centre [nm], photon wavelength [nm]); FISM2 bins are
# 0.1 nm wide and centred on x.x5 nm.
BANDS = {
    "oi": (130.15, 130.65, 130.4),
    "lyb": (102.45, 102.75, 102.572),
}
SNAPSHOT_QUERIES = {
    "fism2_daily_hr_1300_1308A_20251001_20251031.csv": ("2025-10-01", "2025-11-01", 129.9, 130.9),
    "fism2_daily_hr_1300_1308A_20090301.csv": ("2009-03-01", "2009-03-02", 129.9, 130.9),
    "fism2_daily_hr_1024_1028A_20251001_20251031.csv": ("2025-10-01", "2025-11-01", 102.4, 102.8),
    "fism2_daily_hr_1024_1028A_20090301.csv": ("2009-03-01", "2009-03-02", 102.4, 102.8),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def download() -> None:
    SNAPSHOT.mkdir(parents=True, exist_ok=True)
    for name, (start, stop, low, high) in SNAPSHOT_QUERIES.items():
        query = (
            f"time>={start}&time<{stop}&wavelength>={low}&wavelength<={high}"
        )
        request = Request(
            LISIRD + "?" + quote(query, safe="=&<>"),
            headers={"User-Agent": "AMPS-reference-data/1.0"},
        )
        with urlopen(request, timeout=120) as response:
            data = response.read()
        if not data.startswith(b"time (yyyyDDD),wavelength (nm),irradiance (W/m^2/nm)"):
            raise RuntimeError(f"unexpected LISIRD response for {name}")
        (SNAPSHOT / name).write_bytes(data)


def load_snapshot() -> dict[str, dict[float, float]]:
    spectra: dict[str, dict[float, float]] = {}
    for name in SNAPSHOT_QUERIES:
        with (SNAPSHOT / name).open(newline="") as stream:
            reader = csv.reader(stream)
            next(reader)
            for day, wavelength, irradiance, _ in reader:
                spectra.setdefault(day, {})[round(float(wavelength), 2)] = float(irradiance)
    return spectra


def band_photon_flux(spectrum: dict[float, float], band: str) -> float:
    """Photons cm^-2 s^-1 at 1 AU in the FISM2 bins of ``band``."""
    first, last, photon_nm = BANDS[band]
    centres = [round(first + 0.1 * index, 2) for index in range(round((last - first) / 0.1) + 1)]
    missing = [centre for centre in centres if centre not in spectrum]
    if missing:
        raise RuntimeError(f"FISM2 spectrum lacks bins {missing}")
    watts = sum(spectrum[centre] for centre in centres) * 0.1
    return watts / (PLANCK_J_S * LIGHT_SPEED_M_S / (photon_nm * 1.0e-9)) * 1.0e-4


def utc(text: str) -> dt.datetime:
    return dt.datetime.fromisoformat(text.replace("Z", "+00:00"))


def interval_mean(spectra, band: str, start: str, end: str) -> float:
    begin, finish = utc(start), utc(end)
    day = begin.replace(hour=0, minute=0, second=0, microsecond=0)
    weighted = seconds = 0.0
    while day < finish:
        following = day + dt.timedelta(days=1)
        overlap = (min(following, finish) - max(day, begin)).total_seconds()
        if overlap > 0.0:
            key = day.strftime("%Y%j")
            if key not in spectra:
                raise RuntimeError(f"FISM2 snapshot lacks day {key}")
            weighted += overlap * band_photon_flux(spectra[key], band)
            seconds += overlap
        day = following
    return weighted / seconds


def derive() -> list[dict[str, str]]:
    spectra = load_snapshot()
    manifest = json.loads(MANIFEST.read_text())
    reference_triplet = manifest["models"]["O_I_1304_KVB2022_DIRECT_J2"][
        "reference_triplet_flux_1au_photon_cm-2_s-1"
    ]
    oi_reference = band_photon_flux(spectra[REFERENCE_DAY], "oi")
    rows = []
    with CASES.open(newline="") as stream:
        for case in csv.DictReader(stream):
            start = case["earth_equivalent_start_utc"]
            end = case["earth_equivalent_end_utc"]
            oi = interval_mean(spectra, "oi", start, end)
            lyb = interval_mean(spectra, "lyb", start, end)
            rows.append(
                {
                    "case_id": case["case_id"],
                    "earth_equivalent_start_utc": start,
                    "earth_equivalent_end_utc": end,
                    "fism2_130p1_130p7nm_photon_flux_1au_ph_cm-2_s-1": f"{oi:.10e}",
                    "fism2_130p1_130p7nm_ratio_to_20090301": f"{oi / oi_reference:.10f}",
                    "o_i_1304_triplet_photon_flux_1au_ph_cm-2_s-1": (
                        f"{reference_triplet * oi / oi_reference:.10e}"
                    ),
                    "h_lyb_102p4_102p8nm_photon_flux_1au_ph_cm-2_s-1": f"{lyb:.10e}",
                }
            )
    return rows


def render(rows: list[dict[str, str]]) -> str:
    fields = list(rows[0])
    lines = [",".join(fields)]
    lines.extend(",".join(row[field] for field in fields) for row in rows)
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true",
                        help="refresh the FISM2 snapshot from LISIRD first")
    parser.add_argument("--check", action="store_true",
                        help="verify that the tracked CSV is current")
    args = parser.parse_args()
    if args.download:
        download()
    for name in SNAPSHOT_QUERIES:
        print(f"{sha256(SNAPSHOT / name)}  reference/fism2/{name}")
    content = render(derive())
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != content:
            raise SystemExit(f"{OUTPUT} is stale")
        return 0
    OUTPUT.write_text(content)
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
