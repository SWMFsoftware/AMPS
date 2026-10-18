#!/usr/bin/env python3
"""Download and canonicalize the official NIST Stage-4 EM reference tables.

The downloader stores the unmodified response bytes as well as parsed CSVs and
a SHA-256 manifest.  Validation uses only the pinned local files; a normal test
run never depends on network availability.  Re-running this script is an
explicit fixture-update operation and any upstream change is visible as a raw
and canonical checksum change in review.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
from html.parser import HTMLParser
import json
import pathlib
import re
import urllib.parse
import urllib.request


STAR_URL = "https://physics.nist.gov/cgi-bin/Star/apdata.pl"
XRAY_ROOT = "https://physics.nist.gov/PhysRefData/XrayMassCoef"
XCOM_ROOT = "https://physics.nist.gov/cgi-bin/Xcom"
MATERIALS = {
    "aluminum": {"star_number": "013", "xray": "ElemTab/z13.html", "atomic_number": "13"},
    "silicon": {"star_number": "014", "xray": "ElemTab/z14.html", "atomic_number": "14"},
    # Germanium is an exact elemental match to the absorber-library record.
    # It is present in both PSTAR/ASTAR and the elemental photon tables.
    "germanium": {"star_number": "032", "xray": "ElemTab/z32.html", "atomic_number": "32"},
    # NIST PSTAR/ASTAR material 245 is the stoichiometric silicon-dioxide
    # reference.  The X-Ray Mass Attenuation compilation does not publish a
    # matching pre-composed SiO2 page, so only its STAR and formula-driven XCOM
    # records are requested.  We do not substitute glass or elemental silicon.
    "silicon_dioxide": {
        "star_number": "245", "formula": "SiO2", "xcom_name": "silica"
    },
    "water": {"star_number": "276", "xray": "ComTab/water.html", "formula": "H2O"},
}


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(url: str, payload: dict[str, str] | None = None,
          referer: str | None = None) -> bytes:
    data = urllib.parse.urlencode(payload).encode("ascii") if payload else None
    headers = {
            "User-Agent": "shieldSim-stage4-reference-fixture/1.0",
            "Content-Type": "application/x-www-form-urlencoded",
        }
    if referer:
        headers["Referer"] = referer
    request = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def parse_star(raw: bytes) -> list[list[float]]:
    rows: list[list[float]] = []
    number = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:E[+-]?\d+)?"
    pattern = re.compile(rf"^\s*({number})(?:\s+({number})){{6}}\s*$", re.I)
    for line in raw.decode("ascii").splitlines():
        if not pattern.match(line):
            continue
        values = [float(value) for value in line.split()]
        if len(values) == 7:
            rows.append(values)
    if len(rows) < 100:
        raise RuntimeError("NIST STAR response did not contain the expected numeric table")
    return rows


class XrayTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_row = False
        self.in_cell = False
        self.cell_text: list[str] = []
        self.row: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() == "tr":
            self.in_row = True
            self.row = []
        elif self.in_row and tag.lower() == "td":
            self.in_cell = True
            self.cell_text = []

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.cell_text.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "td" and self.in_cell:
            self.row.append("".join(self.cell_text).strip())
            self.in_cell = False
        elif tag.lower() == "tr" and self.in_row:
            if self.row:
                self.rows.append(self.row)
            self.in_row = False


def parse_xray(raw: bytes) -> list[tuple[str, float, float, float]]:
    parser = XrayTableParser()
    parser.feed(raw.decode("latin-1"))
    result: list[tuple[str, float, float, float]] = []
    for row in parser.rows:
        compact = [cell.replace("\xa0", " ").strip() for cell in row if cell.strip()]
        edge = ""
        if len(compact) == 4 and compact[0] in {"K", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5"}:
            edge, values = compact[0], compact[1:]
        elif len(compact) == 3:
            values = compact
        else:
            continue
        try:
            energy, attenuation, energy_absorption = map(float, values)
        except ValueError:
            continue
        result.append((edge, energy, attenuation, energy_absorption))
    if len(result) < 30:
        raise RuntimeError("NIST X-Ray Mass Attenuation response did not contain the expected table")
    return result


def parse_xcom(raw: bytes) -> list[tuple[float, float, float]]:
    """Return energy, total-with-coherent, total-without-coherent rows."""

    parser = XrayTableParser()
    parser.feed(raw.decode("latin-1"))
    result: list[tuple[float, float, float]] = []
    for row in parser.rows:
        compact = [cell.replace("\xa0", " ").strip() for cell in row if cell.strip()]
        # XCOM PIC output has energy + five partial coefficients + two totals.
        if len(compact) != 8:
            continue
        try:
            values = [float(value) for value in compact]
        except ValueError:
            continue
        result.append((values[0], values[6], values[7]))
    if len(result) != 3:
        raise RuntimeError("NIST XCOM response did not contain the three requested energies")
    return result


def write_bytes(path: pathlib.Path, data: bytes, manifest: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    manifest[path.as_posix()] = {"sha256": digest(data), "size_bytes": len(data)}


def write_csv(path: pathlib.Path, header: list[str], rows: list[list[object]],
              manifest: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    data = path.read_bytes()
    manifest[path.as_posix()] = {"sha256": digest(data), "size_bytes": len(data)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("tests/reference_tid/nist"),
        help="fixture directory relative to the current working directory",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    manifest: dict[str, dict] = {}

    for material, metadata in MATERIALS.items():
        for program in ("PSTAR", "ASTAR"):
            payload = {
                "prog": program,
                "matno": metadata["star_number"],
                "ShowDefault": "on",
                "character": "space",
                "electronic": "on",
                "nuclear": "on",
                "total": "on",
                "csda": "on",
                "project": "on",
                "detour": "on",
            }
            raw = fetch(STAR_URL, payload)
            raw_path = output / "raw" / f"{program.lower()}_{material}.txt"
            write_bytes(raw_path, raw, manifest)
            rows = parse_star(raw)
            write_csv(
                output / f"{program.lower()}_{material}.csv",
                [
                    "energy_MeV",
                    "electronic_stopping_MeV_cm2_g",
                    "nuclear_stopping_MeV_cm2_g",
                    "total_stopping_MeV_cm2_g",
                    "csda_range_g_cm2",
                    "projected_range_g_cm2",
                    "detour_factor",
                ],
                rows,
                manifest,
            )

        # Not every PSTAR/ASTAR material has a corresponding tabulated X-Ray
        # Mass Attenuation page.  Absence is explicit in MATERIALS; using a
        # nearby composition as a silent proxy would invalidate the reference.
        if "xray" in metadata:
            xray_url = f"{XRAY_ROOT}/{metadata['xray']}"
            raw = fetch(xray_url)
            write_bytes(output / "raw" / f"xray_{material}.html", raw, manifest)
            xray_rows = parse_xray(raw)
            write_csv(
                output / f"xray_{material}.csv",
                ["absorption_edge", "energy_MeV", "mass_attenuation_cm2_g", "mass_energy_absorption_cm2_g"],
                [list(row) for row in xray_rows],
                manifest,
            )

        xcom_common = {
            "OutOpt": "PIC",
            "Graph6": "on",
            "NumAdd": "1",
            "Energies": "0.1\n1.0\n10.0",
            "WindowXmin": "0.001",
            "WindowXmax": "100000",
            "ResizeFlag": "on",
        }
        if "atomic_number" in metadata:
            xcom_url = f"{XCOM_ROOT}/xcom3_1"
            xcom_payload = {**xcom_common, "ZNum": metadata["atomic_number"]}
            xcom_referer = f"{XCOM_ROOT}/xcom2?Method=Elem&Output2=Hand"
        else:
            xcom_url = f"{XCOM_ROOT}/xcom3_2"
            # XCOM restricts its optional title to a simple token; repository
            # identifiers contain underscores and therefore cannot be sent as
            # the title.  Keep the physics formula separate from that label.
            xcom_payload = {
                **xcom_common,
                "Formula": metadata["formula"],
                "Name": metadata.get("xcom_name", material),
            }
            xcom_referer = f"{XCOM_ROOT}/xcom2?Method=Comp&Output2=Hand"
        xcom_raw = fetch(xcom_url, xcom_payload, xcom_referer)
        write_bytes(output / "raw" / f"xcom_{material}.html", xcom_raw, manifest)
        write_csv(
            output / f"xcom_{material}.csv",
            ["energy_MeV", "mass_attenuation_with_coherent_cm2_g",
             "mass_attenuation_without_coherent_cm2_g"],
            [list(row) for row in parse_xcom(xcom_raw)],
            manifest,
        )

    # Paths in the manifest are fixture-relative so the checksum record remains
    # identical when the repository is checked out elsewhere.
    relative = {
        str(pathlib.Path(path).relative_to(output)): value
        for path, value in sorted(manifest.items())
    }
    document = {
        "schema_version": "shieldSim.nist_stage4_manifest.v1",
        "retrieval_utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        "sources": {
            "PSTAR_ASTAR": STAR_URL,
            "XCOM": XCOM_ROOT,
            "xray_mass_attenuation": XRAY_ROOT,
        },
        "files": relative,
    }
    (output / "manifest.sha256.json").write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {len(relative)} pinned files under {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
