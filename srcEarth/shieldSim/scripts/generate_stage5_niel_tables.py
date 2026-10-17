#!/usr/bin/env python3
"""Generate the canonical Stage-5 NIEL CSVs from saved upstream responses.

This is an explicit fixture-maintenance tool, never a runtime dependency. The
SR-NIEL calculators are session-backed web forms, so the operator first saves
the three result pages produced with the manifest parameters. The script parses
only four-column numeric result rows, obtains the exact 640-group E192 edges
from a saved SAND2019-15194 PDF via ``pdftotext -layout``, and emits stable CSV
bytes. It refuses unexpected row counts, energy order, representative energies,
or hashes instead of guessing how a changed upstream page should be read.

Normal builds/tests consume the committed CSVs offline through the strict C++
loader. Updating a fixture requires review of the raw hashes, generated diff,
manifest provenance/version, upstream terms, and the entire Stage-5 test suite.
"""

from __future__ import annotations

import argparse
import hashlib
from html.parser import HTMLParser
import json
import math
import pathlib
import re
import subprocess


REFERENCE_NIEL = 95.0e-27 * 6.02214076e23 / 28.0855


class NumericRows(HTMLParser):
    """Collect HTML table rows and retain only four finite numeric cells."""

    def __init__(self) -> None:
        super().__init__()
        self.in_row = False
        self.in_cell = False
        self.cell: list[str] = []
        self.row: list[str] = []
        self.rows: list[tuple[float, float, float, float]] = []

    def handle_starttag(self, tag: str,
                        attrs: list[tuple[str, str | None]]) -> None:
        del attrs
        if tag.lower() == "tr":
            self.in_row = True
            self.row = []
        elif self.in_row and tag.lower() == "td":
            self.in_cell = True
            self.cell = []

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.cell.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag.lower() == "td" and self.in_cell:
            self.row.append("".join(self.cell).strip())
            self.in_cell = False
        elif tag.lower() == "tr" and self.in_row:
            if len(self.row) == 4:
                try:
                    values = tuple(float(value) for value in self.row)
                except ValueError:
                    values = ()
                if len(values) == 4 and all(math.isfinite(value) for value in values):
                    self.rows.append(values)  # type: ignore[arg-type]
            self.in_row = False


def parse_response(path: pathlib.Path) -> list[tuple[float, float]]:
    parser = NumericRows()
    parser.feed(path.read_text(encoding="latin-1"))
    # Columns 3 and 4 repeat unit-fluence dose and its Gy conversion; the
    # response table needs only energy and NIEL. Sorting is explicit because the
    # neutron result is returned in descending energy.
    rows = sorted((row[0], row[1]) for row in parser.rows)
    if not rows or any(energy <= 0.0 or niel < 0.0 for energy, niel in rows):
        raise RuntimeError(f"{path} contains no valid positive-energy response table")
    if any(rows[index][0] >= rows[index + 1][0] for index in range(len(rows) - 1)):
        raise RuntimeError(f"{path} result energies are duplicated or nonmonotonic")
    return rows


def parse_e192_edges(pdf: pathlib.Path) -> list[float]:
    completed = subprocess.run(
        ["pdftotext", "-layout", str(pdf), "-"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"pdftotext failed: {completed.stderr.strip()}")
    text = completed.stdout
    start = text.find("C 640 Groups")
    if start < 0:
        raise RuntimeError("SAND report lacks the expected 'C 640 Groups' marker")
    start = text.find("E192", start)
    end = text.find("C   Multipliers converted from E722-19", start)
    if start < 0 or end < 0:
        raise RuntimeError("could not isolate the Appendix-B E192 edge block")
    # Requiring six mantissa decimals avoids accidentally interpreting PDF page
    # numbers or prose as group edges.
    tokens = re.findall(r"\b\d\.\d{6}E[+-]\d{2}\b", text[start:end], re.I)
    edges = [float(token) for token in tokens]
    if len(edges) != 641 or edges[0] != 1.0e-10 or edges[-1] != 20.0:
        raise RuntimeError(f"expected 641 E192 edges from 1e-10 to 20 MeV, found {len(edges)}")
    if any(edges[index] >= edges[index + 1] for index in range(640)):
        raise RuntimeError("E192 edges are not strictly increasing")
    return edges


def write_charged(path: pathlib.Path, rows: list[tuple[float, float]],
                  expected_count: int) -> None:
    positive = [(energy, niel) for energy, niel in rows if niel > 0.0]
    if len(positive) != expected_count:
        raise RuntimeError(
            f"{path.name}: expected {expected_count} positive nodes, found {len(positive)}"
        )
    lines = ["energy_MeV,niel_MeV_cm2_g"]
    lines.extend(f"{energy:.12e},{niel:.12e}" for energy, niel in positive)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_neutron(path: pathlib.Path, rows: list[tuple[float, float]],
                  edges: list[float]) -> None:
    if len(rows) != 640:
        raise RuntimeError(f"expected 640 neutron response values, found {len(rows)}")
    lines = [
        "lower_MeV,upper_MeV,representative_MeV,niel_MeV_cm2_g,"
        "hardness_factor_to_1MeV_neutron"
    ]
    for index, ((representative, niel), lower, upper) in enumerate(
            zip(rows, edges[:-1], edges[1:])):
        midpoint = 0.5 * (lower + upper)
        if not math.isclose(representative, midpoint, rel_tol=5.0e-12, abs_tol=1.0e-300):
            raise RuntimeError(
                f"neutron row {index} representative {representative} != edge midpoint {midpoint}"
            )
        lines.append(
            f"{lower:.12e},{upper:.12e},{representative:.12e},{niel:.12e},"
            f"{niel / REFERENCE_NIEL:.12e}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def sha256(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proton-html", type=pathlib.Path, required=True)
    parser.add_argument("--alpha-html", type=pathlib.Path, required=True)
    parser.add_argument("--neutron-html", type=pathlib.Path, required=True)
    parser.add_argument("--sandia-pdf", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument(
        "--expected-manifest", type=pathlib.Path,
        help="optional manifest whose table hashes must match the generated bytes",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    proton_path = args.output / "sr_niel_v11_1_si_proton.csv"
    alpha_path = args.output / "sr_niel_v11_1_si_alpha.csv"
    neutron_path = args.output / "astm_e722_19_si_neutron_groups.csv"
    write_charged(proton_path, parse_response(args.proton_html), 143)
    write_charged(alpha_path, parse_response(args.alpha_html), 133)
    write_neutron(neutron_path, parse_response(args.neutron_html),
                  parse_e192_edges(args.sandia_pdf))

    generated = {path.name: sha256(path) for path in (proton_path, alpha_path, neutron_path)}
    if args.expected_manifest:
        manifest = json.loads(args.expected_manifest.read_text(encoding="utf-8"))
        expected = {
            table["path"]: (table["sha256"][7:]
                            if table["sha256"].startswith("sha256:")
                            else table["sha256"])
            for table in manifest["tables"]
        }
        if generated != expected:
            raise RuntimeError(
                "generated table hashes differ from the expected manifest:\n"
                f"generated={generated}\nexpected={expected}"
            )
    print(json.dumps(generated, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
