#!/usr/bin/env python3
"""Derive the two pinned inputs of the I24 total-capture provider.

Inputs (both hash-checked against reference/total_capture_case.json):
  * Shingal et al. 1986, J. Phys. B 19, 309, Table 4 (printed page 318,
    PDF page 11): A34 partial capture cross sections for
    H+ + Na(3s) -> H(nl) + Na+, in A^2, versus laboratory energy in keV/u.
    The column order E, 1s, 2s, 2p, 3s, 3p, 3d was checked visually on the
    page image; the PDF text layer loses the sigma subscripts.
  * The retained NIST ASD H I level response from acquire_hydrogen_levels.py.

Outputs (ASCII, LF, native decimal tokens only; no arithmetic is applied):
  * shingal_table4.csv    - every printed row of Table 4;
  * hydrogen_levels.csv   - the nlj levels with n <= 3 and their statistical
                            weights g = 2J + 1.

The production provider pins the SHA-256 of both outputs. Each test run
regenerates them and requires byte equality with reference/. The approval
record must be APPROVED; otherwise the case is SKIPPED.
"""
import argparse
import csv
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CASE = HERE / "reference/total_capture_case.json"
APPROVAL = HERE / "reference/total_capture_approval_20261125.json"
TABLE_HEADER = ("energy_kev_u", "sigma_1s_A2", "sigma_2s_A2", "sigma_2p_A2",
                "sigma_3s_A2", "sigma_3p_A2", "sigma_3d_A2")
LEVEL_HEADER = ("n", "l", "configuration", "term", "J", "g", "level_eV", "uncertainty_eV", "reference")


class NotApproved(Exception):
    """The user approval record is missing or not APPROVED."""


def sha(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def pdf_text(path, page):
    """Read the existing PDF text layer; no OCR and no figure digitization."""
    result = subprocess.run(["pdftotext", "-layout", "-f", str(page), "-l", str(page), str(path), "-"],
                            capture_output=True, encoding="utf-8", check=True, timeout=30)
    return result.stdout


def shingal_rows(path):
    text = pdf_text(path, 11)
    parts = text.split("Table 4. Computed (A34) partial cross sections", 1)
    if len(parts) != 2:
        raise ValueError("Shingal Table 4 caption not found on PDF page 11")
    table = parts[1].split("Figure 5.", 1)[0]
    rows = re.findall(r"^ *(\d+\.\d+)((?: +\d+\.\d+){6}) *$", table, re.M)
    rows = [(energy, *cells.split()) for energy, cells in rows]
    if [r[0] for r in rows] != ["0.5", "1.0", "4.0", "10.0", "20.0"]:
        raise ValueError("unexpected Shingal Table 4 energy rows: %r" % [r[0] for r in rows])
    return rows


def asd_field(value):
    """Undo the spreadsheet quoting of the NIST CSV export, e.g. ="1s"."""
    match = re.fullmatch(r'="(.*)"', value)
    return match.group(1).strip() if match else value.strip()


def hydrogen_rows(path):
    text = Path(path).read_text(encoding="utf-8")
    reader = csv.DictReader(io.StringIO(text))
    rows = []
    for raw in reader:
        configuration = asd_field(raw["Configuration"])
        match = re.fullmatch(r"([1-3])([spd])", configuration)
        if not match:
            continue
        J = asd_field(raw["J"])
        g = int(raw["g"])
        if Fraction(J) * 2 + 1 != g:
            raise ValueError("NIST statistical weight differs from 2J+1 for " + configuration)
        if asd_field(raw["Prefix"]) or asd_field(raw["Suffix"]):
            raise ValueError("bracketed or estimated NIST level for n <= 3: " + configuration)
        rows.append((match.group(1), "spd".index(match.group(2)), configuration, asd_field(raw["Term"]),
                     J, g, asd_field(raw["Level (eV)"]), asd_field(raw["Uncertainty (eV)"]),
                     asd_field(raw["Reference"])))
    # 1s1/2; 2s1/2, 2p1/2, 2p3/2; 3s1/2, 3p1/2, 3p3/2, 3d3/2, 3d5/2.
    if sorted(r[2] for r in rows) != ["1s", "2p", "2p", "2s", "3d", "3d", "3p", "3p", "3s"]:
        raise ValueError("unexpected NIST H I level set for n <= 3")
    return sorted(rows, key=lambda r: (int(r[0]), r[1], Fraction(r[4])))


def render(header, rows):
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerows([header, *rows])
    return stream.getvalue().encode("ascii")


def generate(case=None):
    """Return {filename: bytes} after authenticating the approval and sources."""
    case = case or json.loads(CASE.read_text())
    approval = json.loads(APPROVAL.read_text()) if APPROVAL.is_file() else {}
    if approval.get("status") != "APPROVED":
        raise NotApproved("I24 total-capture approval record is not APPROVED")
    root = HERE.parents[3]
    sources = {}
    for name in ("shingal1986", "nist_asd_h_i"):
        entry = case["sources"][name]
        path = (root / entry["path"]).resolve()
        if not path.is_file():
            raise FileNotFoundError("required source is unavailable: " + str(path))
        if sha(path) != entry["sha256"]:
            raise ValueError("changed source identity: " + str(path))
        sources[name] = path
    return {"shingal_table4.csv": render(TABLE_HEADER, shingal_rows(sources["shingal1986"])),
            "hydrogen_levels.csv": render(LEVEL_HEADER, hydrogen_rows(sources["nist_asd_h_i"]))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="write the outputs to this new directory")
    parser.add_argument("--check", action="store_true", help="require byte equality with reference/")
    args = parser.parse_args()
    try:
        outputs = generate()
    except NotApproved as error:
        print(str(error), file=sys.stderr)
        return 77
    if args.output_dir:
        args.output_dir.mkdir(parents=True, exist_ok=False)
        for name, data in outputs.items():
            (args.output_dir / name).write_bytes(data)
    report = {name: sha256(data).hexdigest() for name, data in outputs.items()}
    if args.check:
        frozen = {name: sha(HERE / "reference" / name) for name in outputs}
        if frozen != report:
            print(json.dumps({"generated": report, "frozen": frozen}, indent=2), file=sys.stderr)
            return 1
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
