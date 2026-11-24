#!/usr/bin/env python3
"""Bind the sodium ground-capture scope review to the reviewed primary sources."""

import argparse
from decimal import Decimal
import hashlib
import html
import json
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

from native_m3_literature import exactly_one, pdf_text, sha256


PACKET = "sodium-ground-capture-benchmark-20261124"
ARCHIVE = "CollisionDB-2024.1.zip"
ARCHIVE_SHA = "a485306d28edac0f11d0977cd558dfd4d2a0e8851abe508c66ed3f121afea895"
MEMBER_SHA = "418a70144dce1719c38ed4d72d5994ecf370c9bb8f5778722c689475df618537"
CHART_SHA = "939ccf8359d42778986eb34b3119315e3418446198e2ceb3be7949d12e8aa24e"
PAPER_SHA = "58bc33f7f8e7dd960ff2972edf0bd3999ada697e8ec666f04da97bc08e7e8078"
FOLLOWUP_SHA = "f63e74687a27bb09a14771c871acdf00dd06c9e1ccebd044f15808472c64d66e"
LICENCE_SHA = "037950c8409104e74fe94af67fe71956c6ecc02c88f933e7b89ce7ae52c31c71"
CAPTURE_DOI = "10.1103/physreva.70.050702"
FIGURE_SCOPE = ("Knoop et al. Figure 2 groups capture by sodium shell. The outer-shell series "
                "sums over hydrogen final states. The figure has no experimental H(1s)-only "
                "outer-shell curve. The archived ground-product label cannot be used to select a physical rate.")


def verify(path, digest):
    """Reject even a well-formed input if it is not the reviewed edition."""
    if sha256(path) != digest:
        raise ValueError(f"changed source bytes: {path.name}")
    return path


def verify_archive(archive_path, member):
    """Optional full-release audit. The small member is sufficient for offline packet replay."""
    # Hash incrementally because the official release exceeds 500 MiB.
    digest = hashlib.sha256()
    with archive_path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    if digest.hexdigest() != ARCHIVE_SHA:
        raise ValueError("unexpected CollisionDB release")
    with ZipFile(archive_path) as archive:
        if archive.read("119990.txt") != member.read_bytes():
            raise ValueError("CollisionDB member differs from the local native record")


def chart_field(markup, label):
    """Read one rendered field. Empty values retain their native meaning."""
    return html.unescape(exactly_one(
        r"<b>\s*" + re.escape(label) + r"\s*=\s*</b>\s*([^<>]*?)\s*<br>", markup)).strip()


def compare_archives(member, chart):
    """Verify the paper and figure, process, units and every knot before pairing these records."""
    lines = member.read_text().splitlines()
    header = json.loads(lines[0])
    if (header["qid"] != "D119990" or header["refs"] != {"B176": {"doi": CAPTURE_DOI}} or
            header["reaction"] != "H+ + Na → H 1s + Na+ 2p6" or
            header["method"] != "experiment" or
            header["columns"] != [{"name": "E", "units": "eV.u-1"}, {"name": "sigma", "units": "cm2"}] or
            not re.fullmatch(r"-{5,}", lines[1])):
        raise ValueError("unexpected v2024.1 header or reference")
    markup = chart.read_text()
    if (chart_field(markup, "Id Number") != "38917" or chart_field(markup, "Reference Number") != "2004E503" or
            chart_field(markup, "Comment") != "FIG. 2. , IAEA Collision DB file name = 119990.txt" or
            chart_field(markup, "Initial Excited B") or chart_field(markup, "Y Axis Error Data")):
        raise ValueError("unexpected CHART cross-reference or uncertainty metadata")
    iaea_rows = [tuple(map(Decimal, line.split())) for line in lines[2:]]
    energies = list(map(Decimal, chart_field(markup, "X Axis Data").split(",")))
    cross_sections = list(map(Decimal, chart_field(markup, "Y Axis Data").split(",")))
    if len(energies) != 12 or len(cross_sections) != 12 or list(zip(energies, cross_sections)) != iaea_rows:
        raise ValueError("native CHART and CollisionDB arrays differ")
    return header


def describe_theory(paper):
    """Recover only the isolated published fact used to assess the conflicting integral."""
    text = pdf_text(paper, 11, 11)
    table = text.split("Table 4. Computed (A34) partial cross sections", 1)[1].split("Figure 5.", 1)[0]
    # The visually reviewed printed header order is E, 1s, 2s, 2p, 3s, 3p, 3d.
    # Poppler loses the small sigma subscripts in the PDF's text layer.
    # Read the first channel only, accepting exactly one seven-cell 4 keV/u row.
    area = exactly_one(r"^ *4\.0 +(\d+\.\d+) (?: +\d+\.\d+){5} *$", table)
    return dict(
        doi="10.1088/0022-3700/19/3/014", locator="Table 4, printed page 318, 4.0 keV/u, sigma_1s",
        attribution="R. Shingal, B. H. Bransden, A. M. Ermolaev, D. R. Flower, C. W. Newby and C. J. Noble (1986)",
        reaction="H+ + Na(3s) → H(1s) + Na+ (frozen closed-shell core)",
        method="A34 two-centre travelling atomic orbital calculation, one active Na valence electron",
        # Preserve the printed precision and axis in fields accessible to
        # any application. No interpolatable cross-section table is issued.
        benchmark_point=dict(projectile_energy="4.0", energy_unit="keV/u", energy_frame="target-rest laboratory",
                             cross_section=area, area_unit="Å²", cross_section_cm2=f"{Decimal(area) * Decimal('1e-16'):.1E}",
                             uncertainty=None, method="A34 theory"),
        assessment=(f"At projectile laboratory energy 4.0 keV/u the A34 ground-state partial is {area} Å², "
                    f"equivalent to {Decimal(area) * Decimal('1e-16'):.1E} cm². This theoretical result "
                    "independently supports the scale of the state-resolved measurement in Knoop's thesis."),
        uncertainty="No 4 keV/u H(1s) uncertainty reported. The authors describe basis sensitivity for small partials.",
        scope="One integral benchmark. It supplies no differential distribution or transport error bound.")


def describe_followup(paper):
    """Read the next basis calculation on its own projectile-energy axis."""
    text = pdf_text(paper, 5, 5)
    table = text.split("Table 1. Our computed partial and total cross sections", 1)[1].split(
        "3.2. Excitation cross section", 1)[0]
    # Visually reviewed columns are E_lab (keV), speed (au), H(1s), H(2s),
    # H(2p), and total. The text layer splits the last token into "64.71 1".
    # The native H(1s) cell is directly readable without repairing that column.
    area = exactly_one(r"^ *4\.00 +0\.400 +(\d+\.\d+) +\d+\.\d+ +\d+\.\d+ +[\d. ]+ *$", table)
    return dict(
        doi="10.1088/0022-3700/20/18/019", locator="Table 1, printed page 4818, E_lab=4.00 keV, sigma_H(1s)",
        attribution="R. Shingal and B. H. Bransden (1987)",
        reaction="H+ + Na(3s) → H(1s) + Na+ (valence capture)",
        method="63-state two-centre expansion in travelling atomic orbitals with rectilinear projectile trajectories",
        benchmark_point=dict(projectile_energy="4.00", energy_unit="keV", energy_frame="target-rest laboratory",
                             cross_section=area, area_unit="Å²", cross_section_cm2=f"{Decimal(area) * Decimal('1e-16'):.3E}",
                             uncertainty=None, method="63-state theory"),
        assessment=(f"At 4.00 keV projectile laboratory energy the calculated partial is {area} Å² "
                    f"({Decimal(area) * Decimal('1e-16'):.3E} cm²). The published lab energy axis is keV. "
                    "The A34 paper labels its axis in keV/u. Both give much smaller ground partials "
                    "than the shell-summed scale of Knoop Figure 2."),
        uncertainty="No experimental or model error bar for this point. The two basis calculations overlap in authorship.")


def prepare(acquisition, output, archive=None):
    member = verify(acquisition / "119990.txt", MEMBER_SHA)
    chart = verify(acquisition / "nifs-chart-38917.html", CHART_SHA)
    paper = verify(acquisition / "shingal1986.pdf", PAPER_SHA)
    followup = verify(acquisition / "shingal1987.pdf", FOLLOWUP_SHA)
    notice = verify(acquisition / "iaea-licence.html", LICENCE_SHA)
    if archive is not None:
        verify_archive(archive, member)
    native_header = compare_archives(member, chart)
    theory = describe_theory(paper)
    followup_theory = describe_followup(followup)
    receipts = json.loads((acquisition / "retrieval.json").read_text())
    output.mkdir(parents=True, exist_ok=False)
    sources = []
    shutil.copyfile(notice, output / notice.name)

    def source(identifier, native, receipt_name, url, quantity, iaea=False):
        receipt = receipts[receipt_name]
        if receipt["source_uri"] != url or receipt["http_status"] != 200:
            raise ValueError("unexpected acquisition for " + identifier)
        # The packet builder verifies these private local inputs. Only IAEA's
        # licensed release member is copied to the distributed object store.
        shutil.copyfile(native, output / native.name)
        sources.append(dict(
            id=identifier, logical_name=native.name, path=native.name, sha256=sha256(native),
            source_uri=url, final_url=url, retrieved_utc=receipt["retrieved_utc"], http_status=200,
            role="source-native", policy="iaea-cc-by-4.0" if iaea else "metadata-only",
            rights_evidence=[dict(path=notice.name, sha256=LICENCE_SHA,
                                  source_uri="https://db-amdis.org/collisiondb/licence/")] if iaea else [],
            physical_activation=False, quantity_metadata=quantity))

    source("iaea-d119990:release-2024.1", member, ARCHIVE,
           "https://db-amdis.org/collisiondb/data/releases/" + ARCHIVE,
           dict(archive_member=member.name, archive_sha256=ARCHIVE_SHA, archive_bytes=540735640,
                attribution="IAEA CollisionDB, JAEA/Kyoto, Knoop, Morgenstern and Hoekstra (2004). CC BY 4.0.",
                refs=native_header["refs"], native_reaction=native_header["reaction"],
                units=native_header["columns"], rows=12, assessment=FIGURE_SCOPE,
                energy_convention="eV/u; target frame is documented on the IAEA dataset web page.",
                uncertainty="No global or per-knot error in release member."), iaea=True)
    source("nifs-chart:38917-knoop", chart, chart.name, "https://dbshino.nifs.ac.jp/nifsdb/chart/list",
           dict(locator="NIFS CHART 38917, bibliography 2004E503, comment identifies Figure 2 and 119990.txt.",
                selected_by="Standard public search, H projectile, Na target, author Knoop%, custom display.",
                scope=FIGURE_SCOPE, same_archive_knots=12, uncertainty="Native Y Axis Error Data is empty."))
    source("shingal1986-na:ground-partial-benchmark", paper, paper.name,
           "https://iopscience.iop.org/article/10.1088/0022-3700/19/3/014/pdf", theory)
    source("shingal1987-na:ground-partial-benchmark", followup, followup.name,
           "https://iopscience.iop.org/article/10.1088/0022-3700/20/18/019/pdf", followup_theory)
    packet = dict(schema="amps-chemical-source-packet-v1", id=PACKET,
                  title="Proton–sodium capture cross-section provenance and ground-state benchmark",
                  provider="IAEA CollisionDB, NIFS CHART, Shingal et al.", version="20261124",
                  description="Audit of the D119990 figure lineage and a separate frozen-core calculation.",
                  limitations=[FIGURE_SCOPE, "Physical activation and a ground-product angular law remain unavailable."],
                  physical_activation=False, sources=sources)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2, ensure_ascii=False) + "\n")
    return output / "source_packet.json"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-archive", type=Path, help="Optionally verify the entire v2024.1 ZIP against its digest.")
    args = parser.parse_args()
    print(prepare(args.acquisition, args.output_dir, args.verify_archive))
