#!/usr/bin/env python3
"""Prepare the reviewed November 24 native chemistry source packet offline."""

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import shutil
from zipfile import ZipFile

import native_iaea


ZENODO_POLICY = "zenodo-4555716-cc-by-4.0"
ZIP_NAME = "zenodo-he-diff-cross-sections.zip"
ZIP_HASH = "16027a414c603d84a2a500964a8a877a2bf2e8d5573f53827330f75e5ab7e0e8"
RECORD_NAME = "zenodo-swa-record-4555716.json"
RECORD_HASH = "11538b55d67c1b89aa2ad69fd0cfbf6682580028c96c15ca0027b2ba1b1e3ef1"
IAEA_NAME = "iaea-d119990-current.html"
IAEA_HASH = "1015ec7a843e02ef7fc1ddf97f3589ac4dc4823f12bb28e9434313876dd90fd6"
IAEA_LICENCE_HASH = "037950c8409104e74fe94af67fe71956c6ecc02c88f933e7b89ce7ae52c31c71"
ENERGY_COUNT = 661
ANGLE_COUNT = 768
AREA_UNIT_CM2 = "2.8002852056E-17"
# This compares printed header integrals with printed quadrature inputs.
# It is a source-format tolerance, unrelated to physical model accuracy.
INTEGRAL_RELATIVE_TOLERANCE = 1e-9
INDEX_COLUMNS = ("log10_e_cm_x1000", "native_member", "energy_cm_eV_native",
                 "relative_speed_cm_s_native", "reduced_mass_eV_c2_native",
                 "integrated_sigma_au_native", "angular_rows")
ANGLE_COLUMNS = ("log10_e_cm_x1000", "angle_rad_native",
                 "weighted_dsigma_au_native", "quadrature_weight_native")
MEMBER = re.compile(r"He0He1ct/dcsHe0He1ct_1000log10e_([+-][0-9]{4})\.dat")


def sha256(path):
    with path.open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def view_tables(zip_path, output):
    """Index the CEX branch and retain the source's exact angular tokens.

    Native filenames encode 1000 log10(E_CM/eV), and table columns describe
    the CM angle of the neutral. Column 2 already contains 2*pi*sin(theta).
    The recipe never inserts missing angles, rescales, or treats the native
    quadrature weights as midpoint-bin edges.
    """
    index_path, angles_path = output / "he-ct-index.csv", output / "he-ct-angles.csv"
    max_error = 0.0
    with ZipFile(zip_path) as archive:
        members = {}
        for name in archive.namelist():
            # The producer ZIP also preserves the native directory and Finder
            # metadata. Neither contains angle rows, and both stay in the ZIP.
            if not name.startswith("He0He1ct/") or name in ("He0He1ct/", "He0He1ct/.DS_Store"):
                continue
            match = MEMBER.fullmatch(name)
            if not match or int(match[1]) in members:
                raise ValueError("unexpected or duplicate Zenodo He0He1ct member: " + name)
            members[int(match[1])] = match[1], name
        if len(members) != ENERGY_COUNT or min(members) != -4000 or max(members) != 4000:
            raise ValueError("native charge-transfer energy grid changed")
        with index_path.open("x", encoding="ascii", newline="") as index_file, \
                angles_path.open("x", encoding="ascii", newline="") as angles_file:
            index, angles = csv.writer(index_file, lineterminator="\n"), csv.writer(angles_file, lineterminator="\n")
            index.writerow(INDEX_COLUMNS)
            angles.writerow(ANGLE_COLUMNS)
            for label, name in (members[code] for code in sorted(members)):
                lines = archive.read(name).decode("ascii").splitlines()
                header_lines, data_lines = lines[:18], lines[18:]
                if ("# Charge transfer collision: He^0 + He^+ -> He^+ + He^0" not in
                        [line.strip() for line in header_lines] or
                        tuple(header_lines[-1].split()) != ("#", "the", "dcs", "wei") or
                        not any("dcs <- dsigma/dOmega*2*pi*sin(the)" in line and
                                AREA_UNIT_CM2 in line for line in header_lines)):
                    raise ValueError("native branch, columns or units changed: " + name)
                header = {}
                for field, unit in (("CM energy", "eV"), ("Rel speed", "cm/s"),
                                    ("Red. mass", "eV/c^2"), ("Int sigma", "a.u.")):
                    matches = re.findall(r"^# " + re.escape(field) + r":\s*(\S+)\s+" +
                                         re.escape(unit) + r"\s*$", "\n".join(header_lines), re.M)
                    if len(matches) != 1 or not math.isfinite(float(matches[0])) or float(matches[0]) <= 0:
                        raise ValueError("invalid native header " + field + " in " + name)
                    header[field] = matches[0]
                if not math.isclose(float(header["CM energy"]), 10 ** (int(label)/1000), rel_tol=1e-9):
                    raise ValueError("native CM energy disagrees with its filename: " + name)
                if len(data_lines) != ANGLE_COUNT:
                    raise ValueError("native angular row count changed: " + name)
                previous_angle, terms = 0.0, []
                for line in data_lines:
                    tokens = line.split()
                    angle, differential, weight = map(float, tokens)
                    if not (previous_angle < angle < math.pi and math.isfinite(differential) and
                            differential >= 0 and math.isfinite(weight) and weight > 0):
                        raise ValueError("invalid native angular data: " + name)
                    previous_angle = angle
                    terms.append(differential * weight)
                    angles.writerow((label, *tokens))
                relative_error = abs(math.fsum(terms) / float(header["Int sigma"]) - 1.0)
                if relative_error > INTEGRAL_RELATIVE_TOLERANCE:
                    raise ValueError("native quadrature integral disagrees with its header: " + name)
                max_error = max(max_error, relative_error)
                index.writerow((label, name, header["CM energy"], header["Rel speed"],
                                header["Red. mass"], header["Int sigma"], ANGLE_COUNT))
    return index_path, angles_path, max_error


def iaea_projection_hash(path):
    """Authenticate the current HTML and compare with the frozen CSV digest."""
    # The live page has a new Bokeh container identity. This audit can use its
    # separately reviewed container hash without changing the r5 recipe pin.
    specification = dict(native_iaea._SOURCES[119990], sha256=IAEA_HASH)
    series, = native_iaea._decode(path.read_bytes(), 119990, specification)
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\r\n")
    writer.writerow(native_iaea._HEADER)
    writer.writerows(series["rows"])
    digest = hashlib.sha256(buffer.getvalue().encode("utf-8")).hexdigest()
    if digest != specification["outputs"][0][2]:
        raise ValueError("live D119990 numerical projection differs from the frozen data")
    return digest


def prepare(evidence, iaea_licence, output):
    """Create a fresh, reproducible packet from explicit retrieval receipts."""
    indexes = {
        folder: json.loads((evidence / folder / "retrieval.json").read_text())
        for folder in ("", "followup", "additional", "confirmed")
    }
    output.mkdir(parents=True, exist_ok=False)
    native = output / "native"
    native.mkdir()
    recipe = output / Path(__file__).name
    shutil.copyfile(__file__, recipe)
    dependency = output / "native_iaea.py"
    shutil.copyfile(native_iaea.__file__, dependency)

    def acquired(folder, filename, identifier, quantity, *, sha=None,
                 policy="metadata-only", role="source-native", rights=()):
        receipt = indexes[folder][filename]
        source = evidence / folder / receipt["file"]
        digest = sha256(source)
        if (receipt["status"] != "RETRIEVED" or receipt["http_code"] != "200"
                or digest != receipt["sha256"] or (sha and digest != sha)):
            raise ValueError("unavailable or changed reviewed source: " + filename)
        target = native / filename
        shutil.copyfile(source, target)
        return dict(id=identifier, logical_name=filename, path=str(target.relative_to(output)),
                    sha256=digest, source_uri=receipt["requested_url"],
                    retrieved_utc=receipt["retrieved_utc"], http_status=200,
                    role=role, policy=policy, rights_evidence=list(rights),
                    quantity_metadata=quantity, physical_activation=False)

    if sha256(iaea_licence) != IAEA_LICENCE_HASH:
        raise ValueError("IAEA licence differs from the reviewed rights evidence")
    licence = native / "iaea-licence.html"
    shutil.copyfile(iaea_licence, licence)
    iaea_rights = dict(path=str(licence.relative_to(output)), sha256=IAEA_LICENCE_HASH,
                       source_uri="https://db-amdis.org/collisiondb/licence/")
    archive_rights = dict(path="native/" + RECORD_NAME, sha256=RECORD_HASH,
                          source_uri="https://zenodo.org/api/records/4555716")
    sources = [acquired("", RECORD_NAME, "swaczyna-4555716:record",
                        dict(doi="10.5281/zenodo.4555716", creator="Pawel Swaczyna", licence="CC BY 4.0"),
                        sha=RECORD_HASH, policy=ZENODO_POLICY, role="source-metadata", rights=[archive_rights])]
    source_scope = dict(
        method="author-calculated differential cross sections",
        process="He^0 + He^+ -> He^+ + He^0 (He0He1ct branch)",
        angle_convention="CM deflection of neutral, Swaczyna et al. Appendix A equation (4)",
        energy_frame="center-of-mass collision energy; eV, not eV/u",
        per_energy_angle_rows=ANGLE_COUNT, energy_support_eV=[0.0001, 10000],
        dcs_token="(d sigma/d Omega) * 2*pi*sin(theta); native area units, no second Jacobian",
        area_unit_cm2=AREA_UNIT_CM2, source_uncertainties="not provided by the ZIP",
        licence="CC BY 4.0", licence_url="https://creativecommons.org/licenses/by/4.0/",
        creator="Pawel Swaczyna", doi="10.5281/zenodo.4555716",
        paper_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC8861972/",
    )
    original = acquired("confirmed", ZIP_NAME, "swaczyna-4555716:native-zip",
                        dict(source_scope, native_archive_scope="all six original collision branches"),
                        sha=ZIP_HASH, policy=ZENODO_POLICY, rights=[archive_rights])
    sources.append(original)
    index_path, angles_path, max_error = view_tables(output / original["path"], output)
    for path, name, columns, count in (
        (index_path, "he-ct-index", INDEX_COLUMNS, ENERGY_COUNT),
        (angles_path, "he-ct-angles", ANGLE_COLUMNS, ENERGY_COUNT * ANGLE_COUNT),
    ):
        identifier = "swaczyna-4555716:" + name
        sources.append(dict(id=identifier + ":csv", logical_name=path.name, path=path.name,
                            sha256=sha256(path), source_uri=original["source_uri"],
                            role="derived-table", policy=ZENODO_POLICY, rights_evidence=[archive_rights],
                            derived_from=[original["id"]],
                            recipe=dict(path=recipe.name, sha256=sha256(recipe)),
                            quantity_metadata=dict(
                                source_scope, operation="lossless header/token CSV view of He0He1ct",
                                recipe_dependencies=[dict(path=dependency.name, sha256=sha256(dependency))]),
                            table=dict(id=identifier, format="csv", row_count=count, columns=list(columns),
                                       schema_id="swaczyna-he0he1ct-native-v1",
                                       description="Swaczyna He0He1ct native-token " + name),
                            physical_activation=False))
    live = acquired("", IAEA_NAME, "iaea-d119990-20261124:live-html",
                    dict(reaction="H+ + Na -> H(1s) + Na+(2p^6)",
                         energy_frame="target rest; projectile eV/u", role="corroborating capture subset",
                         limitations="Original citation, initial Na state, uncertainties and differential data absent."),
                    sha=IAEA_HASH, policy="iaea-cc-by-4.0", rights=[iaea_rights])
    comparison_hash = iaea_projection_hash(output / live["path"])
    live["quantity_metadata"]["frozen_numeric_projection_sha256"] = comparison_hash
    sources.append(live)
    # Register documentary identities for acquired material whose licence
    # or scope does not clear general shared numerical access. The publisher
    # never copies these source bytes into the installed object store.
    for folder, filename, identifier, observable, limitation in (
        ("", "norad-he-partial.txt", "norad-he-ground-residual", "He I photoionization to ground He II",
         "OSU CC BY-NC 4.0. Separate native energy axis and unresolved broadband spectral/state crosswalk."),
        ("", "norad-he-ground-total.txt", "norad-he-ground-total", "Total ground-parent He I photoionization",
         "OSU CC BY-NC 4.0. Total and ground-residual files have unequal native grid membership."),
        ("", "curtin-e-He_sig.zip", "curtin-he-ccc", "TICS.s1S ground-parent one-electron ionization",
         "No managed-storage grant located. Unpublished CCC evaluation gives total cross section only."),
        ("additional", "nifs-ralchenko-index.pdf", "nifs-059", "Ground He integral ionization fit, Eq. 9/Table IV",
         "Publisher directs reproduction inquiries to NIFS. A fit alone does not bound the conditional mean."),
        ("additional", "uk-verner-table1.dat", "verner-yakovlev-subshells", "Atomic subshell ionization fits",
         "No managed-data grant located. Shell labels do not establish residual-ion or satellite populations."),
        ("", "nist-gallagher-jpcrd339.pdf", "gallagher-water-fragments", "H2O fragment ion marginal yields",
         "No source-specific managed-copy grant established. Joint electron/neutral/ion observations are missing."),
        ("confirmed", "groningen-knoop-thesis-c5.pdf", "knoop-h-na", "Na(3s) outer-shell capture to H(1s)",
         "Peak-normalized differential figures from 7 keV/u. Native angular bins/uncertainty and IAEA lineage missing."),
        ("followup", "upm-courbin.pdf", "courbin-h-na", "Na(3s)/Na(3p) capture differential theory for H(n=2)",
         "Figure-level angular data for an excited channel. No managed-copy grant or H(1s) law established."),
        ("confirmed", "heii-nist-levels.html", "nist-heii-levels", "He II levels and first excited state energy",
         "No reviewed grant for these numerical energy levels in the shared release."),
        ("", "iaea-nagata-record.html", "nagata-h-na", "Inferred eventual H(1s) capture, native 0.4–5 keV",
         "Bibliographic inspection only. Different state/time definition from prompt ground-state capture."),
    ):
        sources.append(acquired(folder, filename, "chemistry-archive-20261124:" + identifier,
                                dict(observable=observable, limitation=limitation)))

    packet = dict(schema="amps-chemical-source-packet-v1", id="chemistry-archive-20261124",
                  title="November 24 atomic and molecular collision archive inputs",
                  provider="Swaczyna/Zenodo, IAEA, NORAD, Curtin and primary research archives",
                  version="20261124",
                  description="Calculated native He resonant charge-transfer tables and corroborating sodium evidence.",
                  limitations=["Source data only. No lunar or application reaction has been activated.",
                               "No joint charged-water product distribution or complete error envelope admitted.",
                               "NILS contact-surface counts and Bodi threshold study remain bibliography.",
                               "No transfer of helium angles to proton–sodium or proton–helium."],
                  sources=sources, physical_activation=False)
    (output / "source_packet.json").write_text(json.dumps(packet, indent=2) + "\n", encoding="utf-8")
    return dict(source_count=len(sources), energy_count=ENERGY_COUNT, angular_rows=ENERGY_COUNT*ANGLE_COUNT,
                maximum_relative_quadrature_discrepancy=max_error, iaea_projection_sha256=comparison_hash)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence-root", required=True, type=Path)
    parser.add_argument("--iaea-licence", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.evidence_root, args.iaea_licence, args.output_dir), indent=2))


if __name__ == "__main__":
    main()
