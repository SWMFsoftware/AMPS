#!/usr/bin/env python3
"""Compare completed bounded numerical experiments, not campaign qualification.

Keep exact epoch/channel results and artifact fingerprints. These experiments
usually contain only two difficult epochs and cannot satisfy SA04 population
coverage. Direction-level accounting is available only for the newer all-trace
diagnostic logger; older incomplete logs are explicitly marked unavailable.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re
import sys

from sa04_common import CHANNELS, atomic_json, fractional_change, read_point_table, read_spectrum_zones, sha256_file, utc_now


def completed(root: Path) -> dict:
    """Reject failed/incomplete runs and decks changed after recorded execution."""
    status = json.loads((root / "diagnostic_status.json").read_text())
    if status["status"] != "COMPLETE" or status["return_code"] != 0 or not status["complete_products"]:
        raise ValueError(f"{root}: incomplete diagnostic")
    if sha256_file(root / "input/AMPS_PARAM.in") != status["input_sha256"]:
        raise ValueError(f"{root}: changed diagnostic input")
    # Older preserved diagnostics predate these extra hashes. Keep them readable,
    # but verify every hash that was actually recorded rather than inventing one.
    for key in ("referenced_input_sha256", "output_sha256"):
        for path, expected in status.get(key, {}).items():
            if sha256_file(Path(path)) != expected:
                raise ValueError(f"{root}: changed {key} artifact {path}")
    if status.get("log_sha256") and sha256_file(root / "amps.log") != status["log_sha256"]:
        raise ValueError(f"{root}: changed diagnostic trace log")
    return status


def final_directions(root: Path) -> dict:
    """Select the last bounded attempt for each exact epoch/rigidity/direction seed."""
    results = {}
    for line in (root / "amps.log").read_text().splitlines():
        if not line.startswith("[density-trace]"):
            continue
        fields = dict(field.split("=", 1) for field in line.split()[1:])
        key = (fields["epoch"], fields["rigidity_GV"], fields["direction"])
        if key not in results or int(fields["pass"]) >= int(results[key]["pass"]):
            results[key] = fields
    return results


def complete_direction_population(root: Path, directions: dict) -> bool:
    """Check log coverage against native spectrum sample counts at every energy.

    Equal key sets in two truncated logs do not establish completeness. Native
    spectrum rows provide an independent expected population. The logger converts
    rigidity back to kinetic energy (DensityMode3D.cpp), introducing a few ulps of
    round-trip error relative to the original spectrum grid. Match uniquely within
    eight machine epsilons, never with the 2% scientific convergence tolerance.
    """
    expected = {}
    for path in root.glob("mode3d_points_spectrum_snapshot_*.dat"):
        for zone in read_spectrum_zones(path):
            energy_col, count_col = zone.columns.index("E_MeV"), zone.columns.index("N_sampled")
            for row in zone.rows:
                count = float(row[count_col])
                if not count.is_integer() or count <= 0:
                    raise ValueError(f"{path}: invalid native sampled count")
                expected[(zone.utc, float(row[energy_col]))] = set(range(int(count)))
    actual = {}
    for row in directions.values():
        if "particle_energy_MeV" not in row:
            return False
        key = (row["epoch"].rstrip("Z") + "Z", float(row["particle_energy_MeV"]))
        indices = actual.setdefault(key, set())
        direction = int(row["direction"])
        if direction in indices:
            return False
        indices.add(direction)
    if not expected or len(actual) != len(expected):
        return False
    matched = set()
    for (epoch, energy), indices in actual.items():
        candidates = [key for key in expected if key[0] == epoch and
                      math.isclose(key[1], energy, rel_tol=8 * sys.float_info.epsilon, abs_tol=0)]
        if len(candidates) != 1 or candidates[0] in matched or expected[candidates[0]] != indices:
            return False
        matched.add(candidates[0])
    return matched == set(expected)


def field_audits(root: Path) -> dict:
    """Extract actual field evidence, not driver-derived snapshot identifiers.

    A driver's snapshot ID can be identical even when parallel construction writes
    different sampled fields. The opt-in content fingerprint covers every packed B/E
    value and geometry; the two sampled audits independently test initialization and
    interpolation. Missing evidence remains missing and cannot be read as zero errors.
    """
    epochs = {}
    for line in (root / "amps.log").read_text().splitlines():
        if "Prepared compact global cell-centered fields:" in line:
            fields = dict(re.findall(r"(content|meshRevision|epoch)=([^,\s]+)", line))
            if len(fields) == 3 and fields["content"].startswith("standalone-field-v1-"):
                record = epochs.setdefault(fields.pop("epoch").rstrip("Z") + "Z", {})
                record.update(fields)
        match = re.search(r"\[field-audit\] (owner|interpolation) epoch=(\S+) sampled=(\d+) differing=(\d+)", line)
        if match:
            kind, epoch, sampled, differing = match.groups()
            record = epochs.setdefault(epoch.rstrip("Z") + "Z", {})
            if kind in record:
                raise ValueError(f"{root}: duplicate {kind} field audit for {epoch}")
            record[kind] = {"sampled": int(sampled), "differing": int(differing)}
    return epochs


def repeatability_evidence(first: dict, second: dict, epochs: list[str]) -> dict:
    """Require complete, identical field identities and zero sampled discrepancies.

    This is an exact-repeatability check, deliberately separate from the 2% physical
    refinement gate. It never permits that scientific tolerance to hide a race.
    """
    expected = {epoch.rstrip("Z") + "Z" for epoch in epochs}
    for records in (first, second):
        if set(records) != expected or any(
            set(row) != {"content", "meshRevision", "owner", "interpolation"} or
            any(row[kind]["sampled"] <= 0 for kind in ("owner", "interpolation"))
            for row in records.values()
        ):
            raise ValueError("missing/incomplete field-audit evidence for repeatability")
    matching = all(first[epoch][key] == second[epoch][key] for epoch in expected
                   for key in ("content", "meshRevision"))
    clean = all(row[kind]["differing"] == 0 for records in (first, second)
                for row in records.values() for kind in ("owner", "interpolation"))
    return {"identical_field_content": matching, "sampled_audits_clean": clean,
            "first": first, "second": second, "status": "PASS" if matching and clean else "FAIL"}


def compare(coarse: Path, fine: Path, output: Path, check_repeatability: bool = False) -> dict:
    """Compare all nine integrated quantities without claiming campaign qualification."""
    a_status, b_status = completed(coarse), completed(fine)
    for key in ("epochs", "source_deck_sha256"):
        if a_status[key] != b_status[key]:
            raise ValueError(f"diagnostics differ in {key}")
    rows, artifacts = [], {}
    for kind, columns in (("flux", ["F_tot_m2s1"] + [f"F_{name}_m2s1" for name, _, _ in CHANNELS]),
                          ("density", ["N_cm^-3"])):
        tables = []
        for root in (coarse, fine):
            paths = sorted(root.glob(f"mode3d_points_{kind}_snapshot_*.dat"))
            records = [row for path in paths for row in read_point_table(path)]
            table = {row["TimeUTC"]: row for row in records}
            if len(table) != len(records) or set(table) != {t.rstrip("Z") + "Z" for t in a_status["epochs"]}:
                raise ValueError(f"{root}: inconsistent diagnostic product epochs")
            tables.append(table)
            artifacts.update({str(path.resolve()): sha256_file(path) for path in paths})
        for epoch in sorted(tables[0]):
            for column in columns:
                a, b = float(tables[0][epoch][column]), float(tables[1][epoch][column])
                if not all(math.isfinite(value) and value >= 0 for value in (a, b)):
                    raise ValueError(f"invalid native integrated quantity {epoch}/{column}")
                rows.append(dict(epoch=epoch, column=column, coarse=a, fine=b,
                                 fractional_change=fractional_change(a, b)))
    first, second = final_directions(coarse), final_directions(fine)
    changes = Counter()
    direction_value_changes = 0
    # Only compare complete direction populations, not the unresolved-only logs
    # emitted by the initial forensic build.
    complete_logs = (bool(first) and set(first) == set(second) and
                     complete_direction_population(coarse, first) and
                     complete_direction_population(fine, second))
    if complete_logs:
        for key in first:
            direction_value_changes += first[key] != second[key]
            old, new = first[key]["termination"], second[key]["termination"]
            if old != new:
                changes[old + " -> " + new] += 1
    maximum = max(row["fractional_change"] for row in rows)
    audit = {"coarse": field_audits(coarse), "fine": field_audits(fine)}
    threshold = 0.0 if check_repeatability else 0.02
    passed = maximum <= threshold
    difference_gate = "PASS" if passed else "FAIL"
    field_identity = "NOT_ESTABLISHED"
    if check_repeatability:
        # Rank/thread layout may change; physics, mover, timestep, trace budgets and
        # executable may not. Those are independent numerical experiments instead.
        for key in ("executable_sha256", "mover_override", "dt_override", "extension_passes_override"):
            if a_status.get(key) != b_status.get(key):
                raise ValueError(f"repeatability controls differ in {key}")
        audit = repeatability_evidence(audit["coarse"], audit["fine"], a_status["epochs"])
        field_identity = audit["status"]
        passed = passed and complete_logs and direction_value_changes == 0 and audit["status"] == "PASS"
    elif audit["coarse"] and audit["fine"]:
        # This diagnostic changes only mover/timestep/trace budgets, never mesh or
        # physical inputs. Different constructed fields invalidate an isolated
        # integration-refinement interpretation even if products differ by <2%.
        audit = repeatability_evidence(audit["coarse"], audit["fine"], a_status["epochs"])
        field_identity = audit["status"]
        passed = passed and field_identity == "PASS"
    status = "PASS" if passed else "FAIL"
    if field_identity == "NOT_ESTABLISHED":
        status = "INCOMPLETE"
    result = dict(created_utc=utc_now(), coarse=str(coarse.resolve()), fine=str(fine.resolve()),
                  interpretation="bounded diagnostic only; not full-population SA04 qualification",
                  check="exact_repeatability" if check_repeatability else "bounded_numerical_difference",
                  threshold=threshold, status=status, field_audit=audit,
                  integrated_difference_gate=difference_gate, controlled_field_identity=field_identity,
                  maximum_fractional_change=maximum, rows=rows, artifact_sha256=artifacts,
                  numerical_controls={"coarse": a_status, "fine": b_status},
                  direction_comparison={"available": complete_logs, "population": len(first) if complete_logs else None,
                                        "changed_direction_records": direction_value_changes if complete_logs else None,
                                        "termination_changes": dict(changes)})
    atomic_json(output, result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coarse", required=True, type=Path)
    parser.add_argument("--fine", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--check-repeatability", action="store_true",
                        help="Require exact products/trajectory records and complete clean field audits, not the 2%% convergence tolerance")
    args = parser.parse_args()
    result = compare(args.coarse, args.fine, args.output, args.check_repeatability)
    print(f"[{result['status']}] bounded diagnostic maximum change={result['maximum_fractional_change']:.6g}; {args.output}")
    raise SystemExit(0 if result["status"] == "PASS" else 2)
