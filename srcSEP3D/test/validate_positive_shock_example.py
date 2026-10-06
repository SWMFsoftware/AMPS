#!/usr/bin/env python3
"""Run and audit the production-style reduced-front positive example.

This is an evidence runner, not the example's runtime interface.  It creates
relocated copies of the maintained deck so one-/four-rank and output-cadence
runs cannot overwrite one another, captures complete console logs, and checks
native receipts/products.  It never edits the source deck or event.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import pathlib
import re
import subprocess
import sys
import time


ROOT = pathlib.Path(__file__).resolve().parents[2]
DECK = ROOT / "srcSEP3D/examples/shock-front/positive_1au.in"
EVENT = ROOT / "srcSEP3D/examples/shock-front/positive_1au.event"
EXPECTED_BASE_TICKS = {0, 1, 2, 3, 16, 17, 18, 60, 120, 180, 206, 207, 208}
EXPECTED_ARRIVAL_S = 124200.0
ROOT_TIME_TOLERANCE_S = 1.0e-6


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_header(path: pathlib.Path, line_count: int) -> str:
    """Read a bounded Tecplot header without loading a ~1 GiB volume file.

    AMPS volume files put VARIABLES on line one and ZONE on line two, while
    the reduced-front surface adds TITLE and therefore puts ZONE on line
    three.  Keeping this helper line-count based makes the audit inexpensive
    and, unlike repeated ``path.open().readline()`` expressions, advances one
    stream so every requested header line is actually inspected.
    """
    if not path.exists():
        return ""
    lines: list[str] = []
    with path.open(errors="replace") as stream:
        for _ in range(line_count):
            line = stream.readline()
            if not line:
                break
            lines.append(line)
    return "".join(lines)


def tecplot_variables(line: str) -> list[str]:
    """Return quoted Tecplot variable names without guessing column offsets."""
    return re.findall(r'"([^"]+)"', line)


def audit_ambient_rows(path: pathlib.Path, receipt: dict) -> tuple[bool, str]:
    """Check a physical AMPS node rather than accepting a header-only file.

    AMPS writes interpolation-padding vertices with ``background_valid=0``;
    zeros at those vertices are visualization topology, not plasma.  Streaming
    until the first valid vertex avoids loading a roughly 170-MB file and
    proves that the appended fields are the values installed in native center
    storage at the committed epoch.
    """
    if not path.exists():
        return False, f"missing {path}"
    with path.open(errors="replace") as stream:
        variable_line = stream.readline()
        names = tecplot_variables(variable_line)
        required = ("B_x_T", "B_y_T", "B_z_T", "U_x_m_per_s",
                    "U_y_m_per_s", "U_z_m_per_s", "number_density_m-3",
                    "mass_density_kg_m-3", "pressure_Pa", "temperature_K",
                    "simulation_time_s", "background_generation",
                    "front_generation", "background_valid")
        if not names or any(name not in names for name in required):
            return False, "required ambient columns are absent"
        index = {name: names.index(name) for name in required}
        for line in stream:
            if not line or line.startswith(("ZONE", "VARIABLES")):
                continue
            words = line.split()
            if len(words) < len(names):
                continue
            try:
                values = [float(word) for word in words[:len(names)]]
            except ValueError:
                continue
            if values[index["background_valid"]] != 1.0:
                continue
            physical = tuple(index[name] for name in required[:10])
            finite = all(math.isfinite(values[i]) for i in physical)
            positive = all(values[index[name]] > 0.0 for name in (
                "number_density_m-3", "mass_density_kg_m-3", "pressure_Pa",
                "temperature_K"))
            epoch = abs(values[index["simulation_time_s"]] -
                        float(receipt["time_s"])) <= 1.0e-9 and \
                    int(values[index["background_generation"]]) == \
                        int(receipt["background_generation"]) and \
                    int(values[index["front_generation"]]) == \
                        int(receipt["front_generation"])
            return finite and positive and epoch, (
                f"valid native row finite={finite} positive={positive} "
                f"time={values[index['simulation_time_s']]} "
                f"background_generation={values[index['background_generation']]} "
                f"front_generation={values[index['front_generation']]}")
    return False, "no background_valid=1 native vertex was found"


def audit_front_rows(path: pathlib.Path, receipt: dict) -> tuple[bool, str]:
    """Verify serialized geometry and one-sided RH limits row by row.

    Tecplot-facing rows are finite by construction.  Non-shock patches must
    carry the declared neutral display encoding (Mf=0, X=CB=1 and downstream
    columns copied from upstream) while both acceptance/downstream flags stay
    false.  Conversely, every accepted patch must carry finite, compressive
    upstream/downstream states.  This prevents a syntactically valid surface
    file from passing because only its header named the requested quantities.
    """
    if not path.exists():
        return False, f"missing {path}"
    with path.open(errors="replace") as stream:
        title = stream.readline()
        names = tecplot_variables(stream.readline())
        zone = stream.readline()
        match = re.search(r"\bN=(\d+)", zone)
        required = ("x_m", "y_m", "z_m", "normal_x", "normal_y",
                    "normal_z", "normal_speed_m_s", "time_s", "generation",
                    "status_code", "shock_accepted", "upstream_valid",
                    "downstream_valid", "fast_mach", "theta_Bn_rad",
                    "theta_Bn_valid", "density_compression",
                    "magnetic_compression", "magnetic_compression_valid",
                    "rho1_kg_m3", "p1_Pa", "u1x_m_s", "u1y_m_s",
                    "u1z_m_s", "b1x_T", "b1y_T", "b1z_T",
                    "rho2_kg_m3", "p2_Pa", "u2x_m_s", "u2y_m_s",
                    "u2z_m_s", "b2x_T", "b2y_T", "b2z_T")
        if "RH limits valid only when shock_accepted=1" not in title or \
                not match or \
                any(name not in names for name in required):
            return False, "surface title, node count, or required columns are absent"
        node_count = int(match.group(1))
        index = {name: names.index(name) for name in required}
        # Skip AUXDATA lines, then consume exactly N point records; connectivity
        # consists of four integer indices and is intentionally not plasma data.
        rows: list[list[float]] = []
        for line in stream:
            if line.startswith("AUXDATA"):
                continue
            words = line.split()
            if len(words) != len(names):
                if rows:
                    break
                continue
            try:
                rows.append([float(word) for word in words])
            except ValueError:
                return False, "non-numeric surface point record"
            if len(rows) == node_count:
                break
    accepted = 0
    rejected_with_ambient = 0
    good = len(rows) == node_count and node_count > 0
    for values in rows:
        good = good and all(math.isfinite(value) for value in values)
        normal = math.sqrt(sum(values[index[name]] ** 2 for name in (
            "normal_x", "normal_y", "normal_z")))
        good = good and math.isfinite(normal) and abs(normal - 1.0) <= 1.0e-12
        good = good and abs(values[index["time_s"]] -
                            float(receipt["time_s"])) <= 1.0e-9
        good = good and int(values[index["generation"]]) == \
            int(receipt["front_generation"])
        is_accepted = values[index["shock_accepted"]] == 1.0
        good = good and is_accepted == \
            (int(values[index["status_code"]]) == 5)
        if is_accepted:
            accepted += 1
            finite_names = required
            good = good and values[index["upstream_valid"]] == 1.0 and \
                values[index["downstream_valid"]] == 1.0 and \
                all(math.isfinite(values[index[name]]) for name in finite_names) and \
                values[index["fast_mach"]] > 1.0 and \
                0.0 <= values[index["theta_Bn_rad"]] <= 0.5 * math.pi and \
                values[index["theta_Bn_valid"]] == 1.0 and \
                values[index["density_compression"]] > 1.0 and \
                values[index["magnetic_compression"]] > 0.0 and \
                values[index["magnetic_compression_valid"]] == 1.0 and \
                values[index["rho1_kg_m3"]] > 0.0 and \
                values[index["p1_Pa"]] > 0.0 and \
                values[index["rho2_kg_m3"]] > 0.0 and \
                values[index["p2_Pa"]] > 0.0
        else:
            upstream_downstream = (
                ("rho1_kg_m3", "rho2_kg_m3"), ("p1_Pa", "p2_Pa"),
                ("u1x_m_s", "u2x_m_s"), ("u1y_m_s", "u2y_m_s"),
                ("u1z_m_s", "u2z_m_s"), ("b1x_T", "b2x_T"),
                ("b1y_T", "b2y_T"), ("b1z_T", "b2z_T"))
            ambient_fill = values[index["upstream_valid"]] == 1.0 and \
                values[index["downstream_valid"]] == 0.0 and \
                values[index["fast_mach"]] == 0.0 and \
                values[index["density_compression"]] == 1.0 and \
                values[index["magnetic_compression"]] == 1.0 and \
                values[index["magnetic_compression_valid"]] == 0.0 and \
                all(values[index[upstream]] == values[index[downstream]]
                    for upstream, downstream in upstream_downstream)
            good = good and ambient_fill
            rejected_with_ambient += int(ambient_fill)
    good = good and accepted > 0
    return good, (f"nodes={len(rows)}/{node_count} accepted={accepted} "
                  f"no_shock_ambient_fill={rejected_with_ambient}")


def replace_once(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise RuntimeError(f"expected exactly one deck token {old!r}")
    return text.replace(old, new, 1)


def make_deck(run_root: pathlib.Path, cadence: int) -> pathlib.Path:
    products = run_root / "products"
    text = DECK.read_text()
    text = replace_once(text, "model_asset = positive_1au.event",
                        f"model_asset = {EVENT}")
    text = replace_once(text, "cadence_steps = 60",
                        f"cadence_steps = {cadence}")
    text = replace_once(text,
        "directory = test_output/reduced-front/positive-1au",
        f"directory = {products}")
    for leaf in ("initialization-mesh.dat", "initialization-parker-line.dat",
                 "initialization-data.dat"):
        text = replace_once(text,
            f"test_output/reduced-front/positive-1au/{leaf}",
            str(run_root / leaf))
    text = replace_once(text,
        "output_path = test_output/reduced-front/positive-1au/restart.chk",
        f"output_path = {run_root / 'restart.chk'}")
    run_root.mkdir(parents=True, exist_ok=False)
    path = run_root / "resolved.in"
    path.write_text(text)
    return path


def execute(label: str, ranks: int, cadence: int, root: pathlib.Path,
            executable: pathlib.Path) -> dict:
    run_root = root / label
    deck = make_deck(run_root, cadence)
    command = ["mpiexec", "-n", str(ranks), str(executable),
               "--input", str(deck)]
    log_path = run_root / "execution.log"
    started = time.monotonic()
    with log_path.open("w") as log:
        log.write("command=" + " ".join(command) + "\n")
        log.flush()
        completed = subprocess.run(command, cwd=ROOT, stdout=log,
                                   stderr=subprocess.STDOUT, text=True)
    return {"label": label, "ranks": ranks, "cadence_steps": cadence,
            "returncode": completed.returncode,
            "elapsed_s": time.monotonic() - started,
            "deck": str(deck), "log": str(log_path),
            "products": str(run_root / "products")}


def audit(run: dict) -> tuple[list[dict], dict]:
    checks: list[dict] = []
    products = pathlib.Path(run["products"])

    def check(name: str, condition: bool, detail: str) -> None:
        checks.append({"name": name, "status": "PASS" if condition else "FAIL",
                       "detail": detail, "log": run["log"]})

    check("process-exit", run["returncode"] == 0,
          f"returncode={run['returncode']}")
    receipts = {}
    for path in products.glob("positive-1au-tick-*-receipt.json"):
        data = json.loads(path.read_text())
        receipts[int(data["tick"])] = (path, data)
    expected = set(EXPECTED_BASE_TICKS)
    if run["cadence_steps"] == 30:
        expected |= {30, 90, 150}
    check("snapshot-cadence", set(receipts) == expected,
          f"actual={sorted(receipts)} expected={sorted(expected)}")
    arrival = receipts.get(207, (None, {}))[1]
    exact_arrival = float(arrival.get("exact_observer_arrival_time_s", -1))
    check("accepted-observer-arrival",
          arrival.get("geometric_observer_arrival") is True and
          arrival.get("observer_shock_accepted") is True and
          arrival.get("observer_shock_status") == "solved-fast-shock" and
          abs(exact_arrival - EXPECTED_ARRIVAL_S) <= ROOT_TIME_TOLERANCE_S and
          float(arrival.get("observer_fast_mach", 0)) > 1,
          f"receipt={arrival}")
    all_receipts = [item[1] for item in receipts.values()]
    check("native-owner-epochs-zero-particles", bool(all_receipts) and all(
          r.get("owner_readback_match") is True and
          r.get("received_ghost_readback_match") is True and
          r.get("provider_epoch_match") is True and
          r.get("background_generation") == r.get("front_generation") and
          r.get("particle_count") == 0 and
          r.get("injected_particle_count") == 0 for r in all_receipts),
          f"receipts={len(all_receipts)}")
    if run["ranks"] > 1:
        check("received-mpi-ghost-evidence", bool(all_receipts) and all(
              int(r.get("received_ghost_cells_checked", 0)) > 0
              for r in all_receipts),
              "all four-rank receipts must check received ghosts")

    history_path = products / "shock-history.csv"
    rows = list(csv.DictReader(history_path.open())) if history_path.exists() else []
    check("contiguous-native-history", len(rows) == 209 and
          rows[0].get("tick") == "0" and rows[-1].get("tick") == "208" and
          float(rows[-1].get("time_s", -1)) == 124800 and all(
          int(r["particle_count"]) == 0 and
          int(r["injected_particle_count"]) == 0 for r in rows),
          f"rows={len(rows)} final={rows[-1] if rows else None}")

    ambient = products / "positive-1au-tick-00000207-ambient.dat"
    front = products / "positive-1au-tick-00000207-front.dat"
    ambient_head = read_header(ambient, 2)
    # TITLE, VARIABLES, ZONE and four AUXDATA records precede the first point.
    # Read all seven metadata lines so the encoding policy itself, rather than
    # only the numerical rows, is part of the visualization contract.
    front_head = read_header(front, 7)
    check("visit-tecplot-products", ambient.exists() and front.exists() and
          "mass_density_kg_m-3" in ambient_head and
          "background_generation" in ambient_head and
          "ZONETYPE=FEBRICK" in ambient_head and
          "ZONETYPE=FEQUADRILATERAL" in front_head and
          "theta_Bn_valid" in front_head and
          "magnetic_compression_valid" in front_head and
          "no_shock_fill=\"upstream-ambient-visualization-placeholder\"" in
              front_head and "rho2_kg_m3" in front_head,
          f"ambient={ambient} front={front}")
    ambient_ok, ambient_detail = audit_ambient_rows(ambient, arrival)
    check("native-ambient-state-rows", ambient_ok, ambient_detail)
    front_ok, front_detail = audit_front_rows(front, arrival)
    check("front-geometry-rh-state-rows", front_ok, front_detail)

    facts = {"arrival": arrival,
             "arrival_front_sha256": sha256(front) if front.exists() else "",
             "arrival_ambient_bytes": ambient.stat().st_size if ambient.exists() else 0,
             "event_identity": arrival.get("event_identity", ""),
             "owner_fingerprint_xor": arrival.get("owner_fingerprint_xor"),
             "owner_fingerprint_sum": arrival.get("owner_fingerprint_sum")}
    return checks, facts


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--amps", default=str(ROOT / "amps"))
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--skip-one-rank", action="store_true")
    parser.add_argument("--skip-cadence", action="store_true")
    args = parser.parse_args()
    output_root = pathlib.Path(args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=False)
    executable = pathlib.Path(args.amps).resolve()
    runs = []
    if not args.skip_one_rank:
        runs.append(execute("rank-1-cadence-60", 1, 60, output_root, executable))
    runs.append(execute("rank-4-cadence-60", 4, 60, output_root, executable))
    if not args.skip_cadence:
        runs.append(execute("rank-4-cadence-30", 4, 30, output_root, executable))

    checks, facts = [], {}
    for run in runs:
        run_checks, run_facts = audit(run)
        for item in run_checks:
            item["run"] = run["label"]
        checks.extend(run_checks)
        facts[run["label"]] = run_facts

    base = facts.get("rank-4-cadence-60", {})
    one = facts.get("rank-1-cadence-60")
    if one is not None:
        same = all(one.get(k) == base.get(k) for k in (
            "event_identity", "owner_fingerprint_xor",
            "owner_fingerprint_sum", "arrival_front_sha256"))
        checks.append({"run": "cross-run", "name": "one-four-rank-agreement",
                       "status": "PASS" if same else "FAIL",
                       "detail": f"one={one} four={base}", "log": ""})
    dense = facts.get("rank-4-cadence-30")
    if dense is not None:
        dense_arrival = float(dense.get("arrival", {}).get(
            "exact_observer_arrival_time_s", -1))
        base_arrival = float(base.get("arrival", {}).get(
            "exact_observer_arrival_time_s", -2))
        same = all(dense.get(k) == base.get(k) for k in (
            "event_identity", "owner_fingerprint_xor",
            "owner_fingerprint_sum", "arrival_front_sha256")) and \
            dense_arrival == base_arrival and \
            abs(dense_arrival - EXPECTED_ARRIVAL_S) <= ROOT_TIME_TOLERANCE_S
        checks.append({"run": "cross-run", "name": "output-cadence-convergence",
                       "status": "PASS" if same else "FAIL",
                       "detail": f"cadence60={base} cadence30={dense}", "log": ""})

    counts = {name: sum(c["status"] == name for c in checks)
              for name in ("PASS", "FAIL", "SKIP", "ERROR")}
    summary = {"schema": "srcsep3d-positive-shock-validation-v1",
               "root": str(ROOT), "binary": str(executable),
               "binary_sha256": sha256(executable),
               "input_sha256": sha256(DECK), "event_sha256": sha256(EVENT),
               "runs": runs, "checks": checks, "counts": counts}
    (output_root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (output_root / "summary.txt").open("w") as stream:
        stream.write("Summary: " + " ".join(f"{k}={v}" for k, v in counts.items()) + "\n")
        for check in checks:
            stream.write(f"[{check['status']}] {check['run']}/{check['name']}: "
                         f"{check['detail']} log={check['log']}\n")
    print((output_root / "summary.txt").read_text(), end="")
    return 0 if counts["FAIL"] == counts["ERROR"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
