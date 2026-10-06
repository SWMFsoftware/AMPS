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
# The maintained positive event fixes numerics.azimuth_cells=48.  Its finite
# SSE cap is a topological disk, so exactly these 48 rim edges—and no apex or
# seam edges—must occur once in the triangle incidence map.
EXPECTED_FRONT_BOUNDARY_EDGES = 48


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
    """Verify triangular geometry and one-sided RH limits face by face.

    Production coordinates are nodal while shock state is cell-centred BLOCK
    data. Non-shock patches must
    carry the declared neutral display encoding (Mf=0, X=CB=1 and downstream
    columns copied from upstream) while both acceptance/downstream flags stay
    false.  Conversely, every accepted patch must carry finite, compressive
    upstream/downstream states.  This prevents a syntactically valid surface
    file from passing because only its header named the requested quantities.
    """
    if not path.exists():
        return False, f"missing {path}"
    lines = path.read_text(errors="replace").splitlines()
    if len(lines) < 4:
        return False, "surface is truncated before its data"
    title, names, zone = lines[0], tecplot_variables(lines[1]), lines[2]
    node_match = re.search(r"\bN=(\d+)", zone)
    face_match = re.search(r"\bE=(\d+)", zone)
    required = ("normal_x", "normal_y",
                    "normal_z", "normal_speed_m_s", "time_s", "generation",
                    "quadrature_area_m2", "planar_chord_area_m2",
                    "triangle_stable_id",
                    "status_code", "shock_accepted", "upstream_valid",
                    "downstream_valid", "fast_mach", "theta_Bn_rad",
                    "theta_Bn_valid", "density_compression",
                    "magnetic_compression", "magnetic_compression_valid",
                    "rho1_kg_m3", "p1_Pa", "u1x_m_s", "u1y_m_s",
                    "u1z_m_s", "b1x_T", "b1y_T", "b1z_T",
                    "rho2_kg_m3", "p2_Pa", "u2x_m_s", "u2y_m_s",
                    "u2z_m_s", "b2x_T", "b2y_T", "b2z_T")
    if "RH limits valid only when shock_accepted=1" not in title or \
            not node_match or not face_match or \
            "DATAPACKING=BLOCK" not in zone or \
            "ZONETYPE=FETRIANGLE" not in zone or \
            f"VARLOCATION=([4-{len(names)}]=CELLCENTERED)" not in zone or \
            any(name not in names for name in ("x_m", "y_m", "z_m") + required):
        return False, "surface title, triangular BLOCK contract, or columns are absent"
    node_count, face_count = int(node_match.group(1)), int(face_match.group(1))
    line_index = 3
    while line_index < len(lines) and lines[line_index].startswith("AUXDATA"):
        line_index += 1
    tokens = " ".join(lines[line_index:]).split()
    blocks: dict[str, list[float]] = {}
    cursor = 0
    try:
        for variable, name in enumerate(names):
            count = node_count if variable < 3 else face_count
            words = tokens[cursor:cursor + count]
            if len(words) != count:
                return False, f"truncated block {name}"
            blocks[name] = [float(word) for word in words]
            cursor += count
    except ValueError:
        return False, "non-numeric surface variable block"
    connectivity: list[tuple[int, int, int]] = []
    try:
        for _face in range(face_count):
            words = tokens[cursor:cursor + 3]
            if len(words) != 3:
                return False, "truncated triangle connectivity"
            connectivity.append(tuple(int(word) - 1 for word in words))
            cursor += 3
    except ValueError:
        return False, "non-integer triangle connectivity"
    if cursor != len(tokens):
        return False, "unexpected tokens after triangle connectivity"
    finite_blocks = all(math.isfinite(value) for values in blocks.values()
                        for value in values)
    coordinates = list(zip(blocks["x_m"], blocks["y_m"], blocks["z_m"]))
    edges: dict[tuple[int, int], int] = {}
    total_planar_area = 0.0
    total_curved_area = 0.0
    topology = finite_blocks and node_count > 0 and face_count > 0
    topology_problem = "" if topology else "non-finite data or empty topology"
    for face, vertices in enumerate(connectivity):
        if len(set(vertices)) != 3 or not all(
                0 <= vertex < node_count for vertex in vertices):
            topology = False
            topology_problem = f"face {face + 1} has invalid indices {vertices}"
            break
        a, b, c = (coordinates[vertex] for vertex in vertices)
        ab = tuple(b[i] - a[i] for i in range(3))
        ac = tuple(c[i] - a[i] for i in range(3))
        cross = (ab[1] * ac[2] - ab[2] * ac[1],
                 ab[2] * ac[0] - ab[0] * ac[2],
                 ab[0] * ac[1] - ab[1] * ac[0])
        planar = 0.5 * math.sqrt(sum(value * value for value in cross))
        normal = tuple(blocks[name][face] for name in
                       ("normal_x", "normal_y", "normal_z"))
        published_planar = blocks["planar_chord_area_m2"][face]
        orientation = sum(cross[i] * normal[i] for i in range(3))
        area_error = abs(planar - published_planar)
        area_bound = 1.0e-12 * max(planar, 1.0)
        curved = blocks["quadrature_area_m2"][face]
        # A chart-diagonal is curved on the sphere whereas the Tecplot edge is
        # straight.  The flat triangle is the minimum surface for its straight
        # boundary, not for that curved chart boundary, so its area need not be
        # smaller than the assigned half-cell curved measure face by face.
        # Requiring that false local ordering rejected valid checkerboard
        # triangles.  Reproduction of the stored chord area is checked here;
        # positivity and the global inscribed-mesh area deficit are checked
        # independently without changing any physics tolerance.
        if orientation <= 0.0 or area_error > area_bound or not (
                published_planar > 0.0 and curved > 0.0):
            topology = False
            topology_problem = (
                f"face {face + 1} geometry mismatch: orientation={orientation} "
                f"recomputed_planar={planar} published_planar={published_planar} "
                f"area_error={area_error} area_bound={area_bound} curved={curved}")
            break
        total_planar_area += published_planar
        total_curved_area += curved
        for left, right in ((vertices[0], vertices[1]),
                            (vertices[1], vertices[2]),
                            (vertices[2], vertices[0])):
            edge = tuple(sorted((left, right)))
            edges[edge] = edges.get(edge, 0) + 1
    boundary = [edge for edge, count in edges.items() if count == 1]
    if topology:
        incidence_valid = all(count in (1, 2) for count in edges.values())
        euler = node_count - len(edges) + face_count
        topology = incidence_valid and euler == 1 and \
            len(boundary) == EXPECTED_FRONT_BOUNDARY_EDGES and \
            0.0 < total_planar_area < total_curved_area
        if not topology:
            topology_problem = (
                f"global topology mismatch: incidence_valid={incidence_valid} "
                f"euler={euler} boundary_edges={len(boundary)} "
                f"planar_area={total_planar_area} curved_area={total_curved_area}")

    accepted = 0
    rejected_with_ambient = 0
    good = topology
    for face in range(face_count):
        values = {name: blocks[name][face] for name in required}
        normal = math.sqrt(sum(values[name] ** 2 for name in (
            "normal_x", "normal_y", "normal_z")))
        good = good and math.isfinite(normal) and abs(normal - 1.0) <= 1.0e-12
        good = good and abs(values["time_s"] -
                            float(receipt["time_s"])) <= 1.0e-9
        good = good and int(values["generation"]) == \
            int(receipt["front_generation"])
        good = good and values["triangle_stable_id"] == face + 1
        is_accepted = values["shock_accepted"] == 1.0
        good = good and is_accepted == \
            (int(values["status_code"]) == 5)
        if is_accepted:
            accepted += 1
            good = good and values["upstream_valid"] == 1.0 and \
                values["downstream_valid"] == 1.0 and \
                values["fast_mach"] > 1.0 and \
                0.0 <= values["theta_Bn_rad"] <= 0.5 * math.pi and \
                values["theta_Bn_valid"] == 1.0 and \
                values["density_compression"] > 1.0 and \
                values["magnetic_compression"] > 0.0 and \
                values["magnetic_compression_valid"] == 1.0 and \
                values["rho1_kg_m3"] > 0.0 and values["p1_Pa"] > 0.0 and \
                values["rho2_kg_m3"] > 0.0 and values["p2_Pa"] > 0.0
        else:
            upstream_downstream = (
                ("rho1_kg_m3", "rho2_kg_m3"), ("p1_Pa", "p2_Pa"),
                ("u1x_m_s", "u2x_m_s"), ("u1y_m_s", "u2y_m_s"),
                ("u1z_m_s", "u2z_m_s"), ("b1x_T", "b2x_T"),
                ("b1y_T", "b2y_T"), ("b1z_T", "b2z_T"))
            ambient_fill = values["upstream_valid"] == 1.0 and \
                values["downstream_valid"] == 0.0 and \
                values["fast_mach"] == 0.0 and \
                values["density_compression"] == 1.0 and \
                values["magnetic_compression"] == 1.0 and \
                values["magnetic_compression_valid"] == 0.0 and \
                all(values[upstream] == values[downstream]
                    for upstream, downstream in upstream_downstream)
            good = good and ambient_fill
            rejected_with_ambient += int(ambient_fill)
    good = good and accepted > 0
    return good, (f"nodes={node_count} triangles={face_count} "
                  f"boundary_edges={len(boundary)} accepted={accepted} "
                  f"no_shock_ambient_fill={rejected_with_ambient} "
                  f"topology_problem={topology_problem or 'none'}")


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
    # TITLE, VARIABLES, ZONE and five AUXDATA records precede the first block.
    # Read the complete metadata contract without loading the surface arrays.
    front_head = read_header(front, 8)
    check("visit-tecplot-products", ambient.exists() and front.exists() and
          "mass_density_kg_m-3" in ambient_head and
          "background_generation" in ambient_head and
          "ZONETYPE=FEBRICK" in ambient_head and
          "ZONETYPE=FETRIANGLE" in front_head and
          "DATAPACKING=BLOCK" in front_head and
          "surface_topology=\"triangular-sse-cap-v1\"" in front_head and
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
