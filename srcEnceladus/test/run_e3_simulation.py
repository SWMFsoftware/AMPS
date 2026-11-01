#!/usr/bin/env python3
"""Run and validate the one-rank E3 H2O+CO2 simulation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import plot_flyby_density


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=REPO / "amps")
    parser.add_argument("--input", type=Path,
                        default=REPO / "input" / "enceladus.in")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "test_output" / "srcEnceladus" /
                        "e3_10000_run")
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--launcher", choices=("mpiexec", "direct"),
                        default="mpiexec",
                        help="use a normal MPI launcher or MPI singleton mode")
    parser.add_argument("--timeout-s", type=float, default=7200.0)
    parser.add_argument("--expected-iterations", type=int, default=10000)
    parser.add_argument("--expected-sample-length", type=int, default=1000)
    parser.add_argument("--particles-per-species-step", type=int, default=100)
    parser.add_argument("--maximum-flyby-cell-edge-m", type=float,
                        default=20000.0,
                        help="maximum allowed cell edge at every in-domain flyby row")
    parser.add_argument("--skip-flyby-density-check", action="store_true",
                        help="do not require/plot the E3 flyby Tecplot product")
    args = parser.parse_args()
    if (args.ranks <= 0 or args.timeout_s <= 0 or
            args.expected_iterations <= 0 or args.expected_sample_length <= 0 or
            args.particles_per_species_step <= 0 or
            args.maximum_flyby_cell_edge_m <= 0):
        parser.error("rank, timeout, iteration, sample, and particle values must be positive")
    if args.expected_iterations % args.expected_sample_length != 0:
        parser.error("--expected-iterations must be divisible by --expected-sample-length")
    if args.launcher == "direct" and args.ranks != 1:
        parser.error("--launcher=direct supports exactly one MPI rank")

    executable = args.executable.resolve()
    input_file = args.input.resolve()
    output = args.output.resolve()
    log = output / "run.log"
    receipt_path = output / "test_receipt.json"
    if log.exists() or receipt_path.exists():
        print(f"ERROR: refusing to overwrite prior run evidence in {output}",
              file=sys.stderr)
        return 2
    (output / "PT" / "plots").mkdir(parents=True, exist_ok=True)

    if args.launcher == "mpiexec":
        command = ["mpiexec", "-n", str(args.ranks), str(executable),
                   "--input", str(input_file)]
    else:
        # Open MPI singleton mode still executes the complete MPI-native main
        # path, but it needs no out-of-band network interface for one rank.
        # This is useful in restricted validation containers where mpiexec
        # cannot open even the loopback transport.
        command = [str(executable), "--input", str(input_file)]
    started = time.time()
    timed_out = False
    with log.open("x", encoding="utf-8") as stream:
        process = subprocess.Popen(command, cwd=output, text=True,
                                   stdout=stream, stderr=subprocess.STDOUT)
        next_update = started + 30.0
        while process.poll() is None:
            now = time.time()
            if now - started > args.timeout_s:
                timed_out = True
                process.terminate()
                try:
                    process.wait(timeout=15.0)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
                break
            if now >= next_update:
                print(f"RUNNING elapsed_s={now-started:.1f}", flush=True)
                next_update = now + 30.0
            time.sleep(1.0)
    elapsed = time.time() - started
    returncode = process.returncode
    text = log.read_text(encoding="utf-8", errors="replace")

    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, detail: str) -> None:
        checks.append({"test": name, "status": "PASS" if passed else "FAIL",
                       "detail": detail})

    check("E3R00_process", not timed_out and returncode == 0,
          f"returncode={returncode} timed_out={timed_out}")
    check("E3R01_configuration",
          "source_case=tenishev2014_hybrid_8" in text and
          "plasma_background_mode=none" in text and
          "reaction_products=discard" in text and
          f"sample_length_steps={args.expected_sample_length}" in text and
          f"total_iterations={args.expected_iterations}" in text and
          (args.skip_flyby_density_check or
           ("flyby_density_output=on" in text and
            "flyby_mesh_refinement=on" in text)),
          "expected E3 neutral/loss-only run receipt")
    dt_match = re.search(r"global_dt_s=([0-9.eE+-]+)", text)
    dt = float(dt_match.group(1)) if dt_match else None
    check("E3R02_global_time_step", dt is not None and dt > 0.0,
          f"global_dt_s={dt!r}")
    final_match = re.search(
        r"final_receipt completed_iterations=(\d+).*?"
        r"injection_step_min=(\d+) injection_step_max=(\d+).*?"
        r"resident_macroparticles=(\d+)", text)
    final_values = tuple(int(value) for value in final_match.groups()) \
        if final_match else None
    check("E3R03_completed_window", final_values is not None and
          final_values[:3] == (args.expected_iterations,
                               args.expected_iterations,
                               args.expected_iterations),
          f"final_receipt={final_values!r}")

    species_values: dict[str, dict[str, int]] = {}
    pattern = re.compile(
        r"final_species=(\w+) injected_macroparticles=(\d+) "
        r"expected_injected_macroparticles=(\d+) "
        r"photolytic_parent_losses=(\d+) surface_losses=(\d+) "
        r"domain_escapes=(\d+) resident_macroparticles=(\d+) "
        r"ledger_residual=(-?\d+)")
    for match in pattern.finditer(text):
        species_values[match.group(1)] = {
            "injected": int(match.group(2)),
            "expected": int(match.group(3)),
            "photolytic_losses": int(match.group(4)),
            "surface_losses": int(match.group(5)),
            "domain_escapes": int(match.group(6)),
            "resident": int(match.group(7)),
            "ledger_residual": int(match.group(8)),
        }
    expected_injected = (args.expected_iterations *
                         args.particles_per_species_step)
    injection_ok = (set(species_values) == {"H2O", "CO2"} and
                    all(value["injected"] == expected_injected and
                        value["expected"] == expected_injected
                        for value in species_values.values()))
    check("E3R04_exact_injection", injection_ok,
          json.dumps(species_values, sort_keys=True))
    check("E3R05_co2_transport_only",
          species_values.get("CO2", {}).get("photolytic_losses") == 0,
          json.dumps(species_values.get("CO2", {}), sort_keys=True))
    conservation_ok = (set(species_values) == {"H2O", "CO2"} and
                       all(value["ledger_residual"] == 0
                           for value in species_values.values()) and
                       sum(value["resident"]
                           for value in species_values.values()) ==
                       (final_values[3] if final_values is not None else -1))
    check("E3R06_particle_conservation", conservation_ok,
          json.dumps(species_values, sort_keys=True))
    check("E3R07_no_plasma_output", not (output / "ices.data.dat").exists(),
          "ices.data.dat must not be created in plasma-free mode")
    field_species = {
        match.group(1)
        for path in (output / "PT" / "plots").glob("pic.*.out=*.dat")
        if (match := re.fullmatch(r"pic\.([^.]+)\.s=\d+\.out=\d+\.dat",
                                 path.name)) is not None
    }
    check("E3R08_active_species_output", field_species == {"H2O", "CO2"},
          f"generic_field_species={sorted(field_species)}")
    check("E3R09_clean_termination", "End of the run:" in text,
          "native main reached MPI finalization")

    flyby_tecplot = output / "PT" / "plots" / "flyby_density.E3.dat"
    flyby_png = output / "PT" / "plots" / "flyby_density.E3.png"
    flyby_summary: dict[str, object] | None = None
    if not args.skip_flyby_density_check:
        try:
            variables, zones = plot_flyby_density.read_tecplot(flyby_tecplot)
            expected_zones = (args.expected_iterations //
                              args.expected_sample_length)
            output_sequence = [zone.output for zone in zones]
            row_counts = [len(zone.rows) for zone in zones]
            sample_lengths = [zone.sample_steps for zone in zones]
            domain_index = variables.index("in_domain")
            cell_edge_index = variables.index("cell_max_edge_m")
            domain_counts = [sum(row[domain_index] == 1.0
                                 for row in zone.rows) for zone in zones]
            in_domain_cell_edges = [
                row[cell_edge_index]
                for zone in zones for row in zone.rows
                if row[domain_index] == 1.0]
            tecplot_ok = (
                len(zones) == expected_zones and
                output_sequence == list(range(expected_zones)) and
                all(count == 3201 for count in row_counts) and
                sample_lengths == [args.expected_sample_length]*expected_zones and
                bool(domain_counts) and len(set(domain_counts)) == 1 and
                domain_counts[0] > 0 and bool(in_domain_cell_edges))
            check("E3R10_flyby_tecplot", tecplot_ok,
                  f"zones={len(zones)} outputs={output_sequence} "
                  f"rows={row_counts} in_domain={domain_counts}")
            resolution_ok = (
                bool(in_domain_cell_edges) and
                min(in_domain_cell_edges) > 0.0 and
                max(in_domain_cell_edges) <=
                    args.maximum_flyby_cell_edge_m*(1.0+1.0e-12))
            check("E3R11_flyby_mesh_resolution", resolution_ok,
                  f"min_cell_edge_m={min(in_domain_cell_edges) if in_domain_cell_edges else None} "
                  f"max_cell_edge_m={max(in_domain_cell_edges) if in_domain_cell_edges else None} "
                  f"limit_m={args.maximum_flyby_cell_edge_m}")
            plotted_zones, plotted_domain = plot_flyby_density.plot_density(
                flyby_tecplot, flyby_png)
            png_ok = (flyby_png.is_file() and flyby_png.stat().st_size > 0 and
                      plotted_zones == expected_zones and
                      plotted_domain == domain_counts[-1])
            check("E3R12_flyby_png", png_ok,
                  f"zones={plotted_zones} in_domain={plotted_domain} "
                  f"bytes={flyby_png.stat().st_size if flyby_png.exists() else 0}")
            flyby_summary = {
                "tecplot": str(flyby_tecplot),
                "png": str(flyby_png),
                "zones": len(zones),
                "rows_per_zone": row_counts,
                "in_domain_points_per_zone": domain_counts,
                "minimum_in_domain_cell_edge_m": min(in_domain_cell_edges),
                "maximum_in_domain_cell_edge_m": max(in_domain_cell_edges),
                "variables": variables,
            }
        except Exception as exception:
            check("E3R10_flyby_tecplot", False, repr(exception))
            check("E3R11_flyby_mesh_resolution", False,
                  "Tecplot resolution column not validated")
            check("E3R12_flyby_png", False, "plot not generated")

    passed = all(item["status"] == "PASS" for item in checks)
    receipt = {
        "runner": str(Path(__file__).resolve()),
        "command": command,
        "working_directory": str(output),
        "input": str(input_file),
        "executable": str(executable),
        "elapsed_s": elapsed,
        "returncode": returncode,
        "timed_out": timed_out,
        "status": "PASS" if passed else "FAIL",
        "global_time_step_s": dt,
        "final_receipt": final_values,
        "species": species_values,
        "expected_iterations": args.expected_iterations,
        "expected_sample_length": args.expected_sample_length,
        "flyby_density": flyby_summary,
        "log": str(log),
        "checks": checks,
    }
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    for item in checks:
        print(f"{item['status']} {item['test']} {item['detail']}")
    print(f"SUMMARY status={receipt['status']} elapsed_s={elapsed:.1f} "
          f"receipt={receipt_path}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
