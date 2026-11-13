#!/usr/bin/env python3
"""Run I11's linked event-ledger / Na chemistry prerequisite.

See BUDGET_CONTRACT.md. These fixture checks exercise real production dispatch
and actual accepted events. They are not I11 convergence or I19 qualification.
"""
import argparse
import csv
from decimal import Decimal
import json
import math
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from integration_common import (REPO, derive_i11_input, manifest, parse_receipt,
                                run_monitored, sha256)

POPULATION = 4096
ALPHA = 0.05  # Predeclared two-sided Hoeffding coverage, not an agreement fit.
GAS_SIGNS = {"emission": 1, "photolytic_product": 1, "adsorption": -1,
             "photolytic_loss": -1, "escape": -1, "numerical_removal": -1}


def verify_ledger(path):
    """Independently reconstruct entries-minus-exits from recorded event rows.

    Decimal preserves the full printed long-double weights rather than losing
    their small residual when Python converts ~1e27 atom totals to doubles.
    Initial fixture inventory is a declared initial state, never an emission.
    """
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or rows[0]["checkpoint"] != "initial":
        raise ValueError(f"missing initial inventory checkpoint: {path}")
    previous_count = int(rows[0]["gas_model_particles"])
    previous_atoms = Decimal(rows[0]["gas_atoms"])
    maximum_count_residual = 0
    for row in rows:
        current_count = int(row["gas_model_particles"])
        current_atoms = Decimal(row["gas_atoms"])
        entries = sum(sign * int(row[name + "_count"]) for name, sign in GAS_SIGNS.items())
        atoms = sum((sign * Decimal(row[name + "_atoms"])
                     for name, sign in GAS_SIGNS.items()), Decimal(0))
        residual = current_count - previous_count - entries
        maximum_count_residual = max(maximum_count_residual, abs(residual))
        if residual != int(row["model_count_residual"]):
            raise ValueError(f"reported integer budget disagrees with events: {path}")
        bound = Decimal(row["weighted_roundoff_bound_atoms"])
        if (residual != 0 or abs(Decimal(row["weighted_residual_atoms"])) > bound or
                abs(current_atoms - previous_atoms - atoms) > bound or
                int(row["numerical_removal_count"]) != 0):
            return {"status": "FAIL", "checkpoint": row["checkpoint"],
                    "maximum_count_residual": maximum_count_residual}
        previous_count, previous_atoms = current_count, current_atoms
    return {"status": "PASS", "checkpoints": len(rows),
            "maximum_count_residual": maximum_count_residual}


def run(executable, template, output, ranks):
    output = output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError(f"refusing to overwrite budget evidence: {output}")
    output.mkdir(parents=True, exist_ok=True)
    runtime = output / "runtime"
    runtime.mkdir()
    result = {"id": "I11", "scope": "production-event-budget-prerequisite",
              "formal_i11_status": "NOT EVALUATED", "status": "ERROR"}
    try:
        input_path = output / "application.in"
        derive_i11_input(template.resolve(), input_path, 480)
        receipt = output / "receipt.txt"
        command = ["mpiexec", "-n", str(ranks), str(executable.resolve()),
                   "--moon-integration-test", "I11", "--i11-budget-check",
                   "--integration-output", str(receipt), "-input", str(input_path)]
        identity = manifest("I11", executable, command, 731993, input_path, ranks)
        identity.update(purpose=result["scope"], contract_sha256=sha256(
            Path(__file__).with_name("BUDGET_CONTRACT.md")))
        (output / "run_manifest.json").write_text(json.dumps(identity, indent=2, sort_keys=True) + "\n")
        code = run_monitored(command, runtime, dict(os.environ, OMP_NUM_THREADS="1"),
            output / "stdout.log", output / "stderr.log",
            output / "process_resource_metrics.csv")
        if code != 0:
            raise RuntimeError(f"linked budget fixture exited {code}")
        data = parse_receipt(receipt)
        if data.get("test_id") != "I11" or data.get("mpi_ranks") != ranks:
            raise ValueError("budget fixture receipt ID/MPI layout mismatch")
        half_width = math.sqrt(math.log(2.0 / ALPHA) / (2 * POPULATION))
        chemistry_ok = (data["generic_processor_control_survivors"] == 2 * POPULATION and
                        data["initial_population_per_site"] == POPULATION and
                        data["shadow_survivors"] == POPULATION and
                        abs(float(data["sunlit_survivors"]) / POPULATION - 0.5) <= half_width)
        ledgers = {name: verify_ledger(output / name / "particle_budget.csv")
                   for name in ("chemistry", "transport")}
        conditions = {"chemistry_survival": chemistry_ok,
            "chemistry_bookkeeping": data["chemistry_budget_valid"] == 1,
            "transport_bookkeeping": data["transport_budget_valid"] == 1,
            "accepted_emissions_observed": data["emission_count"] > 0,
            "surface_adsorption_observed": data["adsorption_count"] > 0,
            "outer_domain_escape_observed": data["escape_count"] > 0,
            "no_numerical_removals": data["numerical_removal_count"] == 0,
            "independent_ledger_reconstruction": all(x["status"] == "PASS" for x in ledgers.values())}
        result.update(status="PASS" if all(conditions.values()) else "FAIL",
                      conditions=conditions, metrics=data, ledger_checks=ledgers,
                      sunlit_survival_half_width=half_width)
        identity["artifacts"] = {str(p.relative_to(output)): sha256(p)
            for p in output.rglob("*") if p.is_file() and p.name in (
                "particle_budget.csv", "particle_sources.csv", "receipt.txt")}
        (output / "run_manifest.json").write_text(json.dumps(identity, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError, KeyError, RuntimeError) as error:
        result["reason"] = str(error)
    (output / "result.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(f"{result['status']} I11 production budget/Na chemistry: {output}", flush=True)
    return {"PASS": 0, "FAIL": 1, "ERROR": 2}[result["status"]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amps", type=Path, required=True)
    parser.add_argument("--input", type=Path, default=REPO / "srcMoon/examples/m1_memory200.in")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=8)
    args = parser.parse_args()
    if args.mpi_ranks < 1:
        parser.error("--mpi-ranks must be positive")
    raise SystemExit(run(args.amps, args.input, args.output_dir, args.mpi_ranks))
