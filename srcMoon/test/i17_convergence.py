"""Numerical error contract for June I17, approved before residual inspection.

Only differences between production solutions enter this test. Observational
one-sigma uncertainties set the approved numerical error budget; observed
temperatures never enter the convergence formula. This gate does not establish
observational agreement or spatial convergence.
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

# Approved 2026-10-14 before any June residual inspection. This controls
# numerical truncation error only and is not an observational fit tolerance.
OBSERVATIONAL_SIGMA_FRACTION = 0.1


def compare_time_steps(directories: list[Path], output: Path) -> tuple[bool, dict]:
    """Compare dt, dt/2, dt/4 at exactly identical June site/time bins.

    Every bin must have |T(dt/4)-T(dt/2)| < 0.1*sigma_obs, and this difference
    must be smaller than |T(dt/2)-T(dt)|. Exact equality of all three computed
    temperatures is allowed as zero measured sensitivity. No tolerance is
    added for noisy/nonmonotone convergence. Missing, duplicate, mismatched,
    or nonfinite rows are infrastructure errors, not convergence failures.
    """
    if len(directories) != 3:
        raise ValueError("I17 convergence requires three time steps")
    tables = []
    for directory in directories:
        with (directory / "receipt.txt.I17_training_model.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        table = {row["site_time_bin_id"]: row for row in rows}
        if not table or len(table) != len(rows):
            raise ValueError("I17 convergence has empty or duplicate bin identities")
        tables.append(table)
    if any(set(table) != set(tables[0]) for table in tables[1:]):
        raise ValueError("I17 time-step runs sampled different site/time bins")
    metrics = []
    for key, coarse in tables[0].items():
        rows = [table[key] for table in tables]
        # Comparing different positions, UTCs, or observations could mimic
        # time-step sensitivity. Require byte-identical exported metadata.
        for field in ("utc", "ephemeris_time_tdb_s", "latitude_deg",
                      "longitude_east_deg", "observed_uncertainty_K",
                      "observed_bolometric_brightness_temperature_K"):
            if any(row[field] != coarse[field] for row in rows[1:]):
                raise ValueError(f"I17 {key}: changed {field} across time steps")
        values = [float(row["model_bolometric_brightness_temperature_K"])
                  for row in rows]
        sigma = float(coarse["observed_uncertainty_K"])
        if not all(math.isfinite(value) and value > 0 for value in values + [sigma]):
            raise ValueError(f"I17 {key}: invalid temperature or uncertainty")
        coarse_difference = abs(values[1] - values[0])
        fine_difference = abs(values[2] - values[1])
        decreasing = (fine_difference < coarse_difference or
                      fine_difference == coarse_difference == 0.0)
        passed = fine_difference < OBSERVATIONAL_SIGMA_FRACTION * sigma and decreasing
        metrics.append({"site_time_bin_id": key,
                        "coarse_difference_K": coarse_difference,
                        "fine_difference_K": fine_difference,
                        "numerical_budget_K": OBSERVATIONAL_SIGMA_FRACTION * sigma,
                        "decreasing": decreasing, "pass": passed})
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]))
        writer.writeheader()
        writer.writerows(metrics)
    failed = sum(not row["pass"] for row in metrics)
    return failed == 0, {
        "compared_bins": len(metrics), "failed_bins": failed,
        "maximum_fine_difference_K": max(row["fine_difference_K"] for row in metrics),
        "convergence_csv": str(output),
        "observational_sigma_fraction": OBSERVATIONAL_SIGMA_FRACTION,
    }
