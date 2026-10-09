#!/usr/bin/env python3
"""Fixture-driven tests for the revision-1.4 advanced C++ backends.

Expected values are read from benchmark_points.json at runtime.  They are not
duplicated in this test, which preserves the companion bundle as the single
machine-readable reference and makes a changed fixture digest visible to the
separate reference verifier.
"""

import argparse
import json
import math
import pathlib
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parent
FIXTURES = ROOT / "parallel_diffusion_model_data" / "benchmark_points.json"


def relative(actual, expected):
    return abs(actual / expected - 1.0) if expected else abs(actual)


def evaluate(driver, *arguments):
    completed = subprocess.run(
        [str(driver), *(str(value) for value in arguments)],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())
    fields = completed.stdout.split()
    if len(fields) != 4 or int(fields[0]) != 0:
        raise RuntimeError(f"unexpected driver output: {completed.stdout!r}")
    return tuple(float(value) for value in fields[1:])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", type=pathlib.Path, required=True)
    parser.add_argument("--json", type=pathlib.Path, required=True)
    args = parser.parse_args()
    fixtures = json.loads(FIXTURES.read_text(encoding="utf-8"))
    records = []

    def check(name, model, ok, detail):
        records.append({"name": name, "model": model, "passed": bool(ok),
                        "detail": detail})
        print(("PASS " if ok else "FAIL ") + name + " [" + model + "]: " + detail)

    try:
        value, _, _ = evaluate(args.driver, "selfcheck")
        check("PD04-08-10-API", "shared", value == 1.0,
              "spectrum, pitch-angle, adapter, table, batch, and concurrency identities")
    except Exception as error:
        check("PD04-08-10-API", "shared", False, str(error))

    qlt_tolerance = fixtures["tolerances"]["qlt_eq30"]["rel"]
    for key, point in fixtures["qlt_eq30"].items():
        try:
            parallel, _, _ = evaluate(args.driver, "qlt",
                                      point["r_L_over_ell_s"], point["eps_s2"])
            error = relative(parallel, point["lambda_par_over_ell_s_Eq30"])
            check("PD04-QLT-" + key, "qlt_slab_spectrum",
                  error <= qlt_tolerance, f"relative error {error:.3e}")
        except Exception as error:
            check("PD04-QLT-" + key, "qlt_slab_spectrum", False, str(error))

    polynomial_tolerance = fixtures["tolerances"]["nlgce_f_2014_values"]["rel"]
    nonlinear_tolerance = fixtures["tolerances"]["nlgce_n"]["rel"]
    nlgc_tolerance = fixtures["tolerances"]["nlgc_e"]["rel"]
    for label, point in fixtures["nlgce_f_2014"].items():
        state = point["inputs"]
        arguments = (state["r_L_over_ell_s"], state["f_s"],
                     state["eps2_total"], state["ell_s_over_ell_2"])
        for operation, group, tolerance, model in (
            ("nlgce_f", "nlgce_f_2014", polynomial_tolerance, "nlgce_f_2014"),
            ("nlgce_n", "nlgce_n", nonlinear_tolerance, "nlgce_n"),
            ("nlgc_e", "nlgc_e", nlgc_tolerance, "nlgc_e"),
        ):
            expected = fixtures[group][label]
            try:
                parallel, perpendicular, residual = evaluate(
                    args.driver, operation, *arguments)
                discrepancy = max(
                    relative(parallel, expected["lambda_par_over_ell_s"]),
                    relative(perpendicular, expected["lambda_perp_over_ell_s"]),
                )
                residual_ok = math.isnan(residual) or residual < 1.0e-8
                check("PD05-06-" + operation.upper() + "-" + label, model,
                      discrepancy <= tolerance and residual_ok,
                      f"max relative {discrepancy:.3e}, residual {residual:.3e}")
            except Exception as error:
                check("PD05-06-" + operation.upper() + "-" + label,
                      model, False, str(error))

    nlpa_tolerance = fixtures["tolerances"]["nlpa_given_perp"]["rel"]
    for key, point in fixtures["nlpa_given_perp"].items():
        # The four non-redundant kx=1e-3 cases cover multiple turbulence
        # regimes while keeping the default verification runtime bounded.
        if not key.endswith("kx_1e-3"):
            continue
        state = fixtures["nlgce_f_2014"][point["state"]]["inputs"]
        try:
            parallel, perpendicular, residual = evaluate(
                args.driver, "nlpa", state["r_L_over_ell_s"], state["f_s"],
                state["eps2_total"], state["ell_s_over_ell_2"],
                point["kappa_perp_supplied_over_v_ell_s"])
            discrepancy = relative(parallel, point["lambda_par_over_ell_s"])
            check("PD06-NLPA-" + key, "nlpa_given_perp",
                  discrepancy <= nlpa_tolerance and residual < 1.0e-8,
                  f"relative {discrepancy:.3e}, residual {residual:.3e}")
        except Exception as error:
            check("PD06-NLPA-" + key, "nlpa_given_perp", False, str(error))

    broadened_tolerance = fixtures["tolerances"]["broadened_slab"]["rel"]
    # One finite-width case from each kernel family exercises the nested
    # resonance/pitch-angle quadrature. The companion verifier retains the
    # full width-convergence sequence, including the narrow-width QLT limit.
    for key in ("gaussian:0.1", "lorentzian:0.1"):
        point = fixtures["broadened_slab"]["cases"][key]
        try:
            parallel, _, _ = evaluate(
                args.driver, "broadened", point["r_L_over_ell_s"],
                point["eps_s2"], point["kernel"], point["width_over_Omega"])
            discrepancy = relative(parallel, point["lambda_par_over_ell_s"])
            check("PD07-BROADENED-" + key, "broadened_slab",
                  discrepancy <= broadened_tolerance,
                  f"relative error {discrepancy:.3e}")
        except Exception as error:
            check("PD07-BROADENED-" + key, "broadened_slab", False, str(error))

    report = {
        "schema": "parallel_diffusion_advanced_test_report_v1",
        "fixture": str(FIXTURES.relative_to(ROOT)),
        "checks": len(records),
        "failed": sum(not record["passed"] for record in records),
        "records": records,
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"parallel_diffusion advanced: {'PASS' if not report['failed'] else 'FAIL'} "
          f"({report['checks']} checks)")
    return 1 if report["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())
