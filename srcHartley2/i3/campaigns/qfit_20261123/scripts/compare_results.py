#!/usr/bin/env python3
"""Compare a rerun of the I3 campaign qfit_20261123 with its reference results (../expected_results.json).

Make the rerun's numbers first (last output of each run, as for the reference):
    python3 srcHartley2/i3/tools/make_qfit_figures.py --output final --campaign RERUN_DIR --report-dir REPORT_DIR
then:
    python3 srcHartley2/i3/campaigns/qfit_20261123/scripts/compare_results.py REPORT_DIR/data/qfit_numbers_final.json

Prints rerun / reference for OH and H model/obs (geometric mean over the MAVEN/IUVS points), the O I 130.4-nm
brightness at 1e4 and 1e5 km and the case-5 IUVS/model geometric mean, and marks values that differ by more
than --tolerance (default 0.05). A rerun is stochastic (another machine, compiler or rank count gives other
random numbers), so small differences are expected; the exit status is 1 if any value is outside the tolerance.
"""
import argparse
import json
import sys
from pathlib import Path

REFERENCE = Path(__file__).resolve().parents[1] / "expected_results.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("rerun", type=Path, help="qfit_numbers_final.json of the rerun")
    parser.add_argument("--reference", type=Path, default=REFERENCE, help="default: %(default)s")
    parser.add_argument("--tolerance", type=float, default=0.05, help="relative tolerance (default 0.05)")
    args = parser.parse_args()
    ref = json.loads(args.reference.read_text())["cases"]
    new = json.loads(args.rerun.read_text())["cases"]
    worst, rows = 0.0, []
    for c in sorted(ref):
        r, n = ref[c], new.get(c)
        if n is None:
            print(f"case {c}: missing in the rerun")
            worst = float("inf")
            continue
        pairs = [("OH model/obs", r["OH_model_over_obs"]["geomean"], n["OH_model_over_obs"]["geomean"]),
                 ("H model/obs", r["H_model_over_obs"]["geomean"], n["H_model_over_obs"]["geomean"]),
                 ("O at 1e4 km (R)", r["O_KVB2022_bowen_edge_corrected_R"]["1e+04 km"],
                  n["O_KVB2022_bowen_edge_corrected_R"]["1e+04 km"]),
                 ("O at 1e5 km (R)", r["O_KVB2022_bowen_edge_corrected_R"]["1e+05 km"],
                  n["O_KVB2022_bowen_edge_corrected_R"]["1e+05 km"])]
        if "O_case5_iuvs_over_model_binmean" in r:
            pairs.append(("case-5 IUVS/model", r["O_case5_iuvs_over_model_binmean"]["geomean"],
                          n["o_case5_obs_over_model_binmean"]["geomean"]))
        for name, a, b in pairs:
            dev = abs(b / a - 1.0)
            worst = max(worst, dev)
            rows.append((c, name, a, b, dev))
        if not n.get("converged", False):
            print(f"case {c}: the rerun's last output did not converge")
    for c, name, a, b, dev in rows:
        flag = "  <-- outside tolerance" if dev > args.tolerance else ""
        print(f"case {c}  {name:18s} reference {a:10.4g}  rerun {b:10.4g}  rerun/reference {b / a:6.3f}{flag}")
    print(f"largest relative difference: {worst:.3f} (tolerance {args.tolerance})")
    return 0 if worst <= args.tolerance else 1


if __name__ == "__main__":
    sys.exit(main())
