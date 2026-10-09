#!/usr/bin/env python3
"""Generate/check specification tables against the verified numeric JSON.

This script checks printed rounding against stored equation fixtures. Run the
reference_verification.py script to independently recompute those fixtures.
"""
import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def table(headers,rows):
    return "\n".join(["| "+" | ".join(headers)+" |",
                      "| "+" | ".join(["---"]*len(headers))+" |"]+
                     ["| "+" | ".join(str(x) for x in row)+" |" for row in rows])


def tables(f):
    fmt = lambda x:format(x,".9g")
    rows = f["flpd_states"]
    values = {
        "TABLE_LENGTHS":table(["$q$","$L_U/\\ell_2$","$L_\\perp/\\ell_2$","$L_\\perp^2/(4L_U^2)$"],
                              [[fmt(r[k]) for k in ["q","L_U_over_ell_2","L_perp_over_ell_2","CLRR_shape_prefactor"]]
                               for r in f["lengths"]]),
        "TABLE_FIELDLINES":table(["State","$q$","$f_s$","$\\delta B^2/B_0^2$","$\\kappa_s$","$\\kappa_2$","$\\kappa_{\\rm FL}$"],
                                 [[r["state"]]+[fmt(r[k]) for k in ["q","slab_fraction","total_variance_over_B0_squared",
                                                                   "kappa_s_over_ell","kappa_2_over_ell","kappa_FL_over_ell"]]
                                  for r in rows]),
        "TABLE_FLPD":table(["State"]+["$\\lambda_\\parallel/\\ell_2="+fmt(lp)+"$" for lp in rows[0]["parallel_lengths"]],
                           [[r["state"]]+[fmt(x) for x in r["ratio"]] for r in rows]),
        "TABLE_ENLGC":table(["$\\lambda_\\parallel/\\ell_2$","$\\lambda_\\perp/\\ell_2$"],
                            [[fmt(r[k]) for k in ["lambda_parallel","lambda_perp"]] for r in f["enlgc"]]),
        "TABLE_RBD":table(["$\\delta B_2^2/B_0^2$","$\\lambda_\\parallel/\\lambda_2$","$\\kappa_\\perp/(v\\lambda_2)$"],
                          [[fmt(r[k]) for k in ["variance","lambda_parallel","kappa_perp_over_v_lambda_2"]] for r in f["rbd"]]
                          +[[fmt(var),"$\\infty$",fmt(f["rbd_long_limit_coefficient"]*var**0.5)] for var in [0.1,0.4,0.8]]),
        "TABLE_CLOSED":table(["$\\lambda_\\parallel/\\ell$","$\\lambda_\\perp/\\ell$"],
                             [[fmt(r[k]) for k in ["lambda_parallel_over_ell","lambda_perp_over_ell"]] for r in f["closed_form"]]),
        "TABLE_IMPLICIT_KERNEL":table(["$\\xi$","$\\mathcal K(\\xi)$","Rational relative error"],
                                      [[fmt(r[k]) for k in ["xi","K_exact","rational_relative_error"]]
                                       for r in f["implicit_slab_kernel"]]),
    }
    review = json.loads((ROOT/"review_additions_fixtures.json").read_text())
    paired = json.loads((ROOT/"paired_benchmark_points.json").read_text())
    high = lambda x:format(x,".12g")
    values.update({
        "TABLE_FLPD_SHORT":table(
            ["$q$","$\\alpha_{\\rm fluid}^{(2)}$","$\\alpha_{\\rm D4}^{(2)}$","$\\zeta$","Root $\\alpha_{\\rm FLPD}$"],
            [[high(r[k]) for k in ["q","alpha_fluid","alpha_D4","zeta","alpha_root"]]
             for r in review["short_flpd"]]),
        "TABLE_CLOSED_LENGTHS":table(
            ["$q$","$\\ell_{\\perp,\\rm int}/\\ell_2$"],
            [[high(r[k]) for k in ["q","integral_length_over_ell_2"]] for r in review["short_flpd"]]),
        "TABLE_FLPD_TAIL":table(
            ["$A$","$J(0,\\alpha)-J(A,\\alpha)$","Difference / $A^{5/6}$","Local log exponent"],
            [[high(r["A"]),high(r["J_zero_minus_J_A"]),high(r["scaled_difference"]),
              high(r["local_log_exponent"]) if r["local_log_exponent"] is not None else "unavailable"]
             for r in review["nonanalytic_tail"]]),
        "TABLE_RBD_COMPARE":table(
            ["$\\delta B_2^2/B_0^2$","$\\lambda_\\parallel/\\lambda_2$","RBD/BC","Plain RBD"],
            [[high(r[k]) for k in ["variance","lambda_parallel","corrected","plain"]]
             for r in review["rbd_plain_vs_bc"]]),
        "TABLE_NLGCE_STATES":table(
            ["State","$r$","$f_s$","$\\epsilon^2$","$\\ell_s/\\ell_2$"],
            [[r["label"]]+[format(r[k],".17g") for k in ["r","slab_fraction","epsilon_squared","slab_to_2d_length_ratio"]]
             for r in paired["nlgce"]])
    })
    return values


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--document",type=Path,default=ROOT.parent/"PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md")
    p.add_argument("--write",action="store_true")
    args = p.parse_args()
    fixtures = json.loads((ROOT/"benchmark_points.json").read_text())
    document = args.document.read_text()
    for name,content in tables(fixtures).items():
        block = f"<!-- BEGIN GENERATED {name} -->\n{content}\n<!-- END GENERATED {name} -->"
        pattern = rf"<!-- BEGIN GENERATED {name} -->.*?<!-- END GENERATED {name} -->"
        existing = re.findall(pattern,document,re.S)
        marker = f"<!-- {name} -->"
        if args.write:
            if len(existing)==1:
                document = document.replace(existing[0],block)
            elif not existing and document.count(marker)==1:
                document = document.replace(marker,block)
            else:
                raise AssertionError(f"Missing/duplicate table {name}")
        elif existing != [block]:
            raise AssertionError(f"Printed table mismatch {name}")
    if args.write:
        args.document.write_text(document)
    print(f"{'Wrote' if args.write else 'Verified'} {len(tables(fixtures))} specification tables against numeric fixtures.")


if __name__ == "__main__":
    main()
