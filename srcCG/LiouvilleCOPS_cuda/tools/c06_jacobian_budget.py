#!/usr/bin/env python3
"""Evidence for the C06 Jacobian budget: per-entry error in units of eps·M_J (CUDA env).

VALIDATION.md left the Jacobian budget provisional (1e-11 × column max). The
planned replacement is a per-entry bound like the pressure one,
``|ΔJ_ij| <= K eps M_J,ij``, where ``M_J,ij`` is the sum of the magnitudes of every
contribution entering entry (i, j), computed with the reference's own
per-contribution expressions (``fit_cops.Model.evaluate``, ``jac=True``):

* activity/template column B_j:
  NG ``|cNG bm| Σ|wn u [act>0] B_j|``; RG ``|cRG| Σ|(wp + bm wv u) [act>0] B_j|``;
* temperatures: ``Σ|w dl| + KT Σ|w dl x|`` (Tmax), ``KT Σ|w dl x|`` (T75), ``Σ|w dn|`` (Tmin),
  for the NG weight gN, and for RG ``|cRG_H| (H2O sums) + |κ cRG_C| (CO2 sums)``;
* κ: the CO2 pressure magnitudes.

The GPU values are those of the tree mode on a CPU-loadable subset. Point-set
invariance (bitwise; ``c06_full_run.py``) makes them the full-set values at
those points. Subsets: the C00b 47-point subset, and the 23-point verification
subset of the full run. All C00b fixture states for B29/N30/C30/NC31/R30.

Reports per configuration the maximum of ``|ΔJ| / (eps M_J)`` and, for
comparison, the column-max-relative error. Entries with ``M_J == 0`` must be
exactly zero on both sides. Writes
``$COPS_OUTPUT_ROOT/cuda_c06/jacobian_budget_<UTC>/report.json``.
"""

import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))
sys.path.insert(0, str(APP_ROOT / "tools"))

DEFAULT_OUTPUT_ROOT = "/data/vtenishe/rosetta_data_end_of_mission_cops/cops_output"
EPS = float(np.finfo(float).eps)


def log(message):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), message), flush=True)


def jacobian_magnitudes(P, M, x, request, fit_cops):
    """{key: (M_NG, M_RG)} per requested column, from the reference's contribution terms."""
    a, b, T, kappa = P.unpack(x)
    act_lin = M.Y @ a
    if P.nb:
        act_lin = act_lin + P.enhancement_scale * (P.enhancement_basis @ b)
    active = act_lin > 0.0
    af = np.maximum(act_lin, 0.0)[M.face]
    Tmax, T75, Tmin = T
    Tl = Tmax + fit_cops.KT * (Tmax - T75) * M.xlit
    use_lit = M.lit & (Tl > Tmin)
    Te = np.where(use_lit, Tl, Tmin)
    u = 1.0 / np.sqrt(Te)
    bc = lambda w: np.bincount(M.pt, weights=np.abs(w), minlength=M.npt)
    cNG, cRG, bm, KT = M.cNG, M.cRG, M.bm, fit_cops.KT
    out = {}
    wN = M.wn["H2O"] * u
    wR = M.wp["H2O"] + bm["H2O"] * M.wv["H2O"] * u
    for key in request.keys():
        kind, j = key
        if kind in ("a", "b"):
            column = (M.Y[:, j] if kind == "a" else
                      P.enhancement_scale * P.enhancement_basis[:, j])
            mb = (active * column)[M.face]
            out[key] = (abs(cNG["H2O"] * bm["H2O"]) * bc(wN * mb),
                        np.abs(cRG["H2O"]) * bc(wR * mb))
    du = -0.5 * u / Te
    dl, dn = np.where(use_lit, du, 0.0), np.where(use_lit, 0.0, du)
    dlx = dl * M.xlit

    def mags(w):
        # KT = 1/(1 - 1/cos75) < 0: magnitudes take |KT| (the reference forms
        # KT*s_lx + s_l and -KT*s_lx from the two separately folded sums).
        s_l, s_lx, s_n = bc(w * dl), bc(w * dlx), bc(w * dn)
        return s_l + abs(KT) * s_lx, abs(KT) * s_lx, s_n
    if request.temperature:
        # gN as the reference forms it; its magnitude bounds both of its parts.
        gN = cNG["H2O"] * bm["H2O"] * M.wn["H2O"] * af + kappa * cNG["CO2"] * bm["CO2"] * M.wn["CO2"]
        tN = mags(gN)
        tRH = mags(bm["H2O"] * M.wv["H2O"] * af)
        tRC = mags(bm["CO2"] * M.wv["CO2"])
        for k in range(3):
            out[("T", k)] = (tN[k], np.abs(cRG["H2O"]) * tRH[k] + abs(kappa) * np.abs(cRG["CO2"]) * tRC[k])
    if request.kappa:
        mC = abs(cNG["CO2"] * bm["CO2"]) * bc(M.wn["CO2"] * u)
        mRC = np.abs(cRG["CO2"]) * (bc(M.Q["CO2"] * M.ACO2[M.face] * M.KPhi) +
                                    bm["CO2"] * bc(M.wv["CO2"] * u))
        out[("kappa", 0)] = (mC, mRC)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    import c00_capture_reference as c00
    from cops_cuda import backend as bk, evaluator, manifest, provenance, reference
    from cops_cuda.memory import DeviceBudget
    root = reference.bind()
    capture = Path(args.capture).resolve()
    fixtures = json.loads((capture / "joint" / "joint_fixtures.json").read_text())
    templates = np.load(str(capture / "joint" / "templates.npz"))
    modules = c00.load_reference_modules()
    fit_cops = modules["fit_cops"]
    parse_utc = modules["cops_selection"].parse_utc
    run = c00.base_directory(modules, c00.DEFAULT_BASE)
    header, names, rows, columns, _, _, _ = c00.fit_selection(modules, run, c00.DEFAULT_PROFILE)
    coef = fit_cops.read_activity_coefficients(str(root), header)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    library = {n: templates[n].astype(float)
               for n in ("late_ng_not_rg", "late_rg_not_ng", "syn_cliff_standin")}
    subsets = {"c00b_47": c00.choose_subset(rows, columns, 40, c00.LATE_INTERVAL, parse_utc),
               "full_run_verification_23": c00.choose_subset(rows, columns, 16, c00.LATE_INTERVAL,
                                                            parse_utc)}
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c06" /
                                        ("jacobian_budget_" + provenance.utc_stamp()))
    budget = DeviceBudget(args.device)
    report = {}
    try:
        for sname, subset in subsets.items():
            M = fit_cops.Model(str(run), header, columns, [rows[k] for k in subset], elements, coef)
            M.load_kernels(M.sel)
            per_config = {}
            for label, template_names, role in c00.JOINT_CONFIGS:
                record = fixtures["configs"][label]
                arrays = np.load(str(capture / "joint" / ("joint_%s.npz" % label)))
                config = {"name": label, "degree": 4, "enhancement_degree": 4,
                          "fix_legacy_activity": False, "fix_temperature": False,
                          "fix_kappa": False, "ordered_temperature": False,
                          "templates": bool(template_names)}
                tpl = (list(template_names),
                       np.column_stack([library[n] for n in template_names]) if template_names
                       else np.empty((len(elements["area"]), 0)), {})
                P = c00.make_problem(fit_cops, M, elements, coef, config, tpl)
                request = bk.ColumnRequest.from_problem(P)
                ref = bk.ReferenceBackend.from_problem(P)
                worst, worst_key, colmax_rel, zero_ok = 0.0, None, 0.0, True
                by_kind = {}
                with evaluator.from_problem(P, run, device_id=args.device, mode="tree",
                                            budget=budget) as tree:
                    for v in record["vectors"]:
                        x = arrays["%s__%s__x" % (label, v)]
                        state = bk.ParameterState.from_problem(P, x)
                        got, want = tree.evaluate(state, request), ref.evaluate(state, request)
                        mags = jacobian_magnitudes(P, M, x, request, fit_cops)
                        if any(np.any(m < 0.0) for pair in mags.values() for m in pair):
                            raise SystemExit("negative Jacobian magnitude (formula error)")
                        for key in request.keys():
                            for side, (g, w, m) in enumerate(zip(got.columns[key], want.columns[key],
                                                                 mags[key])):
                                d = np.abs(g - w)
                                zero = m == 0.0
                                zero_ok = zero_ok and bool(np.all(g[zero] == 0.0) and np.all(w[zero] == 0.0))
                                r = float(np.max(d[~zero] / m[~zero]) / EPS) if np.any(~zero) else 0.0
                                kind = key[0]
                                by_kind[kind] = max(by_kind.get(kind, 0.0), r)
                                if r > worst:
                                    worst, worst_key = r, [v, list(key), ("NG", "RG")[side]]
                                scale = max(float(np.max(np.abs(w))), 1e-300)
                                colmax_rel = max(colmax_rel, float(np.max(d)) / scale)
                per_config[label] = {"max_eps_MJ": worst, "worst_entry": worst_key,
                                     "max_eps_MJ_by_column_kind": by_kind,
                                     "max_rel_to_column_max": colmax_rel,
                                     "zero_magnitude_entries_exactly_zero": zero_ok,
                                     "states": len(record["vectors"])}
                log("%s %s: max %.1f eps*M_J (%s), by kind %s, col-max rel %.2e, zeros %s" %
                    (sname, label, worst, worst_key, {k: round(v, 1) for k, v in by_kind.items()},
                     colmax_rel, zero_ok))
            report[sname] = {"points": int(M.npt), "contributions": int(M.indptr[-1]),
                             "configs": per_config}
            del M
    finally:
        budget.close()
    overall = max(c["max_eps_MJ"] for s in report.values() for c in s["configs"].values())
    doc = {"schema": "liouville_cops_cuda.c06_jacobian_budget", "schema_version": 1,
           "created_utc": provenance.utc_now(), "capture": str(capture), "subsets": report,
           "max_eps_MJ_overall": overall,
           "manifest": manifest.backend_manifest("cuda_fp64_v1", args.device)}
    provenance.write_json_exclusive(out / "report.json", doc)
    log("overall max %.1f eps*M_J; report %s" % (overall, out / "report.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
