#!/usr/bin/env python3
"""C06 full-set run: the CUDA evaluator on every fit point of the 788k base (CUDA env).

The C06 numerical gate (``c06_evaluator_check.py``) ran on the 47-point C00b
subset, because the CPU reference ``fit_cops.Model`` holds every kernel
record on the host. This tool runs the production problem on the full fit
selection:

* the problem is built by the unmodified ``fit_cops.Problem`` on a
  **header-only model** (``cops_cuda.host_model``). Kernel records exist only on
  the device, so host memory stays small;
* the backend is ``CudaBackend`` in ``tree`` mode, installed with
  ``backend.install``.

Verification against the CPU uses per-point independence. Each point's
pressures and Jacobian rows depend only on that point's records (the
reference's ``np.bincount`` folds each bin in record order, and the tree
blocks are fixed per point). A CPU reference ``Model`` on a small
**verification subset** (``c00.choose_subset``, evenly spaced plus late rows) is
therefore an exact reference for those points of the full problem. Per
configuration and fixture state:

``subset_sequential_bitwise``
    a GPU sequential backend on the subset equals the CPU reference bitwise;
``point_set_invariance_bitwise``
    a tree backend on the subset equals the full tree backend bitwise at the
    subset points (the full set adds no cross-point effect);
``full_tree_pressure_eps_M``
    the full tree pressures at the subset points are within 512 eps·M of the CPU
    (M = sum of contribution magnitudes, from the CPU subset model);
``full_tree_jacobian_rel``
    the max |ΔJ| / column max over the subset rows, for every requested column.
    This is the evidence for the Jacobian budget proposal, against the
    provisional 1e-11.

Timing (warm, full set): pressures only, the full Jacobian request, and
``Problem.data_residuals`` / ``data_jacobian`` through the facade. The CPU
cost is estimated from the CPU subset per record and labelled as an estimate.

The parameter states are the C00b fixture vectors: parameters do not depend
on the point set. Templates are the C00b ``joint/templates.npz`` columns
(full mesh). Writes ``$COPS_OUTPUT_ROOT/cuda_c06/full_<UTC>/report.json``.
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


def magnitudes(P, M, x, fit_cops):
    """Per-point sums of |contribution| for NG and RG (the M of the budget), CPU model M."""
    a, b, T, kappa = P.unpack(x)
    act = M.Y @ a
    if P.nb:
        act = act + P.enhancement_scale * (P.enhancement_basis @ b)
    af = np.maximum(act, 0.0)[M.face]
    Tmax, T75, Tmin = T
    Tl = Tmax + fit_cops.KT * (Tmax - T75) * M.xlit
    u = 1.0 / np.sqrt(np.where(M.lit & (Tl > Tmin), Tl, Tmin))
    bc = lambda w: np.bincount(M.pt, weights=np.abs(w), minlength=M.npt)
    c, r, bm = M.cNG, M.cRG, M.bm
    mNG = (abs(c["H2O"] * bm["H2O"]) * bc(M.wn["H2O"] * af * u) +
           abs(kappa * c["CO2"] * bm["CO2"]) * bc(M.wn["CO2"] * u))
    mRG = (np.abs(r["H2O"]) * (bc(M.wp["H2O"] * af) + bm["H2O"] * bc(M.wv["H2O"] * af * u)) +
           np.abs(kappa * r["CO2"]) * (bc(M.Q["CO2"] * M.ACO2[M.face] * M.KPhi) +
                                       bm["CO2"] * bc(M.wv["CO2"] * u)))
    return mNG, mRG


def bitwise(a, b, request):
    return (a.pNG.tobytes() == b.pNG.tobytes() and a.pRG.tobytes() == b.pRG.tobytes() and
            all(u.tobytes() == v.tobytes() for key in request.keys()
                for u, v in zip(a.columns[key], b.columns[key])))


def restrict(result, index, request):
    """The rows ``index`` of an EvaluationResult, as plain arrays."""
    cols = {key: (result.columns[key][0][index], result.columns[key][1][index])
            for key in request.keys()}
    return result.pNG[index], result.pRG[index], cols


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True, help="C00b capture (joint fixtures, templates)")
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--configs", default="B29,N30,NC31,R30")
    parser.add_argument("--subset-count", type=int, default=16)
    parser.add_argument("--vectors", default="prior,legacy_strong_clip,archived_b29_state,"
                                            "enhancement_strong_clip,template_on_bound")
    parser.add_argument("--timing-repeats", type=int, default=3)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    import c00_capture_reference as c00
    from cops_cuda import backend as bk, evaluator, host_model, manifest, provenance, reference
    from cops_cuda.memory import DeviceBudget
    root = reference.bind()
    capture = Path(args.capture).resolve()
    fixtures = json.loads((capture / "joint" / "joint_fixtures.json").read_text())
    templates = np.load(str(capture / "joint" / "templates.npz"))
    modules = c00.load_reference_modules()
    fit_cops = modules["fit_cops"]
    parse_utc = modules["cops_selection"].parse_utc
    run = c00.base_directory(modules, c00.DEFAULT_BASE)
    header, names, rows, columns, profile_identity, _, _ = c00.fit_selection(
        modules, run, c00.DEFAULT_PROFILE)
    coef = fit_cops.read_activity_coefficients(str(root), header)
    elements = fit_cops.read_elements(str(run / "elements.bin"))

    # Full problem on the header-only model; verification subset on the CPU model.
    t0 = time.perf_counter()
    MH = host_model.header_model(run, header, columns, rows, elements, coef)
    MH.load_kernels(MH.sel)
    header_s = time.perf_counter() - t0
    subset = c00.choose_subset(rows, columns, args.subset_count, c00.LATE_INTERVAL, parse_utc)
    subset_rows = [rows[k] for k in subset]
    t0 = time.perf_counter()
    MS = fit_cops.Model(str(run), header, columns, subset_rows, elements, coef)
    MS.load_kernels(MS.sel)
    cpu_load_s = time.perf_counter() - t0
    full_ids = [int(p) for p in MH.ipoint[MH.rows]]
    sub_ids = [int(p) for p in MS.ipoint[MS.rows]]
    position = {p: k for k, p in enumerate(full_ids)}
    index = np.array([position[p] for p in sub_ids])       # subset rows within the full problem
    library = {n: templates[n].astype(float)
               for n in ("late_ng_not_rg", "late_rg_not_ng", "syn_cliff_standin")}
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c06" /
                                        ("full_" + provenance.utc_stamp()))
    log("report %s; full %d points / %d records (headers %.1f s); subset %d points / %d records"
        % (out, MH.npt, MH.indptr[-1], header_s, MS.npt, MS.indptr[-1]))

    wanted = [c for c in c00.JOINT_CONFIGS if c[0] in args.configs.split(",")]
    vector_names = args.vectors.split(",")
    budget = DeviceBudget(args.device)
    results, timing, all_ok = {}, {}, True
    try:
        for label, template_names, role in wanted:
            record = fixtures["configs"][label]
            arrays = np.load(str(capture / "joint" / ("joint_%s.npz" % label)))
            config = {"name": label, "degree": 4, "enhancement_degree": 4,
                      "fix_legacy_activity": False, "fix_temperature": False, "fix_kappa": False,
                      "ordered_temperature": False, "templates": bool(template_names)}
            tpl = (list(template_names),
                   np.column_stack([library[n] for n in template_names]) if template_names
                   else np.empty((len(elements["area"]), 0)), {})
            PF = c00.make_problem(fit_cops, MH, elements, coef, config, tpl)
            PS = c00.make_problem(fit_cops, MS, elements, coef, config, tpl)
            if list(PF.names) != record["names"] or list(PS.names) != record["names"]:
                raise SystemExit("parameter names differ from the C00b fixture for %s" % label)
            request = bk.ColumnRequest.from_problem(PF)
            ref = bk.ReferenceBackend.from_problem(PS)
            vectors = [v for v in vector_names if v in record["vectors"]]
            states = {v: arrays["%s__%s__x" % (label, v)] for v in vectors}
            per = {v: {} for v in vectors}
            cpu_t = []
            # 1. Subset: GPU sequential vs CPU (bitwise); keep CPU and subset-tree results.
            cpu_res, sub_tree = {}, {}
            with evaluator.from_problem(PS, run, device_id=args.device, mode="sequential",
                                        budget=budget) as seq:
                for v, x in states.items():
                    state = bk.ParameterState.from_problem(PS, x)
                    t = time.perf_counter()
                    cpu_res[v] = ref.evaluate(state, request)
                    cpu_t.append(time.perf_counter() - t)
                    per[v]["subset_sequential_bitwise"] = bitwise(seq.evaluate(state, request),
                                                                  cpu_res[v], request)
            with evaluator.from_problem(PS, run, device_id=args.device, mode="tree",
                                        budget=budget) as st:
                for v, x in states.items():
                    sub_tree[v] = st.evaluate(bk.ParameterState.from_problem(PS, x), request)
            # 2. Full set, tree mode.
            t = time.perf_counter()
            with evaluator.from_problem(PF, run, device_id=args.device, mode="tree",
                                        budget=budget) as full:
                bind_s = time.perf_counter() - t
                log("%s: bound %s in %.1f s" % (label, full.summary(), bind_s))
                for v, x in states.items():
                    got = full.evaluate(bk.ParameterState.from_problem(PF, x), request)
                    gNG, gRG, gcols = restrict(got, index, request)
                    s = sub_tree[v]
                    per[v]["point_set_invariance_bitwise"] = bool(
                        gNG.tobytes() == s.pNG.tobytes() and gRG.tobytes() == s.pRG.tobytes() and
                        all(gcols[k][0].tobytes() == s.columns[k][0].tobytes() and
                            gcols[k][1].tobytes() == s.columns[k][1].tobytes()
                            for k in request.keys()))
                    w = cpu_res[v]
                    mNG, mRG = magnitudes(PS, MS, x, fit_cops)
                    floor = 1e-300
                    eps_m = max(float(np.max(np.abs(gNG - w.pNG) / np.maximum(mNG, floor))),
                                float(np.max(np.abs(gRG - w.pRG) / np.maximum(mRG, floor)))) / EPS
                    jac = 0.0
                    for key in request.keys():
                        for gpu_col, cpu_col in zip(gcols[key], w.columns[key]):
                            scale = max(float(np.max(np.abs(cpu_col))), floor)
                            jac = max(jac, float(np.max(np.abs(gpu_col - cpu_col))) / scale)
                    per[v].update({"full_tree_pressure_eps_M": eps_m,
                                   "full_tree_pressure_within_512_eps_M": eps_m <= 512.0,
                                   "full_tree_jacobian_rel": jac,
                                   "full_tree_jacobian_within_provisional_1e-11": jac <= 1e-11,
                                   "all_finite": bool(np.isfinite(got.pNG).all() and
                                                      np.isfinite(got.pRG).all())})
                    ok = (per[v]["subset_sequential_bitwise"] and
                          per[v]["point_set_invariance_bitwise"] and eps_m <= 512.0 and
                          jac <= 1e-11 and per[v]["all_finite"])
                    all_ok = all_ok and ok
                # 3. Warm timing on the full set, and the unmodified Problem through the facade.
                x = states[vectors[0]]
                state = bk.ParameterState.from_problem(PF, x)
                full.evaluate(state, request)
                t_p, t_j = [], []
                for _ in range(args.timing_repeats):
                    t = time.perf_counter(); full.evaluate(state); t_p.append(time.perf_counter() - t)
                    t = time.perf_counter(); full.evaluate(state, request); t_j.append(time.perf_counter() - t)
                facade = bk.install(PF, full)
                PF.lam = 0.0
                t = time.perf_counter(); r = PF.data_residuals(x); t_r = time.perf_counter() - t
                t = time.perf_counter(); J = PF.data_jacobian(x); t_dj = time.perf_counter() - t
                PF.M = facade._model
                n_full, n_sub = int(MH.indptr[-1]), int(MS.indptr[-1])
                timing[label] = {
                    "bind_s": bind_s, "dataset_load_s": full.timing.get("dataset_load_s"),
                    "pressures_s_median": float(np.median(t_p)),
                    "full_jacobian_s_median": float(np.median(t_j)),
                    "problem_data_residuals_s": t_r, "problem_data_jacobian_s": t_dj,
                    "columns": len(request.keys()), "records": n_full,
                    "residuals": int(r.size), "jacobian_shape": list(J.shape),
                    "cpu_full_jacobian_s_estimate": float(np.median(cpu_t)) * n_full / max(n_sub, 1),
                    "cpu_estimate_note": "CPU subset time per record x full records; an estimate"}
            results[label] = {"role": role, "columns": len(request.keys()), "vectors": per}
            log("%s: %s; timing %s" % (label, {v: (d["subset_sequential_bitwise"],
                                                   d["point_set_invariance_bitwise"],
                                                   round(d["full_tree_pressure_eps_M"], 2),
                                                   "%.1e" % d["full_tree_jacobian_rel"])
                                               for v, d in per.items()},
                                       {k: (round(v, 3) if isinstance(v, float) else v)
                                        for k, v in timing[label].items() if k.endswith("_s_median")
                                        or k.endswith("_estimate")}))
    finally:
        budget.close()

    report = {"schema": "liouville_cops_cuda.c06_full_run", "schema_version": 1,
              "created_utc": provenance.utc_now(), "capture": str(capture), "run": str(run),
              "profile": profile_identity, "full_points": int(MH.npt),
              "full_records": int(MH.indptr[-1]), "subset_point_ids": sub_ids,
              "subset_records": int(MS.indptr[-1]), "cpu_subset_load_s": cpu_load_s,
              "header_model_s": header_s, "results": results, "timing": timing,
              "all_pass": bool(all_ok),
              "manifest": manifest.backend_manifest("cuda_fp64_v1", args.device)}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; report %s" % ("PASS" if all_ok else "FAIL", out / "report.json"))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
