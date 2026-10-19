#!/usr/bin/env python3
"""C06 evidence: the CUDA evaluator against the C00b joint fixtures (CUDA env).

It rebuilds, in the CUDA environment, exactly the problems of the C00b
``joint`` stage: the same 47-point subset of the 788k base, the same templates
(``joint/templates.npz``) and the configurations B29, N30, C30, NC31 and R30.
Then, for every fixture state:

``sequential``
    The GPU verification mode against the CPU reference evaluated in this
    interpreter. The pressures and all requested Jacobian columns must be
    **bitwise** identical.
``tree``
    The GPU production mode against the C00b fixture, computed by the CPU
    reference interpreter:
    * pressures: ≤ 512 eps·M (M is the fixture's contribution magnitude);
    * data residuals: the derived bound w·512 eps·M/p;
    * data Jacobian: the provisional ≤ 1e-11 × column max, computed through
      the unmodified ``fit_cops.Problem`` via the C02 facade.
``fixture_reproducibility``
    Whether the CPU reference in this interpreter reproduces the fixture
    bitwise. This is informational: two interpreters with the same NumPy are
    expected to agree.
``timing``
    Warm tree-mode times for pressures only and for the full Jacobian request,
    on this subset.

Writes ``$COPS_OUTPUT_ROOT/cuda_c06/check_<UTC>/report.json``.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--capture", required=True)
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--timing-repeats", type=int, default=5)
    parser.add_argument("--output-root", default=os.environ.get("COPS_OUTPUT_ROOT", DEFAULT_OUTPUT_ROOT))
    args = parser.parse_args()

    import c00_capture_reference as c00
    from cops_cuda import backend as bk, evaluator, manifest, provenance, reference
    reference.bind()
    capture = Path(args.capture).resolve()
    fixtures = json.loads((capture / "joint" / "joint_fixtures.json").read_text())
    templates = np.load(str(capture / "joint" / "templates.npz"))
    modules = c00.load_reference_modules()
    fit_cops = modules["fit_cops"]
    parse_utc = modules["cops_selection"].parse_utc
    root = reference.bind()
    run = c00.base_directory(modules, c00.DEFAULT_BASE)
    # The C00b stage used the default profile; the subset identity check below
    # proves the same rows were selected.
    header, names, rows, columns, _, _, _ = c00.fit_selection(modules, run, c00.DEFAULT_PROFILE)
    subset = c00.choose_subset(rows, columns, 40, c00.LATE_INTERVAL, parse_utc)
    subset_rows = [rows[k] for k in subset]
    point_ids = [int(r[columns["iPoint"]]) for r in subset_rows]
    if point_ids != fixtures["subset_point_ids"]:
        raise SystemExit("subset differs from the C00b fixture")
    coef = fit_cops.read_activity_coefficients(str(root), header)
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    M = fit_cops.Model(str(run), header, columns, subset_rows, elements, coef)
    M.load_kernels(M.sel)
    library = {"late_ng_not_rg": templates["late_ng_not_rg"].astype(float),
               "late_rg_not_ng": templates["late_rg_not_ng"].astype(float),
               "syn_cliff_standin": templates["syn_cliff_standin"].astype(float)}
    out = provenance.make_new_directory(Path(args.output_root) / "cuda_c06" /
                                        ("check_" + provenance.utc_stamp()))
    log("report directory %s; %d points, %d contributions" % (out, len(point_ids), M.indptr[-1]))

    results, all_ok, timing = {}, True, {}
    for label, template_names, role in c00.JOINT_CONFIGS:
        record = fixtures["configs"][label]
        config = {"name": label, "degree": 4, "enhancement_degree": 4,
                  "fix_legacy_activity": False, "fix_temperature": False, "fix_kappa": False,
                  "ordered_temperature": False, "templates": bool(template_names)}
        tpl = (list(template_names),
               np.column_stack([library[n] for n in template_names]) if template_names
               else np.empty((len(elements["area"]), 0)), {})
        P = c00.make_problem(fit_cops, M, elements, coef, config, tpl)
        if list(P.names) != record["names"]:
            raise SystemExit("parameter names differ from the C00b fixture for %s" % label)
        arrays = np.load(str(capture / "joint" / ("joint_%s.npz" % label)))
        request = bk.ColumnRequest.from_problem(P)
        ref = bk.ReferenceBackend.from_problem(P)
        per_vector = {}
        with evaluator.from_problem(P, run, device_id=args.device, mode="sequential") as seq:
            for vlabel in record["vectors"]:
                x = arrays["%s__%s__x" % (label, vlabel)]
                state = bk.ParameterState.from_problem(P, x)
                want, got = ref.evaluate(state, request), seq.evaluate(state, request)
                bitwise = (got.pNG.tobytes() == want.pNG.tobytes() and
                           got.pRG.tobytes() == want.pRG.tobytes() and
                           all(g.tobytes() == w.tobytes() for key in request.keys()
                               for g, w in zip(got.columns[key], want.columns[key])))
                fixture_bitwise = (want.pNG.tobytes() == arrays["%s__%s__png" % (label, vlabel)].tobytes()
                                   and want.pRG.tobytes() == arrays["%s__%s__prg" % (label, vlabel)].tobytes())
                per_vector[vlabel] = {"sequential_bitwise_vs_cpu": bool(bitwise),
                                      "cpu_here_bitwise_vs_fixture": bool(fixture_bitwise)}
        with evaluator.from_problem(P, run, device_id=args.device, mode="tree") as tree:
            facade = bk.install(P, tree)
            for vlabel in record["vectors"]:
                key = "%s__%s" % (label, vlabel)
                x = arrays[key + "__x"]
                state = bk.ParameterState.from_problem(P, x)
                got = tree.evaluate(state)
                mNG, mRG = arrays[key + "__png_magnitude"], arrays[key + "__prg_magnitude"]
                pNG_ok = bool(np.all(np.abs(got.pNG - arrays[key + "__png"]) <= 512 * EPS * mNG))
                pRG_ok = bool(np.all(np.abs(got.pRG - arrays[key + "__prg"]) <= 512 * EPS * mRG))
                P.lam = 0.0
                r_gpu, J_gpu = P.data_residuals(x), P.data_jacobian(x)
                r_ref, J_ref = arrays[key + "__data_residuals"], arrays[key + "__data_jacobian"]
                floor = 1e-300
                delta = np.concatenate([
                    P.wNG * 512 * EPS * mNG[P.useNG] / np.maximum(arrays[key + "__png"][P.useNG], floor),
                    P.wRG * 512 * EPS * mRG[P.useRG] / np.maximum(arrays[key + "__prg"][P.useRG], floor)])
                r_ok = bool(np.all(np.abs(r_gpu - r_ref) <= delta))
                scale = np.maximum(np.max(np.abs(J_ref), axis=0), 1e-300)
                jac_rel = float(np.max(np.abs(J_gpu - J_ref) / scale))
                cost = 0.5 * float(np.sum(r_ref ** 2))
                cost_bound = float(np.sum(np.abs(r_ref) * delta) + 0.5 * np.sum(delta ** 2))
                cost_ok = abs(0.5 * float(np.sum(r_gpu ** 2)) - cost) <= cost_bound
                pmax = float(max(np.max(np.abs(got.pNG - arrays[key + "__png"]) / np.maximum(mNG, floor)),
                                 np.max(np.abs(got.pRG - arrays[key + "__prg"]) / np.maximum(mRG, floor))) / EPS)
                per_vector[vlabel].update({
                    "tree_pressure_within_512_eps_M": pNG_ok and pRG_ok,
                    "tree_pressure_max_eps_M": pmax,
                    "tree_residuals_within_derived_bound": r_ok,
                    "tree_cost_within_derived_bound": bool(cost_ok),
                    "tree_jacobian_rel_to_column_max": jac_rel,
                    "tree_jacobian_within_provisional_1e-11": jac_rel <= 1e-11})
                vector_ok = (per_vector[vlabel]["sequential_bitwise_vs_cpu"] and pNG_ok and pRG_ok
                             and r_ok and cost_ok and jac_rel <= 1e-11)
                all_ok = all_ok and vector_ok
            # Timing (warm): pressures only, and the full request.
            x = arrays["%s__prior__x" % label]
            state = bk.ParameterState.from_problem(P, x)
            tree.evaluate(state, request)
            t_p, t_j = [], []
            for _ in range(args.timing_repeats):
                t0 = time.perf_counter(); tree.evaluate(state); t_p.append(time.perf_counter() - t0)
                t0 = time.perf_counter(); tree.evaluate(state, request); t_j.append(time.perf_counter() - t0)
            timing[label] = {"pressures_s_median": float(np.median(t_p)),
                             "full_jacobian_s_median": float(np.median(t_j)),
                             "columns": len(request.keys()),
                             "records": int(M.indptr[-1])}
            P.M = facade._model
        results[label] = {"role": role, "columns": len(request.keys()), "vectors": per_vector}
        log("%s: %s; timing %s" % (label, {v: (d["sequential_bitwise_vs_cpu"], round(d["tree_pressure_max_eps_M"], 2),
                                              "%.1e" % d["tree_jacobian_rel_to_column_max"])
                                          for v, d in per_vector.items()}, timing[label]))

    report = {"schema": "liouville_cops_cuda.c06_check", "schema_version": 1,
              "created_utc": provenance.utc_now(), "capture": str(capture),
              "subset_points": len(point_ids), "contributions": int(M.indptr[-1]),
              "results": results, "timing": timing, "all_pass": bool(all_ok),
              "manifest": manifest.backend_manifest("cuda_fp64_v1", args.device)}
    provenance.write_json_exclusive(out / "report.json", report)
    log("ALL %s; report %s" % ("PASS" if all_ok else "FAIL", out / "report.json"))
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
