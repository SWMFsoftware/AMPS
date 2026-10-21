#!/usr/bin/env python3
"""Late-interval activity fine tuning of a completed fit, kernel responses on the GPU (C13).

    .venv/bin/python tools/fine_tune_gpu.py BASE --restart-from FIT/fit_restart.npz \\
        --minutes 5 --intervals 5 [--method cpu|multiplier] [--third-group all|mixed] \\
        [--dominance-ratio 1e-3] [--lambda 0] [--max-nfev 1000] [--prepare-only] \\
        [--devices 0,1] [--headroom-fraction F] [--precision mixed] --out NEW_DIR

A second fitting stage that changes only the activity of facets first seen in the last
``--minutes``. The window, intervals, nearest-sample weights, first-interval NG/RG/mixed
dominance and the log-residual objective (each gauge's pressure on its own, equal total weight
per gauge; no NG/RG ratio term) are those of ``../LiouvilleCOPS/FINE_TUNING_STEP.md``.

``--method cpu`` (default, user decision 2026-10-21): the CPU branch's method exactly
(``tools/fine_tune_activity.py``). Strict support (no Kn, KPhi, KV != 0 at any earlier row), so
every earlier pressure is bitwise unchanged; additive correction A_final = A0 + A_scale (B c),
c >= 0, on the (interval, NG/RG/third) groups of ``--third-group``; the CPU's own
``fine_tuning.group_basis`` and ``fit_corrections`` (bounded TRF, lambda Gram penalty). Only the
per-atom kernel responses are computed here, on the GPU (``restart.StateEvaluator``
block sums; equal to the CPU's ``pressure_matrices`` to rounding, ``tests/test_fine_tuning_gpu.py``).
Products as the CPU tool's.

``--method multiplier`` (variant): A_final = f_g A0 + c_g A_scale on disjoint groups, f_g >= 0,
c_g >= 0; ``--confine total --epsilon E`` keeps the summed change before the window within E
(support by the ``restart_tool enhance`` prefix rule, constrained fit); ``--gauges joint`` (one
fit to NG and RG) or ``independent`` (NG only, RG only). See ``cops_cuda/fine_tuning.py``.

Checks before anything is written: the restart's mesh, base elements/table/input hashes, row
order and prior coefficients against the base; static kernels and analytical temperatures; the
loaded kernel set against ``run_manifest.json``; the restart's stored predictions reproduced
(relative 1e-6); exact zero response of every group before its interval.

Outputs in ``--out`` (must not exist):
* ``fine_tuning.json``: identities, settings, groups, parameters, optimizer, diagnostics;
* ``fine_tuning_prepared.npz``: first_interval, atom, scores, bins, weights, edges, baseline,
  WA and W1 (the group responses);
* ``fine_tuning_state.npz``: fitted f and c, frozen and final activity, predictions;
* ``pressure_comparison.txt``: per row before/after pressures, references, interval, flags;
* ``figures/pressures_{time,altitude}`` (.eps, .png; the fine-tuning window shaded) and
  ``figures/groups_{ng,rg,mixed}`` lon/lat maps of the groups (``--no-maps`` skips them);
* ``ActivityFit.txt``, ``ActivityEnhancementFit.txt``, ``cops_fine_tuning.in``: the final field
  for an independent C++ replay (``./liouville_cops <out>/cops_fine_tuning.in``).
"""

import argparse
import datetime as dt
import json
import os
from pathlib import Path
import re
import sys
import time

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def save_arrays(path, arrays, meta):
    from cops_cuda import provenance
    meta = dict(meta, array_sha256={k: provenance.array_sha256(np.ascontiguousarray(v)) for k, v in arrays.items()})
    with Path(path).open("xb") as stream:
        np.savez_compressed(stream, meta_json=np.array(json.dumps(meta, sort_keys=True)), **arrays)


def export_cpp(out, run, elements, r, final, body_frame, lam, header_value,
               regularization="late_interval_group_multiplier_floor_change_relative_frozen_power_v1"):
    """The final field in the strict per-facet format of ../LiouvilleCOPS (as its fine-tuning
    export): legacy coefficients plus the explicit residual final - Y a, one clip in C++."""
    from cops_cuda import provenance, reference
    fit_cops = reference.load("fit_cops")
    coefficient = out / "ActivityFit.txt"
    coef = np.array(r.coef_prior, dtype=float).copy()
    coef[0] = r.a
    coef[1] *= float(np.asarray(r.kappa).reshape(-1)[0])
    np.savetxt(coefficient, coef, fmt="%.17e")
    residual = final - elements["Y"] @ r.a
    enhancement = out / "ActivityEnhancementFit.txt"
    with enhancement.open("x") as stream:
        stream.write("# Frozen fitted state with late-interval group multipliers and floors; explicit facet field\n")
        fields = dict(Format=fit_cops.FACET_ENHANCEMENT_FORMAT, FormatVersion=1,
                      BasisID=fit_cops.COMPOSITE_BASIS_ID, Species="H2O", BodyFrame=body_frame,
                      Ordering="increasing_face_index", Clipping="max(0,legacy_linear+enhancement)",
                      HarmonicMinimumDegree=5, HarmonicMaximumDegree=4, Scale=1,
                      ScaleInterpretation="explicit_dimensionless_residual_per_facet",
                      LegacyActivityIdentity="sha256:" + provenance.file_identity(coefficient)["sha256"],
                      Regularization=regularization,
                      Lambda="%.17g" % lam, FacetCount=len(residual))
        for key, value in fields.items():
            stream.write("%s = %s\n" % (key, value))
        for i, (xyz, area, delta) in enumerate(zip(elements["centroid"], elements["area"], residual)):
            stream.write("Facet = %d %.17e %.17e %.17e %.17e %.17e\n" % (i, *xyz, area, delta))
    text = (run / "input.snapshot.in").read_text()
    T = np.asarray(r.T, float)
    replacement = dict(ActivityModel="file", ActivityFile=str(coefficient),
                       ActivityFileHeliocentricModel=header_value("HeliocentricModel"),
                       AnalyticalTmax="%.17g" % T[0], AnalyticalT75="%.17g" % T[1], AnalyticalTmin="%.17g" % T[2],
                       ActivityEnhancementModel=fit_cops.FACET_ENHANCEMENT_MODEL,
                       ActivityEnhancementDegree="4", ActivityEnhancementScale="1",
                       ActivityEnhancementFile=str(enhancement), OutputDirectory=str(out) + "_verify")
    for key, value in replacement.items():
        pattern = r"^\s*" + re.escape(key) + r"\s*=.*$"
        line = key + " = " + value
        text = (re.sub(pattern, lambda _: line, text, flags=re.M) if re.search(pattern, text, re.M)
                else text + "\n" + line + "\n")
    (out / "cops_fine_tuning.in").write_text(text)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("run", help="the original unenhanced fitting base of the fit")
    ap.add_argument("--restart-from", required=True, help="the completed fit's fit_restart.npz (frozen)")
    ap.add_argument("--minutes", type=float, default=5.0, help="window length m [min] (default 5)")
    ap.add_argument("--intervals", type=int, required=True, help="number n of equal intervals")
    ap.add_argument("--end-utc", help="window end; default the last modeled row (recorded)")
    ap.add_argument("--dominance-ratio", type=float, default=1e-3)
    ap.add_argument("--method", choices=("cpu", "multiplier"), default="cpu",
                    help="cpu (default, user decision 2026-10-21): the CPU branch's method exactly "
                         "(../LiouvilleCOPS tools/fine_tune_activity.py): strict support, additive c >= 0, "
                         "--third-group groups, its own group_basis and fit_corrections; "
                         "multiplier: f*A0 + c*A_scale on disjoint groups with --confine/--epsilon/--gauges")
    ap.add_argument("--third-group", choices=("all", "mixed"), default="all",
                    help="--method cpu: the CPU's third column per interval (default all)")
    ap.add_argument("--prepare-only", action="store_true",
                    help="--method cpu: write the support and responses (fine_tuning_prepared.npz) without fitting")
    ap.add_argument("--confine", choices=("total", "strict"), default="total",
                    help="total (default, user decision 2026-10-20): earlier pressures change by at most "
                         "--epsilon in total; strict: the CPU definition (no earlier structural response, "
                         "earlier pressures bitwise unchanged)")
    ap.add_argument("--epsilon", type=float, default=1e-3,
                    help="relative limit for --confine total (default 1e-3, user decision 2026-10-20)")
    ap.add_argument("--gauges", choices=("joint", "independent"), default="joint",
                    help="joint (default, user decision 2026-10-20): one fit to the NG and the RG pressures, "
                         "each pressure on its own (no NG/RG ratio term), written to fit_joint/; "
                         "independent: one fit to NG only and one to RG only (fit_ng/, fit_rg/)")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.0)
    ap.add_argument("--max-nfev", type=int, default=1000)
    ap.add_argument("--devices", default="0", help="comma-separated GPUs; rows are split over them")
    ap.add_argument("--headroom-fraction", type=float,
                    help="device memory headroom of the budget (default 0.20; 3M on two A6000: 0.15)")
    ap.add_argument("--precision", choices=("fp64", "mixed"), default="fp64",
                    help="fp64 (default) or mixed (C14 experiment): FP32 kernel records, FP64 arithmetic; "
                         "needed for SHAP8 12M (with --devices 0,1 --headroom-fraction 0.06)")
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--no-maps", action="store_true", help="skip the group lon/lat maps")
    ap.add_argument("--out", required=True, help="new output directory")
    args = ap.parse_args(argv)
    from cops_cuda import memory
    if args.headroom_fraction is not None:
        lo, hi = memory.HEADROOM_RANGE
        if not lo <= args.headroom_fraction <= hi:
            ap.error("--headroom-fraction must be in [%g, %g]" % (lo, hi))
        os.environ[memory.HEADROOM_ENVIRONMENT] = repr(args.headroom_fraction)
    from cops_cuda import precision as cuda_precision
    os.environ[cuda_precision.RECORD_PRECISION_ENVIRONMENT] = "fp32" if args.precision == "mixed" else "fp64"
    from cops_cuda import fine_tuning as ftg
    from cops_cuda import provenance, reference
    from cops_cuda import restart as rs
    fit_cops = reference.load("fit_cops")

    run = Path(args.run).resolve()
    out = Path(args.out).resolve()
    if out.exists():
        ap.error("output already exists: %s" % out)
    if not np.isfinite(args.lam) or args.lam < 0 or args.max_nfev < 1:
        ap.error("lambda must be finite/nonnegative and max-nfev positive")
    if not 0 < args.dominance_ratio < 1:
        ap.error("dominance ratio must be between zero and one")
    devices = [int(d) for d in args.devices.split(",")]
    r = rs.read(args.restart_from)
    t_start = time.perf_counter()

    # ---- input checks (the CPU tool's, against the restart's recorded base)
    for key, name in (("base_elements", "elements.bin"), ("base_model_table", "COPS_model.dat"),
                      ("base_input", "input.snapshot.in")):
        expected = r.meta.get("base", {}).get(key)
        if not expected or provenance.file_identity(run / name)["sha256"] != expected.get("sha256"):
            ap.error("restart/base identity mismatch or missing provenance: %s" % key)
    header, names, rows = fit_cops.read_model_file(str(run / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    hv = lambda key: fit_cops.header_value(header, key)
    if hv("FlightTimeRotation") != "none" or hv("TemperatureModel") != "analytical":
        ap.error("fine tuning requires static kernels and analytical temperatures")
    if fit_cops.optional_header_value(header, "ActivityEnhancementModel") not in (None, "none"):
        ap.error("give the original unenhanced fitting base, not a verification run")
    if not np.array_equal(r.point_ids, [int(x[cols["iPoint"]]) for x in rows]) or \
            not np.array_equal(r.utc, [x[cols["UTC"]] for x in rows]):
        ap.error("restart observation rows differ from the base table")
    coef = fit_cops.read_activity_coefficients(str(reference.bind()), header)
    if not np.array_equal(coef, r.coef_prior):
        ap.error("restart prior coefficients differ from the base")
    manifest = json.loads((run / "run_manifest.json").read_text())
    body_frames = manifest.get("effective_settings", {}).get("BodyFrame", [])
    if len(body_frames) != 1:
        ap.error("base must declare one BodyFrame")
    bins, weights, edges = ftg.cpu_intervals(r.utc, args.minutes, args.intervals, args.end_utc)
    frozen = np.maximum(np.asarray(r.activity_linear, float), 0.0)
    area = np.asarray(r.area, float)
    power = float(np.dot(area, frozen * frozen))
    if not np.isfinite(power) or power <= 0:
        ap.error("frozen activity has zero/nonfinite area-weighted power")
    scale = float(r.enhancement_scale)
    inside = np.array([int(x[cols["inside"]]) for x in rows], dtype=bool)
    used = np.column_stack([r.use_ng, r.use_rg]).astype(bool) & ~inside[:, None]
    refs = np.column_stack([r.ref_ng, r.ref_rg]).astype(float)
    late = bins >= 0

    # ---- GPU preparation
    ev = ftg.SplitEvaluator(run, devices)
    timing = {}
    try:
        if ev.mesh_identity != r.mesh_identity():
            raise SystemExit("restart mesh differs from the base mesh")
        t0 = time.perf_counter()
        png, prg = ev.pressures(r)
        timing["bind_and_pressures_s"] = time.perf_counter() - t0
        kernels = ev.kernel_identity()
        if not kernels["verified_against_run_manifest"]:
            raise SystemExit("the base's run_manifest.json records no kernel set digest")
        baseline = np.column_stack([png, prg])
        stored = np.column_stack([r.pred_ng, r.pred_rg])
        replay = float(np.max(np.abs(baseline - stored) / np.maximum(np.abs(stored), 1e-300)))
        if not np.isfinite(stored).all() or not np.allclose(baseline, stored, rtol=1e-6, atol=1e-30):
            raise SystemExit("the restart's stored pressures do not reproduce on this base (max rel %.3g)" % replay)
        t0 = time.perf_counter()
        if args.method == "cpu":
            args.confine = "strict"          # the CPU definition of the support
        confinement = {"mode": args.confine}
        if args.confine == "strict":
            first = ftg.first_intervals(ev.first_structural_row(), bins)
        else:
            delta_ref = ftg.reference_change(frozen, scale)
            tau, eligible, confinement = ftg.confine_total(ev, r, delta_ref, baseline, bins, args.epsilon)
            confinement["mode"] = "total"
            first = np.full(len(frozen), ftg.BEFORE, dtype=np.int32)
            if eligible.any():
                sig = ev.first_significant_row(r, delta_ref, 0.0, row_mask=late)   # first window row seen
                if np.any(sig[eligible] < 0) or np.any(bins[sig[eligible]] < 0):
                    raise SystemExit("internal: an eligible facet is not first significant in the window")
                first[eligible] = bins[sig[eligible]]
            print("confine total: epsilon %g, tau %.4g, %d of %d window-seen facets eligible "
                  "(reference change outside %.4g)" % (args.epsilon, confinement["tau"], confinement["eligible"],
                  confinement["candidates_seen_in_window"], confinement["reference_outside_max_rel"]), flush=True)
        ng_score, rg_score = ev.first_interval_scores(r, first, bins, weights)
        atom = ftg.classify(first, ng_score, rg_score, args.dominance_ratio)
        near = ftg.threshold_near(first, ng_score, rg_score, args.dominance_ratio)
        groups = ftg.group_names(atom, args.intervals)
        timing["support_s"] = time.perf_counter() - t0
        counts = {"first_seen_before_window": int((first == -1).sum()), "never_seen": int((first == -2).sum()),
                  "first_seen_per_interval": [int((first == j).sum()) for j in range(args.intervals)],
                  "grouped": int((atom >= 0).sum())}
        print("facets: %s" % json.dumps(counts), flush=True)
        if args.method == "cpu":
            cpu = cpu_method_fit(ev, r, args, atom, groups, frozen, area, power, scale, baseline, refs, used,
                                 bins, late)
            timing["responses_and_fit_s"] = time.perf_counter() - t0
        else:
            if not groups:
                raise SystemExit("no facet is first seen inside the window with a nonzero response: nothing to tune")
            t0 = time.perf_counter()
            WA, W1 = ftg.response_matrices(ev, r, atom, groups, frozen)
            timing["responses_s"] = time.perf_counter() - t0
            if args.confine == "strict":
                for k, (code, name) in enumerate(groups):
                    if np.any(WA[bins < code // 3, :, k] != 0) or np.any(W1[bins < code // 3, :, k] != 0):
                        raise SystemExit("nonzero response of group %s before its interval" % name)
            if np.any(WA < 0) or np.any(W1 < 0):
                raise SystemExit("negative group response: review the base pressure formulation")
            pnames, D, xref, kinds, index = ftg.parameters(groups, WA, W1, scale)
            G = ftg.gram_matrix(groups, atom, area, frozen, scale, power, kinds, index)
            solutions = []
            for label, gauges in MODES[args.gauges]:
                t0 = time.perf_counter()
                x, prediction, optimizer = ftg.fit(
                    baseline, D, xref, refs, used, late, G, args.lam, args.max_nfev, gauges=gauges,
                    confine=None if args.confine == "strict" else (~late, args.epsilon))
                fit_s = time.perf_counter() - t0
                final = ftg.final_activity(frozen, atom, groups, x, kinds, index, scale)
                # independent GPU check of the fitted field: p0 + (exact change of the final field)
                dn, dr = ev.point_changes(r, final - frozen)
                check = np.column_stack([baseline[:, 0] + dn, baseline[:, 1] + dr])
                check_rel = float(np.max(np.abs(check - prediction) / np.abs(prediction)))
                outside_rel = float(np.max(np.abs(check[~late] - baseline[~late]) / np.abs(baseline[~late])))
                if args.confine == "strict":
                    if not np.array_equal(prediction[~late], baseline[~late]) or np.any(dn[~late] != 0) or \
                            np.any(dr[~late] != 0):
                        raise SystemExit("earlier predictions changed")
                elif outside_rel > args.epsilon * (1 + 1e-6):
                    raise SystemExit("earlier predictions changed by %.4g > epsilon %g" % (outside_rel, args.epsilon))
                solutions.append(dict(label=label, gauges=gauges, x=x, prediction=prediction, optimizer=optimizer,
                                      final=final, check_rel=check_rel, outside_rel=outside_rel, fit_s=fit_s))
    finally:
        ev.release()
    if args.method == "cpu":
        return cpu_method_outputs(cpu, out, run, ev, r, args, first, atom, ng_score, rg_score, near, bins,
                                  weights, edges, baseline, refs, used, frozen, area, scale, kernels, replay,
                                  devices, timing, t_start, body_frames[0], hv)

    common = {
        "schema": ftg.SCHEMA, "created_utc": now(),
        "correction_law": "A_final = f_g*A0 + c_g*A_scale on group g (f_g >= 0 multiplier, c_g >= 0 floor); "
                          "A0 = max(activity_linear, 0); facets outside the groups unchanged",
        "objective": "log(p_model/p_data) of the selected late points of the fitted gauge(s), each pressure on "
                     "its own (never the NG/RG ratio); gauges=independent: one fit to NG, one to RG",
        "support_definition": ("--confine total: window-seen facets kept by the enhance prefix rule "
                               "(summed reference change max(A0, A_scale) within epsilon before the window); "
                               "first window interval with a nonzero response; first-interval dominance"
                               if args.confine == "total" else
                               "../LiouvilleCOPS/FINE_TUNING_STEP.md sections 2-4 (first structural response)"),
        "base_run": str(run), "source_restart": provenance.file_identity(args.restart_from),
        "mesh_identity": r.mesh_identity(), "kernels": kernels, "devices": devices,
        "row_blocks": [list(map(int, b)) for b in ev.blocks],
        "memory_headroom_fraction": memory.headroom_fraction_default(),
        "record_precision": "fp32" if args.precision == "mixed" else "fp64",
        "settings": {"minutes": args.minutes, "intervals": args.intervals, "dominance_ratio": args.dominance_ratio,
                     "confine": args.confine, "epsilon": args.epsilon if args.confine == "total" else None,
                     "gauges": args.gauges, "lambda": args.lam, "max_nfev": args.max_nfev},
        "end_utc": args.end_utc or str(r.utc[-1]), "end_source": "explicit" if args.end_utc else "last modeled row",
        "edges_epoch_s": edges.tolist(),
        "window_rows": int(late.sum()), "selected_late_points": {"NG": int((late & used[:, 0]).sum()),
                                                                  "RG": int((late & used[:, 1]).sum())},
        "frozen": ["H2O activity field (all fitted harmonics and templates)", "Tmax", "T75", "Tmin", "kappa",
                   "all base physics"],
        "enhancement_scale": scale, "frozen_replay_max_relative_difference": replay,
        "confinement": confinement,
        "facets": {"not_eligible": int((first == -1).sum()), "never_seen": int((first == -2).sum()),
                   "per_interval": [int((first == j).sum()) for j in range(args.intervals)],
                   "in_window_without_response": int(((first >= 0) & (atom < 0)).sum()),
                   "grouped": int((atom >= 0).sum()), "threshold_near_rtol_1e-9": near.tolist()[:100],
                   "threshold_near_count": int(len(near))},
        "parameter_names": pnames, "timing_s": timing, "cpp_verification": "not run", "scientific_review": "pending"}

    out.mkdir(parents=True)
    save_arrays(out / "fine_tuning_prepared.npz",
                dict(first_interval=first, atom=atom, ng_score=ng_score, rg_score=rg_score, bins=bins,
                     weights=weights, edges=edges, baseline=baseline, WA=WA, W1=W1),
                {"schema": ftg.SCHEMA, "group_names": [n for _, n in groups], "group_codes": [c for c, _ in groups]})
    summary = dict(common, solutions={})
    for sol in solutions:
        sub = out / ("fit_" + sol["label"])
        sub.mkdir()
        report = solution_report(sol, common, groups, atom, area, frozen, ng_score, rg_score, baseline,
                                 refs, used, late, power)
        save_arrays(sub / "fine_tuning_state.npz",
                    dict(parameters=sol["x"], parameter_reference=xref, design=D, gram=G, activity_frozen=frozen,
                         activity_final=sol["final"], prediction=sol["prediction"], refs=refs, used=used,
                         point_ids=np.asarray(r.point_ids)),
                    {"schema": ftg.SCHEMA, "parameter_names": pnames, "fitted_gauges": report["fitted_gauges"],
                     "source_restart": common["source_restart"]})
        with (sub / "pressure_comparison.txt").open("x") as stream:
            stream.write("# iPoint UTC interval refNG refRG beforeNG beforeRG afterNG afterRG useNG useRG saturatedNG\n")
            for k, point in enumerate(r.point_ids):
                stream.write("%d %s %d %s %d %d %d\n" % (point, r.utc[k], bins[k], " ".join(
                    "%.17e" % v for v in np.r_[refs[k], baseline[k], sol["prediction"][k]]), used[k, 0], used[k, 1],
                    int(r.saturated[k])))
        export_cpp(sub, run, ev.elements, r, sol["final"], body_frames[0], args.lam, hv)
        if not args.no_figures:
            report["figures"] = figures_for(sub, run, r, baseline, sol["prediction"], bins, atom, groups,
                                            devices[0], False, args, "fitted to %s" % report["fitted_gauges"])
        report["outputs"] = {p.name: provenance.file_identity(p) for p in sorted(sub.iterdir()) if p.is_file()}
        (sub / "fine_tuning.json").write_text(json.dumps(report, indent=1))
        summary["solutions"][sol["label"]] = {k: report[k] for k in (
            "fitted_gauges", "late_log_rms", "earlier_pressures_max_relative_change", "groups", "optimizer",
            "identifiability_warning", "activity_change")}
        rms = report["late_log_rms"]
        print("[%s] groups: %s" % (sol["label"], ", ".join("%s (%d facets): f=%.4g c=%.4g" % (
            g["name"], g["face_count"], g["f"], g["c"]) for g in report["groups"])))
        print("[%s] late log rms NG %.4f -> %.4f, RG %.4f -> %.4f; earlier change %.3g; optimizer %s" % (
            sol["label"], rms["NG"]["before"], rms["NG"]["after"], rms["RG"]["before"], rms["RG"]["after"],
            sol["outside_rel"], sol["optimizer"]["message"]))
        if report["identifiability_warning"]:
            print("[%s] WARNING: %s" % (sol["label"], report["identifiability_warning"]))
    if not args.no_figures and not args.no_maps:
        summary["maps"] = figures_for(out, run, r, None, None, bins, atom, groups, devices[0], True, args, None)
    summary["timing_s"]["total_s"] = time.perf_counter() - t_start
    (out / "fine_tuning.json").write_text(json.dumps(summary, indent=1))
    print("wrote %s" % out)
    return 0


MODES = {"independent": (("ng", (0,)), ("rg", (1,))), "joint": (("joint", (0, 1)),)}
CPU_SCHEMA = "liouville_cops.late_interval_fine_tuning.v1"
CPU_REGULARIZATION = "late_interval_total_added_activity_relative_frozen_power_v1"


def cpu_method_fit(ev, r, args, atom, groups, frozen, area, power, scale, baseline, refs, used, bins, late):
    """--method cpu: ../LiouvilleCOPS tools/fine_tune_activity.py on GPU-prepared responses.

    The atom responses are the CPU's ``pressure_matrices`` (sum of the per-unit-activity pressure
    response over each atom's facets, all rows), computed exactly with the evaluator; the group
    basis, design, Gram matrix and the fit are the CPU's own functions."""
    from cops_cuda import reference
    ft = reference.load("fine_tuning")
    n = args.intervals
    atom_response = np.zeros((len(baseline), 2, 3 * n))
    for code, _ in groups:
        atom_response[:, 0, code], atom_response[:, 1, code] = ev.point_changes(r, (atom == code).astype(float))
    for j in range(n):
        if np.any(atom_response[bins < j, :, 3 * j:3 * j + 3] != 0):
            raise SystemExit("nonzero early response in interval group %d" % j)
    basis, labels, omitted, sizes, areas = ft.group_basis(atom, area, n, args.third_group)
    design = scale * (atom_response @ basis)
    gram = scale ** 2 * (basis.T @ (areas[:, None] * basis)) / power
    result = dict(atom_response=atom_response, basis=basis, labels=labels, omitted=omitted, sizes=sizes,
                  areas=areas, design=design, gram=gram)
    if args.prepare_only:
        return result
    coefficients, prediction, optimizer = ft.fit_corrections(baseline, design, refs, used, late, gram,
                                                             args.lam, args.max_nfev)
    change = np.zeros(len(frozen))
    valid = atom >= 0
    change[valid] = scale * (basis @ coefficients)[atom[valid]]
    if not np.array_equal(prediction[~late], baseline[~late]):
        raise SystemExit("earlier predictions changed")
    dn, dr = ev.point_changes(r, change)              # independent GPU check of the final field
    check = np.column_stack([baseline[:, 0] + dn, baseline[:, 1] + dr])
    if np.any(dn[~late] != 0) or np.any(dr[~late] != 0):
        raise SystemExit("earlier predictions changed (GPU check)")
    result.update(coefficients=coefficients, prediction=prediction, optimizer=optimizer, change=change,
                  final=frozen + change,
                  check_rel=float(np.max(np.abs(check - prediction) / np.abs(prediction))))
    return result


def cpu_method_outputs(cpu, out, run, ev, r, args, first, atom, ng_score, rg_score, near, bins, weights, edges,
                       baseline, refs, used, frozen, area, scale, kernels, replay, devices, timing, t_start,
                       body_frame, hv):
    """The CPU tool's products (names, arrays and report fields), plus GPU identities."""
    from cops_cuda import memory, provenance
    labels, basis, sizes, areas = cpu["labels"], cpu["basis"], cpu["sizes"], cpu["areas"]
    n = args.intervals
    late = bins >= 0
    report = dict(
        schema=CPU_SCHEMA, implementation="LiouvilleCOPS_cuda tools/fine_tune_gpu.py --method cpu (GPU responses; "
        "CPU fine_tuning.group_basis and fit_corrections)", created_utc=now(), state="prepared",
        source_restart=provenance.file_identity(args.restart_from), base_run=str(run),
        mesh_identity=r.mesh_identity(), kernels=kernels, devices=devices,
        row_blocks=[list(map(int, b)) for b in ev.blocks], memory_headroom_fraction=memory.headroom_fraction_default(),
        record_precision="fp32" if args.precision == "mixed" else "fp64",
        settings={"minutes": args.minutes, "intervals": n, "dominance_ratio": args.dominance_ratio,
                  "third_group": args.third_group, "lambda": args.lam, "max_nfev": args.max_nfev},
        end_utc=str(args.end_utc or r.utc[-1]), end_source="explicit" if args.end_utc else "last modeled row",
        edges_epoch_s=edges.tolist(), third_group=args.third_group, enhancement_only=True, correction_scale=scale,
        frozen_parameters=["complete fitted H2O field", "Tmax", "T75", "Tmin", "kappa", "all base physics"],
        parameter_names=labels, omitted=cpu["omitted"], atom_face_counts=sizes.tolist(),
        atom_areas_m2=areas.tolist(), requested_parameters=3 * n, support_parameters=len(labels),
        facets={"first_seen_before_window": int((first == -1).sum()), "never_seen": int((first == -2).sum()),
                "first_seen_per_interval": [int((first == j).sum()) for j in range(n)],
                "late_only": int((atom >= 0).sum()), "threshold_near_count": int(len(near))},
        selected_saturated_ng_late=int(np.count_nonzero(np.asarray(r.saturated, bool) & used[:, 0] & late)),
        frozen_replay_max_relative_difference=replay, earlier_response_exactly_zero=True,
        kernel_time_scope="all supplied base rows, including excluded observations; not continuous-time proof")
    valid = atom >= 0
    atom_ng = np.bincount(atom[valid], weights=ng_score[valid], minlength=3 * n)
    atom_rg = np.bincount(atom[valid], weights=rg_score[valid], minlength=3 * n)
    report["groups"] = [dict(name=name, face_count=int(sizes @ basis[:, k]), area_m2=float(areas @ basis[:, k]),
                             first_interval_ng_response_pa_seconds=float(atom_ng @ basis[:, k]),
                             first_interval_rg_response_pa_seconds=float(atom_rg @ basis[:, k]))
                        for k, name in enumerate(labels)]
    out.mkdir(parents=True)
    save_arrays(out / "fine_tuning_prepared.npz",
                dict(first_interval=first, atom=atom, ng_score=ng_score, rg_score=rg_score, baseline=baseline,
                     atom_response=cpu["atom_response"], bins=bins, weights=weights, edges=edges),
                {"schema": CPU_SCHEMA})
    print("Late-only facets: %d; supported coefficients: %d/%d" % (int(valid.sum()), len(labels), 3 * n), flush=True)
    if args.prepare_only:
        report["outputs"] = {"fine_tuning_prepared.npz": provenance.file_identity(out / "fine_tuning_prepared.npz")}
        (out / "fine_tuning.json").write_text(json.dumps(report, indent=1))
        print("Prepared diagnostics: %s" % out)
        return 0
    coefficients, prediction, optimizer, change = cpu["coefficients"], cpu["prediction"], cpu["optimizer"], cpu["change"]
    report.update(state="complete", optimizer=optimizer, regularization_lambda=args.lam,
                  earlier_pressure_max_absolute_change=float(np.max(np.abs(prediction[~late] - baseline[~late]))),
                  added_activity_area_integral=float(np.dot(area, change)),
                  added_activity_area_rms=float(np.sqrt(np.dot(area, change * change) / area.sum())),
                  gpu_check_of_final_field_max_relative_difference=cpu["check_rel"],
                  cpp_verification="not run", scientific_review="pending")
    report["identifiability_warning"] = (
        "Selected late responses are rank deficient; coefficients are not individually resolved."
        if optimizer["data_rank"] < len(optimizer["active_columns"]) else None)
    report["late_log_rms"] = {}
    for g, name in enumerate(("NG", "RG")):
        sel = late & used[:, g] & (refs[:, g] > 0)
        report["late_log_rms"][name] = {lab: float(np.sqrt(np.mean(np.log(v[sel, g] / refs[sel, g]) ** 2)))
                                        for lab, v in (("before", baseline), ("after", prediction))}
    save_arrays(out / "fine_tuning_state.npz",
                dict(coefficients=coefficients, basis=basis, design=cpu["design"], gram=cpu["gram"],
                     activity_frozen=frozen, activity_final=cpu["final"], delta_activity=change, prediction=prediction,
                     refs=refs, used=used, saturated=np.asarray(r.saturated), point_ids=np.asarray(r.point_ids)),
                {"schema": CPU_SCHEMA, "parameter_names": labels, "source_restart": report["source_restart"]})
    with (out / "pressure_comparison.txt").open("x") as stream:
        stream.write("# iPoint UTC interval refNG refRG beforeNG beforeRG afterNG afterRG useNG useRG saturatedNG\n")
        for k, point in enumerate(r.point_ids):
            stream.write("%d %s %d %s %d %d %d\n" % (point, r.utc[k], bins[k], " ".join(
                "%.17e" % v for v in np.r_[refs[k], baseline[k], prediction[k]]), used[k, 0], used[k, 1],
                int(r.saturated[k])))
    export_cpp(out, run, ev.elements, r, cpu["final"], body_frame, args.lam, hv, regularization=CPU_REGULARIZATION)
    if not args.no_figures:
        report["figures"] = figures_for(out, run, r, baseline, prediction, bins, atom, [], devices[0],
                                        not args.no_maps and bool(valid.any()), args, "CPU method")
    report["outputs"] = {p.name: provenance.file_identity(p) for p in sorted(out.iterdir()) if p.is_file()}
    report["timing_s"] = dict(timing, total_s=time.perf_counter() - t_start)
    (out / "fine_tuning.json").write_text(json.dumps(report, indent=1))
    rms = report["late_log_rms"]
    print("coefficients: %s" % ", ".join("%s=%.4g" % (nm, c) for nm, c in zip(labels, coefficients)))
    print("late log rms NG %.4f -> %.4f, RG %.4f -> %.4f; optimizer %s" % (
        rms["NG"]["before"], rms["NG"]["after"], rms["RG"]["before"], rms["RG"]["after"], optimizer["message"]))
    if report["identifiability_warning"]:
        print("WARNING: " + report["identifiability_warning"])
    print("Fine tuning complete: %s; optimizer success=%s" % (out, optimizer["success"]))
    return 0


def solution_report(sol, common, groups, atom, area, frozen, ng_score, rg_score, baseline, refs, used, late, power):
    x, optimizer, final, prediction = sol["x"], sol["optimizer"], sol["final"], sol["prediction"]
    group_info = []
    for k, (code, name) in enumerate(groups):
        m = atom == code
        group_info.append({
            "name": name, "interval": code // 3, "kind": ("ng", "rg", "mixed")[code % 3], "face_count": int(m.sum()),
            "area_m2": float(area[m].sum()),
            "area_fraction_clipped": float(area[m & (frozen <= 0)].sum() / area[m].sum()),
            "frozen_activity_area_mean": float(np.dot(area[m], frozen[m]) / area[m].sum()),
            "f": float(x[2 * k]), "c": float(x[2 * k + 1]), "f_fixed": 2 * k not in optimizer["active_columns"],
            "c_fixed": 2 * k + 1 not in optimizer["active_columns"],
            "final_activity_area_mean": float(np.dot(area[m], final[m]) / area[m].sum()),
            "first_interval_ng_response_pa_s": float(ng_score[m].sum()),
            "first_interval_rg_response_pa_s": float(rg_score[m].sum())})
    late_rms = {}
    for g, name in enumerate(("NG", "RG")):
        sel = late & used[:, g] & (refs[:, g] > 0)
        late_rms[name] = {lab: float(np.sqrt(np.mean(np.log(v[sel, g] / refs[sel, g]) ** 2)))
                          for lab, v in (("before", baseline), ("after", prediction))}
    change = final - frozen
    return dict(common, state="complete", solution=sol["label"],
                fitted_gauges="+".join(("NG", "RG")[g] for g in sol["gauges"]),
                groups=group_info, parameters=x.tolist(), optimizer=optimizer,
                identifiability_warning=("selected late responses are rank deficient (%d of %d); f and c are not "
                                         "individually resolved" % (optimizer["data_rank"],
                                                                     len(optimizer["active_columns"]))
                                         if optimizer["data_rank"] < len(optimizer["active_columns"]) else None),
                late_log_rms=late_rms, earlier_pressures_bitwise_unchanged=bool(sol["outside_rel"] == 0.0),
                earlier_pressures_max_relative_change=sol["outside_rel"],
                gpu_check_of_final_field_max_relative_difference=sol["check_rel"], fit_s=sol["fit_s"],
                activity_change={"area_integral": float(np.dot(area, change)),
                                 "area_rms": float(np.sqrt(np.dot(area, change * change) / area.sum())),
                                 "relative_to_frozen_power": float(np.dot(area, change * change) / power)})


def figures_for(out, run, r, baseline, prediction, bins, atom, groups, device, maps, args, label):
    """Pressure figures with the fine-tuning window shaded (if ``prediction``); lon/lat maps of
    the groups (if ``maps``)."""
    from cops_cuda import figures
    written, errors = [], []
    (out / "figures").mkdir(exist_ok=True)
    try:
        if prediction is None:
            raise StopIteration
        table = {"utc": list(r.utc), "time": np.array([figures.parse_utc(u) for u in r.utc]),
                 "altitude": np.asarray(r.altitude_m, float), "refNG": np.asarray(r.ref_ng, float),
                 "refRG": np.asarray(r.ref_rg, float), "useNG": np.asarray(r.use_ng, float),
                 "useRG": np.asarray(r.use_rg, float), "iPoint": np.asarray(r.point_ids, float)}
        rows = np.flatnonzero(bins >= 0)
        window = figures.late_window(table, start=str(r.utc[rows[0]]), stop=str(r.utc[rows[-1]]))
        window["name"] = "fine-tuning window"
        series = [("fitted state (frozen)", dict(table, pNG=baseline[:, 0], pRG=baseline[:, 1]), (0, (3, 1.2))),
                  ("fine tuned, %s" % label, dict(table, pNG=prediction[:, 0], pRG=prediction[:, 1]), "-")]
        written += figures.pressure_figures(None, series, table, window, str(out / "figures" / "pressures"),
                                            title="%s; fine tuning of the last %g min, %s" % (
                                                Path(args.restart_from).parent.name, args.minutes, label),
                                            zoom_minutes=max(10.0, 2 * args.minutes))
    except StopIteration:
        pass
    except Exception as exc:                      # a figure never fails the run
        errors.append("pressures: %s" % exc)
    if maps:
        try:
            gt = figures.groundtrack_data(run, device_id=device)
            for kind in range(3):
                mask = (atom >= 0) & (atom % 3 == kind)
                if mask.any():
                    written += figures.region_map(gt, mask, str(out / "figures" / ("groups_" + figures_kind(kind))),
                                                  "fine-tuning groups, %s-dominated" % ("NG", "RG")[kind]
                                                  if kind < 2 else "fine-tuning groups, mixed NG/RG",
                                                  device_id=device)
        except Exception as exc:
            errors.append("maps: %s" % exc)
    for e in errors:
        print("figure warning: %s" % e)
    return {"written": [str(p) for p in written], "errors": errors}


def figures_kind(kind):
    return ("ng", "rg", "mixed")[kind]


if __name__ == "__main__":
    sys.exit(main())
