#!/usr/bin/env python3
"""Inspect, plot, edit and translate fit restart files (ROADMAP C12; CUDA env).

Subcommands
-----------
``show RESTART``
    Summary: mesh, base, settings, optimizer outcome, physical parameters, templates,
    misfit of the stored predictions.

``plot RESTART [RESTART ...] --out PREFIX [--labels ...] [--x time|altitude|both]``
    Comparison plots from the stored predictions, with no recomputation, in the
    publication format of ``cops_cuda/figures.py``: NG pressure, RG pressure and NG/RG
    against the measurements, with the late-visibility window shaded. The window is the
    restart's template interval, else the c05_write_bundle default.

    Writes ``PREFIX_time`` and ``PREFIX_altitude`` (.eps and .png; both regardless of
    ``--x``) and ``PREFIX_summary.json``.

``enhance RESTART --run BASE --region R (--alt-min/--alt-max and/or --utc-start/--utc-stop)
          --mode add|factor --value V --out PREFIX [--save NEW_RESTART] [--name NAME]``
    Raises the H2O activity on a surface region defined by visibility from the
    observation points of an altitude (and/or time) range, then compares the model with
    and without the change. Regions (``masks.visibility_masks``, union over the selected
    points of the same-point predicates):
    * ``ng_not_rg``: seen by the NG and not by the RG;
    * ``rg_not_ng``: seen by the RG and not by the NG;
    * ``ng_seen``, ``rg_seen``, ``any_seen``.

    The region is then **confined to the window** by a significance threshold (user decision
    2026-10-20). Strict exclusion of every facet visible from outside empties the region:
    facets seen during the descent are also seen, weakly, from higher up. A candidate facet is
    ranked by its own influence: the largest relative change it alone causes in the NG or RG
    pressure of any point outside the selected set. The change is evaluated at the restart
    state, including the clip (``kernels/restart.cu``).
    * ``--confine total`` (default): the kept facets are the largest prefix of that ranking whose
      **summed** change keeps every outside point within ``--epsilon`` (relative; default 1e-3).
      The pressures are linear in each facet's clipped activity, so the sum is exact.
    * ``--confine facet``: each facet's own change is ≤ ε. The sum can exceed ε; on a large
      region it reached 16%.

    The total outside change, from the full GPU evaluation, is reported and agrees with the
    prediction.

    The change is a new template column (nonnegative field), so it enters the model exactly as
    the fitter's templates do: linearly, before the single clip.
    * ``add``: column = region mask, amplitude V ≥ 0. This adds V·A_scale to the linear activity;
      it also acts on facets whose activity is clipped.
    * ``factor``: column = mask · max(act_lin, 0)/A_scale, amplitude V − 1, for any V > 0. This
      multiplies the positive activity on the region by V; V < 1 reduces it. Clipped facets
      are unchanged.
    A saved reduction (amplitude < 0) evaluates and plots correctly. The fitter keeps template
    amplitudes ≥ 0, so a refit with ``--restart-from`` clips it to 0 and reports that.

    The pressures of both states are computed on the GPU (``restart.StateEvaluator``).
    ``--save`` writes the modified restart and a template bundle with the new column
    (``<save>.templates.npz``, written by the CPU's own bundle writer), so that
    ``gpu_fit.py --restart-from`` can refit with the new region's amplitude free.

``translate RESTART --to-run BASE_B --save NEW_RESTART [--target-bundle BUNDLE_B] [--evaluate]``
    Carries a restart from its mesh (A) to the mesh of base run B.
    * Legacy coefficients, harmonic amplitudes, T and κ are surface-independent functions or
      scalars and are carried unchanged.
    * Harmonic and template amplitudes are rescaled by A_scale(A)/A_scale(B), so that the
      absolute activity they add is preserved.
    * A template whose name exists in ``--target-bundle`` (masks recomputed physically on mesh B,
      e.g. by ``c05_write_bundle.py``; preferred) takes that column. Any other template is
      mapped by nearest facet centroid and renamed ``xl_<name>``; the mapping statistics are
      recorded. Cliff masks mapped this way carry no gravity provenance, so recompute them on B
      instead.
    * ``--evaluate`` computes the predictions on base B (GPU).
"""

import argparse
import datetime as dt
import json
from pathlib import Path
import sys

import numpy as np

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))

REGION_BITS = {"ng_seen": ("ng_seen",), "rg_seen": ("rg_seen",), "ng_not_rg": ("late_ng_not_rg",),
               "rg_not_ng": ("late_rg_not_ng",), "any_seen": ("ng_seen", "rg_seen")}


def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()


# ----------------------------------------------------------------------------- misfit
def misfit(pred_ng, pred_rg, ref_ng, ref_rg, use_ng, use_rg):
    """rms log residuals and the fitter's data cost (equal weight per gauge)."""
    floor = 1e-30
    sel_n = (use_ng > 0) & (ref_ng > 0) & np.isfinite(pred_ng)
    sel_r = (use_rg > 0) & (ref_rg > 0) & np.isfinite(pred_rg)
    rn = np.log(np.maximum(pred_ng[sel_n], floor)) - np.log(ref_ng[sel_n])
    rr = np.log(np.maximum(pred_rg[sel_r], floor)) - np.log(ref_rg[sel_r])
    if not rn.size or not rr.size:
        return {"points_ng": int(rn.size), "points_rg": int(rr.size)}
    return {"rms_log_ng": float(np.sqrt(np.mean(rn ** 2))), "rms_log_rg": float(np.sqrt(np.mean(rr ** 2))),
            "data_cost": float(0.5 * (np.mean(rn ** 2) + np.mean(rr ** 2))),
            "points_ng": int(rn.size), "points_rg": int(rr.size)}


# ------------------------------------------------------------------------------- plots
def parse_utc(s):
    day, clock = s.split("T")
    h, m, sec = clock.split(":")
    return dt.datetime.fromisoformat(day) + dt.timedelta(hours=int(h), minutes=int(m), seconds=float(sec))


def late_window_of(r, table):
    """The late-visibility window of a restart: its templates' interval, else the default."""
    from cops_cuda import figures
    for src in r.meta.get("templates", {}).get("sources", []) or []:
        interval = ((src or {}).get("definition") or {}).get("interval_utc")
        if interval:
            return figures.late_window(table, start=interval[0], stop=interval[1])
    return figures.late_window(table)


def plot_models(series, data, out_prefix, axes=("time", "altitude"), title="", restart=None):
    """Publication pressure figures (cops_cuda/figures.py) of restart predictions.

    series: [(label, pred_ng, pred_rg)]; data: dict with utc, altitude_m, ref/use arrays.
    Both ``<out_prefix>_time`` and ``_altitude`` are written, as .eps and .png; ``axes`` is
    kept for the command line's ``--x`` and selects nothing.
    """
    from cops_cuda import figures
    del axes
    table = {"utc": list(data["utc"]), "time": np.array([figures.parse_utc(u) for u in data["utc"]]),
             "altitude": np.asarray(data["altitude_m"], float), "refNG": np.asarray(data["ref_ng"], float),
             "refRG": np.asarray(data["ref_rg"], float), "useNG": np.asarray(data["use_ng"], float),
             "useRG": np.asarray(data["use_rg"], float), "iPoint": np.arange(len(data["utc"]), dtype=float)}
    if "point_ids" in data:
        table["iPoint"] = np.asarray(data["point_ids"], float)
    window = late_window_of(restart, table) if restart is not None else figures.late_window(table)
    styles = ["-", (0, (3, 1.2)), (0, (1, 1.2)), (0, (4, 1.2, 1, 1.2))]
    series_tables = [(label, dict(table, pNG=np.asarray(png, float), pRG=np.asarray(prg, float)), styles[k % 4])
                     for k, (label, png, prg) in enumerate(series)]
    return figures.pressure_figures(None, series_tables, table, window, out_prefix, title=title or None)


def plot_region(prefix, run, mask, title, device=0):
    """Longitude/latitude map of a facet region with the ground track (cops_cuda/figures.py)."""
    from cops_cuda import figures
    gt = figures.groundtrack_data(run, device_id=device)
    return figures.region_map(gt, mask, prefix, title, device_id=device)


def not_fitted(meta, note):
    """An edited or translated state keeps its source's optimizer record, but not as its own."""
    if "optimizer" in meta:
        meta["source_optimizer"] = meta.pop("optimizer")
    meta["state_note"] = note


def stored_data(r):
    return {k: np.asarray(getattr(r, k)) for k in ("utc", "altitude_m", "ref_ng", "ref_rg", "use_ng", "use_rg",
                                                   "point_ids")}


# --------------------------------------------------------------------------- commands
def cmd_show(args):
    from cops_cuda import restart as rs
    r = rs.read(args.restart)
    m = r.meta
    print("restart %s (%s)" % (args.restart, m.get("kind")))
    print("  mesh %s, %d facets; base elements %s" % (m["mesh"]["identity"][:16], m["mesh"]["faces"],
                                                      (m.get("base") or {}).get("base_elements", {}).get("path")))
    print("  settings %s" % m.get("settings"))
    print("  optimizer %s" % m.get("optimizer", "none (%s)" % m.get("state_note", "")))
    print("  T = %s K, kappa = %.6g, A_scale = %.6g" % (np.round(r.T, 4).tolist(), float(r.kappa[0]),
                                                         r.enhancement_scale))
    for n, v in zip(r.template_names, r.template_values):
        print("  template %-28s c = %.6g  (%d facets)" % (n, v, int(np.count_nonzero(
            r.template_columns[:, r.template_names.index(n)]))))
    if len(r.harmonic_values):
        print("  harmonic amplitudes: %d (|b| max %.3g)" % (len(r.harmonic_values), np.max(np.abs(r.harmonic_values))))
    print("  stored predictions: %s" % misfit(r.pred_ng, r.pred_rg, r.ref_ng, r.ref_rg, r.use_ng, r.use_rg))
    print("  history: %s" % [h.get("action") for h in m.get("history", [])])
    return 0


def cmd_plot(args):
    from cops_cuda import restart as rs
    states = [rs.read(p) for p in args.restarts]
    labels = args.labels or [Path(p).parent.name or Path(p).stem for p in args.restarts]
    if len(labels) != len(states):
        raise SystemExit("give one label per restart")
    data = stored_data(states[0])
    for s in states[1:]:
        if not np.array_equal(s.point_ids, states[0].point_ids):
            raise SystemExit("restarts cover different observation rows; plot them separately")
    series = [(lab, np.asarray(s.pred_ng), np.asarray(s.pred_rg)) for lab, s in zip(labels, states)]
    axes = ("time", "altitude") if args.x == "both" else (args.x,)
    written = plot_models(series, data, args.out, axes, title="restart comparison", restart=states[0])
    summary = {lab: misfit(s.pred_ng, s.pred_rg, s.ref_ng, s.ref_rg, s.use_ng, s.use_rg)
               for lab, s in zip(labels, states)}
    with open(args.out + "_summary.json", "x") as f:
        json.dump({"restarts": args.restarts, "labels": labels, "misfit": summary}, f, indent=2)
    print("wrote %s and %s_summary.json" % (written, args.out))
    for lab, v in summary.items():
        print("  %-30s %s" % (lab, v))
    return 0


def select_points(ev, args):
    """Rows of the base table in the altitude/time window, outside the nucleus."""
    from cops_cuda import reference
    parse = reference.load("cops_selection").parse_utc
    c = ev.cols
    pts = []
    for row in ev.rows:
        alt = float(row[c["altitude"]])
        if int(row[c["inside"]]):
            continue
        if args.alt_min is not None and alt < args.alt_min:
            continue
        if args.alt_max is not None and alt > args.alt_max:
            continue
        if args.utc_start and parse(row[c["UTC"]]) < parse(args.utc_start):
            continue
        if args.utc_stop and parse(row[c["UTC"]]) > parse(args.utc_stop):
            continue
        pts.append(int(row[c["iPoint"]]))
    return pts


def cmd_enhance(args):
    from cops_cuda import masks, provenance, reference
    from cops_cuda import restart as rs
    templates = reference.load("activity_templates")
    r = rs.read(args.restart)
    if args.alt_min is None and args.alt_max is None and not (args.utc_start or args.utc_stop):
        raise SystemExit("give an altitude range (--alt-min/--alt-max) and/or a time range")
    if args.mode == "factor" and not args.value > 0.0:
        raise SystemExit("--mode factor needs V > 0 (V < 1 reduces the local activity)")
    if args.mode == "add" and args.value < 0.0:
        raise SystemExit("--mode add needs V >= 0")
    ev = rs.StateEvaluator(args.run, args.device)
    try:
        points = select_points(ev, args)
        vis, vmeta = masks.visibility_masks(args.run, points, device_id=args.device)
        candidate = np.zeros(len(r.area), bool)
        for key in REGION_BITS[args.region]:
            candidate |= vis[key]
        if not candidate.any():
            raise SystemExit("no facet satisfies %s for the selected points" % args.region)
        act = np.asarray(r.activity_linear, float)
        scale = r.enhancement_scale
        # Confine the effect to the window (user decision 2026-10-20: significance threshold).
        # A candidate facet is kept only if its change alters no OUTSIDE point's NG or RG pressure
        # by more than epsilon (relative), at the restart state.
        if args.mode == "add":
            changed = np.where(candidate, act + args.value * scale, act)
        else:
            changed = np.where(candidate, act + (args.value - 1.0) * np.maximum(act, 0.0), act)
        delta_af = np.maximum(changed, 0.0) - np.maximum(act, 0.0)
        chosen = set(points)
        outside = [int(p) for p in ev.point_ids if int(p) not in chosen]
        share_ng, share_rg, state_ng, state_rg = ev.facet_point_shares(r, delta_af, outside)
        share = np.maximum(share_ng, share_rg)
        out_flag = np.isin(ev.point_ids, np.asarray(outside, dtype=int))

        def outside_total(threshold):
            """Max over outside points of the relative total pressure change (NG, RG)."""
            keep = candidate & (share <= threshold)
            dng, drg = ev.point_changes(r, np.where(keep, delta_af, 0.0))
            with np.errstate(divide="ignore", invalid="ignore"):
                rn = np.where(state_ng[out_flag] > 0, np.abs(dng[out_flag]) / state_ng[out_flag], 0.0)
                rr = np.where(np.abs(state_rg[out_flag]) > 0, np.abs(drg[out_flag]) / np.abs(state_rg[out_flag]), 0.0)
            return float(max(rn.max(initial=0.0), rr.max(initial=0.0)))

        if args.confine == "facet":
            threshold = args.epsilon
        else:
            # Facets in order of increasing share; the largest prefix whose SUMMED change keeps every
            # outside point within epsilon (monotone: all changes of a region have one sign).
            levels = np.unique(share[candidate])
            lo, hi = -1, len(levels) - 1
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if outside_total(levels[mid]) <= args.epsilon:
                    lo = mid
                else:
                    hi = mid - 1
            threshold = levels[lo] if lo >= 0 else -1.0
        mask = candidate & (share <= threshold)
        predicted_outside = outside_total(threshold) if mask.any() else 0.0
        if not mask.any():
            raise SystemExit("the region is empty: no subset of the %d candidate facets keeps the outside "
                             "pressures within epsilon=%g (confine=%s; smallest facet share %.3g)"
                             % (int(candidate.sum()), args.epsilon, args.confine, float(share[candidate].min())))
        if args.mode == "add":
            column, amplitude = mask.astype(float), float(args.value)
        else:
            # column >= 0 (a valid template field); amplitude V-1 (negative for V < 1: a reduction)
            column, amplitude = mask * np.maximum(act, 0.0) / scale, float(args.value) - 1.0
            if not np.any(column > 0):
                raise SystemExit("factor mode: the activity is clipped (zero) on the whole region; use --mode add")
        name = args.name or ("enh_%s_%s" % (args.region, "_".join(
            "%s%g" % (k, v) for k, v in (("a", args.alt_min), ("b", args.alt_max)) if v is not None))
                             ).replace(".", "p").replace("-", "m")[:64]
        templates.validate_name(name)
        if name in r.template_names:
            raise SystemExit("template %s already exists in the restart" % name)
        # baseline and modified pressures on the GPU (same backend basis, new amplitude 0 / value)
        basis = np.column_stack((rs.enhancement_basis(r), column))
        b0 = np.concatenate((rs.b_vector(r), [0.0]))
        b1 = np.concatenate((rs.b_vector(r), [amplitude]))
        p0 = ev.pressures(r, basis=basis, b=b0)
        p1 = ev.pressures(r, basis=basis, b=b1)
    finally:
        ev.release()
    data = stored_data(r)
    clipped = float(r.area[mask & (act <= 0)].sum() / r.area[mask].sum())
    label = "%s %s %g on %s" % (name, args.mode, args.value, args.region)
    written = plot_models([("restart state", p0[0], p0[1]), (label, p1[0], p1[1])], data, args.out,
                          title="%s; %d points selected; region %d facets (%.1f%% of its area clipped)"
                                % (Path(args.restart).name, len(points), int(mask.sum()), 100 * clipped),
                          restart=r)
    written += plot_region(args.out + "_region", args.run, mask, "%s, points %d" % (args.region, len(points)),
                           args.device)
    out_rows = ~np.isin(np.asarray(r.point_ids), np.asarray(points))
    outside_change = {"NG": float(np.max(np.abs(p1[0][out_rows] / p0[0][out_rows] - 1))) if out_rows.any() else 0.0,
                      "RG": float(np.max(np.abs(p1[1][out_rows] / p0[1][out_rows] - 1))) if out_rows.any() else 0.0,
                      "bitwise_unchanged": bool(np.array_equal(p0[0][out_rows], p1[0][out_rows]) and
                                                np.array_equal(p0[1][out_rows], p1[1][out_rows]))}
    in_rows = ~out_rows
    inside_change = {"NG": float(np.max(np.abs(p1[0][in_rows] / p0[0][in_rows] - 1))) if in_rows.any() else 0.0,
                     "RG": float(np.max(np.abs(p1[1][in_rows] / p0[1][in_rows] - 1))) if in_rows.any() else 0.0}
    summary = {"region": args.region, "points": len(points), "alt_range_m": [args.alt_min, args.alt_max],
               "utc_range": [args.utc_start, args.utc_stop], "mode": args.mode, "value": args.value,
               "template_name": name, "facets": int(mask.sum()),
               "candidate_facets": int(candidate.sum()),
               "epsilon": args.epsilon, "confine": args.confine,
               "excluded_facets": int((candidate & ~mask).sum()),
               "max_facet_share_kept": float(share[mask].max()),
               "predicted_outside_change": predicted_outside,
               "amplitude": amplitude,
               "outside_window_change": outside_change,
               "inside_window_change": inside_change,
               "area_fraction": float(r.area[mask].sum() / r.area.sum()),
               "region_area_fraction_clipped": clipped,
               "misfit_restart_state": misfit(p0[0], p0[1], r.ref_ng, r.ref_rg, r.use_ng, r.use_rg),
               "misfit_enhanced": misfit(p1[0], p1[1], r.ref_ng, r.ref_rg, r.use_ng, r.use_rg),
               "max_rel_change": {"NG": float(np.max(np.abs(p1[0] / p0[0] - 1))),
                                  "RG": float(np.max(np.abs(p1[1] / p0[1] - 1)))}}
    if args.save:
        save = Path(args.save)
        bundle_path = save.with_name(save.stem + ".templates.npz")
        definition = {"definition": "union over points of the same-point predicate %s" % args.region,
                      "points": points, "alt_range_m": [args.alt_min, args.alt_max],
                      "utc_range": [args.utc_start, args.utc_stop], "mode": args.mode, "value": args.value,
                      "computed_by": vmeta["backend"], "kernel_set_sha256": vmeta["kernel_set_sha256"],
                      "interpretation": "user-specified activity enhancement region; not a fitted result"}
        templates.write_bundle(bundle_path, [name], column[:, None], r.centroid, r.area,
                               {"created_utc": now(), "base_run": str(Path(args.run).resolve()),
                                "templates": {name: definition},
                                "source_restart": provenance.file_identity(args.restart)})
        new = r.copy()
        new.arrays["template_columns"] = np.column_stack((r.template_columns, column))
        new.arrays["template_values"] = np.concatenate((r.template_values, [amplitude]))
        new.arrays["activity_linear"] = act + scale * amplitude * column
        new.arrays["pred_ng"], new.arrays["pred_rg"] = p1[0], p1[1]
        new.meta["templates"]["names"].append(name)
        new.meta["templates"]["sources"].append({"name": name, "bundle": provenance.file_identity(bundle_path),
                                                 "mesh_identity": r.mesh_identity(), "definition": definition})
        new.meta["kind"] = "edited"
        not_fitted(new.meta, "enhanced by restart_tool.py; the new amplitude is user-set, not fitted")
        new.meta["history"].append({"utc": now(), "action": "enhance", "parameters": {
            k: getattr(args, k) for k in ("region", "alt_min", "alt_max", "utc_start", "utc_stop", "mode", "value")},
            "template": name, "source": provenance.file_identity(args.restart)})
        summary["saved"] = rs.write(save, new)
        if amplitude < 0:
            print("note: the saved amplitude of %s is %.6g < 0 (a reduction). Templates are nonnegative "
                  "in the fitter, so a refit with --restart-from clips it to 0 (reported in the start "
                  "record); the restart state itself evaluates and plots the reduction." % (name, amplitude))
        summary["bundle"] = provenance.file_identity(bundle_path)
    with open(args.out + "_summary.json", "x") as f:
        json.dump(summary, f, indent=2)
    print("wrote %s and %s_summary.json" % (written, args.out))
    print(json.dumps({k: summary[k] for k in ("facets", "candidate_facets", "epsilon", "confine", "excluded_facets",
                                              "max_facet_share_kept", "predicted_outside_change",
                                              "amplitude", "inside_window_change", "outside_window_change",
                                              "area_fraction", "region_area_fraction_clipped",
                                              "misfit_restart_state", "misfit_enhanced", "max_rel_change")},
                     indent=1))
    return 0


def cmd_translate(args):
    from scipy.spatial import cKDTree
    from cops_cuda import provenance, reference
    from cops_cuda import restart as rs
    templates = reference.load("activity_templates")
    r = rs.read(args.restart)
    ev = rs.StateEvaluator(args.to_run, args.device)
    try:
        el = ev.elements
        cB, aB = np.asarray(el["centroid"], float), np.asarray(el["area"], float)
        dist, idx = cKDTree(r.centroid).query(cB)
        scale_a, scale_b = r.enhancement_scale, ev.enhancement_scale()
        ratio = scale_a / scale_b
        target = templates.load_bundle(args.target_bundle, cB, aB) if args.target_bundle else None
        names, columns, sources, values = [], [], [], []
        mapped = []
        for k, (n, v) in enumerate(zip(r.template_names, r.template_values)):
            if target is not None and n in target["names"]:
                col = templates.select_templates(target, [n])[:, 0]
                names.append(n)
                sources.append({"name": n, "bundle": provenance.file_identity(args.target_bundle),
                                "mesh_identity": target["mesh_identity"],
                                "definition": target["metadata"].get("templates", {}).get(n)})
            else:
                col = r.template_columns[idx, k]
                new_name = ("xl_" + n)[:64]
                names.append(new_name)
                mapped.append((new_name, n, col))
                sources.append(None)
            columns.append(col)
            values.append(float(v) * ratio)
        stats = {"method": "nearest facet centroid", "distance_m": {
            "median": float(np.median(dist)), "p95": float(np.percentile(dist, 95)), "max": float(dist.max())},
            "enhancement_scale_source": scale_a, "enhancement_scale_target": scale_b, "amplitude_ratio": ratio}
        save = Path(args.save)
        if mapped:
            bundle_path = save.with_name(save.stem + ".templates.npz")
            defs = {}
            for new_name, old, col in mapped:
                src = next(s for s in r.meta["templates"]["sources"] if s["name"] == old)
                defs[new_name] = {"definition": "translated column (nearest facet centroid) of template %s" % old,
                                  "translation": stats, "source_restart": provenance.file_identity(args.restart),
                                  "source_bundle": src.get("bundle"), "source_mesh_identity": r.mesh_identity(),
                                  "interpretation": "geometric transfer, not recomputed on this mesh"}
            templates.write_bundle(bundle_path, [m[0] for m in mapped], np.column_stack([m[2] for m in mapped]),
                                   cB, aB, {"created_utc": now(), "base_run": str(Path(args.to_run).resolve()),
                                            "templates": defs})
            ident = provenance.file_identity(bundle_path)
            for i, s in enumerate(sources):
                if s is None:
                    sources[i] = {"name": names[i], "bundle": ident, "mesh_identity": templates.mesh_identity(cB, aB),
                                  "definition": defs[names[i]]}
        new = rs.Restart(json.loads(json.dumps(r.meta)), {})
        nf = len(aB)
        cols = np.column_stack(columns) if columns else np.empty((nf, 0))
        harmonic_b = np.asarray(r.harmonic_values, float) * ratio
        tmp = rs.Restart(new.meta, {"centroid": cB, "template_columns": cols, "harmonic_values": harmonic_b,
                                    "template_values": np.asarray(values, float), "a": r.a})
        new.meta["templates"] = {"names": names, "sources": sources}
        new.meta["enhancement_scale"] = scale_b
        new.meta["mesh"] = {"identity": templates.mesh_identity(cB, aB), "faces": int(nf)}
        new.meta["base"] = {"base_elements": provenance.file_identity(Path(args.to_run) / "elements.bin"),
                            "base_model_table": provenance.file_identity(Path(args.to_run) / "COPS_model.dat")}
        new.meta["kind"] = "translated"
        not_fitted(new.meta, "translated to another mesh; not a fit result on this mesh. Selection flags "
                             "and observation rows are the target base table's")
        new.meta["translation"] = stats
        new.meta["history"] = list(r.meta.get("history", [])) + [{
            "utc": now(), "action": "translate", "source": provenance.file_identity(args.restart),
            "to_run": str(Path(args.to_run).resolve()), "statistics": stats}]
        act = el["Y"] @ np.asarray(r.a, float)
        b = np.concatenate((harmonic_b, values))
        if b.size:
            act = act + scale_b * (rs.enhancement_basis(tmp) @ b)
        table = ev.table_arrays()
        new.arrays = {"a": np.asarray(r.a, float), "coef_prior": np.asarray(ev.coef, float),
                      "T": np.asarray(r.T, float), "kappa": np.asarray(r.kappa, float),
                      "harmonic_values": harmonic_b, "template_values": np.asarray(values, float),
                      "template_columns": cols, "activity_linear": act, "centroid": cB, "area": aB,
                      "saturated": np.zeros(len(table["point_ids"]), bool),
                      "x": np.empty(0), "x0": np.empty(0), "lo": np.empty(0), "hi": np.empty(0), **table,
                      "pred_ng": np.full(len(table["point_ids"]), np.nan),
                      "pred_rg": np.full(len(table["point_ids"]), np.nan)}
        if args.evaluate:
            png, prg = ev.pressures(new)
            new.arrays["pred_ng"], new.arrays["pred_rg"] = png, prg
    finally:
        ev.release()
    ident = rs.write(save, new)
    print("wrote %s; templates %s; mapping %s" % (save, names, stats))
    if args.evaluate:
        print("  predictions on the target base: %s" % misfit(new.pred_ng, new.pred_rg, new.ref_ng, new.ref_rg,
                                                               new.use_ng, new.use_rg))
    del ident
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("show")
    s.add_argument("restart")
    s = sub.add_parser("plot")
    s.add_argument("restarts", nargs="+")
    s.add_argument("--out", required=True, help="output prefix (new files)")
    s.add_argument("--labels", nargs="+")
    s.add_argument("--x", choices=("time", "altitude", "both"), default="both")
    s = sub.add_parser("enhance")
    s.add_argument("restart")
    s.add_argument("--run", required=True, help="base run of the restart's mesh")
    s.add_argument("--region", choices=sorted(REGION_BITS), required=True)
    s.add_argument("--alt-min", type=float)
    s.add_argument("--alt-max", type=float)
    s.add_argument("--utc-start")
    s.add_argument("--utc-stop")
    s.add_argument("--mode", choices=("add", "factor"), required=True)
    s.add_argument("--epsilon", type=float, default=1e-3,
                   help="bound on the relative change of any outside point's NG or RG pressure (default 1e-3)")
    s.add_argument("--confine", choices=("total", "facet"), default="total",
                   help="total: the summed change of the kept facets is <= epsilon at every outside point "
                        "(default); facet: each facet's own change is <= epsilon (the sum may exceed it)")
    s.add_argument("--value", type=float, required=True)
    s.add_argument("--name")
    s.add_argument("--out", required=True, help="output prefix for plots and summary (new files)")
    s.add_argument("--save", help="write the modified restart here (new file)")
    s.add_argument("--device", type=int, default=0)
    s = sub.add_parser("translate")
    s.add_argument("restart")
    s.add_argument("--to-run", required=True, help="base run of the target mesh")
    s.add_argument("--save", required=True, help="new restart for the target mesh")
    s.add_argument("--target-bundle", help="template bundle computed on the target mesh (preferred source)")
    s.add_argument("--evaluate", action="store_true", help="compute the predictions on the target base (GPU)")
    s.add_argument("--device", type=int, default=0)
    args = p.parse_args()
    return {"show": cmd_show, "plot": cmd_plot, "enhance": cmd_enhance, "translate": cmd_translate}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
