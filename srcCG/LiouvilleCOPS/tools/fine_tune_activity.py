#!/usr/bin/env python3
"""Refine a completed CPU/CUDA fit using first-visible late-interval groups.

No individual-facet optimization and no new harmonics. Read FINE_TUNING_STEP.md
before using the result scientifically. Inputs are immutable; all outputs go
to a new directory outside LiouvilleCOPS. The CUDA branch is not imported.
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np

import fit_cops
import fit_restart as rs
import fine_tuning as ft
from cops_io import file_identity, validate_fit_base, write_json_atomic, utc_now
from cops_paths import output_path, product_input

SCHEMA = "liouville_cops.late_interval_fine_tuning.v1"


def response_iterator(run, header, cols, rows, elements, state):
    for k, row in enumerate(rows):
        yield ft.row_response(run, header, cols, row, elements, state.coef_prior, state.T)
        if (k+1) % 100 == 0 or k+1 == len(rows):
            print("kernel response rows: %d/%d" % (k+1, len(rows)), flush=True)


def save_arrays(path, arrays, metadata):
    """Portable, pickle-free interchange with byte-level numeric identities."""
    metadata = dict(metadata, array_sha256={k: rs.array_sha256(v) for k, v in arrays.items()})
    with Path(path).open("xb") as stream:
        np.savez_compressed(stream, meta_json=np.array(json.dumps(metadata, sort_keys=True)), **arrays)


def read_prepared(path, expected):
    with np.load(path, allow_pickle=False) as archive:
        meta = json.loads(str(archive["meta_json"]))
        arrays = {k: np.array(archive[k]) for k in archive.files if k != "meta_json"}
    if meta.get("definition") != expected:
        raise ValueError("prepared response identity/settings differ from this request")
    if set(meta.get("array_sha256", {})) != set(arrays):
        raise ValueError("prepared response array inventory differs from its hashes")
    for key, value in arrays.items():
        if rs.array_sha256(value) != meta["array_sha256"][key]:
            raise ValueError("corrupt prepared response array: "+key)
    return arrays


def export_cpp(out, run, el, state, final_activity, body_frame, lam):
    """Replay the final physical field using the existing strict facet format.

    Fold all frozen old enhancements and the new correction into the explicit
    residual relative to the saved legacy coefficients. The C++ solver still
    applies one final clip. Do not label this as new spherical harmonics.
    """
    coefficient = out/"ActivityFit.txt"
    coef = state.coef_prior.copy()
    coef[0] = state.a
    coef[1] *= float(state.kappa[0])
    np.savetxt(coefficient, coef, fmt="%.17e")
    residual = final_activity-el["Y"]@state.a
    enhancement = out/"ActivityEnhancementFit.txt"
    with enhancement.open("x") as stream:
        stream.write("# Frozen fitted state plus late-interval enhancement; explicit facet field\n")
        fields = dict(Format=fit_cops.FACET_ENHANCEMENT_FORMAT, FormatVersion=1,
                      BasisID=fit_cops.COMPOSITE_BASIS_ID, Species="H2O", BodyFrame=body_frame,
                      Ordering="increasing_face_index", Clipping="max(0,legacy_linear+enhancement)",
                      HarmonicMinimumDegree=5, HarmonicMaximumDegree=4, Scale=1,
                      ScaleInterpretation="explicit_dimensionless_residual_per_facet",
                      LegacyActivityIdentity="sha256:"+file_identity(coefficient)["sha256"],
                      Regularization="late_interval_total_added_activity_relative_frozen_power_v1",
                      Lambda="%.17g" % lam, FacetCount=len(residual))
        for key, value in fields.items():
            stream.write("%s = %s\n" % (key, value))
        for i, (xyz, area, delta) in enumerate(zip(el["centroid"], el["area"], residual)):
            stream.write("Facet = %d %.17e %.17e %.17e %.17e %.17e\n" %
                         (i, *xyz, area, delta))
    source = run/"input.snapshot.in"
    text = source.read_text()
    replacement = dict(ActivityModel="file", ActivityFile=str(coefficient),
                       ActivityFileHeliocentricModel=state.meta["fine_tuning_helio"],
                       AnalyticalTmax="%.17g" % state.T[0], AnalyticalT75="%.17g" % state.T[1],
                       AnalyticalTmin="%.17g" % state.T[2],
                       ActivityEnhancementModel=fit_cops.FACET_ENHANCEMENT_MODEL,
                       ActivityEnhancementDegree="4", ActivityEnhancementScale="1",
                       ActivityEnhancementFile=str(enhancement), OutputDirectory=str(out)+"_verify")
    # Indented and repeated assignments are all replaced; later definitions
    # must not silently override the fitted temperature or field.
    for key, value in replacement.items():
        pattern = r"^\s*"+re.escape(key)+r"\s*=.*$"
        line = key+" = "+value
        text = re.sub(pattern, lambda _: line, text, flags=re.M) if re.search(pattern, text, re.M) else text+"\n"+line+"\n"
    (out/"cops_fine_tuning.in").write_text(text)


def plot_comparison(path, utc, refs, used, before, after, edges):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    times = np.array([ft.parse_utc(str(u)) for u in utc])
    x = (times-edges[-1])/60
    fig, axes = plt.subplots(2, 2, figsize=(12, 7))
    for g, name in enumerate(("NG", "RG")):
        for k, ax in enumerate(axes[g]):
            ax.semilogy(x, before[:, g], label="Frozen fit", color="C0")
            ax.semilogy(x, after[:, g], label="Fine tuning", color="C1")
            valid = (refs[:, g] > 0) & np.isfinite(refs[:, g])
            selected = valid & used[:, g]
            ax.scatter(x[selected], refs[selected, g], s=8, c="black", label="Selected reference")
            excluded = valid & ~used[:, g]
            ax.scatter(x[excluded], refs[excluded, g], s=14, facecolors="none", edgecolors="gray", label="Excluded reference")
            for boundary in edges:
                ax.axvline((boundary-edges[-1])/60, color="gray", alpha=.25)
            if k:
                ax.set_xlim((edges[0]-edges[-1])/60, 0)
            ax.set_ylabel(name+" pressure [Pa]")
            ax.set_xlabel("Minutes relative to declared end")
            ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", help="original unenhanced static-kernel fitting base")
    parser.add_argument("--restart-from", required=True, help="completed CPU/CUDA fit_restart.npz to freeze")
    parser.add_argument("--minutes", type=float, default=5., help="final window length m [minutes], default 5")
    parser.add_argument("--intervals", type=int, required=True, help="number n of equal-duration intervals")
    parser.add_argument("--end-utc", help="explicit end UTC; default last modeled UTC, recorded in report")
    parser.add_argument("--dominance-ratio", type=float, default=1e-3)
    parser.add_argument("--third-group", choices=("all", "mixed"), default="all")
    parser.add_argument("--lambda", dest="lam", type=float, default=0.)
    parser.add_argument("--max-nfev", type=int, default=1000)
    parser.add_argument("--prepared", help="reuse a matching fine_tuning_prepared.npz (still audit input kernels)")
    parser.add_argument("--prepare-only", action="store_true", help="save group/response diagnostics without fitting")
    parser.add_argument("--out", required=True, help="new output directory outside LiouvilleCOPS")
    args = parser.parse_args(argv)
    run = Path(product_input(args.run))
    restart_path = Path(product_input(args.restart_from))
    out = output_path(args.out).resolve()
    if out.exists():
        parser.error("output already exists: "+str(out))
    if not np.isfinite(args.lam) or args.lam < 0 or args.max_nfev < 1:
        parser.error("lambda must be finite/nonnegative and max-nfev positive")
    if not 0 < args.dominance_ratio < 1:
        parser.error("dominance ratio must be between zero and one")
    audit = validate_fit_base(run)
    if audit["legacy_manifest"]:
        parser.error("fine tuning requires a manifested static base")
    header, names, rows, cols = (audit[k] for k in ("headers", "names", "rows", "columns"))
    state = rs.read(restart_path)
    el = fit_cops.read_elements(run/"elements.bin")
    rs.check_mesh(state, el)
    if (fit_cops.header_value(header, "FlightTimeRotation") != "none" or
            fit_cops.header_value(header, "TemperatureModel") != "analytical"):
        parser.error("fine tuning currently requires static kernels and analytical temperatures")
    # This stage does not silently move a state onto different physical inputs.
    # The mesh identity alone cannot detect changed gas laws or gauge settings.
    for key, path in (("base_elements", run/"elements.bin"),
                      ("base_model_table", run/"COPS_model.dat"),
                      ("base_input", run/"input.snapshot.in")):
        expected = state.meta.get("base", {}).get(key)
        if not expected or file_identity(path)["sha256"] != expected.get("sha256"):
            parser.error("restart/base identity mismatch or missing provenance: "+key)
    if not np.array_equal(state.point_ids, [int(r[cols["iPoint"]]) for r in rows]) or not np.array_equal(
            state.utc, [r[cols["UTC"]] for r in rows]):
        parser.error("restart observation row order differs from base")
    if fit_cops.optional_header_value(header, "ActivityEnhancementModel") not in (None, "none"):
        parser.error("supply the original unenhanced fitting base, not a verification run")
    coef = fit_cops.read_activity_coefficients(str(Path(__file__).resolve().parents[1]), header)
    if not np.array_equal(coef, state.coef_prior):
        parser.error("restart prior coefficients differ from base")
    body_frames = audit["manifest"].get("effective_settings", {}).get("BodyFrame", [])
    if len(body_frames) != 1:
        parser.error("base must declare one BodyFrame")
    bins, weights, edges = ft.intervals(state.utc, args.minutes, args.intervals, args.end_utc)
    frozen = np.maximum(state.activity_linear, 0)
    if np.any(el["area"] <= 0):
        parser.error("nonpositive facet area")
    power = float(np.dot(el["area"], frozen*frozen))
    if not np.isfinite(power) or power <= 0:
        parser.error("frozen activity has zero/nonfinite area-weighted power")
    identity = file_identity(restart_path)
    definition = dict(schema=SCHEMA, restart_sha256=identity["sha256"],
                      elements_sha256=audit["elements"]["sha256"],
                      table_sha256=audit["model_table"]["sha256"],
                      kernels_sha256=audit["kernels"]["set_sha256"],
                      minutes=args.minutes, intervals=args.intervals, edges=edges.tolist(),
                      dominance_ratio=args.dominance_ratio,
                      implementation_sha256=file_identity(ft.__file__)["sha256"],
                      forward_implementation_sha256=file_identity(fit_cops.__file__)["sha256"])
    responses = lambda: response_iterator(run, header, cols, rows, el, state)
    if args.prepared:
        arrays = read_prepared(product_input(args.prepared), definition)
    else:
        first, atom, ng, rg = ft.prepare_support(len(frozen), bins, weights, responses(), args.dominance_ratio)
        baseline, atom_response = ft.pressure_matrices(atom, 3*args.intervals, frozen, float(state.kappa[0]), responses())
        arrays = dict(first_interval=first, atom=atom, ng_score=ng, rg_score=rg,
                      baseline=baseline, atom_response=atom_response, bins=bins, weights=weights, edges=edges)
    baseline, atom_response, atom = arrays["baseline"], arrays["atom_response"], arrays["atom"]
    old_prediction = np.column_stack([state.pred_ng, state.pred_rg])
    if not np.isfinite(old_prediction).all() or not np.allclose(baseline, old_prediction, rtol=1e-6, atol=1e-30):
        parser.error("frozen restart pressures do not reproduce on this base within rtol=1e-6")
    # Structural protection is checked for EACH interval, not merely before
    # the final window. This also catches preparation-cache implementation bugs.
    for j in range(args.intervals):
        if np.any(atom_response[bins < j, :, 3*j:3*j+3] != 0):
            parser.error("nonzero early response in interval group %d" % j)
    basis, labels, omitted, sizes, areas = ft.group_basis(atom, el["area"], args.intervals, args.third_group)
    print("Late-only facets: %d; supported coefficients: %d/%d" %
          (np.count_nonzero(atom >= 0), len(labels), 3*args.intervals), flush=True)
    scale = state.enhancement_scale
    design = scale*(atom_response@basis)
    gram = scale**2*(basis.T@(areas[:, None]*basis))/power
    refs = np.column_stack([state.ref_ng, state.ref_rg])
    used = np.column_stack([state.use_ng, state.use_rg]).astype(bool)
    used &= np.array([not int(r[cols["inside"]]) for r in rows])[:, None]
    report = dict(schema=SCHEMA, created_utc=utc_now(), state="prepared", definition=definition,
                  source_restart=identity, mesh_identity=state.mesh_identity(),
                  end_utc=str(args.end_utc or state.utc[-1]), end_source="explicit" if args.end_utc else "last modeled row",
                  third_group=args.third_group, enhancement_only=True, correction_scale=scale,
                  frozen_parameters=["complete fitted H2O field", "Tmax", "T75", "Tmin", "kappa", "all base physics"],
                  parameter_names=labels, omitted=omitted, atom_face_counts=sizes.tolist(), atom_areas_m2=areas.tolist(),
                  requested_parameters=3*args.intervals, support_parameters=len(labels),
                  selected_saturated_ng_late=int(np.count_nonzero(state.saturated & used[:, 0] & (bins >= 0))),
                  frozen_replay_max_relative_difference=float(np.max(np.abs(baseline-old_prediction)/np.maximum(np.abs(old_prediction), 1e-300))),
                  earlier_response_exactly_zero=True,
                  kernel_time_scope="all supplied base rows, including excluded observations; not continuous-time proof")
    # Compact group inventory avoids storing a dense facet mask for each of
    # the 3n columns. Atom membership plus this small basis is sufficient to
    # recover every group exactly on either CPU or GPU.
    valid_atom = atom >= 0
    atom_ng = np.bincount(atom[valid_atom], weights=arrays["ng_score"][valid_atom], minlength=3*args.intervals)
    atom_rg = np.bincount(atom[valid_atom], weights=arrays["rg_score"][valid_atom], minlength=3*args.intervals)
    report["groups"] = [dict(name=name, face_count=int(sizes@basis[:, k]),
                             area_m2=float(areas@basis[:, k]),
                             first_interval_ng_response_pa_seconds=float(atom_ng@basis[:, k]),
                             first_interval_rg_response_pa_seconds=float(atom_rg@basis[:, k]))
                        for k, name in enumerate(labels)]
    out.mkdir(parents=True, exist_ok=False)
    manifest = out/"fine_tuning.json"
    write_json_atomic(manifest, report)
    try:
        save_arrays(out/"fine_tuning_prepared.npz", arrays, {"definition": definition})
        if args.prepare_only:
            report["outputs"] = {"fine_tuning_prepared.npz": file_identity(out/"fine_tuning_prepared.npz")}
            write_json_atomic(manifest, report)
            print("Prepared diagnostics: %s" % out, flush=True)
            return 0
        coefficients, prediction, optimizer = ft.fit_corrections(
            baseline, design, refs, used, bins >= 0, gram, args.lam, args.max_nfev)
        change = np.zeros(len(frozen))
        valid = atom >= 0
        change[valid] = scale*(basis@coefficients)[atom[valid]]
        final = frozen+change
        if not np.array_equal(prediction[bins < 0], baseline[bins < 0]):
            raise ValueError("earlier predictions changed")
        report.update(state="complete", optimizer=optimizer, regularization_lambda=args.lam,
                      earlier_pressure_max_absolute_change=float(np.max(np.abs(prediction[bins < 0]-baseline[bins < 0]))),
                      added_activity_area_integral=float(np.dot(el["area"], change)),
                      added_activity_area_rms=float(np.sqrt(np.dot(el["area"], change*change)/el["area"].sum())),
                      cpp_verification="not run", scientific_review="pending")
        report["identifiability_warning"] = (
            "Selected late responses are rank deficient; coefficients are not individually resolved."
            if optimizer["data_rank"] < len(optimizer["active_columns"]) else None)
        if report["identifiability_warning"]:
            print("WARNING: "+report["identifiability_warning"], flush=True)
        report["late_log_rms"] = {}
        for g, name in enumerate(("NG", "RG")):
            select = (bins >= 0) & used[:, g] & (refs[:, g] > 0)
            report["late_log_rms"][name] = {label: float(np.sqrt(np.mean(np.log(values[select, g]/refs[select, g])**2)))
                                           for label, values in (("before", baseline), ("after", prediction))}
        save_arrays(out/"fine_tuning_state.npz", dict(coefficients=coefficients, basis=basis, design=design,
                    gram=gram, activity_frozen=frozen, activity_final=final, delta_activity=change,
                    prediction=prediction, refs=refs, used=used, saturated=state.saturated,
                    point_ids=state.point_ids), {"schema": SCHEMA, "parameter_names": labels, "source_restart": identity})
        with (out/"pressure_comparison.txt").open("x") as stream:
            stream.write("# iPoint UTC interval refNG refRG beforeNG beforeRG afterNG afterRG useNG useRG saturatedNG\n")
            for k, point in enumerate(state.point_ids):
                stream.write("%d %s %d %s %d %d %d\n" % (point, state.utc[k], bins[k],
                    " ".join("%.17e" % v for v in np.r_[refs[k], baseline[k], prediction[k]]),
                    used[k, 0], used[k, 1], state.saturated[k]))
        plot_comparison(out/"pressure_comparison.png", state.utc, refs, used, baseline, prediction, edges)
        # Use the same evaluated rows as the optimizer, not the rounded C++
        # replay, and explicitly show the frozen fit when no groups survive.
        from cops_publication import (read_table, evaluated_table, fine_interval,
                                      write_result_figures)
        base_table = read_table(run)
        before_table = evaluated_table(base_table, state.utc, state.point_ids, baseline, refs, used)
        after_table = evaluated_table(base_table, state.utc, state.point_ids, prediction, refs, used)
        figure_products = write_result_figures(
            out, base_table, after_table,
            series=[("frozen fit", before_table, "--"), ("fine tuning", after_table, "-")],
            interval=fine_interval(report["end_utc"], args.minutes),
            window_name="fine-tuning window")
        state.meta["fine_tuning_helio"] = fit_cops.header_value(header, "HeliocentricModel")
        export_cpp(out, run, el, state, final, body_frames[0], args.lam)
        report["outputs"] = {p.name: file_identity(p) for p in out.iterdir() if p.is_file() and p != manifest}
        report["outputs"].update({name: file_identity(out/name) for name in figure_products})
        write_json_atomic(manifest, report)
        print("Fine tuning complete: %s; optimizer success=%s" % (out, optimizer["success"]), flush=True)
    except BaseException as exc:
        report.update(state="failed", error=str(exc)); write_json_atomic(manifest, report)
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
