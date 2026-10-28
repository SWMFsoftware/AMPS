#!/usr/bin/env python3
"""tools/fine_tune_gpu.py (--method cpu) with the dust-grain term refitted (user decisions 2026-10-27).

    .venv/bin/python <this> <fine_tune_gpu.py arguments, --method cpu>

Local, untracked extension (model_ext/); fine_tune_gpu.py, cops_cuda and the CPU reference are not edited.
The fit being tuned must come from model_ext/grain_fit.py (its directory holds grain_term.json).

What changes against the unchanged fine tuning (--method cpu):
* the frozen baseline is the fit's full model: surface pressures (GPU, from the restart) + the fit's grain on NG,
  so the restart's stored predictions (which include the grain) reproduce;
* the optimizer is ``fit_corrections_with_grain``: the CPU's ``fit_corrections`` (nonnegative additive group
  amplitudes c, equal total weight per gauge, log residuals, lambda Gram penalty on c, bounded TRF, x_scale jac,
  tolerances 1e-10, zero-response columns fixed at 0) plus the grain amplitude A and slope B (>= 0, start = the
  fit's values) refitted jointly (user 2026-10-27: "Refit A, B with the groups");
* the selected data are the CPU's late points (both gauges, last --minutes) plus the NG points from the grain
  onset t_end to the window start (user 2026-10-27: "Window + grain-affected NG"), equal total weight per gauge;
* earlier rows: before t_end the predictions must stay exactly the baseline (asserted); between t_end and the
  window only the grain may change them (asserted: the change equals the grain change).
Products: those of fine_tune_gpu.py, plus grain_finetune.json (A, B before/after, data, crossing check).
"""

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(APP), str(APP / "tools"), str(HERE)]

import grain_crossing as gc  # noqa: E402


def main():
    argv = sys.argv[1:]
    restart = Path(argv[argv.index("--restart-from") + 1])
    out = Path(argv[argv.index("--out") + 1])
    run = Path(argv[0])
    term = json.loads((restart.parent / "grain_term.json").read_text())
    A0, B0 = term["result"]["A_Pa"], term["result"]["B_Pa_per_s"]
    t_end = float(term["t_end_et"])
    ets = gc.et_by_point(run)
    state = {"A0_Pa": A0, "B0_Pa_per_s": B0, "t_end_et": t_end, "source_grain_term": str(restart.parent / "grain_term.json"),
             "module_sha256": {p: hashlib.sha256((HERE / p).read_bytes()).hexdigest()
                               for p in ("grain_crossing.py", "grain_finetune.py")},
             "policy": "user 2026-10-27: grain A, B refitted with the late groups; data = the CPU late points "
                       "+ NG points from the grain onset to the window start; equal total weight per gauge"}

    import fine_tune_gpu as ftg_tool
    from cops_cuda import fine_tuning as ftg
    from cops_cuda import reference
    reference.bind()
    ft = reference.load("fine_tuning")

    def grain(A, B, et):
        s = et - t_end
        g = np.where(s >= 0, A - B * np.maximum(s, 0), 0.0)
        dA = np.where((s >= 0) & (g >= 0), 1.0, 0.0)
        dB = np.where((s >= 0) & (g > 0), -np.maximum(s, 0), 0.0)
        return np.maximum(g, 0.0), dA, dB

    # ---- baseline = surface (GPU) + the fit's grain on NG
    original_pressures = ftg.SplitEvaluator.pressures

    def pressures(self, r):
        png, prg = original_pressures(self, r)
        et = np.array([ets[int(p)] for p in r.point_ids])
        state["et"] = et
        return png + grain(A0, B0, et)[0], prg
    ftg.SplitEvaluator.pressures = pressures

    # ---- the optimizer
    def fit_corrections_with_grain(baseline, design, refs, used, late, gram, lam=0., max_nfev=1000):
        from scipy.optimize import least_squares
        et = state["et"]
        g0 = grain(A0, B0, et)[0]
        surface = baseline.copy()
        surface[:, 0] -= g0
        onset = (et >= t_end) & ~late
        selected = used & late[:, None] & (refs > 0)
        selected[:, 0] |= used[:, 0] & onset & (refs[:, 0] > 0)
        counts = selected.sum(axis=0)
        if np.any(counts == 0):
            raise ValueError("late refinement requires selected positive NG and RG observations")
        row, gauge = np.nonzero(selected)
        d, p, target = design[row, gauge], surface[row, gauge], refs[row, gauge]
        isng = gauge == 0
        weights = 1 / np.sqrt(counts[gauge])
        eig, vec = np.linalg.eigh(gram)
        penalty = np.sqrt(lam) * np.sqrt(np.maximum(eig, 0))[:, None] * vec.T
        nc = design.shape[-1]
        sA, sB = max(A0, 1e-9), max(B0, 1e-12) if B0 > 0 else max(A0, 1e-9) / 400.0   # parameter units

        def model(z):
            c, A, B = z[:nc], z[nc] * sA, z[nc + 1] * sB
            gr, dA, dB = grain(A, B, et[row])
            return p + d @ c + np.where(isng, gr, 0.0), dA, dB

        def residual(z):
            m = model(z)[0]
            return np.r_[weights * np.log(m / target), penalty @ z[:nc]]

        def jacobian(z):
            m, dA, dB = model(z)
            J = np.zeros((len(m) + penalty.shape[0], nc + 2))
            J[:len(m), :nc] = weights[:, None] * d / m[:, None]
            J[:len(m), nc] = np.where(isng, weights * dA * sA / m, 0.0)
            J[:len(m), nc + 1] = np.where(isng, weights * dB * sB / m, 0.0)
            J[len(m):, :nc] = penalty
            return J

        active = np.r_[np.flatnonzero(np.any(d != 0, axis=0)), [nc, nc + 1]]
        z = np.r_[np.zeros(nc), A0 / sA, B0 / sB]

        def expand(v):
            w = z.copy()
            w[active] = v
            return w
        result = least_squares(lambda v: residual(expand(v)), z[active], jac=lambda v: jacobian(expand(v))[:, active],
                               bounds=(0, np.inf), x_scale="jac", max_nfev=max_nfev, ftol=1e-10, xtol=1e-10, gtol=1e-10)
        z = expand(result.x)
        c, A, B = z[:nc], z[nc] * sA, z[nc + 1] * sB
        status = dict(success=bool(result.success), status=int(result.status), message=str(result.message),
                      nfev=int(result.nfev), optimality=float(result.optimality))
        singular = np.linalg.svd(jacobian(z)[:len(p)][:, active], compute_uv=False)
        tolerance = (np.finfo(float).eps * max(len(p), len(active)) * singular[0]) if len(singular) else 0.
        z0 = np.r_[np.zeros(nc), A0 / sA, B0 / sB]
        status.update(coefficients=c.tolist(), active_columns=[int(k) for k in active if k < nc],
                      zero_response_columns=np.setdiff1d(np.arange(nc), active).tolist(),
                      data_singular_values=singular.tolist(), data_rank=int(np.count_nonzero(singular > tolerance)),
                      selected_counts=counts.tolist(), cost=float(np.sum(residual(z) ** 2) / 2),
                      initial_cost=float(np.sum(residual(z0) ** 2) / 2), grain_A_Pa=A, grain_B_Pa_per_s=B)
        g1 = grain(A, B, et)[0]
        prediction = surface + np.einsum("igk,k->ig", design, c)
        prediction[:, 0] += g1
        state.update(A_Pa=A, B_Pa_per_s=B, grain_new=g1, grain_old=g0, onset=onset,
                     selected_onset_ng_points=int((used[:, 0] & onset & (refs[:, 0] > 0)).sum()),
                     selected_counts=counts.tolist(), optimizer=status)
        return c, prediction, status
    ft.fit_corrections = fit_corrections_with_grain

    # ---- cpu_method_fit with the earlier-row check that allows the grain change (copy of the tool's function)
    def cpu_method_fit(ev, r, args, atom, groups, frozen, area, power, scale, baseline, refs, used, bins, late):
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
        dgrain = state["grain_new"] - state["grain_old"]
        early = ~late
        if not np.array_equal(prediction[early, 1], baseline[early, 1]) or \
                not np.allclose(prediction[early, 0] - baseline[early, 0], dgrain[early], rtol=1e-12, atol=1e-18) or \
                not np.array_equal(prediction[early & ~state["onset"], 0], baseline[early & ~state["onset"], 0]):
            raise SystemExit("earlier predictions changed other than by the grain")
        dn, dr = ev.point_changes(r, change)
        if np.any(dn[early] != 0) or np.any(dr[early] != 0):
            raise SystemExit("earlier predictions changed (GPU check)")
        check = np.column_stack([baseline[:, 0] + dn + dgrain, baseline[:, 1] + dr])
        result.update(coefficients=coefficients, prediction=prediction, optimizer=optimizer, change=change,
                      final=frozen + change,
                      check_rel=float(np.max(np.abs(check - prediction) / np.abs(prediction))))
        state["check_rel"] = result["check_rel"]
        # crossing check (not a penalty here): sign of ln(NG/RG) of the tuned model where both gauges are used
        both = used[:, 0] & used[:, 1]
        et = state["et"]
        wrong_after = both & (et >= t_end) & (prediction[:, 0] < prediction[:, 1])
        wrong_before = both & (et < t_end) & (prediction[:, 0] > prediction[:, 1])
        state["crossing_check"] = {"wrong_side_after_onset": int(wrong_after.sum()),
                                   "wrong_side_before_onset": int(wrong_before.sum())}
        state["earlier_ng_max_relative_change_from_grain"] = float(np.max(np.abs(dgrain[early]) / baseline[early, 0]))
        return result
    ftg_tool.cpu_method_fit = cpu_method_fit

    sys.argv = [str(APP / "tools" / "fine_tune_gpu.py")] + argv
    rc = ftg_tool.main()
    if out.exists():
        keep = {k: v for k, v in state.items() if not isinstance(v, np.ndarray)}
        (out / "grain_finetune.json").write_text(json.dumps(keep, indent=1, default=float))
        print("grain: A %.4e -> %.4e Pa, B %.4e -> %.4e Pa/s; NG points %d (+%d grain-onset); crossing %s" % (
            A0, state.get("A_Pa", float("nan")), B0, state.get("B_Pa_per_s", float("nan")),
            state.get("selected_counts", [0])[0], state.get("selected_onset_ng_points", 0),
            state.get("crossing_check")), flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
