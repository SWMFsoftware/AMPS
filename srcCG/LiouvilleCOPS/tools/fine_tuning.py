"""Late-first-response activity refinement; no new spherical harmonics.

Geometry preparation streams the existing kernels. Optimization acts only on
small response matrices, never on all facet/observation pairs. Structural zero
support, not a small-response threshold, protects every earlier cached row.
See FINE_TUNING_STEP.md for the CPU/GPU interchange contract and equations.
"""
import contextlib
import io

import numpy as np
from scipy.optimize import least_squares

import fit_cops
from cops_selection import parse_utc


def intervals(utc, minutes, count, end=None):
    """Half-open bins, last bin closed; duration weights within each bin.

    Each sample represents its nearest-time cell bounded by adjacent sample
    midpoints and the bin edges. A one-sample bin receives its full duration.
    These weights classify gauge dominance only, not the fitting objective.
    """
    t = np.array([parse_utc(str(u)) for u in utc], dtype=float)
    if len(t) < 2 or np.any(np.diff(t) <= 0):
        raise ValueError("base UTCs must be strictly increasing with at least two rows")
    if not np.isfinite(minutes) or minutes <= 0 or count < 1:
        raise ValueError("minutes and interval count must be positive")
    if count > len(t):
        raise ValueError("at least one interval has no modeled samples; reduce --intervals")
    end = t[-1] if end is None else parse_utc(end)
    if end < t[-1]:
        raise ValueError("end UTC precedes a modeled row; use a final-mission window")
    edges = np.linspace(end-60*minutes, end, count+1)
    if not t[0] < edges[0]:
        raise ValueError("window needs earlier modeled rows for the no-change audit")
    bins = np.searchsorted(edges, t, side="right")-1
    bins[t == edges[-1]] = count-1
    bins[t < edges[0]] = -1
    weights = np.zeros(len(t))
    for j in range(count):
        rows = np.flatnonzero(bins == j)
        if not len(rows):
            raise ValueError("interval %d contains no modeled samples; reduce --intervals" % j)
        times = t[rows]
        cell_edges = np.r_[edges[j], (times[:-1]+times[1:])/2, edges[j+1]]
        weights[rows] = np.diff(cell_edges)
    return bins, weights, edges


def row_response(run, header, cols, row, elements, coef, temperature):
    """Actual NG/RG pressure response per unit final H2O facet activity [Pa].

    Reuse the production Model's source, SZA, shadow, gauge, velocity and
    temperature factors. Constructing a one-row Model reads only that kernel;
    its global geometry arrays are references, not a full transport-cache copy.
    Unlike the activity-coefficient Jacobian, these weights have no derivative
    of max(0,L): a currently inactive facet can still acquire an enhancement.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        model = fit_cops.Model(str(run), header, cols, [row], elements, coef)
    tmax, t75, tmin = temperature
    tl = tmax + fit_cops.KT*(tmax-t75)*model.xlit
    te = np.where(model.lit & (tl > tmin), tl, tmin)
    if np.any(te <= 0) or not np.isfinite(te).all():
        raise ValueError("nonpositive/nonfinite frozen source temperature")
    inv = 1/np.sqrt(te)
    ng = model.cNG["H2O"]*model.bm["H2O"]*model.wn["H2O"]*inv
    rg = model.cRG["H2O"][0]*(model.wp["H2O"]+model.bm["H2O"]*model.wv["H2O"]*inv)
    if not np.isfinite(ng).all() or not np.isfinite(rg).all():
        raise ValueError("nonfinite per-unit-activity pressure response")
    # Nonnegative source additions must not yield negative gauge pressure.
    # Do not silently take abs/clamp in the forward response to hide a model
    # incompatibility. Absolute values below are only dominance diagnostics.
    if np.any(ng < 0) or np.any(rg < 0):
        raise ValueError("negative gauge response: review the base pressure formulation")
    structural = (model.Kn != 0) | (model.KPhi != 0) | (model.KV != 0)
    co2 = np.array([
        model.cNG["CO2"]*model.bm["CO2"]*np.sum(model.wn["CO2"]*inv),
        model.cRG["CO2"][0]*(model.PhiCO2[0]+model.bm["CO2"]*np.sum(model.wv["CO2"]*inv))])
    return model.face, structural, ng, rg, co2


def prepare_support(nfaces, bins, weights, response_rows, epsilon):
    """Partition first-seen facets into disjoint atomic (interval,type) groups.

    Atom type 0=NG dominated, 1=RG dominated, 2=mixed. -1 means ineligible.
    Structural first response is assigned across ALL rows, irrespective of
    observation-selection flags; facets first seen before the window stay out.
    """
    if not np.isfinite(epsilon) or not 0 < epsilon < 1:
        raise ValueError("dominance ratio must be strictly between zero and one")
    first = np.full(nfaces, -2, dtype=np.int32)  # -2 never seen; -1 seen before window
    ng_score, rg_score = np.zeros(nfaces), np.zeros(nfaces)
    for j, weight, response in zip(bins, weights, response_rows):
        face, structural, ng, rg, _ = response
        visible = face[structural]
        new = visible[first[visible] == -2]
        first[new] = j
        if j >= 0:
            eligible = first[face] == j
            np.add.at(ng_score, face[eligible], weight*np.abs(ng[eligible]))
            np.add.at(rg_score, face[eligible], weight*np.abs(rg[eligible]))
    has_signal = (ng_score > 0) | (rg_score > 0)
    kind = np.full(nfaces, 2, dtype=np.int8)
    kind[(ng_score > 0) & (rg_score <= epsilon*ng_score)] = 0
    kind[(rg_score > 0) & (ng_score <= epsilon*rg_score)] = 1
    atom = np.where((first >= 0) & has_signal, 3*first+kind, -1).astype(np.int32)
    return first, atom, ng_score, rg_score


def group_basis(atom, area, count, third="all"):
    """Small atom-to-coefficient matrix, never a dense faces-by-3n array.

    Empty columns are omitted. With overlapping 'all' masks and no mixed
    facets, all=NG+RG; remove all deterministically. For nonnegative amplitudes
    this leaves exactly the same attainable physical corrections. Further
    response degeneracy is diagnosed later, not silently interpreted as data.
    """
    if third not in ("all", "mixed"):
        raise ValueError("third group must be all or mixed")
    valid = atom >= 0
    sizes = np.bincount(atom[valid], minlength=3*count)
    areas = np.bincount(atom[valid], weights=area[valid], minlength=3*count)
    columns, names, omitted = [], [], []
    for j in range(count):
        for k, label in enumerate(("ng", "rg", third)):
            name = "late_%03d_%s" % (j, label)
            column = np.zeros(3*count)
            if k == 2 and third == "all":
                column[3*j:3*j+3] = 1
            else:
                column[3*j+k] = 1
            if not np.any(sizes[column != 0]):
                omitted.append({"name": name, "reason": "empty support"})
            elif k == 2 and third == "all" and sizes[3*j+2] == 0:
                omitted.append({"name": name, "reason": "all is a sum of existing nonempty NG/RG columns"})
            else:
                columns.append(column)
                names.append(name)
    return (np.column_stack(columns) if columns else np.empty((3*count, 0))), names, omitted, sizes, areas


def pressure_matrices(atom, n_atoms, activity, kappa, response_rows):
    """Stream baseline pressures and persistent atom responses over ALL rows.

    Membership is spatial, not a time switch. Every later occurrence of an
    assigned facet contributes to its column, including subsequent intervals.
    """
    baseline, matrices = [], []
    for face, _, ng, rg, co2 in response_rows:
        baseline.append([np.dot(ng, activity[face])+kappa*co2[0],
                         np.dot(rg, activity[face])+kappa*co2[1]])
        valid = atom[face] >= 0
        matrices.append(np.stack([
            np.bincount(atom[face[valid]], weights=w[valid], minlength=n_atoms)
            for w in (ng, rg)]))
    return np.asarray(baseline), np.asarray(matrices)


def fit_corrections(baseline, design, refs, used, late, gram, lam=0., max_nfev=1000):
    """Nonnegative additive activity amplitudes, equal-gauge log residuals.

    The field is A_final=max(0,L_fitted)+positive correction. Thus the reduced
    pressure model baseline+design*c is EXACT at fixed physics, including for
    initially inactive facets. We never optimize temperatures or old amplitudes.
    gram is the area-weighted TOTAL correction norm, including overlaps.
    """
    if not np.isfinite(lam) or lam < 0 or max_nfev < 1:
        raise ValueError("lambda must be finite/nonnegative and max-nfev positive")
    selected = used & late[:, None] & (refs > 0)
    if not np.isfinite(refs[selected]).all():
        raise ValueError("nonfinite selected reference pressure")
    counts = selected.sum(axis=0)
    if np.any(counts == 0):
        raise ValueError("late refinement requires selected positive NG and RG observations")
    row, gauge = np.nonzero(selected)
    d, p, target = design[row, gauge], baseline[row, gauge], refs[row, gauge]
    if np.any(p <= 0) or not np.isfinite(p).all():
        raise ValueError("selected baseline pressure must be finite and positive")
    weights = 1/np.sqrt(counts[gauge])
    eig, vec = np.linalg.eigh(gram)
    penalty = np.sqrt(lam)*np.sqrt(np.maximum(eig, 0))[:, None]*vec.T
    def residual(c):
        return np.r_[weights*np.log((p+d@c)/target), penalty@c]
    def jacobian(c):
        return np.vstack([weights[:, None]*d/(p+d@c)[:, None], penalty])
    # Exact zero-sensitivity columns stay at zero rather than acquiring a
    # meaningless fitted amplitude. Every nominal/support column is retained
    # in the report so this reduction remains visible to users and GPU ports.
    active = np.flatnonzero(np.any(d != 0, axis=0))
    x = np.zeros(design.shape[-1])
    if len(active):
        def expand(z):
            c = np.zeros_like(x); c[active] = z
            return c
        result = least_squares(lambda z: residual(expand(z)), x[active],
                               jac=lambda z: jacobian(expand(z))[:, active],
                               bounds=(0, np.inf), x_scale="jac", max_nfev=max_nfev,
                               ftol=1e-10, xtol=1e-10, gtol=1e-10)
        x[active] = result.x
        status = dict(success=bool(result.success), status=int(result.status),
                      message=str(result.message), nfev=int(result.nfev),
                      optimality=float(result.optimality))
    else:
        status = dict(success=True, status=0, message="no identifiable late-response columns", nfev=0)
    singular = np.linalg.svd(jacobian(x)[:len(p)], compute_uv=False)
    tolerance = (np.finfo(float).eps*max(d.shape)*singular[0]) if len(singular) else 0.
    status.update(coefficients=x.tolist(), active_columns=active.tolist(),
                  zero_response_columns=np.setdiff1d(np.arange(len(x)), active).tolist(),
                  data_singular_values=singular.tolist(),
                  data_rank=int(np.count_nonzero(singular > tolerance)),
                  selected_counts=counts.tolist(), cost=float(np.sum(residual(x)**2)/2),
                  initial_cost=float(np.sum(residual(np.zeros_like(x))**2)/2))
    return x, baseline+np.einsum("igk,k->ig", design, x), status
