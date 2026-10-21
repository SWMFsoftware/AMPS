"""GPU late-interval fine tuning (cops_cuda/fine_tuning.py) against the CPU definitions.

* support: first structural row, first interval, first-interval dominance scores and groups equal
  the CPU tool's ``fine_tuning.prepare_support`` (exactly for the integers, to rounding for the
  scores); the group responses equal the CPU ``pressure_matrices`` atom responses;
* a two-part row split gives bitwise the one-part results;
* the multiplier/floor fit recovers known factors, a factor < 1 lowers the pressures, and the
  rows before the window are bitwise unchanged;
* the fitted final activity, evaluated on the GPU, reproduces the reduced-problem prediction.
"""

import numpy as np
import pytest

from cops_cuda import fine_tuning as ftg
from cops_cuda import reference
import test_restart
import test_backend_contract as contract

cpu_model = contract.cpu_model

cp = pytest.importorskip("cupy")
pytestmark = pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() == 0, reason="no CUDA device")

MINUTES, INTERVALS, EPS = 4.0, 2, 1e-3


def setup(cpu_model, tmp_path):
    _, _, r = test_restart.make_restart(cpu_model, tmp_path, test_restart.template_config())
    run = cpu_model["synthetic"]["run"]
    bins, weights, edges = ftg.cpu_intervals(r.utc, MINUTES, INTERVALS)
    return r, run, bins, weights, edges


def gpu_support(ev, r, bins, weights):
    p0 = ev.pressures(r)
    first = ftg.first_intervals(ev.first_structural_row(), bins)
    ng, rg = ev.first_interval_scores(r, first, bins, weights)
    return p0, first, ng, rg


def test_support_and_responses_equal_the_cpu(cpu_model, tmp_path, clean_pool):
    r, run, bins, weights, edges = setup(cpu_model, tmp_path)
    ft = reference.load("fine_tuning")
    fit_cops = cpu_model["fit_cops"]
    header, names, rows = fit_cops.read_model_file(str(run / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    responses = lambda: (ft.row_response(run, header, cols, row, cpu_model["elements"], r.coef_prior, r.T)
                         for row in rows)
    frozen = np.maximum(r.activity_linear, 0.0)
    c_first, c_atom, c_ng, c_rg = ft.prepare_support(len(frozen), bins, weights, responses(), EPS)
    c_base, c_resp = ft.pressure_matrices(c_atom, 3 * INTERVALS, frozen, float(r.kappa[0]), responses())

    ev = ftg.SplitEvaluator(run, [0])
    try:
        p0, first, ng, rg = gpu_support(ev, r, bins, weights)
        atom = ftg.classify(first, ng, rg, EPS)
        groups = ftg.group_names(atom, INTERVALS)
        WA, W1 = ftg.response_matrices(ev, r, atom, groups, frozen)
    finally:
        ev.release()
    assert np.array_equal(first, c_first)
    assert (first >= 0).sum() > 0 and (first == -1).sum() > 0          # both kinds present
    assert np.allclose(ng, c_ng, rtol=1e-12, atol=0) and np.allclose(rg, c_rg, rtol=1e-12, atol=0)
    assert np.array_equal(atom, c_atom)
    assert len(groups) >= 2
    for k, (code, _) in enumerate(groups):
        ref = c_resp[:, :, code]
        assert np.all(np.abs(W1[:, :, k] - ref) <= 1e-12 * np.maximum(np.abs(ref), 1e-300))
        member = atom == code
        assert np.all(W1[bins < code // 3, :, k] == 0.0) and np.all(WA[bins < code // 3, :, k] == 0.0)
        assert np.any(W1[:, :, k] > 0) == bool(member.any())
    base = np.column_stack(p0)
    assert np.all(np.abs(base - c_base) <= 1e-12 * np.abs(c_base))


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs two CUDA devices")
def test_split_is_bitwise_one_part(cpu_model, tmp_path, clean_pool):
    """One DeviceBudget owns each device's pool, so the parts must be on different devices."""
    r, run, bins, weights, edges = setup(cpu_model, tmp_path)
    frozen = np.maximum(r.activity_linear, 0.0)
    results = []
    for devices in ([0], [0, 1]):
        ev = ftg.SplitEvaluator(run, devices)
        try:
            p0, first, ng, rg = gpu_support(ev, r, bins, weights)
            atom = ftg.classify(first, ng, rg, EPS)
            WA, W1 = ftg.response_matrices(ev, r, atom, ftg.group_names(atom, INTERVALS), frozen)
            ident = ev.kernel_identity()
            results.append((p0[0], p0[1], first, ng, rg, WA, W1))
            assert ident["verified_against_run_manifest"]
            assert len(ev.parts) == len(devices)
        finally:
            ev.release()
    for a, b in zip(*results):
        assert a.tobytes() == b.tobytes()


def test_fit_recovers_factors_and_keeps_earlier_rows(cpu_model, tmp_path, clean_pool):
    r, run, bins, weights, edges = setup(cpu_model, tmp_path)
    frozen = np.maximum(r.activity_linear, 0.0)
    ev = ftg.SplitEvaluator(run, [0])
    try:
        p0, first, ng, rg = gpu_support(ev, r, bins, weights)
        atom = ftg.classify(first, ng, rg, EPS)
        groups = ftg.group_names(atom, INTERVALS)
        WA, W1 = ftg.response_matrices(ev, r, atom, groups, frozen)
        scale = r.enhancement_scale
        names, D, xref, kinds, index = ftg.parameters(groups, WA, W1, scale)
        baseline = np.column_stack(p0)
        late = bins >= 0
        used = np.ones_like(baseline, dtype=bool)
        # truth: every multiplier and floor of the identifiable columns perturbed
        x_true = xref.copy()
        rng = np.random.default_rng(4)
        for j, kind in enumerate(kinds):
            x_true[j] = rng.uniform(0.4, 0.8) if kind == "f" else rng.uniform(0.0, 0.3)
        refs = baseline + np.einsum("igk,k->ig", D, x_true - xref)
        power = float(np.dot(r.area, frozen ** 2))
        G = ftg.gram_matrix(groups, atom, np.asarray(r.area), frozen, scale, power, kinds, index)
        x, pred, status = ftg.fit(baseline, D, xref, refs, used, late, G, 0.0, 2000)
        assert pred[~late].tobytes() == baseline[~late].tobytes()
        active = status["active_columns"]
        assert status["cost"] < 1e-10 * status["initial_cost"]
        # f and c of a group are degenerate when WA is proportional to W1 (for example a
        # one-facet group); the data fit is then exact but the split between them is not.
        if status["data_rank"] == len(active):
            assert np.allclose(x[active], x_true[active], rtol=1e-4, atol=1e-7)
        assert np.allclose(pred[late], refs[late], rtol=1e-7, atol=0)
        # factors < 1 lower the activity of positive facets, so the late pressures can fall
        final = ftg.final_activity(frozen, atom, groups, x, kinds, index, scale)
        dn, dr = ev.point_changes(r, final - frozen)
        assert np.all(np.abs(baseline[:, 0] + dn - pred[:, 0]) <= 1e-12 * np.abs(pred[:, 0]))
        assert np.all(np.abs(baseline[:, 1] + dr - pred[:, 1]) <= 1e-12 * np.abs(pred[:, 1]))
        assert np.all(dn[~late] == 0.0) and np.all(dr[~late] == 0.0)
    finally:
        ev.release()


def test_classify_and_threshold_near():
    first = np.array([0, 0, 1, -1, 1, -2], dtype=np.int32)
    ng = np.array([1.0, 1e-3, 1.0, 1.0, 0.0, 0.0])
    rg = np.array([1e-3, 1.0, 1.0, 0.0, 0.0, 0.0])
    atom = ftg.classify(first, ng, rg, 1e-3)
    assert atom.tolist() == [0, 1, 5, -1, -1, -1]
    assert ftg.threshold_near(first, ng, rg, 1e-3).tolist() == [0, 1]


def test_confine_total_and_independent_gauges(cpu_model, tmp_path, clean_pool):
    """--confine total keeps the summed reference change within epsilon before the window, the
    constrained fit keeps the fitted change within epsilon, and a one-gauge fit ignores the
    other gauge's data entirely (pressures fitted on their own, never their ratio)."""
    r, run, bins, weights, edges = setup(cpu_model, tmp_path)
    frozen = np.maximum(r.activity_linear, 0.0)
    scale = r.enhancement_scale
    ev = ftg.SplitEvaluator(run, [0])
    try:
        baseline = np.column_stack(ev.pressures(r))
        ref = ftg.reference_change(frozen, scale)
        eps = 1e-2
        tau, eligible, record = ftg.confine_total(ev, r, ref, baseline, bins, eps)
        assert record["reference_outside_max_rel"] <= eps and eligible.any()
        dn, dr = ev.point_changes(r, np.where(eligible, ref, 0.0))
        before = bins < 0
        assert np.max(np.abs(dn[before]) / baseline[before, 0]) <= eps
        first = np.full(len(frozen), ftg.BEFORE, dtype=np.int32)
        sig = ev.first_significant_row(r, ref, 0.0, row_mask=bins >= 0)
        assert np.all(bins[sig[eligible]] >= 0)
        first[eligible] = bins[sig[eligible]]
        ng, rg = ev.first_interval_scores(r, first, bins, weights)
        atom = ftg.classify(first, ng, rg, EPS)
        groups = ftg.group_names(atom, INTERVALS)
        WA, W1 = ftg.response_matrices(ev, r, atom, groups, frozen)
    finally:
        ev.release()
    names, D, xref, kinds, index = ftg.parameters(groups, WA, W1, scale)
    late, used = bins >= 0, np.ones_like(baseline, dtype=bool)
    G = np.zeros((len(xref), len(xref)))
    refs = baseline * np.where(late[:, None], 3.0, 1.0)          # a late excess on both gauges
    x, pred, status = ftg.fit(baseline, D, xref, refs, used, late, G, 0.0, 500, confine=(~late, eps))
    assert status["outside_max_rel"] <= eps * (1 + 1e-6)
    assert np.max(np.abs(pred[~late] / baseline[~late] - 1)) <= eps * (1 + 1e-6)
    refs2 = refs.copy()
    refs2[late, 1] *= 0.1                                        # change only the RG data
    xa, _, _ = ftg.fit(baseline, D, xref, refs, used, late, G, 0.0, 500, gauges=(0,), confine=(~late, eps))
    xb, _, _ = ftg.fit(baseline, D, xref, refs2, used, late, G, 0.0, 500, gauges=(0,), confine=(~late, eps))
    assert xa.tobytes() == xb.tobytes()


@pytest.mark.parametrize("third", ["all", "mixed"])
def test_cpu_method_equals_the_cpu_pipeline(cpu_model, tmp_path, clean_pool, third):
    """tools/fine_tune_gpu.py --method cpu: the GPU-prepared atom responses with the CPU's own
    group_basis and fit_corrections reproduce the CPU pipeline (prepare_support, pressure_matrices,
    group_basis, fit_corrections) of ../LiouvilleCOPS tools/fine_tune_activity.py."""
    import argparse
    import fine_tune_gpu as tool
    r, run, bins, weights, edges = setup(cpu_model, tmp_path)
    ft = reference.load("fine_tuning")
    fit_cops = cpu_model["fit_cops"]
    header, names, rows = fit_cops.read_model_file(str(run / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    responses = lambda: (ft.row_response(run, header, cols, row, cpu_model["elements"], r.coef_prior, r.T)
                         for row in rows)
    frozen = np.maximum(r.activity_linear, 0.0)
    area = np.asarray(r.area, float)
    power = float(np.dot(area, frozen ** 2))
    scale = r.enhancement_scale
    late = bins >= 0
    used = np.column_stack([r.use_ng, r.use_rg]).astype(bool)
    refs = np.column_stack([r.ref_ng, r.ref_rg]).astype(float)
    refs[late] *= 1.7                                   # a late excess to fit
    # CPU pipeline
    c_first, c_atom, _, _ = ft.prepare_support(len(frozen), bins, weights, responses(), EPS)
    c_base, c_resp = ft.pressure_matrices(c_atom, 3 * INTERVALS, frozen, float(r.kappa[0]), responses())
    used &= c_base > 0                                  # the synthetic empty point has no pressure
    c_basis, c_labels, _, _, c_areas = ft.group_basis(c_atom, area, INTERVALS, third)
    c_design = scale * (c_resp @ c_basis)
    c_gram = scale ** 2 * (c_basis.T @ (c_areas[:, None] * c_basis)) / power
    c_x, c_pred, _ = ft.fit_corrections(c_base, c_design, refs, used, late, c_gram, 0.0, 200)
    # GPU method
    args = argparse.Namespace(intervals=INTERVALS, third_group=third, lam=0.0, max_nfev=200, prepare_only=False)
    ev = ftg.SplitEvaluator(run, [0])
    try:
        p0, first, ng, rg = gpu_support(ev, r, bins, weights)
        atom = ftg.classify(first, ng, rg, EPS)
        baseline = np.column_stack(p0)
        cpu = tool.cpu_method_fit(ev, r, args, atom, ftg.group_names(atom, INTERVALS), frozen, area, power, scale,
                                  baseline, refs, used, bins, late)
    finally:
        ev.release()
    assert np.array_equal(atom, c_atom) and cpu["labels"] == c_labels and np.array_equal(cpu["basis"], c_basis)
    assert np.all(np.abs(cpu["atom_response"] - c_resp) <= 1e-12 * np.maximum(np.abs(c_resp), 1e-300))
    assert len(c_labels) >= 2
    assert np.allclose(cpu["coefficients"], c_x, rtol=1e-6, atol=1e-9)
    assert np.allclose(cpu["prediction"], c_pred, rtol=1e-9, atol=0)
    assert cpu["prediction"][~late].tobytes() == baseline[~late].tobytes()


def test_fp32_records_match_fp64(cpu_model, tmp_path, clean_pool, monkeypatch):
    """C14 experiment: the restart.cu ``_m`` kernels on FP32 records (factors recomputed in FP64)
    reproduce the FP64 kernels: integers (first rows) exactly, since FP32 keeps zero-ness; sums
    and shares within 1e-6 relative (inputs rounded by <= 2^-24)."""
    from cops_cuda import precision
    r, run, bins, weights, edges = setup(cpu_model, tmp_path)
    frozen = np.maximum(r.activity_linear, 0.0)
    delta = np.where(np.arange(len(frozen)) % 3 == 0, 0.5 * frozen + 1e-3, 0.0)
    late = bins >= 0
    out = {}
    for rp in ("fp64", "fp32"):
        monkeypatch.setenv(precision.RECORD_PRECISION_ENVIRONMENT, rp)
        ev = ftg.SplitEvaluator(run, [0])
        try:
            p0, first, ng, rg = gpu_support(ev, r, bins, weights)
            assert ev.parts[0]._backend.mixed == (rp == "fp32")
            out[rp] = dict(p0=np.column_stack(p0), first=first, ng=ng, rg=rg,
                           sig=ev.first_significant_row(r, delta, 0.0, late),
                           sig3=ev.first_significant_row(r, delta, 1e-3),
                           share=ev.max_shares(r, delta, ~late),
                           change=np.column_stack(ev.point_changes(r, delta)))
        finally:
            ev.release()
    a, b = out["fp64"], out["fp32"]
    close = lambda x, y: np.all(np.abs(x - y) <= 1e-6 * np.abs(x))
    assert np.array_equal(a["first"], b["first"]) and np.array_equal(a["sig"], b["sig"])
    assert (a["sig"] >= 0).any() and (a["sig3"] >= 0).any()
    assert close(a["p0"], b["p0"]) and close(a["ng"], b["ng"]) and close(a["rg"], b["rg"])
    assert close(a["share"], b["share"]) and close(a["change"], b["change"])
    # a threshold can only move where a share lies within rounding of tau
    moved = a["sig3"] != b["sig3"]
    assert moved.sum() <= max(1, int(0.001 * len(moved)))
