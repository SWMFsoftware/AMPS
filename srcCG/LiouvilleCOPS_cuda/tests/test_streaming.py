"""Streamed records (SCALING_PLAN.md P2, invariant I1): streamed = resident, BITWISE.

Records beyond a resident prefix stay in pinned host memory and are copied through device
slots each pass, overlapped with the kernels. The kernels see block for block the same records,
so the per-block partials and their block-order reduction are those of the resident path:
* pressures, every Jacobian column (with templates, after a basis rebind) and PhiCO2 are
  bitwise equal for resident caps from 0 to almost everything, with many small chunks (every
  slot reused), with and without face compaction, and from a record store;
* ``seen_faces`` and the device+host record digests are those of the resident load;
* two GPUs streaming equal one resident GPU;
* the restart/fine-tuning kernels run over ``record_chunks`` (P2b) and equal the resident ones
  bitwise; the chunks cover every block once, in order.
"""

import numpy as np
import pytest

from cops_cuda import backend as bk
from cops_cuda import host_model
import test_backend_contract as contract

cp = pytest.importorskip("cupy")
if cp.cuda.runtime.getDeviceCount() == 0:
    pytest.skip("no CUDA device", allow_module_level=True)
from cops_cuda import evaluator, masks  # noqa: E402

cpu_model = contract.cpu_model


def same(u, v):
    return np.asarray(u).tobytes() == np.asarray(v).tobytes()


def full_model(env):
    fit_cops = env["fit_cops"]
    run = str(env["synthetic"]["run"])
    header, names, table = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {n: i for i, n in enumerate(names)}
    M = host_model.header_model(run, header, cols, table, env["elements"], env["coef"])
    M.load_kernels(M.sel)
    return M


def evaluate_all(be, states, request, basis):
    res = [be.evaluate(bk.ParameterState(s.a, s.T, s.kappa)) for s in states]
    be.rebind_basis(basis, 1.3)
    res += [be.evaluate(s, request) for s in states]
    return res, be.PhiCO2, masks.seen_positive(be.dataset), be.dataset.device_digests()


def check_same(ref, got):
    (ra, rphi, rseen, rdig), (ga, gphi, gseen, gdig) = ref, got
    assert same(rphi, gphi) and same(rseen, gseen)
    assert {k: v for k, v in rdig.items() if k.startswith("kernel.")} == \
        {k: v for k, v in gdig.items() if k.startswith("kernel.")}
    for x, y in zip(ra, ga):
        assert same(x.pNG, y.pNG) and same(x.pRG, y.pRG)
        for key in x.columns:
            assert same(x.columns[key][0], y.columns[key][0]), key
            assert same(x.columns[key][1], y.columns[key][1]), key


def setup(env):
    M = full_model(env)
    nf = len(env["elements"]["area"])
    basis = np.random.default_rng(5).normal(size=(nf, 2))
    request = bk.ColumnRequest(range(25), range(2), temperature=True, kappa=True)
    states = [bk.ParameterState(env["coef"][0] * s, M.T0, k, b)
              for s, k, b in ((1.0, 1.0, [0.2, -0.1]), (0.6, 0.3, [1.0, 0.4]))]
    return M, basis, request, states


@pytest.mark.parametrize("compact", [True, False])
def test_streamed_is_bitwise_resident(cpu_model, clean_pool, monkeypatch, compact):
    M, basis, request, states = setup(cpu_model)
    run = cpu_model["synthetic"]["run"]
    nrec = int(M.indptr[-1])
    kw = dict(mode="tree", records_per_block=17, record_precision="fp32", compact_faces=compact)
    with evaluator.CudaBackend(M, run, None, 1.0, **kw) as be:
        assert not be.dataset.streamed
        ref = evaluate_all(be, states, request, basis)
    monkeypatch.setenv(evaluator.STREAM_CHUNK_ENVIRONMENT, "40")
    for cap in (0, 17, 300, nrec // 2, nrec - 1):
        monkeypatch.setenv(evaluator.RESIDENT_ENVIRONMENT, str(cap))
        with evaluator.CudaBackend(M, run, None, 1.0, **kw) as be:
            ds = be.dataset
            assert ds.streamed and ds.nresident <= cap and ds.nresident % 1 == 0
            assert len(be._stream["chunks"]) >= (evaluator.STREAM_SLOTS + 2 if cap == 0 else 1)
            got = evaluate_all(be, states, request, basis)
            assert ds.audit()["face_out_of_range"] == 0
        check_same(ref, got)


def test_streamed_from_record_store(cpu_model, clean_pool, monkeypatch, tmp_path):
    from cops_cuda import record_store
    M, basis, request, states = setup(cpu_model)
    run = cpu_model["synthetic"]["run"]
    record_store.build(run, tmp_path / "store")
    kw = dict(mode="tree", records_per_block=17, record_precision="fp32")
    with evaluator.CudaBackend(M, run, None, 1.0, **kw) as be:
        ref = evaluate_all(be, states, request, basis)
    monkeypatch.setenv(evaluator.STREAM_CHUNK_ENVIRONMENT, "64")
    monkeypatch.setenv(evaluator.RESIDENT_ENVIRONMENT, "100")
    with evaluator.CudaBackend(M, run, None, 1.0, record_store=tmp_path / "store", **kw) as be:
        assert be.dataset.streamed and be.dataset.record_source == "record_store"
        got = evaluate_all(be, states, request, basis)
    check_same(ref, got)


@pytest.mark.skipif(cp.cuda.runtime.getDeviceCount() < 2, reason="needs two CUDA devices")
def test_two_gpus_streaming_equal_one_resident(cpu_model, clean_pool, monkeypatch):
    from cops_cuda.multi_gpu import MultiGpuBackend
    from test_gpu_model import make_gpu_model
    monkeypatch.setenv("COPS_CUDA_RECORD_PRECISION", "fp32")
    G = make_gpu_model(cpu_model)
    run = cpu_model["synthetic"]["run"]
    P = contract.make_problem(dict(cpu_model, model=G), contract.CONFIGS[-1])
    request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 8))
    with evaluator.CudaBackend(G, run, P.enhancement_basis, P.enhancement_scale,
                               records_per_block=17) as one:
        want = one.evaluate(state, request)
    monkeypatch.setenv(evaluator.STREAM_CHUNK_ENVIRONMENT, "50")
    monkeypatch.setenv(evaluator.RESIDENT_ENVIRONMENT, "60")
    with MultiGpuBackend(G._spec, G, P.enhancement_basis, P.enhancement_scale, devices=(0, 1),
                         records_per_block=17) as two:
        got = two.evaluate(state, request)
    assert same(got.pNG, want.pNG) and same(got.pRG, want.pRG)
    for key in request.keys():
        assert same(got.columns[key][0], want.columns[key][0]) and same(got.columns[key][1], want.columns[key][1])


@pytest.mark.parametrize("cap,chunk", [(0, 40), (50, 40), (130, 1000)])
def test_restart_and_fine_tuning_kernels_stream_bitwise(cpu_model, clean_pool, monkeypatch, tmp_path, cap, chunk):
    """P2b: the restart.cu kernels (facet shares, point changes, first structural / significant row,
    first-interval scores) over ``CudaBackend.record_chunks`` equal the resident ones bitwise,
    including points split over chunks and over the resident/streamed boundary."""
    from cops_cuda import fine_tuning as ftg
    import test_fine_tuning_gpu as tft
    r, run, bins, weights, edges = tft.setup(cpu_model, tmp_path)
    frozen = np.maximum(r.activity_linear, 0.0)
    delta = np.where(np.arange(len(frozen)) % 3 == 0, 0.5 * frozen + 1e-3, 0.0)
    late = bins >= 0
    from cops_cuda import restart as rs
    monkeypatch.setenv("COPS_CUDA_RECORD_PRECISION", "fp32")
    monkeypatch.setattr(rs, "RECORDS_PER_BLOCK", 17)
    out = {}
    for streamed in (False, True):
        if streamed:
            monkeypatch.setenv(evaluator.RESIDENT_ENVIRONMENT, str(cap))
            monkeypatch.setenv(evaluator.STREAM_CHUNK_ENVIRONMENT, str(chunk))
        ev = ftg.SplitEvaluator(run, [0])
        try:
            p0, first, ng, rg = tft.gpu_support(ev, r, bins, weights)
            be = ev.parts[0]._backend
            assert be.dataset.streamed == streamed and be.scratch_per_face > 0
            if streamed:
                assert len(be._stream["chunks"]) >= (2 if cap < 130 else 1)
            sn, sr, _, _ = ev.parts[0].facet_point_shares(r, delta, ev.point_ids[~late])
            out[streamed] = dict(p0=np.column_stack(p0), first=first, ng=ng, rg=rg, sn=sn, sr=sr,
                                 sig=ev.first_significant_row(r, delta, 0.0, late),
                                 sig3=ev.first_significant_row(r, delta, 1e-3),
                                 change=np.column_stack(ev.point_changes(r, delta)))
        finally:
            ev.release()
    for key in out[False]:
        assert same(out[False][key], out[True][key]), key
    assert (out[True]["sig"] >= 0).any() and np.any(out[True]["ng"] > 0)


def test_streamed_record_chunks_cover_every_block_once(cpu_model, clean_pool, monkeypatch):
    M, basis, request, states = setup(cpu_model)
    run = cpu_model["synthetic"]["run"]
    monkeypatch.setenv(evaluator.RESIDENT_ENVIRONMENT, "70")
    monkeypatch.setenv(evaluator.STREAM_CHUNK_ENVIRONMENT, "40")
    with evaluator.CudaBackend(M, run, None, 1.0, records_per_block=17, record_precision="fp32") as be:
        table = be.arrays["blocks"].get().reshape(-1, 3)
        seen, faces = [], []
        for g0, g1, t, rec in be.record_chunks():
            t = t.get().reshape(-1, 3)
            assert len(t) == g1 - g0 and np.array_equal(t[:, 0], table[g0:g1, 0])
            assert np.array_equal(t[:, 2] - t[:, 1], table[g0:g1, 2] - table[g0:g1, 1])
            f = rec["face"].get()
            faces += [f[b:e] for _, b, e in t]
            seen += list(range(g0, g1))
        assert seen == list(range(len(table)))
        want = np.concatenate([be.dataset.arrays["kernel.face"].get()[:be.dataset.nresident],
                               be.dataset.host_records["face"]])
        assert np.array_equal(np.concatenate(faces), want)
