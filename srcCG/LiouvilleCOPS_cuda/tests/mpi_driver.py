"""MPI test driver for tests/test_mpi_transport.py (run under mpirun, not collected by pytest).

    mpirun -n N python tests/mpi_driver.py CASE RUN OUT.npz ACTIVITY

Rank 0 builds the synthetic problem of the test fixtures (``test_backend_contract.cpu_model``,
``CONFIGS[-1]`` = NC31), evaluates through ``MultiGpuBackend(transport="mpi")`` and writes the
results to OUT.npz; ranks >= 1 serve (``transport.mpi_worker_main``) until rank 0 releases them.
``sequence`` is shared with the in-process reference of the tests.
"""

import json
import sys
from pathlib import Path

APP = Path(__file__).resolve().parents[1]
for path in (APP, APP / "tests"):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import numpy as np  # noqa: E402


def environment(run, config=-1):
    """The ``cpu_model`` fixture of test_backend_contract.py, with a GPU model."""
    from cops_cuda import reference
    import test_backend_contract as contract
    from test_gpu_model import make_gpu_model
    fit_cops = reference.load("fit_cops")
    run = str(run)
    header, names, rows = fit_cops.read_model_file(run + "/COPS_model.dat")
    cols = {name: index for index, name in enumerate(names)}
    elements = fit_cops.read_elements(run + "/elements.bin")
    rng = np.random.default_rng(3)
    coef = np.zeros((3, 25))
    coef[0] = rng.normal(0.0, 0.4, 25)
    coef[0][0] = 0.8
    coef[1:] = 0.5
    model = fit_cops.Model(run, header, cols, rows, elements, coef)
    model.load_kernels(model.sel)
    env = {"fit_cops": fit_cops, "model": model, "elements": elements, "coef": coef,
           "header": header, "synthetic": {"run": Path(run)}}
    G = make_gpu_model(env)
    P = contract.make_problem(dict(env, model=G), contract.CONFIGS[config])
    return env, G, P


def sequence(be, G, P):
    """Jacobians, pressures, a basis rebind, Jacobians again, seen_faces: [(name, array)]."""
    from cops_cuda import backend as bk
    import test_backend_contract as contract
    request = bk.ColumnRequest(P.fit_j, range(P.nb), temperature=True, kappa=True)
    states = [bk.ParameterState.from_problem(P, contract.perturbed(P, seed)) for seed in (1, 2, 3)]
    results = [be.evaluate(s, request) for s in states] + [be.evaluate(s) for s in states]
    basis2 = np.random.default_rng(9).normal(size=(len(G.Y), 3))
    be.rebind_basis(basis2, 0.7)
    request2 = bk.ColumnRequest(P.fit_j, range(3), temperature=True, kappa=True)
    results += [be.evaluate(bk.ParameterState(s.a, s.T, s.kappa, np.array([0.3, -0.2, 0.5])), request2)
                for s in states]
    out = []
    for i, r in enumerate(results):
        out += [("%d.pNG" % i, r.pNG), ("%d.pRG" % i, r.pRG)]
        for key in sorted(r.columns):
            out += [("%d.%s.%d.NG" % ((i,) + key), r.columns[key][0]),
                    ("%d.%s.%d.RG" % ((i,) + key), r.columns[key][1])]
    seen = be.seen_faces() if hasattr(be, "seen_faces") else None
    if seen is None:
        from cops_cuda import masks
        seen = masks.seen_positive(be.dataset)
    return out + [("seen", seen)]


def case_evaluate(run, activity):
    from cops_cuda.multi_gpu import MultiGpuBackend
    env, G, P = environment(run)
    with MultiGpuBackend(G._spec, G, P.enhancement_basis, P.enhancement_scale, transport="mpi",
                         activity=activity, records_per_block=17) as be:
        out = sequence(be, G, P)
        identity = be.identity()
    info = {"devices": identity["devices"], "transport": identity["transport"],
            "activity": identity["activity"], "ranks": identity["ranks"],
            "partition": identity["partition"],
            "worker_activity": [w["activity"]["mode"] for w in identity["workers"]],
            "streamed_records": [w["streamed_records"] for w in identity["workers_ready"]]}
    return dict(out, info=json.dumps(info))


def case_failure(run, activity):
    """A request failing on the workers is a BackendError on rank 0; the ranks serve the next bind."""
    from cops_cuda import backend as bk
    from cops_cuda.errors import BackendError
    from cops_cuda.multi_gpu import MultiGpuBackend
    import test_backend_contract as contract
    env, G, P = environment(run)
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 4))
    message = None
    be = MultiGpuBackend(G._spec, G, P.enhancement_basis, P.enhancement_scale, transport="mpi",
                         activity="device", records_per_block=17)
    try:
        be.evaluate(state)
        bad = np.full((len(G.Y), 2), "x")                 # every worker fails converting it
        try:
            be.rebind_basis(bad, 1.0)
        except BackendError as exc:
            message = str(exc)
        broken = be.broken
    finally:
        be.close()
    with MultiGpuBackend(G._spec, G, P.enhancement_basis, P.enhancement_scale, transport="mpi",
                         activity="device", records_per_block=17) as again:
        after = again.evaluate(state)
    return {"message": np.array(message or ""), "broken": np.array(broken),
            "pNG": after.pNG, "pRG": after.pRG}


def case_no_basis(run, activity):
    """B29: no template basis at bind (``basis=None`` travels to the workers)."""
    from cops_cuda import backend as bk
    from cops_cuda.multi_gpu import MultiGpuBackend
    import test_backend_contract as contract
    env, G, P = environment(run, 0)
    request = bk.ColumnRequest(P.fit_j, (), temperature=True, kappa=True)
    state = bk.ParameterState.from_problem(P, contract.perturbed(P, 5))
    with MultiGpuBackend(G._spec, G, None, 1.0, transport="mpi", activity=activity,
                         records_per_block=17) as be:
        r = be.evaluate(state, request)
    out = {"pNG": r.pNG, "pRG": r.pRG}
    for key in sorted(r.columns):
        out["%s.%d.NG" % key], out["%s.%d.RG" % key] = r.columns[key]
    return out


CASES = {"evaluate": case_evaluate, "failure": case_failure, "no_basis": case_no_basis}


def main():
    from mpi4py import MPI
    from cops_cuda import transport
    case, run, out, activity = sys.argv[1:5]
    comm = MPI.COMM_WORLD
    if comm.Get_rank() != 0:
        transport.mpi_worker_main(None, comm)
        return 0
    try:
        np.savez(out, **CASES[case](run, activity))
    finally:
        transport.mpi_release_workers(comm)
    return 0


if __name__ == "__main__":
    sys.exit(main())
