# SPDX-License-Identifier: GPL-3.0-only
#
# Parts of this module (``face_table``, the final scaling in ``_finish`` and
# the constants) are translated from polyhedral-gravity-model 3.3,
# Copyright (C) Jonas Schuhmacher and contributors,
# https://github.com/esa/polyhedral-gravity-model, GNU GPL version 3
# (license text: third_party/src_pypi/polyhedral_gravity-3.3/LICENSE).
# Derived from GravityEvaluable.cpp (prepare, evaluate), GravityModelDetail.cpp
# (buildVectorsOfSegments, buildUnitNormalOfPlane, buildUnitNormalOfSegments),
# UtilityContainer.h (cross, euclideanNorm, normal) and UtilityConstants.h.
# Modified 2026-10-19 for LiouvilleCOPS_cuda C03 (NumPy vectorization, device
# driver).  This module is distributed under the GPL v3 like its origin.
"""Exact FP64 polyhedral gravity on one GPU (ROADMAP C03; CUDA environment).

The device code is in ``kernels/gravity.cu``. This module prepares the face
table, plans and enforces device memory, launches the kernels, and applies the
package's final scaling. Usage::

    body = GravityBody.from_bdf(bdf_path, density=533.0)   # audited, no HEAL
    with GravityEvaluator(body, device_id=0) as ev:
        potential, acceleration = ev.evaluate(points)       # (n,), (n, 3), SI

Precomputed face terms (``face_table``)
---------------------------------------
``GravityEvaluable::prepare`` computes, per face and once, the segment
vectors G_q = v_{q+1} - v_q, the plane unit normal N = (G_0 x G_1)/|G_0 x G_1|
and the segment unit normals n_q = (G_q x N)/|G_q x N| from the original
vertices. ``face_table`` evaluates the same expressions in NumPy, with the
operation order of ``util::cross`` and ``util::euclideanNorm``. NumPy rounds
every elementwise operation separately (no fused multiply-add), so the table
is bitwise identical to the package's own. ``tests/test_gravity.py`` checks
this against ``GravityEvaluable.__getstate__``. The three original vertices
are stored with each face, so the kernel needs no indexed gathers. A table row
is 30 doubles (240 bytes): 189 MB for the 788k body.

Evaluation variants
-------------------
``fold``
    One thread per target, faces in index order: the package's summation
    order, so results differ from the package only through libdevice versus
    glibc/xsimd transcendental functions. Targets are processed in chunks.
    Each launch covers a bounded face range and the running sums carry over
    exactly, so any chunking gives bitwise identical results (tested).
``split``
    Faces split into ``chunks`` with partial sums reduced in a fixed chunk
    order. It is deterministic, but its rounding differs from the package
    (bounded by the C00 budget). It keeps the device busy when there are few
    targets (tests, audits); for the full 788k target set, ``fold`` already
    has about 6000 blocks.

Units and scaling
-----------------
Vertices and points are in metres. The package's METER scaling is
``prefix = G * density * orientation`` with G = 6.67430e-11 m^3 kg^-1 s^-2
and orientation +1 (OUTWARDS). It returns ``V = (sum*prefix)/2`` and
``g = sum * (-1.0*prefix)``; the acceleration points into the body. The
same expressions are evaluated here, on the host, in NumPy.
"""

import time

import numpy as np

from . import cuda_sources, precision, provenance
from .errors import CapacityError
from .memory import DeviceBudget, MemoryPlan, human

GRAVITATIONAL_CONSTANT = 6.67430e-11       # util::GRAVITATIONAL_CONSTANT, m^3 kg^-1 s^-2
ORIENTATION_FACTOR = 1.0                   # NormalOrientation::OUTWARDS
FACE_STRIDE = 30
BLOCK = 128                                # must equal GRAVITY_BLOCK in gravity.cu
NSUMS = 5                                  # potential, ax, ay, az, M (gravity.cu)
VARIANTS = ("fold", "split")
PRECISIONS = ("fp64", "mixed")
FACE32_STRIDE = 20                         # gravity.cu: c q0 q1 q2 N wA rnear2 pad pad pad
MIXED_POLICY = "mixed_near_fp64_exact_far_fp32_quad3_v1"
MIXED_TILE = 64                            # must equal MIXED_TILE in gravity.cu
DEFAULT_NEAR_FACTOR = 6.0


def _cross(a, b):
    """util::cross, component by component."""
    return np.stack((a[:, 1] * b[:, 2] - a[:, 2] * b[:, 1],
                     a[:, 2] * b[:, 0] - a[:, 0] * b[:, 2],
                     a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]), axis=1)


def _norm(a):
    """util::euclideanNorm: sqrt((x*x + y*y) + z*z)."""
    return np.sqrt((a[:, 0] * a[:, 0] + a[:, 1] * a[:, 1]) + a[:, 2] * a[:, 2])


def _normal(a, b):
    """util::normal: cross(a, b) / |cross(a, b)|, element-wise division."""
    c = _cross(a, b)
    return c / _norm(c)[:, None]


def face_terms(vertices, faces):
    """(G[nf,3,3], N[nf,3], n[nf,3,3]) exactly as GravityEvaluable::prepare."""
    v0, v1, v2 = vertices[faces[:, 0]], vertices[faces[:, 1]], vertices[faces[:, 2]]
    G = np.stack((v1 - v0, v2 - v1, v0 - v2), axis=1)          # buildVectorsOfSegments
    N = _normal(G[:, 0], G[:, 1])                               # buildUnitNormalOfPlane
    n = np.stack([_normal(G[:, q], N) for q in range(3)], axis=1)   # buildUnitNormalOfSegments
    return G, N, n


def face_table(body):
    """The (nfaces, 30) float64 device table: v0 v1 v2 G0 G1 G2 N n0 n1 n2."""
    V, F = body.vertices, body.faces
    G, N, n = face_terms(V, F)
    table = np.empty((len(F), FACE_STRIDE), dtype=np.float64)
    table[:, 0:3], table[:, 3:6], table[:, 6:9] = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    table[:, 9:18] = G.reshape(len(F), 9)
    table[:, 18:21] = N
    table[:, 21:30] = n.reshape(len(F), 9)
    if not np.isfinite(table).all():
        raise ValueError("non-finite precomputed face terms (degenerate face?)")
    return table


def face_table32(body, near_factor=DEFAULT_NEAR_FACTOR):
    """The (nfaces, 20) float32 far-field table of ``gravity_mixed_fold`` (mixed precision).

    Per face: centroid c, the 3-point quadrature points q_k (barycentric 2/3, 1/6, 1/6),
    the unit normal N of ``face_terms``, the weight A/3 and the near radius squared
    (near_factor * circumradius)^2. Computed in FP64, stored in FP32."""
    if not near_factor > 0:
        raise ValueError("near_factor must be positive")
    V, F = body.vertices, body.faces
    v0, v1, v2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    _, N, _ = face_terms(V, F)
    c = (v0 + v1 + v2) / 3.0
    area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0), axis=1)
    radius = np.max(np.stack([np.linalg.norm(v - c, axis=1) for v in (v0, v1, v2)]), axis=0)
    t = np.zeros((len(F), FACE32_STRIDE))
    t[:, 0:3] = c
    for k, (a, b, d) in enumerate(((v0, v1, v2), (v1, v2, v0), (v2, v0, v1))):
        t[:, 3 + 3 * k:6 + 3 * k] = (2.0 / 3.0) * a + (1.0 / 6.0) * b + (1.0 / 6.0) * d
    t[:, 12:15] = N
    t[:, 15] = area / 3.0
    t[:, 16] = (near_factor * radius) ** 2
    return np.ascontiguousarray(t, dtype=np.float32)


def scaling_prefix(density):
    """Polyhedron::getGravityModelScaling for METER: G * density * orientation."""
    return GRAVITATIONAL_CONSTANT * float(density) * ORIENTATION_FACTOR


class GravityEvaluator:
    """A gravity body resident on one device, evaluating batches of targets.

    Parameters
    ----------
    target_chunk : targets per fold launch group. Each chunk's coordinates
        and running sums are device scratch; 65536 targets use 4.2 MB.
    max_pairs_per_launch : bound on face-target pairs in one kernel launch.
        The default of 1e8 is about 0.5 s on an A6000. Both GPUs of this host
        report the kernel watchdog enabled (``KernelExecTimeout``), which kills
        launches lasting several seconds. Short launches also keep evaluation
        interruptible. In ``fold`` the running sums carry between launches,
        so the bound never changes a result.
    split_chunks, split_target_chunk : the face-split variant's chunk count
        and the target chunk, which together bound its partial-sum buffer.
    """

    def __init__(self, body, device_id=0, budget=None, target_chunk=1 << 16,
                 max_pairs_per_launch=None, split_chunks=128, split_target_chunk=1 << 13,
                 precision="fp64", near_factor=DEFAULT_NEAR_FACTOR):
        """``precision="mixed"`` (experiment): exact FP64 terms of the faces within
        ``near_factor`` circumradii of a target, FP32 quadrature for all others
        (``kernels/gravity.cu: gravity_mixed_fold``); "fp64" is the reference."""
        import cupy
        if precision not in PRECISIONS:
            raise ValueError("precision must be one of %s" % (PRECISIONS,))
        self.cupy = cupy
        self.body = body
        self.precision = precision
        self.near_factor = float(near_factor)
        self.device_id = int(device_id)
        self.target_chunk = int(target_chunk)
        if max_pairs_per_launch is None:
            max_pairs_per_launch = 10 ** 8 if precision == "fp64" else 10 ** 9
        self.max_pairs_per_launch = int(max_pairs_per_launch)
        if self.max_pairs_per_launch < 1:
            raise ValueError("max_pairs_per_launch must be positive")
        self.split_chunks = int(min(split_chunks, max(1, body.nfaces)))
        self.split_target_chunk = int(split_target_chunk)
        self.prefix = scaling_prefix(body.density)
        self.released = False
        t0 = time.perf_counter()
        table = face_table(body)
        self.table_sha256 = provenance.array_sha256(table)
        table32 = face_table32(body, self.near_factor) if precision == "mixed" else None
        self.table32_sha256 = provenance.array_sha256(table32) if table32 is not None else None
        self.timing = {"precompute_s": time.perf_counter() - t0}
        self.near_counts = None

        nf = body.nfaces
        self.plan = MemoryPlan("gravity evaluator")
        self.plan.add_array("gravity.faces", "static", (nf, FACE_STRIDE), np.float64)
        self.plan.add_array("gravity.targets", "scratch", (3, self.target_chunk), np.float64)
        self.plan.add_array("gravity.sums", "scratch", (NSUMS, self.target_chunk), np.float64)
        self.plan.add_array("gravity.split_partial", "partial_sums",
                            (self.split_chunks, NSUMS, self.split_target_chunk), np.float64)
        self.plan.add_array("gravity.split_targets", "scratch",
                            (3, self.split_target_chunk), np.float64)
        self.plan.add_array("gravity.split_sums", "scratch", (NSUMS, self.split_target_chunk),
                            np.float64)
        if precision == "mixed":
            self.plan.add_array("gravity.faces32", "static", (nf, FACE32_STRIDE), np.float32)
            self.plan.add_array("gravity.near_count", "scratch", (self.target_chunk,), np.uint32)
        self.owns_budget = budget is None
        self.budget = DeviceBudget(self.device_id) if budget is None else budget
        self.reservation = None
        self.arrays = {}
        try:
            self.reservation = self.budget.reserve(self.plan)
            with cupy.cuda.Device(self.device_id), self.budget.allocation_guard(self.plan):
                t0 = time.perf_counter()
                self.arrays["faces"] = cupy.asarray(table)
                # Flat scratch buffers: the first k*m entries viewed as (k, m)
                # are C-contiguous, which is the [k][nt] layout the kernels index.
                self.arrays["targets"] = cupy.empty(3 * self.target_chunk, dtype=np.float64)
                self.arrays["sums"] = cupy.empty(NSUMS * self.target_chunk, dtype=np.float64)
                self.arrays["split_partial"] = cupy.empty(
                    self.split_chunks * NSUMS * self.split_target_chunk, dtype=np.float64)
                self.arrays["split_targets"] = cupy.empty(3 * self.split_target_chunk,
                                                          dtype=np.float64)
                self.arrays["split_sums"] = cupy.empty(NSUMS * self.split_target_chunk,
                                                       dtype=np.float64)
                self.fold = cuda_sources.function("gravity", "gravity_fold", self.device_id)
                self.partial = cuda_sources.function("gravity", "gravity_partial", self.device_id)
                self.reduce = cuda_sources.function("gravity", "gravity_reduce", self.device_id)
                if precision == "mixed":
                    self.arrays["faces32"] = cupy.asarray(table32)
                    self.arrays["near_count"] = cupy.empty(self.target_chunk, dtype=np.uint32)
                    self.mixed = cuda_sources.function("gravity", "gravity_mixed_fold", self.device_id)
                self.timing["upload_and_compile_s"] = time.perf_counter() - t0
        except BaseException:
            self.release()
            raise

    # ----------------------------------------------------------------- evaluate
    def evaluate(self, points, variant="fold", progress=None, with_magnitude=False):
        """Potential [m^2/s^2] and acceleration [m/s^2] at ``points`` (n, 3).

        With ``with_magnitude`` also return M [m/s^2]: the scaled sum, over all
        faces, of the magnitudes of the terms that carry a transcendental result
        (see ``kernels/gravity.cu``). ``|g_gpu - g_ref| / M`` is the error
        measure of the C03 budget, and it is independent of cancellation.
        ``progress(done, total, seconds)`` is called after each target chunk.
        """
        if self.released:
            raise RuntimeError("evaluator has been released")
        if variant not in VARIANTS:
            raise ValueError("variant must be one of %s" % (VARIANTS,))
        points = np.ascontiguousarray(points, dtype=np.float64)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be a finite (n, 3) array")
        n = len(points)
        sums = np.empty((NSUMS, n), dtype=np.float64)
        if self.precision == "mixed":
            if variant != "fold":
                raise ValueError("mixed precision supports the fold order only")
            self.near_counts = np.empty(n, dtype=np.int64)
        chunk = self.target_chunk if variant == "fold" else self.split_target_chunk
        started = time.perf_counter()
        cupy = self.cupy
        with cupy.cuda.Device(self.device_id):
            for begin in range(0, n, chunk):
                end = min(begin + chunk, n)
                if self.precision == "mixed":
                    sums[:, begin:end], self.near_counts[begin:end] = self._mixed_chunk(points[begin:end])
                elif variant == "fold":
                    sums[:, begin:end] = self._fold_chunk(points[begin:end])
                else:
                    sums[:, begin:end] = self._split_chunk(points[begin:end])
                if progress is not None:
                    progress(end, n, time.perf_counter() - started)
        self.timing["last_evaluate_s"] = time.perf_counter() - started
        self.timing["last_evaluate_targets"] = n
        potential, acceleration = self._finish(sums)
        if with_magnitude:
            return potential, acceleration, sums[4] * abs(self.prefix)
        return potential, acceleration

    def _finish(self, sums):
        # GravityEvaluable::evaluate step 10, the same IEEE operations:
        # potential = (potential * prefix) / 2.0; acceleration *= (-1.0 * prefix)
        potential = (sums[0] * self.prefix) / 2.0
        acceleration = np.ascontiguousarray(sums[1:4].T) * (-1.0 * self.prefix)
        return potential, acceleration

    def _fold_chunk(self, points):
        cupy, m, nf = self.cupy, len(points), self.body.nfaces
        t_dev = self.arrays["targets"][:3 * m].reshape(3, m)
        t_dev.set(np.ascontiguousarray(points.T))
        s_dev = self.arrays["sums"][:NSUMS * m].reshape(NSUMS, m)
        s_dev.fill(0.0)                 # the package's result{} initial value
        grid = (-(-m // BLOCK),)
        faces_per_launch = max(1, min(nf, self.max_pairs_per_launch // m))
        self.last_faces_per_launch = faces_per_launch
        for face_begin in range(0, nf, faces_per_launch):
            face_end = min(face_begin + faces_per_launch, nf)
            self.fold(grid, (BLOCK,), (self.arrays["faces"], np.int64(face_begin),
                                       np.int64(face_end), t_dev, np.int64(m), s_dev))
        return s_dev.get()

    def _mixed_chunk(self, points):
        m, nf = len(points), self.body.nfaces
        t_dev = self.arrays["targets"][:3 * m].reshape(3, m)
        t_dev.set(np.ascontiguousarray(points.T))
        s_dev = self.arrays["sums"][:NSUMS * m].reshape(NSUMS, m)
        s_dev.fill(0.0)
        near = self.arrays["near_count"][:m]
        near.fill(0)
        grid = (-(-m // BLOCK),)
        # Whole tiles per launch: the FP32 tile sums then never depend on the launch split.
        faces_per_launch = max(MIXED_TILE, (self.max_pairs_per_launch // m) // MIXED_TILE * MIXED_TILE)
        self.last_faces_per_launch = faces_per_launch
        for face_begin in range(0, nf, faces_per_launch):
            face_end = min(face_begin + faces_per_launch, nf)
            self.mixed(grid, (BLOCK,), (self.arrays["faces"], self.arrays["faces32"], np.int64(face_begin),
                                        np.int64(face_end), t_dev, np.int64(m), s_dev, near))
        return s_dev.get(), near.get().astype(np.int64)

    def _split_chunk(self, points):
        cupy, m, nf = self.cupy, len(points), self.body.nfaces
        chunk_faces = -(-nf // self.split_chunks)
        nchunks = -(-nf // chunk_faces)
        t_dev = self.arrays["split_targets"][:3 * m].reshape(3, m)
        t_dev.set(np.ascontiguousarray(points.T))
        partial = self.arrays["split_partial"][:nchunks * NSUMS * m]
        # Several launches of whole chunks, each within the pair bound.
        per_launch = max(1, min(nchunks, self.max_pairs_per_launch // max(1, m * chunk_faces)))
        for first in range(0, nchunks, per_launch):
            count = min(per_launch, nchunks - first)
            self.partial((-(-m // BLOCK), count), (BLOCK,),
                         (self.arrays["faces"], np.int64(nf), np.int64(chunk_faces),
                          np.int64(first), t_dev, np.int64(m), partial))
        out = self.arrays["split_sums"][:NSUMS * m]
        self.reduce((-(-NSUMS * m // 256),), (256,), (partial, np.int64(nchunks), np.int64(m), out))
        self.last_split_chunk_faces = chunk_faces
        return out.reshape(NSUMS, m).get()

    # ------------------------------------------------------------- bookkeeping
    def identity(self):
        return {"backend": "cuda_polyhedral_gravity_v1",
                "translated_from": "polyhedral-gravity-model 3.3 (GPL-3.0)",
                "body": self.body.identity(), "face_table_sha256": self.table_sha256,
                "gravitational_constant": GRAVITATIONAL_CONSTANT, "prefix": self.prefix,
                "precision": (precision.POLICY_NAME if self.precision == "fp64" else
                              {"policy": MIXED_POLICY, "near_factor": self.near_factor,
                               "far_table_sha256": self.table32_sha256,
                               "near_faces_per_target": (None if self.near_counts is None else
                                   {"mean": float(np.mean(self.near_counts)),
                                    "max": int(np.max(self.near_counts)),
                                    "min": int(np.min(self.near_counts))})}),
                "cuda_source": cuda_sources.source_identity("gravity"),
                "launch": {"block": BLOCK, "target_chunk": self.target_chunk,
                           "max_pairs_per_launch": self.max_pairs_per_launch,
                           "split_chunks": self.split_chunks,
                           "split_target_chunk": self.split_target_chunk},
                "memory_plan": self.plan.as_dict(), "timing_s": dict(self.timing)}

    def summary(self):
        return "%d faces, %s planned on device %d" % (self.body.nfaces,
                                                       human(self.plan.total_bytes),
                                                       self.device_id)

    def release(self):
        if self.released:
            return
        self.released = True
        self.arrays.clear()
        if self.reservation is not None:
            self.reservation.release()
            self.reservation = None
        if self.owns_budget and self.budget is not None:
            self.budget.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.release()
        return False


__all__ = ["GravityEvaluator", "face_table", "face_table32", "face_terms", "scaling_prefix",
           "GRAVITATIONAL_CONSTANT", "CapacityError"]
