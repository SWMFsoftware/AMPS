# Gravity reference: polyhedral-gravity 3.3 as implemented

C00 deliverable. Written 2026-10-19 from the source of the pinned release, before any
translation. It records what the CPU reference actually computes, so the CUDA kernel can
reproduce it operation by operation. The paper equations alone are not enough for that.

## 1. Pinned source and binary identity

| Item | Identity |
|---|---|
| Release | `polyhedral-gravity==3.3` (ESA / TU Munich, author Jonas Schuhmacher) |
| Source used for translation | PyPI sdist `polyhedral_gravity-3.3.tar.gz`, sha256 `1476ab9d638538f70b4d31e7cb1efce032c56c95b04264d88d005929769e7715` (matches the PyPI-published digest), unpacked in `third_party/src_pypi/` |
| Campaign binary | `/tmp/liouville_cops_polyhedral_gravity/polyhedral_gravity.cpython-38-x86_64-linux-gnu.so`, sha256 `9ef1df2e…6987`. It is byte-identical to the `.so` in PyPI's `polyhedral_gravity-3.3-cp38-cp38-manylinux_2_17_x86_64` wheel (wheel sha256 `a3b7fbe4…e671`, saved in `third_party/download_polyhedral_gravity_v3.3/`) |
| Build facts reported by the binary | `__parallelization__ = TBB`, `__commit__ = 77157f0b+modified`, logging INFO |
| CUDA environment copy | the same release from the cp311 wheel in `.venv` (`requirements-cuda.lock`) |

The campaign binary sits in `/tmp`, which can be cleared at reboot. If it is lost, it can be
recovered from the saved wheel, since the bytes are identical. The GitHub tag tarball
(`archive/refs/tags/v3.3`) returned HTTP 404 and was not used. The sdist is the canonical
source.

## 2. How the campaign calls it

`../LiouvilleCOPS/tools/build_physical_activity_templates.py::effective_slope` does the following:

* Body: the closed NASTRAN surface read by `read_nastran_surface`. Vertices are indexed from
  zero in sorted GRID-id order, and faces keep CTRIA3 order.
* `Polyhedron(…, density=533, normal_orientation=OUTWARDS, integrity_check=HEAL, metric_unit=METER)`.
* `GravityEvaluable(polyhedron)(points, parallel=True)`, called in chunks of 2048 targets.
* Targets: `centroid + 0.1 m * outward_normal` of every `elements.bin` facet. The
  centrifugal term is evaluated at the **centroid**.
* Only the acceleration (`record[1]`) is kept. The potential and the tensor are discarded.

## 3. Per-face algorithm (GravityEvaluable::evaluateFace, GravityModelDetail.cpp)

All vertices are first translated so the computation point P is the origin
(`Polyhedron::transformIterator(P)`: `v - P` for each vertex, in FP64). Every symbol below is
in those translated coordinates. "eps" means `EPSILON_ZERO_OFFSET = 1e-14`, an **absolute**
value in metres. `sgn(x, eps)` is -1, 0 or +1 with the dead band |x| <= eps.

Precomputed once per face (`prepare()`, independent of P):

1. `G_q = v_{q+1} - v_q` for q = 0,1,2 (cyclic). Tsoulis (2012) eq. 18.
2. `N = normalize(G_0 x G_1)`, the plane unit normal (eq. 19).
3. `n_q = normalize(G_q x N)`, the segment unit normals (eq. 20).

Per face and per point:

4. `sigma_p = sgn(N . v_0, eps)` (eq. 21; the double negation of the paper cancels).
5. Hessian plane: `(a,b,c) = (v_0 - v_1) x (v_0 - v_2)`, `d = sum((0 - v_0) * (a,b,c))`.
   `h_p = |d / sqrt(a^2 + b^2 + c^2)|`.
6. `P' = N * h_p`, with a per-coordinate sign fix from the intercepts `d/a`, `d/b`, `d/c`.
   An intercept is taken as 0 when its denominator is exactly 0. If the intercept is
   negative, the coordinate becomes `|P'_i|`. Otherwise it is negated when `N_i > 0` and
   kept when `N_i <= 0`. Away from h_p ~ eps this equals `sigma_p * h_p * N`. The kernel
   must use the package's literal branches anyway, because sigma_p (cutoff on `N.v_0`) and
   h_p (from the Hessian form) are computed independently and can disagree near zero.
7. `sigma_pq = -sgn(n_q . (P' - v_q), eps)` (eq. 23).
8. `P''_q = P'` when `sigma_pq == 0`. Otherwise it is the solution of the 3x3 system
   (eqs. 24–26) by Cramer's rule. The rows are `G = v_{q+1} - v_q`,
   `r2 = (v_q - P') x G` and `r3 = r2 x G`, with right-hand side
   `(G . P', r2 . P', r3 . v_q)`. The column matrix is the transpose of those rows. `det` is
   the explicit six-term expansion in `UtilityContainer.h`.
9. `h_pq = |P''_q - P'|`.
10. `l1 = |v_q|`, `l2 = |v_{q+1}|`, `s1 = |P''_q - v_q|`, `s2 = |P''_q - v_{q+1}|`. Signs follow
    Tsoulis (2021), checked in this order:
    * Option 4 applies if `|s1 - l1| < eps && |s2 - l2| < eps` (P lies on the segment's line):
      * case 2 (`s2 < s1`): negate all four;
      * case 1 (`|s2 - s1| < eps`): negate s1 and l1;
      * case 3: no change.
    * Otherwise, with `g = |G_q|`:
      * option 1 (`s1 < g && s2 < g`): negate s1;
      * option 2 (`s2 < s1`): negate s1 and s2;
      * option 3: no change.
11. `r_k = |P' - v_k|` for k = 0,1,2.
12. `LN_q = 0` when `(sigma_pq == 0 && (r_{q+1} < eps || r_q < eps)) || (|s1 + s2| < eps && |l1 + l2| < eps)`.
    Otherwise `LN_q = log((s2 + l2) / (s1 + l1))` (eq. 14). Note that `r1` is
    `r_{(q+1)%3}` and `r2` is `r_q`.
13. `AN_q = 0` when `h_p < eps || h_pq < eps`. Otherwise
    `AN_q = atan((s2*h_p)/(l2*h_pq)) + atan((s1*h_p)/(-l1*h_pq))` (eq. 15). The package uses
    xsimd here, with exactly these products and quotients. The second term is the
    subtraction written with a negated denominator.
14. Singularity terms (`computeSingularityTerms`). The first matching case wins:
    1. All `sigma_pq == 1` (P' inside the face): `(-2 pi h_p, -2 pi sigma_p N)`.
    2. Some q has `sigma_pq == 0` and `r_{q+1} < |G_q|` and `r_q < |G_q|` (P' on an edge):
       `(-pi h_p, -pi sigma_p N)`.
    3. Some q has `sigma_pq == 0` and (`r_{q+1} < eps` or `r_q < eps`) (P' at a vertex).
       The selection uses the first such q. If `r_{q+1} == 0.0` **exactly**, then
       `g1 = G_q` and `g2 = G_{q+1}`; otherwise `g1 = G_{q-1}` and `g2 = G_q`.
       `theta = pi/2` if `(-g1).g2 == 0`, else `acos((-g1).g2 / (|g1||g2|))`. The result is
       `(-theta h_p, -theta sigma_p N)`. The selection mixes an exact-zero test with the eps
       test; the translation keeps that.
    4. Otherwise `(0, 0)`.
15. `sum1 = sum_q sigma_pq h_pq LN_q`, `sum2 = sum_q sigma_pq AN_q`,
    `sum1T = sum_q n_q LN_q` (q = 0,1,2, in that order: `std::accumulate`).
16. `S = sum1 + h_p * sum2 + singA`.
17. The face returns:
    * potential `sigma_p h_p S`;
    * acceleration `N * S`;
    * tensor from `sub = (sum1T + N*(sigma_p*sum2)) + singB`, giving
      `(N_x sub_x, N_y sub_y, N_z sub_z, N_x sub_y, N_x sub_z, N_y sub_z)`.

Final scaling (`GravityEvaluable::evaluate`, `Polyhedron::getGravityModelScaling`):
`prefix = G * density * orientation_factor`, with `G = 6.67430e-11` m^3 kg^-1 s^-2
(`UtilityConstants.h`) and `orientation_factor = +1` for OUTWARDS and -1 for INWARDS.
The outputs are `V = prefix * sum / 2`, `g = -prefix * sum` and `Gamma = prefix * sum`.
The acceleration points toward the body. The campaign slope uses `-normal . g_eff`, which
agrees with that sign.

The `isCriticalDifference` warning in `evaluateFace` only logs; it never changes values.

## 4. Summation order (relevant for bitwise comparison)

With `parallel=True` and a list of points, the package parallelizes over **points** (TBB).
Each point calls `evaluate<false>`, a `thrust::host` (sequential CPP backend)
`transform_reduce` over faces 0..n-1. It starts from zero and adds the element-wise
(potential, acceleration, tensor) tuples. So each target's sum is a left fold in face order.

A CUDA kernel that assigns one target per thread and loops over faces in index order
(faces tiled through shared memory) reproduces exactly this fold. The kernel must be
compiled with `--fmad=false`: the manylinux wheel targets baseline x86-64 without FMA, so
its products are rounded separately. With that, the remaining differences come only from
transcendental functions (`log`, `atan`, `acos`) in CUDA libdevice versus glibc/xsimd. A
tiled reduction across threads (C03 benchmarks it) changes the order. Its difference from
the reference must then stay inside the C00 budget, which is set from the measured
reversed-order repeatability (`gravity_fixtures.json`).

## 5. Integrity check and HEAL

`HEAL` (used by the campaign) runs in the constructor:

1. **Degeneracy**: every face must have a surface area > 0, otherwise it throws.
2. **Orientation** (`checkPlaneUnitNormalOrientation`): for each face, a ray starts at
   `centroid + 1e-14 * N` and points along `N`. The tool counts the *distinct* intersection
   points (`std::set<Array3>`) with every face, using Möller–Trumbore with eps = 1e-14 and
   `t > eps`. An odd count marks the face as inward. That is O(n^2): 6.2e11 ray–triangle
   tests for 787,754 faces.
3. Majority rule: if at most n/2 faces are inward, the orientation is OUTWARDS and the
   inward faces are "violating". Otherwise it is INWARDS and the outward faces are
   violating.
4. HEAL sets the polyhedron orientation to the majority and swaps indices 0 and 1 of every
   violating face. The healed table is what gets evaluated.

### 5.1 Defect found in C00: HEAL corrupts kilometre-scale meshes given in metres

The ray origin offset `1e-14 * N` is absolute, in mesh units. For a body of about 2 km given
in metres, one rounding unit of the coordinates is about 4e-13 m, so the offset rounds away
and the origin stays on the face's own plane. Möller–Trumbore then hits the starting face at
`t ~ 1e-13 > 1e-14` about half the time. That makes the count odd, and the face is reported
as inward. Measured on 2026-10-19 (scratch scripts; reproduced by the C00 `gravity` stage):

| Mesh | Independent topology | Package check | HEAL result | Gravity error |
|---|---|---|---|---|
| icosphere L3, r = 1, 10, 100 m | outward | OUTWARDS, 0 flagged | unchanged | 0 |
| icosphere L3, r = 1000 m | outward | INWARDS, 612/1280 flagged | 612 faces flipped, orientation INWARDS | up to 121% |
| icosphere L3, r = 2000 m | outward | INWARDS, 588/1280 flagged | 588 flipped | up to 127% |
| SHAP5 788k BDF (|x| <= 2592 m) | closed 2-manifold, every directed edge once, Euler characteristic 2, V = 1.876e10 m^3 | exact replica: 138/300 sampled faces (46%) flagged, all of them own-face hits | about 362k faces would be flipped | see below |

On 48 evenly spaced 788k facets, 0.1 m offset, flipping a random 46% of faces (the flags
depend on rounding, so they are effectively random) gives these results relative to the
correct winding:
* |g| error: median 92%;
* slope error: median 81°, p95 144°;
* facets at ≥ 60°: 38 of 48, against 2 of 48 for correct gravity.

The correct gravity is plausible: |g| from 1.4e-4 to 2.3e-4 m/s^2, median slope 18°.

**Consequences.**
1. Every slope and cliff product of `build_physical_activity_templates.py` with HEAL and
   METER is invalid. That includes the rejected 100k bundle (its 5° resolution-gate failure
   is probably this artifact) and the 788k build running on 2026-10-19. The late
   NG-not-RG mask does not use gravity and is unaffected.
2. The C00 gravity reference therefore evaluates the input winding with `DISABLE`. It is
   admitted only when the independent topology audit passes: closed, consistently
   oriented, positive volume, no zero-area face, Euler characteristic 2. The package
   HEAL/VERIFY is never used on this body. The replica of the package test on a sample is
   kept only as evidence.
3. **Port policy:** the CUDA gravity backend does not port the ray-cast HEAL. It requires a
   body that passes the topology audit and refuses any other mesh. The evaluation formula
   (§3) is unaffected by this defect.
4. The fix of the CPU reference builder (use DISABLE after the topology audit, or apply
   HEAL to coordinates scaled to km) is the user's decision. It changes the reference
   application and is not made by C00.

## 6. Singular and near-singular cases the tests must cover

| Case | Where in the code | C00 fixture |
|---|---|---|
| P' inside a face (`h_p` arbitrary) | singularity case 1 | `near_face_*`, offset targets |
| P exactly on a face plane (`h_p = 0`, `sigma_p = 0`) | AN = 0 branch, `P' ~ 0` | `on_face_centroid_*`, 788k `on_surface` |
| P' on an edge line | option 4, singularity case 2 | `on_edge_midpoint_*`, `near_edge_*` |
| P' at a vertex | LN = 0 branch, singularity case 3 | `on_vertex_*`, `near_vertex_*` |
| Far field | cancellation in `S` (isCriticalDifference) | `far_5R`, `far_50R` |
| Interior points | sum of the case-1 terms | `interior_*` |

**Measured in C00 (`gravity/gravity_synthetic.npz`):** for points placed *exactly* on a
vertex or an edge midpoint, the package returns NaN potential and acceleration in some
cases and finite values in others. Which cases fail depends on rounding: the same
point on a translated or rotated copy may switch between NaN and finite. Examples are
icosphere L0/L1/L3 `on_vertex_*` and `on_edge_midpoint_*`, and the cube and tetrahedron
vertices after rotation. The exact-coincidence branches (LN = 0, singularity case 3) need
`r_q < 1e-14` or `sigma_pq == 0`, and a rounding residue can defeat both. Face-centroid
points exactly on the surface (32 of them on the 788k body) and all 0.1 m and 0.05 m
offset targets are finite. C03 comparison rule: where the reference is finite, the kernel
must agree within budget. Where the reference is nonfinite, the case is reported and
excluded, never "fixed" silently. Campaign targets never fall in it.

The eps thresholds are absolute (1e-14 m). For 67P coordinates (~2 km) a rounding unit is
about 4.5e-13 m, so these branches trigger only for exact coincidences such as an exactly
on-surface point or a shared vertex, never for the 0.1 m offset targets. The kernel must still
contain every branch.

## 7. License obligations (GPLv3), decided 2026-10-19

The user chose to **translate the GPL source**. The CUDA gravity kernel is therefore a
derivative work of polyhedral-gravity-model, which is under the GNU GPL version 3
(`third_party/src_pypi/polyhedral_gravity-3.3/LICENSE`). Obligations:

* `cops_cuda/kernels/gravity.cu` (C03) carries a GPLv3 notice, the original copyright
  attribution (Jonas Schuhmacher and contributors, github.com/esa/polyhedral-gravity-model;
  the source files carry no per-file copyright lines), the list of
  source files it is derived from (`GravityEvaluable.cpp`, `GravityModelDetail.cpp`,
  `Polyhedron.cpp`, `UtilityContainer.h`, `UtilityConstants.h`), and dated
  modification notes.
* A copy of the GPLv3 text ships with the application (`third_party/…/LICENSE`, referenced
  from `third_party/README.md`).
* Running the code privately carries no release obligation. **Conveying** it (giving
  binaries or sources to anyone else, including collaborators) requires providing the
  corresponding source under GPLv3. If it is conveyed, any program distributed as a single
  combined work with the kernel must be under GPLv3-compatible terms. Keeping the gravity
  kernel and its loader in separate, clearly marked files keeps that boundary visible.
* Cite the package paper (doi:10.21105/joss.06384), Tsoulis (2012), Geophysics 77(2),
  F1–F11 (doi:10.1190/geo2010-0334.1), and Tsoulis & Gavriilidou (2021), Geophysical
  Prospecting 69(8–9), 1745–1760 (doi:10.1111/1365-2478.13134; the sign rules of step 10),
  in any publication using the results.

## 8. The CUDA translation as built (C03)

| Piece | File | License |
|---|---|---|
| Body reader binding, topology audit, identity | `cops_cuda/gravity_body.py` | project code |
| Face precompute (`prepare`), final scaling, device driver | `cops_cuda/gravity.py` | GPLv3 (derived) |
| Per-face evaluation, fold and split kernels | `cops_cuda/kernels/gravity.cu` | GPLv3 (derived) |

**Face table.** One row of 30 doubles per face: v0, v1, v2 (original vertices), G0..G2, N,
n0..n2. `gravity.face_terms` evaluates `prepare()` in NumPy with the C++ operation order. On
every synthetic shape and on the 787,754-face body it is **bitwise identical** to the package's
`GravityEvaluable.__getstate__()`.

**Per face and target.** `face_contribution` in `gravity.cu` follows §3 steps 1-04 to 1-14 and
steps 2, 4, 5 and 7 in the source's evaluation order. It is compiled with `--fmad=false`, so
products and sums are rounded as on the CPU; `sqrt` and division are correctly rounded on
both sides. The only intended differences are libdevice `log`, `atan` and `acos` versus glibc
`log` and `acos` and xsimd `atan`. The tensor and singularity vector are not computed (not
used), and they do not affect potential or acceleration. Every branch is kept, including:
* the exact-zero tests (`sigma_pq == 0.0` and the `r1 == 0.0` vertex selection);
* the −0.0 from `sgn(...) * -1.0`;
* the literal sign rules for P′.

**Summation.**
* `gravity_fold`: one thread per target, faces in index order through a 32-face shared-memory
  tile, sums carried in global memory across launches. This is the package's left fold. The
  result does not depend on chunking (tested bitwise).
* `gravity_partial` + `gravity_reduce`: faces in chunks, partial sums combined in a fixed
  chunk order. It is deterministic but rounded differently; the measured difference from
  fold is ≤ 16 eps·M. It is used when there are few targets.

**Term magnitude M.** This is a fifth accumulator: the sum over faces of the magnitudes of
all terms that carry a transcendental result. The magnitude-based budget is
|Δg| ≤ 16 eps·M (VALIDATION.md, C03 amendment). The probes give the reason: up to 4e-10
relative error at 50 R, but at most 0.33 eps·M. The near-surface 788k targets have
M/|g| ≈ 1e6, and their relative error is still ≤ 2.3e-13.

**NaN parity.** On the exact on-vertex and on-edge probes where the package returns NaN, the
GPU returns NaN at exactly the same points. These points are reported, never repaired.

**Launch bound.** Both A6000s on this host report `KernelExecTimeout` = 1. A launch is
limited to 1e8 face–target pairs (about 0.5 s). In fold mode this changes no result.

**Resources.** `gravity_fold` uses 162 registers and 72 B of local stack per thread with a
7.5 KB shared tile. Throughput is about 2.3e8 face–target pairs/s per A6000 in FP64; see
STATUS.md for measured times.
