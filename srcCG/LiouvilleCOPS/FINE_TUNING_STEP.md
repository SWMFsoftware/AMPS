# Late-interval activity fine tuning

## Purpose and scope

This is an optional **second fitting stage**, applied to a completed CPU or
CUDA fit. It targets a late pressure mismatch without changing the earlier
cached NG/RG pressure curves. It does not increase harmonic order and does not
fit individual facets independently.

The final `m` minutes are split into `n` equal-duration intervals. Facets are
assigned to the first interval in which they contribute, provided they never
contributed before that window. Each interval supplies up to three correction
parameters: NG-dominated, RG-dominated, and all newly contributing facets.
The dominance ratio is adjustable and defaults to `1e-3`.

Implementation is confined to LiouvilleCOPS:

- `tools/fine_tuning.py`: interval assignment, kernel-response construction,
  support classification, compact basis construction, and reduced optimizer.
- `tools/fine_tune_activity.py`: CLI, input/provenance checks, preparation
  reuse, reports, comparison figure, and C++ replay export.
- `tests/test_fine_tuning.py`: independent support/solver tests and an actual
  binary-kernel-to-export synthetic pipeline.

The CUDA sibling is neither imported nor modified. The algorithm and array
contracts below are intended to support a later CUDA implementation.

This first version implements **nonnegative additive enhancements only**.
It cannot lower a group's physical activity below its frozen fitted value.
That is an explicit modeling restriction suitable for testing an additional
source explanation of the NG excess, not a claim that every late residual
can be explained this way. Signed corrections would require a separately
specified positivity-constrained formulation.

## 1. Inputs and the frozen first-stage solution

Required inputs are:

1. The original manifested, unenhanced static-kernel fitting base: model table,
   `elements.bin`, all `kernels/point=NNNN.bin`, and `input.snapshot.in`.
2. The completed fit's CPU/CUDA-compatible `fit_restart.npz`.
3. `m`, `n`, the dominance ratio, and a new output directory.

The stage freezes the **completed fit**, not the original run prior:

- its full H2O activity field, including all previously fitted enhancements;
- its three analytical temperature parameters;
- its CO2 multiplier and original CO2 spatial map;
- all source laws, illumination/shadow geometry, gauge settings, and kernels.

The restart's `activity_linear` is converted to the physical nonnegative
activity `A0 = max(activity_linear, 0)`. This is the frozen baseline. The new
stage does not call `Problem` with `--fix-*`, which would otherwise retain
values from the original base rather than necessarily the completed fit.

Before creating outputs, the runner checks the restart mesh identity, point
IDs/UTC ordering, base elements/table/input hashes, exact prior coefficients,
and the base's archived kernel identities. A manifested BodyFrame is required.
The supported forward mode is `FlightTimeRotation=none` with analytical
temperatures, matching the current version-1 kernel fitter.

It independently recomputes the frozen fit's NG/RG pressures and requires
agreement with the restart predictions within relative tolerance `1e-6`
(absolute tolerance `1e-30 Pa`). The actual maximum relative difference is
reported. The looser serialization/replay tolerance is NOT an allowance for
changing earlier pressures during fine tuning: the reduced-model pre-window
change is checked separately and must be exactly zero.

The restart's reference pressures and selection/saturation flags are retained.
There is no automatic re-enabling of excluded or saturated observations.
Inside-shape rows are never fitted. Geometry/support checks still scan every
supplied row, irrespective of observation selection.

## 2. Window, intervals, and strict first-response membership

Let `t_end` be `--end-utc`, or the **last modeled UTC** if omitted. This default
is recorded explicitly; it is not presented as a separately established
physical impact epoch. The end cannot precede any supplied modeled row.

Define edges

```
e[j] = t_end - 60*m + j*(60*m/n),  j = 0,...,n.
```

Intervals are `[e[j],e[j+1])`, except the last includes its right endpoint.
Rows must have strictly increasing UTCs. Every interval needs at least one
modeled sample, and there must be samples preceding the window. Otherwise the
command rejects the request rather than manufacturing evidence or parameters.

For each facet, inspect its records chronologically. A structural contribution
exists if ANY of `Kn`, `KPhi`, or `KV` is nonzero. All three must be exactly zero
to count as structurally absent. Missing facet records likewise have zero
response. This is conservative: a structural response can exclude a facet
even if a particular frozen source or velocity factor makes its current
pressure contribution vanish.

Assign `first_interval` as:

- `-2`: no structural contribution anywhere in the supplied rows;
- `-1`: first structural contribution before the late window, permanently
  ineligible for this refinement;
- `j >= 0`: first structural contribution in interval `j`.

Crucially, membership excludes response at **any earlier time**, not merely
the preceding interval. A facet that contributes early, disappears, and
reappears late must not become eligible again. The dominance threshold is
never used to turn a small early response into an exact zero.

Every eligible facet has exactly one interval assignment. It keeps its new
activity at every later time: this is a spatial source correction, not a
source that is artificially switched on/off at interval boundaries.

## 3. Compare gauge pressure responses, not raw kernels

For fixed fitted temperature `T_i(t)`, define the actual H2O pressure
responses to unit **final facet activity**:

```
H_NG(t,i) = cNG_H2O * bm_H2O * Q_H2O(t,i) * Kn(t,i) / sqrt(T_i(t))
H_RG(t,i) = cRG_H2O(t) * Q_H2O(t,i) *
            [KPhi(t,i) + bm_H2O * VdotB(t) * KV(t,i) / sqrt(T_i(t))].
```

The existing production `fit_cops.Model` supplies `Q`, `bm`, `cNG`, `cRG`,
SZA/shadow flags, and velocity factors. `row_response` constructs one row's
Model at a time to reuse those definitions without holding the complete
transport cache. The temperature law is the same analytical law used by
`Model.evaluate` but evaluated at the **frozen fitted temperatures**.

Both responses have units Pa per unit dimensionless activity. Raw `Kn` and
`KPhi` are not directly comparable. The response excludes the derivative of
the nonnegative activity clip: zero current activity does not imply zero
potential response to an added source. Nonfinite or negative pressure
responses are rejected rather than silently clipped or made positive in the
forward calculation.

For a facet first assigned to interval `j`, accumulate only within that
interval:

```
R_NG(i) = sum_{t in interval j} w(t) * abs(H_NG(t,i))
R_RG(i) = sum_{t in interval j} w(t) * abs(H_RG(t,i)).
```

Weights are nearest-sample time-cell lengths: cell boundaries are adjacent
sample midpoints and the interval edges. Their sum equals the interval
duration. A one-sample interval receives its full duration. These quadrature
weights classify dominance; they do not change the pressure-fitting weights.

With `epsilon = --dominance-ratio` (`0 < epsilon < 1`, default `1e-3`):

- NG-dominated: `R_NG > 0` and `R_RG <= epsilon*R_NG`.
- RG-dominated: `R_RG > 0` and `R_NG <= epsilon*R_RG`.
- Mixed: neither dominance condition, with at least one nonzero score.
- Both zero: omit the facet from the correction support.

This is **integrated first-interval dominance**, not an assertion that the
facet is invisible to the other gauge at every sample or in later intervals.
All later gauge responses remain in the fitting calculation.

## 4. Three coefficients and a compact basis

The disjoint atomic groups are `(interval, NG/RG/mixed)`. Store facet membership
as an integer `atom[i] = 3*j + type`, where type is 0, 1, or 2. Ineligible facets
have `atom = -1`.

The default `--third-group all` implements the requested overlapping basis:

```
NG atom:     correction = A_scale * (c_j_NG + c_j_all)
RG atom:     correction = A_scale * (c_j_RG + c_j_all)
mixed atom:  correction = A_scale * c_j_all.
```

`--third-group mixed` instead fits three disjoint corrections. This alternative
is explicit, not a silent reinterpretation of "all". For enhancement-only
coefficients, the alternatives can have different feasible fields: overlapping
"all" imposes a common enhancement beneath the two gauge-specific extras.

`A_scale` is the frozen restart's existing enhancement scale. It is not fitted.
All new coefficients are bounded below by zero. No upper bound is imposed in
this version; inspect fitted amplitudes and physical plausibility.

Empty basis columns are omitted. If an interval has no mixed atom, its `all`
column equals the sum of the nonempty NG/RG columns. The implementation omits
that `all` column deterministically, preserving every attainable nonnegative
physical correction. Thus `3*n` is a requested maximum, not a promise of that
many identifiable parameters. Every omission and its reason is reported.

Internally, a small matrix `B[3*n, K]` maps fitted amplitudes to atom amplitudes.
There is **no dense nfaces-by-K mask matrix**. An atom index plus this small
matrix reconstructs every facet correction exactly.

## 5. Precompute the small pressure problem

A second chronological kernel pass accumulates:

- frozen pressure `p0[rows,2]`, including frozen CO2;
- atom response `W[rows,2,3*n]`, summing `H` for the facets of each atom.

No temporal gate is applied in this pass. A group first seen in interval zero
can contribute throughout every later interval and to either gauge.

For EVERY group interval `j`, the implementation asserts:

```
W[bins < j, :, 3*j:3*j+3] == 0 exactly.
```

Define `D = A_scale * W @ B`. At fixed physics,

```
A_final(i) = max(L_fitted(i),0) + A_scale*(B@c)[atom(i)]
p_final    = p0 + D@c.
```

Ineligible facets receive zero correction. Because corrections are
nonnegative, no new clipping boundary complicates the reduced problem.
These pressure equations are exact within the cached model, including for
initially inactive facets; they are not a first-order Jacobian approximation.

Preparation uses two streaming passes, plus the existing byte-level cache
audit. Optimization then costs only small dense operations on rows-by-K
arrays. It never reevaluates all facet contributions each iteration. The
current preparation loop is serial; BLAS threading is library-controlled.
Requesting 36 optimizer processes would not accelerate one small coupled fit.

## 6. Objective, regularization, and diagnostics

All new coefficients are optimized **jointly**, not greedily interval by
interval. Use selected positive late NG and RG reference pressures from the
restart, with log residuals and equal total weight per gauge:

```
r_g(t) = [log(p_final_g(t)) - log(p_reference_g(t))] / sqrt(N_g).
```

The analytic data Jacobian is `D_g(t,k) / [sqrt(N_g)*p_final_g(t)]`.
Both gauges must have selected positive late observations. A coefficient with
exactly zero response across all selected late data is fixed to zero and
reported, not assigned a meaningless fitted value.

Optional regularization (`--lambda`, default zero) penalizes the TOTAL added
activity, including cross terms of overlapping masks:

```
G = A_scale^2 * B.T @ diag(atom_area) @ B / sum_i area_i*A0_i^2
penalty squared = lambda * c.T @ G @ c.
```

The residual square root uses the symmetric eigendecomposition of `G`.
The optimizer is bounded SciPy `least_squares` (TRF), analytic Jacobian,
`x_scale='jac'`, and `ftol=xtol=gtol=1e-10`, with adjustable `--max-nfev`.
Initial correction is zero. SciPy may internally move a zero start slightly
inside its nonnegative bounds.

Report singular values and numerical rank of the fitted weighted data
Jacobian, plus a rank-deficiency warning. Rank is a numerical diagnostic, not
an uncertainty estimate or proof of identifiability under measurement error.
The initial/final objective, late NG/RG log RMS, bound-constrained amplitudes,
group areas/counts, and area-integrated added activity are recorded.
Area-integrated activity is NOT a gas production rate: source production
also depends on illumination and the prescribed source law.

All pre-window reduced pressures must be bitwise unchanged. C++ full-forward
verification remains separate and may differ by floating-point summation.
The guarantee concerns the cached time samples and modeled gauges, not every
continuous time, a different trajectory, a different instrument model, or
non-static/time-of-flight transport absent from the current fitter.

## 7. Commands and outputs

Run from LiouvilleCOPS. Replace the completed-fit path and interval count with
your chosen experiment; no production fit is selected implicitly.

```bash
# Inspect support and build reusable responses without fitting.
python3 tools/fine_tune_activity.py output_fit_shap5_788k_base \
  --restart-from completed_fit/fit_restart.npz \
  --minutes 5 --intervals 5 --dominance-ratio 1e-3 \
  --third-group all --prepare-only --out late5m_prepare

# Fit using that identical preparation and a new output directory.
python3 tools/fine_tune_activity.py output_fit_shap5_788k_base \
  --restart-from completed_fit/fit_restart.npz \
  --minutes 5 --intervals 5 --dominance-ratio 1e-3 \
  --third-group all --lambda 0 --max-nfev 1000 \
  --prepared late5m_prepare/fine_tuning_prepared.npz --out late5m_fit

# Optional independent full C++ replay; 36 threads as requested for production.
env OMP_NUM_THREADS=36 ./liouville_cops /absolute/output/path/late5m_fit/cops_fine_tuning.in
```

The generated replay input preserves the original input settings except for
the fitted source field, temperatures, and new verification output path.
An explicit `Threads` setting in the input can take precedence over the
environment; inspect it when requesting a particular thread count.

Relative products resolve through `COPS_OUTPUT_ROOT` (default
`<AMPS>/cops_output`); absolute output directories must be outside LiouvilleCOPS.
Existing directories are refused. Inputs are not modified.

Outputs:

| File | Contents |
|---|---|
| `fine_tuning.json` | identities, settings, support groups, omissions, diagnostics, optimizer result, output hashes |
| `fine_tuning_prepared.npz` | first membership, atom indices, dominance scores, timing weights/edges, frozen pressures, atom-response matrices |
| `fine_tuning_state.npz` | fitted coefficients, compact basis, reduced design/Gram matrices, frozen/final activity, predicted pressures and selected references |
| `pressure_comparison.txt` | all-row before/after NG/RG pressures, references, interval, selection, saturation |
| `pressure_comparison.png` | full-time and late-window before/after comparisons for both gauges |
| `ActivityFit.txt` | frozen H2O legacy coefficients and CO2 coefficients with fitted kappa folded in |
| `ActivityEnhancementFit.txt` | complete final residual per facet, in the existing strict C++ format |
| `cops_fine_tuning.in` | independent forward replay configuration, output `<fine-tuning-directory>_verify` |

The fine-tuning products use schema `liouville_cops.late_interval_fine_tuning.v1`.
`fine_tuning_state.npz` is deliberately NOT named `fit_restart.npz`: it is not
the C12 optimizer-state representation and cannot be passed to the original
fitter's restart loader. It describes this separate refinement stage. The
original C12 restart remains unchanged, and the C++ export carries the final
physical field. A future GPU port should consume this explicit schema rather
than silently treating its coefficients as harmonics or ordinary templates.

The preparation file verifies all numeric array hashes and a definition tying
it to the input restart, base table/elements/kernel hashes, window, epsilon,
and CPU implementation hashes. Lambda and third-group choice can change
without rebuilding the atom responses; window, n, epsilon, source state,
kernel inputs, or implementation changes invalidate reuse. The base cache is
still audited on reuse. No old gravity or visibility product is reused merely
because its filename or face count matches.

`--prepare-only` saves no optimized activity or replay input. No eligible
columns produces an explicit unchanged control result. Hitting `max-nfev`
can produce a complete output directory with `optimizer.success=false`;
process completion is not convergence. Full-forward verification and
scientific review are recorded as not run/pending.

## 8. CUDA port: contracts and development stages

Preserve these definitions before optimizing execution:

1. **Read-only interchange and frozen state.** Load the same C12 restart,
   validate all base identities, keep the same observation selections, and
   reproduce the saved baseline predictions. Never substitute the original
   base's temperatures or kappa for the completed fit's values.
2. **First-response reduction.** For each kernel record with any nonzero
   geometric kernel, reduce the earliest row index for its facet. Convert that
   row index to `first_interval`, with the same pre-window/never-seen codes.
   Parallel reduction by earliest index avoids chronological GPU races.
3. **Pressure weights.** Reuse the validated GPU physical evaluator to produce
   FP64 `H_NG/H_RG` per contribution at frozen temperatures. Do not use the
   clipped-activity Jacobian or the actual current facet flux as visibility.
4. **First-bin classification.** Reduce weighted responses only where the
   contribution's interval equals the facet's first interval. Apply the exact
   inequalities and zero-denominator cases specified above. GPU summation can
   change threshold-near classifications: compare scores, identify those
   facets explicitly, and reproduce the CPU classification or declare a
   reviewed tolerance policy rather than silently changing masks.
5. **Persistent atom matrix.** Reduce every contribution, including those in
   later intervals, into `W[row,gauge,atom]`. Require exact structural zeros
   in all earlier intervals. Deterministic reductions are preferable for
   reproducibility; never materialize an observation-by-facet dense matrix.
6. **Reuse compact CPU optimization initially.** Transfer only `W`, p0, atom
   areas, and group metadata to the host. Since K is small, the optimizer does
   not need to run on a GPU to benefit from GPU preparation. Verify identical
   residual/Jacobian/regularization definitions before changing solvers.
7. **Export and cross-validation.** Produce the same versioned arrays, group
   metadata, pressure tables, and C++ facet export. Preserve area, facet order,
   source identities, and all frozen components. Mark the GPU implementation
   in the preparation definition rather than pretending it has CPU hashes.

Array layouts (zero-based indexing):

```
first_interval: int32 [Nfacet], codes -2/-1/0..n-1
atom:           int32 [Nfacet], -1 or 3*j+type
ng_score/rg_score: float64 [Nfacet], Pa*s per unit activity
bins:           integer [Nrow], -1 or 0..n-1
weights:        float64 [Nrow], seconds
edges:          float64 [n+1], same parse_utc epoch as CPU
baseline:       float64 [Nrow,2], Pa; gauge order NG,RG
atom_response:  float64 [Nrow,2,3*n], Pa per unit activity
basis:          float64 [3*n,K], 0/1 mapping
design:         float64 [Nrow,2,K], Pa per scaled coefficient
gram:           float64 [K,K], normalized total added activity norm
coefficients:   float64 [K], nonnegative
```

Numeric SHA256 follows `fit_restart.array_sha256`: dtype string, shape tuple,
then contiguous C-order bytes. NPZ files do not require pickle.

## 9. Validation gates and interpretation

```bash
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 \
  python3 -m unittest discover -s tests -p test_fine_tuning.py -v
```

Tests cover historical reappearance, tiny early contributions, dominance
thresholds, interval boundaries and empty bins, overlapping/all redundancy,
later cross-gauge responses, known coefficient recovery, regularization,
zero-sensitivity/empty controls, preparation corruption, cache-identity
mismatch, full binary-base CLI execution, independent `Model.evaluate` replay,
and the actual C++ strict per-facet reader.

Production acceptance additionally requires the chosen real base/restart,
inspection of support maps and threshold-near facets, comparison of both
gauges, a full C++ replay, and sensitivity studies in m, n, epsilon, and lambda.
Inspect effective source strengths: an unconstrained enhancement can fit data
without validating a causal explanation or the collisionless physical regime.
A late-only region can be empty or too weak to explain the jump. Do not weaken
earlier-zero protection or add a hidden time switch to force an improvement.
The partition is trajectory-conditioned, not a claim about intrinsic
geological boundaries. Report that distinction in publications.
