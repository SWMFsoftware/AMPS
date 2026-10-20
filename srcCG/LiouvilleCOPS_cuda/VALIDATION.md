# Numerical acceptance policy

The machine-readable form is `validation/tolerances.json`.

**Status: approved by the user on 2026-10-19**, except for one item that stays
**provisional** until a fixed decision point (marked below). The Jacobian budget was
approved at the end of C06 (2026-10-19). This closes the C00 gate.
Every number below comes from a measurement:
* C00 capture `cops_output/cuda_reference/c00_20261019T064536Z`;
* the C01 doctor reports;
* C03 report `cops_output/cuda_c03/gravity_20261020T022050Z`.

Kernels are validated against these budgets and never tuned toward them. A budget is changed
only by a documented amendment that the user approves.

## Principles

* **Numerical error** is the difference between two implementations of the *same*
  formula. It is budgeted here.
* **Physical-model uncertainty** is separate. Examples: the outward offset (0.25° gate,
  measured 0.005°) and the gravity-shape resolution (5° p95). Those audits are never used as
  GPU round-off allowances.
* Where a quantity is exact by construction (IDs, masks, ordering), the comparison is exact.
* Floating-point budgets scale with the **sum of magnitudes** M of the summed terms, not with
  the result. That keeps them meaningful under cancellation and near zero.

## Measured reference facts behind the numbers

| Quantity | Measurement | Source |
|---|---|---|
| CPU evaluator summation error (`np.bincount`) | ≤ 216 eps·M over 6 configurations × 3–6 vectors, 47 points, 8.78M contributions | `fit/fit_fixtures.json` |
| RG cancellation M/\|p\| | 1.0, so no cancellation for these vectors | same |
| GPU fixed-order segmented sum | ≤ 0.31 eps·M, bitwise repeatable | doctor reports |
| cuSPARSE SpMV/SpMM | not bitwise repeatable | doctor reports |
| C++ replay on the subset | 3.7e-11 (NG), 9.1e-10 (RG) | `fit` stage log |
| Analytic vs central-difference Jacobian | ≤ 6.1e-6 on non-crossing columns | `fit/fit_fixtures.json` |
| Gravity face-order effect (788k) | \|Δg\|/\|g\| ≤ 9.5e-14; Δslope ≤ 5.2e-12° | `gravity/gravity_fixtures.json` |
| Offset-halving audit (788k, 512 samples) | 0.0050° (gate 0.25°) | same |
| GPU vs package gravity, synthetic probes | ≤ 0.33 eps·M; up to 4e-10 relative at 50 R (cancellation) | C03 report |
| GPU vs package gravity, 788k targets | ≤ 2.3e-13 relative; slope ≤ 9e-12°; M/\|g\| up to 1.1e6 | C03 report |
| Facets within 1e-6° of 45° or 60° (full 788k field) | 0 | C03 report |

## Approved budgets

**Exact** (the comparison must be identical):
* point IDs, flags, reasons and row order;
* facet order and geometry identities;
* all visibility masks, bitwise;
* cliff masks, after the threshold band is resolved on the CPU;
* parameter names and bounds;
* regularization rows. They are computed by the shared CPU code under the same NumPy.

**FP64 evaluator:**
* Pressure vs the exact sum: ≤ 64 eps·M_i, where M_i is the per-point sum of the magnitudes
  of the contributions.
* Pressure vs the CPU evaluator: ≤ 512 eps·M_i (CPU ≤ 216 measured, GPU ≤ 64, with about
  1.8× margin). For the C00 vectors that is about 1.1e-13 relative.
* **Log residuals** follow from the pressure budget. The residual of point i in gauge g is
  r_i = w_g (ln p_i − ln p_obs,i), with w_g = 1/√N_g, so
  |Δr_i| ≤ δ_i := w_g · 512 eps · M_i / p_i.
  The pressure floor is applied identically on both sides.
* **Jacobian vs CPU (approved 2026-10-19 at the end of C06; replaces the provisional
  1e-11 × column max):** **|J_gpu[i,j] − J_cpu[i,j]| ≤ 1024 eps · M_J[i,j]** per entry of the
  pressure Jacobian.
  * M_J[i,j] is the sum of the magnitudes of every contribution entering entry (i, j),
    computed from the reference's own per-contribution terms
    (`tools/c06_jacobian_budget.py`, `jacobian_magnitudes`).
  * Entries with M_J = 0 must be exactly 0 on both sides.
  * Evidence: `cops_output/cuda_c06/jacobian_budget_20261020T025421Z/report.json`.
    * Coverage: the 47-point C00b subset and the 23-point verification subset of the full
      run, B29/N30/C30/NC31/R30, all fixture states (7–10 each).
    * Observed maximum 245 eps·M_J (activity column a02, NG). Temperatures reach 236,
      templates 134 and κ 233.
    * Every zero-magnitude entry is exactly zero on both sides.
    * Correction (2026-10-20): the first report (`jacobian_budget_20261019T223008Z`, 329)
      used KT instead of |KT| in the Tmax/T75 magnitudes, where KT = 1/(1 − 1/cos 75°) < 0.
      That left T75 unmeasured and Tmax overestimated. The tool now uses |KT| and refuses
      negative magnitudes. The approved constant 1024 is unchanged; its margin is 4.2×.
  * Why 1024, not 512: the difference is dominated by the CPU's own sequential summation,
    as for pressures. That error grows with a point's record count, and the CPU reference
    can be held only for subsets. 1024 gives a 4.2× margin and is ≈ 2e-13 relative to the
    entry magnitude.
  * The column-max-relative error (≤ 2.9e-14 observed) stays in reports as a diagnostic.
    The analytic vs central-difference comparison (≤ 6.1e-6 away from branch crossings) is
    a sanity check, never a tolerance.
* **Objective (amended at sign-off; replaces the fixed 1e-12 relative bound).** With
  cost = ½ Σ r_i² + the regularization term:
  **|Δcost| ≤ Σ_i |r_i| δ_i + ½ Σ_i δ_i²**.
  The regularization term is computed by the same CPU code from the same x, so it contributes
  nothing. A fixed relative bound would fail a correct evaluator exactly when the fit is good:
  with rms residual ≈ 0.05, 2·max|Δr|/rms(r) ≈ 4e-12.

**Gravity:**
* **Term magnitude.** M = G ρ · Σ_faces [ Σ_q |σ_pq h_pq LN_pq| + h_p Σ_q |σ_pq| (|atan t0_q| +
  |atan t1_q|) + |sing_A| ]. The kernel accumulates it alongside g. It bounds how far the
  few-ulp libdevice vs glibc/xsimd differences of log, atan and acos can propagate,
  whatever the cancellation.
* **Acceleration:** **|Δg| ≤ 16 eps · M** wherever the package result is finite, **and**
  ≤ 1e-12 relative for the campaign targets (788k facets at 0.1 m and 0.05 m offset, and on
  the surface).
  * A budget relative to |g| alone is not used for every target. Far from the body each AN
    term is the sum of two atans near ±π/2, and at the centre of a symmetric body |g| is a
    tiny remainder. There the relative error reaches 4e-10 while |Δg|/M ≤ 0.33 eps.
* **Potential:** |ΔV| ≤ 16 eps · M · (largest target–vertex distance), because h_p is at most
  that distance.
* **Slope:** ≤ 1e-9°.
* **Threshold band:** facets within 1e-6° of a threshold are recomputed by the CPU reference
  and take its classification.
* **NaN parity:** the NaN pattern must equal the package's. Points where the package returns
  NaN (exact vertex or edge positions only) are reported and excluded, never repaired.

**Repeatability:** bitwise for a pinned configuration (device model, launch configuration,
CuPy and NVRTC version). Cross-configuration differences, such as another chunking, variant or
device, are measured and must satisfy the budgets above.

**Mixed precision:** no budget. C08 must propose and justify its own, and that needs a
separate approval.

## GPU-built fitting base (C11; approved by the user 2026-10-20)

A base built by `native/liouville_cops_cuda` is compared with the base the C++ `liouville_cops`
builds from the same input (`tools/c11_compare_bases.py`). Bitwise identity is impossible:
the C++ binary contracts multiply-adds into FMA (GCC `-march=native`), while the GPU uses
FP64 with `--fmad=false`.

| Check | Bound |
|---|---|
| `elements.bin` | bitwise |
| kernel header geometry | bitwise, except the altitude: ≤ 1e-14 relative |
| CosSZA of common records | bitwise |
| shadow flags | every mismatch listed; ≤ 1e-6 of lit facets |
| record set | every flipped facet listed; flipped share of each point's Kn, KPhi, KV totals ≤ 1e-9 |
| common records, per point | Σ\|ΔK\| / Σ\|K\| ≤ 1e-12 |
| single records | reported (a subdivided facet can differ by one partially visible sample) |
| `COPS_model.dat` model columns (11 significant digits printed) | ≤ 2e-10 relative (one unit in the last printed digit) |
| fitter | the unchanged fitter accepts the base; consistency check ≤ 1e-6 |
| fit on the GPU base vs the C++ base | the end-to-end budget (provisional; see below) |

Evidence: `cuda_c11/compare_20261020T063859Z` (788k) and
`cuda_c11/compare100k_20261020T065426Z` (100k). A B29 fit on the GPU-built 788k base is
bitwise identical to the fit on the C++ base.

## Provisional item (not yet approved)

* **End-to-end FP64 fit:** the same termination class, final cost within 1e-9 relative, an
  identical active-bound set, and parameters within 1e-3 posterior σ. Along null or
  near-null directions, predictions must agree within 1e-8 relative instead. *Decision point:*
  after the first full CPU-vs-GPU fit in C06/C07, set from the observed differences. TRF
  amplifies 1e-13 evaluator differences along its path, so the parameter criterion may be
  too strict for weakly determined directions.

## Physical audits (not numerical allowances)

* The 0.25° outward-offset-halving gate and the 5° p95 cross-resolution gate test the
  physical realization. They stay in force and never absorb GPU error.

## Open items

* Fitted production vectors are not in the fixtures yet: the only finished fits have zero
  free parameters. Add them once the L5 scan, or a later fit, has finished and been
  inspected.
* The full cliff-mask equivalence (C05) needs a *correct* CPU reference template bundle. The
  current builder's HEAL step corrupts gravity (`GRAVITY_REFERENCE.md` §5.1).
