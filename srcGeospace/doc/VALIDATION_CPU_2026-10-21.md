# CPU validation, 2026-10-21/22: srcEarth reference vs batched engine (double, float)

> Historical report.  It was written while srcGeospace was still called `srcEarth_gpu`
> (with the application sources in `srcEarth_gpu/app/`) and its tests were run with the
> srcEarth test runner.  Executables, stage names and run roots below are those of that
> time.  The run roots are kept in
> `$GEOSPACE_RUNS/results/2026-10-21_cpu_validation/` (`../AMPS_runs/srcGeospace/...`
> next to the repository), with a README; the executables were deleted.  The tools
> named here are now in `test/tools/`.  The validation of the renamed tree is in
> section 7.

Date: 2026-10-21/22.  Machine: Intel Core i9-10980XE (18 cores / 36 threads), 251 GB RAM,
g++ with `-O3` (the AMPS defaults), OpenMPI.  No GPU was used: the GPU-built executable
was linked but never started.

Test list: `list_cpu`, the 61 tests of `srcEarth/test/list` without SA01–SA03.  Each
test was run by the unchanged runner `srcEarth/test/test_runner.py -j 1`, from the run
roots described below.  `srcEarth/` and the repository build were not
modified.

| root | executable | sources |
|---|---|---|
| `ref/` | `build/ref/amps` | unmodified srcEarth at `a29ec499` (`SOURCE_REVISION`), built in a private stage |
| `double/` | `build/app/amps` | `srcEarth_gpu/app`, batched engine, Real=double, SI units (default build) |
| `float/` | `build/app_float/amps` | same, `-DGEOSPACE_PRECISION_FLOAT` (mixed precision, normalized units; `PRECISION.md`), with the exit-inset fix of 4.5 |
| `float_rerun019/` | `build/app_float/amps` | test #019 (C2 3d MESH) alone: in `float/` it was killed by the out-of-memory killer while an extra C19 run occupied the machine |
| `float_v1/` | float executable before the fix of 4.5 | first float pass; C19 aborted, all other results as in `float/` |
| `fma/` | `build/app_fma/amps` | same as `double/` plus `-mfma -ffp-contract=fast`: a rounding-level baseline, see 4.2.  26 tests only |
| `nondeterminism_checks/` | `build/ref/amps`, `build/app/amps` | `ref_rerun/`: second run of the reference F11 test; `f5rep_old1..3/`: three reruns of F5 with the double executable |
| `double_check/` | `build/app_check/amps` | double build from the sources after the float exit-inset fix (4.5); 5 tests |

## 1. Summary

| root | tests | as expected by the list | not as expected |
|---|---|---|---|
| ref | 61 | 59 | #011 UStep12Release, #053 C10 (both pre-existing, section 2) |
| double | 61 | 59 | the same two |
| float | 61 | 54 | #011, #053 (pre-existing); #015 C1 3d MESH, #019 C2 3d MESH (from `float_rerun019/`), #051 C9 FULL_SCAN, #052 C9 DIRECT_ACCESS, #056 C17 |
| fma (subset) | 26 | 24 | C2 3d MESH, C9 DIRECT_ACCESS: the same criteria fail as in float |

* **Double build**: every test has the same status as the reference.  The outputs are
  identical except for differences that also occur between two runs of the reference
  executable, or that are timings (section 3).  The in-situ per-trajectory check
  (`AMPS_EARTH_ENGINE_VERIFY=1`) found no bitwise mismatch: about 943,000 gridless
  trajectories (C1, F1) and 8,914 Mode3D MESH trajectories (C1).  The Mode3D mesh
  self-check matched AMPS on 139,785 nodes and 26.4M cells.
* **Float build**: 5 additional test failures.  The double+FMA build, which changes only
  floating-point rounding, fails 3 of them (C2 3d MESH, both C9) with the same
  criteria.  C1 3d MESH (Rc vs Stormer) and C17 (one reflection pair over a 1 %
  threshold) fail only in float (section 4).  The first float pass also aborted C19
  with INVALID_FIELD; that was a float-specific defect of the engine, fixed and
  verified (4.5).

## 2. Pre-existing failures (identical in all roots)

* **#011 UStep12Release**: `test_test_list_provenance_and_approved_structure` compares
  the SHA-256 of the normalized `srcEarth/test/list` with an approved digest
  (`83d3ee94…` vs `2cf62089…`).  It checks the list file, not `amps`.
* **#053 C10 GRIDDED DIRECT_ACCESS**: physics acceptance against POES/MEPED fails in the
  reference (`PRIMARY P6/P7 rmse=3.090 bias=2.671 corr=0.761`).  Double gives the same
  numbers; float gives `rmse=3.091 bias=2.742 corr=0.795`.

## 3. Double build vs reference: outputs

`compare_outputs.py ref double` (plots, logs and JSON skipped; `double_vs_ref.csv`):

| status | files | explanation |
|---|---|---|
| IDENTICAL | 1741 | byte-for-byte equal |
| IDENTICAL_AFTER_PATHS | 21 | equal after replacing the run-root path (manifests, generated decks) |
| NONDETERMINISTIC | 198 | AMPS Mode3D mesh cache files (`mesh->msh`, `*.mesh->bin`); they also differ between two runs of the reference executable |
| NUMERIC | 28 | F11_gridless, relative differences <= 4.1e-16 (see below) |
| NUMERIC | 2 | C13_dipole_nt4 `C13_rank_distribution.csv`, `C13_summary.csv`: timing columns only (see below) |
| DIFFERENT / MISSING | 0 | |

* **F11**: the DYNAMIC MPI scheduler hands out chunks in arrival order, so which rank
  sums which tasks, and hence the last bit of the MPI reduction, depends on timing.  A
  second run of the reference executable (`ref_rerun/`) differs from the first in 24
  F11 files with the same magnitude (max relative 4.1e-16 in data files; the
  `abs_error`/`rel_error` columns of `F11_summary.csv` go from 0 to ~3e-16).  The set of
  affected files changes from run to run.  F5 uses the same scheduler: it was identical
  in this pass, but three reruns of the double executable differ from the pass in 2, 4
  and 6 F5 files (max relative 4.4e-16).
* **C13**: only `elapsed_seconds`, `speedup_vs_one_rank`, `parallel_efficiency`,
  `elapsed_ratio_to_best_same_rank`, `task_imbalance_max_over_avg` (summary) and
  `tasks`, `task_fraction` per rank (DYNAMIC scheduler distribution) differ.  The
  physics columns (`numeric_mismatches`, `flag_mismatches`, `max_relative_difference`,
  Stormer fractions) are identical.

## 4. Float build

### 4.1 Test status

| # | test | float | criterion |
|---|---|---|---|
| 015 | C1 3d MESH | F | `max_abs_rel` of Rc vs Stormer 12.2 % at 9000 km/-30 deg and 18.9 % at 9000 km/+30 deg (tolerance 10 %); reference max over the deck 6.0 % |
| 019 | C2 3d MESH | F | longitude and north/south symmetry at +-30 deg: spread/ref 12 % (tolerance 1 %); `float_rerun019/`, same result and outputs as `float_v1/` |
| 051 | C9 FULL_SCAN | F | time of the minimum cutoff latitude 188 min off (tolerance 100 min); RMSE 2.032 deg, bias +1.908 deg, corr 0.964 pass (reference 2.042, +1.910, 0.962) |
| 052 | C9 DIRECT_ACCESS | F | same as 051 |
| 056 | C17 BORIS | F | charge-sign/reflection symmetry: 2 of 130 stable pairs mismatch at point 3 (lon 0, lat 30), 1.5 % against a 1 % tolerance; reference and FMA: 1 of 130 (0.8 %) |

All other tests pass in float as in the reference, including C19 (T05, Mode3D MESH, GOES
channels).  Its channel metrics are within about 2 % of double (ALL P5: RMSE 0.0972 vs
0.1000, corr 0.885 vs 0.888) and every channel has the same pass/fail as the reference.

### 4.2 Rounding-level baseline (double + FMA)

The penumbral cutoff is discontinuous in the trajectory details: any change of
rounding, not only float, moves some cutoffs by a scan step or by a whole forbidden
band (`PRECISION.md` 5.1).  To separate the float effect from that sensitivity,
the double build was recompiled with `-mfma -ffp-contract=fast`.  That executable is
the same double/SI code with different, equally valid rounding (fused multiply-adds;
4943 FMA instructions in the binary).  It also approximates a GPU double build made
with `build_gpu_engine.sh --fma`, i.e. without `nofma`.  It was run on the 24 tests where float outputs differ (`list_rounding_subset`)
and on C9 DIRECT_ACCESS and C17:

* C2 3d MESH fails the same symmetry check (spread/ref 9.2 % at +-30 deg).  The
  reference passes it only because no longitude happens to resolve the forbidden band
  at 1.58–1.62 GV that float and FMA find at some longitudes.
* C9 DIRECT_ACCESS fails the same criterion with the same 188 min error.  The lowest-rigidity
  (0.4236 GV) cutoff latitude has two nearly equal minima, 57.958 deg (23:55) and
  57.765 deg (03:03) in the reference.  FMA gives 57.895/57.956 and float 57.933/58.060,
  which swaps the argmin.  Per-interval changes from the reference are 0.03–0.49 deg (FMA)
  and 0.03–0.40 deg (float), against a 0.19 deg margin between the two minima.
* C1 3d MESH passes with FMA (max 6.0 %, the reference value).
* C2: FMA finds the band at one longitude (1.584 GV at -30 deg), float at several
  (1.583–1.616 GV).
* C17 passes with FMA (1 of 130 stable pairs at point 3, as the reference).

### 4.3 Per-variable differences from the reference

`variable_stats.py ref float` and `variable_stats.py ref fma`: fraction of values
identical as printed (7 significant digits), and the 90th percentile and maximum of
|a-b|/max(|a|,|b|) over the values that differ.  The cutoff, density and flux values of
`float/` are identical to those of `float_v1/`: the fix of 4.5 changes only where the
exit pitch angle is sampled.  C2 3d MESH is from `float_rerun019/`.

Cutoff rigidity `Rc_GV` on the shell grids:

| test | file | values | float: identical | float: p90 / max rel | FMA double: identical | FMA double: p90 / max rel |
|---|---|---|---|---|---|---|
| C1_gridless | cutoff_gridless_shells.dat | 168 | 0.65 | 4.2e-02 / 4.2e-02 | 0.85 | 2.4e-02 / 2.4e-02 |
| C1_3d_analytic | cutoff_3d_shells.dat | 168 | 0.64 | 4.2e-02 / 4.2e-02 | 0.86 | 3.8e-02 / 3.8e-02 |
| C1_3d_mesh | cutoff_3d_shells.dat | 168 | 0.60 | 7.7e-02 / 1.6e-01 | 0.88 | 2.0e-02 / 2.5e-02 |
| C2_gridless | cutoff_gridless_shells.dat | 84 | 0.62 | 5.9e-03 / 5.9e-03 | 1.00 | 0 / 0 |
| C2_3d_analytic | cutoff_3d_shells.dat | 84 | 0.90 | 5.9e-03 / 5.9e-03 | 1.00 | 0 / 0 |
| C2_3d_mesh | cutoff_3d_shells.dat | 84 | 0.51 | 8.3e-02 / 1.0e-01 | 0.89 | 6.5e-03 / 8.4e-02 |
| C3_gridless | cutoff_gridless_shells.dat | 168 | 0.75 | 7.3e-01 / 7.3e-01 | 0.93 | 2.5e-03 / 2.5e-03 |
| C3_3d_analytic | cutoff_3d_shells.dat | 84 | 0.76 | 6.0e-02 / 6.0e-02 | 0.93 | 2.5e-03 / 2.5e-03 |
| C3_3d_mesh | cutoff_3d_shells.dat | 168 | 0.50 | 8.6e-02 / 7.4e-01 | 0.46 | 1.4e-01 / 7.4e-01 |
| C6_initial | cutoff_gridless_shells.dat | 444 | 0.65 | 7.6e-02 / 2.4e-01 | 0.90 | 3.7e-02 / 9.4e-02 |
| C6_complete | cutoff_gridless_shells.dat | 684 | 0.66 | 5.6e-02 / 2.1e-01 | 0.91 | 3.8e-02 / 7.4e-02 |
| C6_modern | cutoff_gridless_shells.dat | 888 | 0.70 | 5.8e-02 / 2.3e-01 | 0.90 | 3.2e-02 / 1.2e-01 |
| C6_initial_gridded | cutoff_3d_shells.dat | 444 | 0.43 | 2.6e-02 / 1.5e-01 | 0.67 | 1.6e-02 / 3.1e-02 |
| C6_complete_gridded | cutoff_3d_shells.dat | 684 | 0.50 | 3.8e-02 / 1.3e-01 | 0.76 | 2.4e-02 / 6.8e-02 |
| C6_modern_gridded | cutoff_3d_shells.dat | 888 | 0.44 | 3.8e-02 / 1.9e-01 | 0.74 | 3.0e-02 / 1.0e-01 |
| C7 | cutoff_3d_shells.dat | 52 | 0.60 | 3.5e-02 / 1.9e-01 | 0.75 | 4.4e-02 / 7.5e-02 |
| C11_mode3d | cutoff_3d_shells.dat | 336 | 0.80 | 1.0e-01 / 1.0e-01 | 0.92 | 6.5e-03 / 6.5e-03 |
| C14_cross_solver | cutoff_gridless_shells.dat | 168 | 0.75 | 2.0e-02 / 5.9e-02 | 0.88 | 6.1e-02 / 6.3e-02 |
| C14_cross_solver | cutoff_3d_shells.dat | 168 | 0.64 | 5.7e-02 / 1.1e-01 | 0.86 | 2.0e-02 / 7.3e-02 |

Density and omnidirectional flux at the F-test points:

| test | variable | values | float: identical | float: p90 / max rel | FMA double: identical | FMA double: p90 / max rel |
|---|---|---|---|---|---|---|
| F3_gridless | N_m^-3 | 1 | 0.00 | 1.5e-02 / 1.5e-02 | 0.00 | 2.9e-02 / 2.9e-02 |
| F3_gridless | F_tot_m2s1 | 1 | 0.00 | 1.4e-02 / 1.4e-02 | 0.00 | 2.4e-02 / 2.4e-02 |
| F4_gridless | N_m^-3 | 5 | 0.00 | 1.5e-01 / 1.5e-01 | 0.00 | 1.1e-01 / 1.1e-01 |
| F4_gridless | F_tot_m2s1 | 5 | 0.00 | 1.5e-01 / 1.5e-01 | 0.00 | 1.1e-01 / 1.1e-01 |
| F5_gridless | N_m^-3 | 9 | 0.22 | 2.7e-05 / 3.4e-05 | 0.11 | 2.7e-05 / 3.4e-05 |
| F5_gridless | F_tot_m2s1 | 9 | 0.22 | 2.7e-05 / 3.4e-05 | 0.00 | 2.7e-05 / 3.4e-05 |
| F11_gridless | N_m^-3 | 32 | 0.12 | 5.5e-04 / 9.5e-04 | 0.12 | 4.6e-09 / 1.2e-08 |
| F11_gridless | F_tot_m2s1 | 32 | 0.12 | 5.5e-04 / 9.5e-04 | 0.22 | 1.0e-08 / 1.2e-08 |

* Cutoffs: in every deck except C3 3d MESH, float changes more grid points than the
  FMA baseline.  The changed fraction is 0.10–0.57 for float and 0–0.33 for FMA (C3 3d
  MESH: 0.50 and 0.54), about 2–3 times as many in most decks.  The largest change per
  deck is 0.6–24 % for float and 0–12 % for FMA (ignoring the C3 outliers below).  Both are dominated by
  penumbral points, where the cutoff jumps by a scan step or by a forbidden band.  C3 gridless max 73 % is one BINARY
  search point (30 deg E, 60 deg S) where the reference returns the scan floor 0.0433 GV
  and float returns 0.1625 GV (Stormer 0.160 GV).  C3 3d MESH has a 74 % point in both float and
  FMA.
* Densities and fluxes from many trajectories: float differs by up to 1e-3 relative
  (F11), against 1e-8 for FMA.  Where few penumbral trajectories dominate the result (F3,
  F4, F5), FMA and float differ from the reference by the same amount.
* `compare_outputs.py ref float` (`float_vs_ref.csv`) also reports 29 DIFFERENT files in
  tests that pass.  They contain classification flags that flip (True/False,
  UNRESOLVED/FORBIDDEN in C12, C14, C3, C16) or lists whose length changes (boundary
  points in C10, pairwise states in C17, directional cutoffs in C19).

### 4.4 Assessment of single precision (task 1)

* Range: with the normalized units (L=Re, v=c, m=amu, q=e) float runs every deck,
  gridless and Mode3D, without overflow, underflow or INVALID_FIELD aborts, after
  the exit-inset fix (4.5).  One related float limitation remains open: the
  outer-box event search (`FindFirstEventR`) uses the 1 m boundary-event tolerance for
  its on-face test.  A step that leaves the box within about one float spacing (24 m
  at 40 Re) of a box edge could miss the event and end with INVALID_FIELD (reason 4)
  on the next iteration.  It did not occur on these decks (`PRECISION.md` 6.1).  In SI
  float, p^2 of cutoff protons and of all electrons is subnormal (`PRECISION.md` 2).
* Accuracy: mixed precision (float state, field and mover; double time, path, trap
  detector, exit state) gives results of the same kind as a rounding-level change of the
  double build.  Its scatter is 2–3 times larger for cutoffs.  For densities and fluxes
  averaged over many trajectories, float differs by up to 1e-3 relative where FMA
  differs by 1e-8 (F11).  On these decks it costs two acceptance
  tests that the FMA baseline passes (C1 3d MESH; C17 by one of 130 pairs); three
  other float failures are reproduced by the baseline.
* The double build remains the reference.  A float build is suitable for production runs whose acceptance
  criteria are statistical (RMSE, bias, correlation, percentiles), not for
  regression comparisons that expect bit-level or 1 %-level agreement at penumbral
  points.

### 4.5 Float C19 abort and the exit-inset fix

The first float pass (`float_v1/`) aborted C19 (T05, Mode3D MESH, 40 Re box): every rank
stopped with INVALID_FIELD at R = 0.168 GV.  With `AMPS_EARTH_DEBUG_INVALID=1` the
check was reason 10.  The trajectory crossed the outer box, but the exit pitch angle
could not be completed because B at the exit sample point was zero, and no field
sample had been flagged unavailable.  srcEarth samples the exit field max(1 m,
boundary tolerance) inside the event point.  Near 40 Re the float spacing of a
position in Re units is 3.8e-6 Re (24 m), so the 1 m inset is lost and the sample can
round onto or past the face, where the mesh has no leaf.

Fix (`gpu/EngineStep.h`, `Tracer::PopulateExit`): in the float build the inset is
at least 64 ulps of the event position (about 2 km at 40 Re, against 1 Re boundary
cells).  The double build keeps the srcEarth expression:

* Compiled before and after the fix, `EngineDevice.o` and `EngineHost.o` of the double
  build contain the same arithmetic instructions.  Only register moves, the order of
  independent instructions and alignment padding differ, which does not change IEEE
  results without FMA.
* The double executable built after the fix (`double_check/`) passes C1 gridless, C4,
  F1, F5 and C16.  Its outputs equal the double pass, except for F5 at <= 4.4e-16, which
  the unchanged double executable also shows from run to run (section 3).

All float results in this document are from the second float pass (`float/`), made
with the fixed executable, unless they are marked `float_v1`.

## 5. Run times

`run_timing.sh` ran `list_timing` (12 tests: gridless cutoff and density, directional
transmission, the C8/C12/C13/C17 campaigns, Mode3D MESH and a gridded C6 deck) for the
three executables in two rounds of opposite order (ref, double, float; float, double,
ref), 2026-10-22 02:17-03:01.  The machine was not idle: two of the user's
LiouvilleCOPS processes used about 2 cores and 53 GB throughout, and the 1-minute load
was 15-53 (the tests themselves start 64 threads on 36 hardware threads;
`timing_machine_state.txt`).  The two reference rounds agree within 7 %; the double
and float rounds differ by up to 33 % (F3) and 23 % (C6 float), so ratios are good to
about 10-30 %.  The table excludes the 0.1 s
`F4/tests/run_self_tests.py`; the total includes it.

| test | ref s (rounds 1/2) | double s | float s | double/ref | float/ref |
|---|---|---|---|---|---|
| `C1/run_C1.py --mode 3d --mode3d-field-eval MESH -np 4 -nt 1...` | 90.6 / 87.1 | 89.7 / 98.6 | 90.3 / 93.8 | 1.06 | 1.04 |
| `C1/run_C1.py --mode gridless -np 4 -nt 16` | 2.0 / 2.1 | 1.7 / 1.7 | 1.7 / 1.9 | 0.84 | 0.89 |
| `C6/run_C6.py --subtest INITIAL --solver GRIDDED --mode3d-me...` | 166.0 / 154.5 | 147.1 / 143.6 | 146.5 / 180.2 | 0.91 | 1.02 |
| `F1/run_F1.py -np 4 -nt 16` | 3.7 / 3.6 | 3.5 / 3.2 | 3.4 / 3.6 | 0.91 | 0.95 |
| `F3/run_F3.py --mover HC4 --lons 0 --lats 45 --scan-n 24 --d...` | 58.4 / 55.5 | 16.5 / 22.0 | 18.2 / 23.7 | 0.34 | 0.37 |
| `F4/run_F4.py -np 4 -nt 16` | 14.3 / 13.2 | 2.7 / 2.7 | 2.6 / 3.3 | 0.20 | 0.21 |
| `F11/run_F11.py -np 4 -nt 16` | 24.9 / 24.3 | 19.8 / 19.6 | 19.7 / 21.3 | 0.80 | 0.83 |
| `C17/run_C17.py --amps ./amps --movers BORIS -np 4 -nt 16` | 18.8 / 18.3 | 5.5 / 5.4 | 5.2 / 5.5 | 0.29 | 0.29 |
| `C13/run_C13.py --profile ROUTINE --models DIPOLE --amps ./a...` | 53.0 / 53.1 | 16.3 / 16.1 | 15.8 / 15.4 | 0.31 | 0.29 |
| `C12/run_C12.py --profile ROUTINE --amps ./amps -np 4 -nt 16` | 35.1 / 34.6 | 14.4 / 14.2 | 17.2 / 17.3 | 0.41 | 0.50 |
| `C8/run_C8.py --profile ROUTINE --amps ./amps -np 4 -nt 16` | 86.1 / 84.6 | 52.7 / 51.3 | 60.8 / 58.4 | 0.61 | 0.70 |
| **total** | 553 / 531 | 370 / 379 | 382 / 424 | 0.69 | 0.74 |

* Both engine builds use the same MPI ranks and host threads (`-np`, `-nt`) as srcEarth.
* Gridless workloads with many directions per point are 2.4-5 times faster on the CPU
  (F3, F4, C13, C17, C12); C8 and F11 1.25-1.6 times; short decks (C1 gridless, F1)
  1.1-1.2 times.  The cause was not profiled per test.
* Mode3D MESH decks (C1 3d MESH, C6 gridded) take the same time as srcEarth: the engine
  does not change mesh generation or field sampling, and C1 3d MESH traces only 8,914
  trajectories.
* Float is not faster than double on the CPU (C12 and C8 are 15-21 % slower).  Its purpose is the
  GPU, where FP32 throughput is much higher (`PRECISION.md` 6); GPU timings have not
  been measured.
* Full-pass totals (61 tests, run under varying load, not a benchmark): reference
  13,598 s, double 11,764 s.  The float passes are not comparable: `float_v1` stopped
  C19 early, and `float/` shared the machine with an extra C19 run.

## 6. Reproduce

Since 2026-10-27 srcGeospace has no copy of the tests: they are srcEarth's, run with
srcEarth's runner through `test/tools/run_srcearth_tests.sh` (`../test/README.md`).  The
equivalent commands are:

```
R=$GEOSPACE_RUNS/results/<run>
T=srcGeospace/test/tools
$T/run_srcearth_tests.sh --run-root $R/double --amps $GEOSPACE_RUNS/builds/cpu_double/amps --exclude 'test/SA0[123]/' -- -j 1
$T/run_srcearth_tests.sh --run-root $R/float  --amps $GEOSPACE_RUNS/builds/cpu_float/amps  --exclude 'test/SA0[123]/' -- -j 1
$T/compare_outputs.py $R/ref $R/double --out $R/double_vs_ref.csv     # ref: a srcEarth build
$T/variable_stats.py $R/ref $R/float
$T/summarize.py --root ref=$R/ref --root double=$R/double --root float=$R/float
srcGeospace/scripts/build_amps.sh -j 8 --stage $GEOSPACE_RUNS/builds/cpu_fma --extra-cxx "-mfma -ffp-contract=fast"
$T/run_srcearth_tests.sh --run-root $R/fma --amps $GEOSPACE_RUNS/builds/cpu_fma/amps --list $T/list_rounding_subset -- -j 1
$T/run_timing.sh $R ref=<srcEarth amps> double=$GEOSPACE_RUNS/builds/cpu_double/amps float=$GEOSPACE_RUNS/builds/cpu_float/amps
```

The reference executable of this report was srcEarth itself, built in a private stage
by the earlier version of `build_amps.sh` (the repository's Earth deck with
`ProjectSourceDirectory` pointing at a copy of srcEarth).  The current script builds
srcGeospace from `srcGeospace/input/`.

## 7. Validation of the renamed tree (srcGeospace)

SRCGEOSPACE_VALIDATION_PENDING
