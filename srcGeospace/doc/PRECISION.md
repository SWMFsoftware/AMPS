# Precision and units of the batched trajectory engine

This note records why the engine's arithmetic type is a compile-time choice, which
unit system each choice implies, and the measurements behind it.  The policy below
is what `gpu/EngineConfig.h` implements now.  The data come from
`tools/precision_experiment` (rerun 2026-10-21).

## 1. Policy (implemented)

| quantity | default build | `-DGEOSPACE_PRECISION_FLOAT` build |
|---|---|---|
| trajectory state `x`, `p` (`EgSlot<R>::x,p`) | double, SI | float, normalized geospace units |
| field sample `B`, mesh field array, row-stencil weights | double, SI | float, normalized units |
| mover arithmetic (`Movers<R,U>`), boundary event search (`FindFirstEventR`) | double | float |
| trace time `EgSlot::t`, path length `EgSlot::s` | double | **double** (normalized units) |
| time step `dt` returned by the selectors | double | **double** (selector arithmetic in float, result converted) |
| trap detector (`EgTrapDetector`, profiles, arena) | double, SI | **double, SI** (state converted on every update) |
| Shue magnetopause classification | double, SI | **double, SI** |
| exit kinematics, `EgResult`, `EgRequest`, `EgContext` | double, SI | **double, SI** |
| dipole-error diagnostic (`RecordDipoleError`) | double, SI | **double, SI** |

* The precision is selected at compile time by `GEOSPACE_PRECISION_FLOAT`
  (`typedef ... Real`, `EG_REAL_NAME`); the unit system follows from it
  (`Units = UnitsSI` or `UnitsGeo`).  There is no run-time switch.
* In the default build all unit conversions are `if constexpr` identities, and every
  expression is a literal port of the legacy tracer, so the engine is designed to be
  bitwise identical to it with the same compiler and flags.  `AMPS_EARTH_ENGINE_VERIFY=1`
  checks this at run time, trajectory by trajectory (`DESIGN.md`).
* The float build is not bitwise comparable with the legacy tracer.  It must be validated
  statistically (cutoffs, transmission functions, trap statistics), as in section 5.
* The trap detector, time and path accumulators stay double because they integrate
  over 1e5..1e7 steps (sections 5.3 and 5.4).
* Host-side quantities passed through the C API (`EgRequest`, `EgContext`,
  `EgResult`, `EgMeshHost`) are always SI double, so the library and the driver
  exchange the same POD types in both builds, and a precision mismatch would not break
  the C ABI.  The g++ side and the nvc++ library must still be compiled with the same
  `GEOSPACE_PRECISION_FLOAT` setting.  Otherwise the host pool and `TraceSync()` (g++
  code, used for host-field models and for requests made outside fibers, including
  NONE/DIPOLE/mesh requests) would integrate in a different precision from the device
  pool.  `scripts/build_amps.sh --gpu --float` passes the define to both sides, and
  `eg_dev_backend_name()` reports the library's precision.

## 2. Why SI single precision fails

Float: FLT_MIN = 1.18e-38 (smallest normal number), smallest subnormal 1.4e-45,
FLT_MAX = 3.4e38, 24-bit significand (relative spacing 6e-8).

| quantity (SI) | value | relative to FLT_MIN |
|---|---|---|
| p^2, 1 MeV proton | 5.4e-40 | 0.046 (subnormal: about 18 significant bits) |
| p^2, 10 MeV proton | 5.4e-39 | 0.46 (subnormal) |
| p^2, 1 GeV proton | 8.2e-37 | 70 |
| (m_p c)^2 | 2.5e-37 | 21 (normal, but with no headroom for further products) |
| (m_e c)^2 | 7.5e-44 | 6e-6 (subnormal, about 6 significant bits) |
| p^2, 1 MeV electron | 5.8e-43 | 5e-5 (subnormal) |
| p^2, 10 keV electron | 2.9e-45 | rounds to 2 x 1.4e-45 (precision lost) |

`GammaFromMomentum` evaluates `p2/(mc*mc)` and every mover evaluates `vdot(p,p)`.
In SI float these are subnormal for typical cutoff protons and for all electrons.
Under flush-to-zero arithmetic (for example nvcc `-ftz=true`, implied by
`--use_fast_math`) they become exactly zero.

## 3. Why the AMPS PIC::Units normalization was rejected

`src/pic/units/units.cpp` builds the AMPS default normalization with
`mSI = _AMU_`, `lSI = 1e-7 e^2/mSI` and `uSI = 1e6 m/s` (the choice that makes
`mass_conv = 1/amu` and `charge_conv = 1/e`); `pic_units_normalization.h` derives
`rho0 = M0/L0^3`, `B0 = sqrt(rho0) U0`, `T0 = L0/U0` (CGS). Evaluated:

| scale | value |
|---|---|
| lSI | 1.546e-18 m |
| B0 | 6.70e15 T |
| T0 | 1.546e-24 s |
| 60 Re | 2.47e26 units, so \|x\|^2 = 6.1e52 > FLT_MAX: **overflow** |
| 1 nT | 1.49e-25 units, so \|B\|^2 = 2.2e-50 < smallest subnormal: **underflow to 0** |

That normalization suits the field solver's plasma scales.  For magnetospheric
trajectories it overflows and underflows single precision itself.

## 4. Normalized geospace units (`UnitsGeo`, used when Real=float)

| unit | definition | value |
|---|---|---|
| length | `_EARTH__RADIUS_` | 6371.0e3 m |
| speed | c | 299792458 m/s (so c = 1) |
| mass | `_AMU_` | 1.66053906660e-27 kg |
| charge | AMPS `ElectronCharge` | 1.602176565e-19 C |
| time | Re/c | 0.021251 s |
| momentum | amu c | 4.978e-19 kg m/s |
| magnetic field | amu c/(e Re) | 4.877e-7 T (about 488 nT) |

Typical magnitudes are then O(1e-4..1e5): positions 1..60, the surface equatorial
field 3.12e-5 T = 64.0, 1 nT = 2.05e-3, a 1 GV proton p = 1.07, a 10 keV electron
p = 1.1e-4, (m_e c)^2 = 3.0e-7, 600 s = 2.8e4, a 1e-5 s step = 4.7e-4.  This fixes the
*range* problem completely.  It does not change the 24-bit significand: a position of
35 Re carries an absolute rounding of 35 Re x 2^-24 = 13 m per operation, and every
step of a long trajectory adds rounding at that level.  Section 5 measures the effect.

Conversion details (as implemented):

* Requests carry SI double; `Tracer::MakeConsts` converts charge, mass, time budget,
  path budget, DT_TRACE, box and inner radius once per launch.
* The initial momentum is computed in SI double
  (`|p| = R_GV 1e9 |q|/c`) and then converted, so it is the same number the legacy tracer uses,
  rounded once.
* The dipole moment is pre-scaled on the host into internal units
  (`cm = (mu0/4pi) M /(B_u L_u^3)`) and the dipole formula is evaluated entirely in
  `R` (no double arithmetic in the float kernel).
* The Mode3D mesh coordinates and B array are converted once when the mesh is uploaded
  (`ConvertCopy` / `HostMesh::Set`).
* Host fields (Geopack/Tsyganenko) are evaluated in SI double and the sample is
  converted to `R`.
* The legacy "no limit" sentinels (1e300/1e299) are replaced by 1e30f/1e29f in float
  (`Sentinel<float>`).
* The experiment below used m_p instead of amu as the mass unit (B_u = 4.9125e-7 T).
  The factor 1.0073 between the two does not affect any conclusion.

## 5. Experiment

Setup (details in `tools/precision_experiment/README.md`): centered dipole,
B_eq = 3.12e-5 T, box +-35 Re, inner sphere 1 Re, 600 s, DT_TRACE 1 s, gyro-angle limit
0.15 rad, per-step adaptive dt; protons.  Variants: **DD** all double; **DDp** all double
with every dt scaled by (1+1e-7), which calibrates how sensitive the answer is to a
perturbation at the double rounding level; **FD** float state and field, double
accumulators (the engine's float build); **XD** float momentum and field but double
position; **FF** everything float, including time and path accumulators.

The rerun on 2026-10-21 reproduced every number of the design-phase run
(`original/` and `2026-10-21_rerun/` in `$GEOSPACE_RUNS/results/precision_experiment/`,
outside srcGeospace) digit for digit in sections A-D.  Section E did not exist in the
earlier build.

### 5.1 Vertical cutoffs (UPPER_SCAN)

Relative difference from DD; points not listed are identical in all variants.

| mover | alt km | lat | Stormer GV | DD GV | DDp | FD | XD | FF |
|---|---|---|---|---|---|---|---|---|
| BORIS | 500 | 60 | 0.8005 | 0.6543 | +15.3 % | +1.95 % | +14.4 % | +1.95 % |
| BORIS | 9000 | 30 | 1.4397 | 1.1107 | +0.79 % | **+20.3 %** | +19.9 % | +20.3 % |
| BORIS | 9000 | 45 | 0.6398 | 0.5809 | 0 | +0.97 % | +2.72 % | +0.97 % |
| BORIS | 9000 | 60 | 0.1600 | 0.1493 | +1.15 % | -1.53 % | 0 | -1.53 % |
| HC4 | 500 | 45 | 3.2021 | 3.6337 | **-11.5 %** | -11.5 % | -11.6 % | -11.5 % |
| HC4 | 500 | 60 | 0.8005 | 0.8145 | -1.23 % | -0.28 % | +0.10 % | -0.28 % |
| HC4 | 9000 | 30 | 1.4397 | 1.4521 | -0.30 % | -0.35 % | -0.15 % | -0.35 % |
| HC4 | 9000 | 45 | 0.6398 | 0.6492 | +0.69 % | +0.98 % | -0.10 % | +0.98 % |
| HC4 | 9000 | 60 | 0.1600 | 0.1672 | 0 | -0.39 % | -1.55 % | -0.39 % |
| RK4 | 500 | 45 | 3.2021 | 3.2125 | -0.09 % | 0 | 0 | 0 |
| RK4 | 500 | 60 | 0.8005 | 0.8091 | +0.10 % | 0 | 0 | 0 |
| RK4 | 9000 | 30 | 1.4397 | 1.4521 | -0.30 % | 0 | 0 | 0 |
| RK4 | 9000 | 45 | 0.6398 | 0.6505 | -0.59 % | 0 | +0.39 % | 0 |
| RK4 | 9000 | 60 | 0.1600 | 0.1640 | 0 | 0 | +0.40 % | 0 |

* The equatorial and 30 deg/500 km points have a sharp cutoff: all variants agree.
* At penumbral points UPPER_SCAN returns the edge of the highest forbidden band, which
  is discontinuous in the trajectory details.  A 1e-7 relative change of dt in double
  (DDp) moves the cutoff by up to 15 % (Boris) and 11.5 % (HC4).  The float variants
  move it by up to 20 % (Boris), 11.5 % (HC4) and 0.4 % (RK4).  Float changes are of
  the same kind and size as a double-rounding-level perturbation.  The penumbral
  cutoff is ill-conditioned, not float-specific.
* At HC4 500 km/45 deg, DD is the outlier: DDp, FD, XD and FF all agree at 3.215 GV.
* FD and FF give identical cutoffs: float time/path accumulators do not matter for
  600 s budgets (they do for longer ones, 5.4).
* Boris truncation error at 0.15 rad per step is larger than any rounding effect
  (equatorial 500 km: 13.51 GV vs Stormer 12.81 GV; HC4 and RK4: 12.809 GV).

### 5.2 Directional transmission (45 deg, 9000 km, 12 rigidities x 1152 directions)

| mover | mismatching directions out of 13824 (DDp / FD / XD) | max \|T - T_DD\| (DDp / FD / XD) |
|---|---|---|
| HC4 | 63 / 103 / 108 | 0.0122 / 0.0078 / 0.0122 |
| RK4 | 63 / 102 / 90 | 0.0044 / 0.0070 / 0.0069 |
| BORIS | 562 / 687 / 672 | 0.0147 / 0.0182 / 0.0113 |

* No direction of the DD run hit the time or step limit.
* Float increases the number of reclassified directions by a factor of 1.2-1.7 over
  the double-perturbation baseline (HC4 103 vs 63, RK4 102 vs 63, Boris 687 vs 562).  Every mismatch is in the penumbral
  rigidities (R/Rc = 0.55..1.30).
* Transmission differences are at most 1.8 percentage points (Boris FD) and are
  comparable to DDp for all movers.  Relative to T_DD (bins with T_DD > 0.1) they are
  up to 4 % (Boris FD), 1.0 % (HC4 FD) and 0.9 % (RK4 FD).

### 5.3 Long trapped orbit (10 MeV proton at 4 Re, 1200 s) and |p| conservation

| mover | DD | FD | FF |
|---|---|---|---|
| HC4 | TimeLimit, 960840 steps, dp/p 1.6e-13 | TimeLimit, 960930 steps, dp/p 1.22e-4 | TimeLimit, 960180 steps, dp/p 1.22e-4 |
| RK4 | TimeLimit, 978826 steps, dp/p 4.96e-2 | TimeLimit, 979911 steps, dp/p 4.96e-2 | TimeLimit, 980342 steps, dp/p 4.96e-2 |
| BORIS | Inner, 93212 steps at 95.8 s, dp/p 2.2e-14 | Inner, 87587 steps at 88.7 s, dp/p 1.2e-5 | Inner, 87587 steps at 88.7 s, dp/p 1.2e-5 |

Maximum relative |p| spread over the resolved transmission traces of 5.2: Boris DD
1.6e-14, FD 9.5e-6, XD 9.4e-6; HC4 DD 7.5e-14, FD 3.9e-5, XD 4.2e-5; RK4 DD 2.3e-3,
FD 2.4e-3, XD 2.4e-3.

* The trap detector's momentum gate `TRAP_ENERGY_REL_TOL` defaults to 1e-4
  (`trapEnergyRelativeTolerance`).  Float HC4 reaches 3.9e-5 on 600 s traces and
  1.2e-4 after 1e6 steps, so a float state would make the detector reject long
  trapped orbits that the double build accepts.  This is one reason the detector
  and its |p| bookkeeping stay double.  The detector still sees float-rounded
  momenta, though, so long HC4 runs in the float build should raise the tolerance or
  be checked against the double build.
* RK4 is not norm-conserving: it drifts by 2.3e-3 in 600 s and 5e-2 in 1200 s even in
  double, far above the 1e-4 gate.  That is a property of the mover, not of precision.
* With Boris at 0.15 rad per step this orbit is not kept trapped even in double: DD
  reaches the inner boundary at 95.8 s, FD at 88.7 s.  HC4 and RK4 stay trapped for
  1200 s.

### 5.4 Float time accumulator

2.4e7 additions of dt = 5e-5 s (2.35e-3 units): double reaches 1200.000 s, float reaches
**1362.634 s**.  Near t = 1200 s (5.6e4 units) the float spacing is 3.9e-3 units; each
`t += dt` rounds to a multiple of it, giving a systematic +14 % error.  For smaller
dt, or more time units (SI), the sum stalls completely once dt < ulp(t)/2.  Therefore
`EgSlot::t` and `EgSlot::s` are double in all builds, and the time-step selectors
return double.

### 5.5 Adaptive per-trajectory dt vs. a global fixed dt (section E)

The 13824 DD trajectories of 5.2, for each mover:

| mover | steps (adaptive) | traced time | smallest gyro-limited dt | steps at that fixed dt | ratio |
|---|---|---|---|---|---|
| HC4 | 4.949e7 | 1.194e5 s | 2.820e-5 s | 4.234e9 | 86x |
| RK4 | 4.893e7 | 1.197e5 s | 2.820e-5 s | 4.243e9 | 87x |
| BORIS | 3.496e7 | 1.155e5 s | 2.758e-5 s | 4.187e9 | 120x |

A batch scheme that advances all trajectories with one global dt would have to use the
smallest gyro-limited step of any trajectory: 86-120x more steps on this workload.  The
engine therefore keeps the legacy per-trajectory adaptive dt.  Each slot carries its own
time and step selection, and a kernel launch advances every slot by up to K of its
own steps (`DESIGN.md`).

### 5.6 In-AMPS validation on the test list

The float build ran the 61 CPU tests of the list (without SA01-SA03); the results are
in `VALIDATION_CPU_2026-10-21.md` section 4.  As the in-AMPS counterpart of DDp, the
double build was also recompiled with `-mfma -ffp-contract=fast` (same code and
precision, different rounding) and run on the 24 tests whose float outputs differ, and
on C9 DIRECT_ACCESS and C17.

* Float runs every deck, gridless and Mode3D (dipole, IGRF/Tsyganenko via the host pool,
  Mode3D MESH and ANALYTIC), without range failures.  The first float pass aborted C19
  (T05 mesh, 40 Re box) because the 1 m exit-field inset was below the float spacing;
  after the fix of 6.1 C19 passes with channel metrics within about 2 % of double.
* Float passes 54 of 61 tests (reference and double: 59).  Of the 5 additional
  failures, the FMA double build fails C2 3d MESH and both C9 criteria in the same way:
  a penumbral forbidden band appears at some longitudes, and two nearly equal minima
  188 min apart swap.  C1 3d MESH (Rc vs Stormer, 18.9 % against a 10 % tolerance) and
  C17 (2 instead of 1 of 130 reflection pairs, against a 1 % limit) fail only in
  float.
* Cutoffs: float changes 10-57 % of the shell-grid points, FMA 0-33 % (C3 3d MESH 54 %);
  the largest change per deck is up to 24 % for float and 12 % for FMA.
* Densities and fluxes over many trajectories: up to 1e-3 relative for float, 1e-8 for
  FMA (F11).  Where a few penumbral trajectories dominate (F3-F5), both differ from the
  reference by the same amount.

## 6. Conclusions and status

1. Default build: double/SI, bitwise identical to the legacy tracer.  It is the reference for
   all validation.
2. Float build: the range problem is solved by `UnitsGeo`.  The precision effects are
   of the same kind as a double-rounding-level perturbation for cutoffs and
   transmission (5.1, 5.2) and smaller than the Boris truncation error.  On the
   test decks their scatter is about 2-3 times that of the rounding-level baseline
   (5.6).  Invariant drift on long
   orbits is larger (5.3) and needs care with tight trap tolerances.
3. Time, path length, the trap detector, exit kinematics and the Shue boundary are
   double in every build (5.3, 5.4).
4. Status: the float build compiles with g++ (`build_amps.sh --float`) and with nvc++
   `-acc=gpu`/`-acc=multicore` (`BUILD.md`).  It was run on the CPU test list inside
   AMPS (5.6, `VALIDATION_CPU_2026-10-21.md`): it runs every deck, and its differences from the
   double build are of the kind of 5.1, with about 2-3 times the scatter of a
   rounding-level change of the double build.
5. Performance motivation: on the RTX A6000 (GA102) FP64 throughput is a small fraction
   of FP32.  The double kernel also uses 255 registers and a 6704-byte stack frame per
   thread, against 168 registers and 5296 bytes for float (ptxas, `DESIGN.md` §5).
   Measured since (`../GPU_STATUS.md` §4): GPU float tracing takes 0.67-0.78 of GPU
   double tracing on the Mode3D mesh decks and 0.27 on the analytic-dipole F11 kernel;
   on the CPU float and double take the same time (ratio 1.00, `../STATUS.md` §5.2).

### 6.1 Absolute lengths in the float build

Lengths that the legacy code gives in metres become tiny in Re units, and near the outer box
they are below the float spacing of a position (3.8e-6 Re = 24 m at 32-64 Re):

| length | where | float handling |
|---|---|---|
| exit-field inset max(1 m, boundary tolerance) | `Tracer::PopulateExit` | **fixed**: at least 64 ulps of the event position (C19 aborted before, `VALIDATION_CPU_2026-10-21.md` 4.5) |
| configured-domain faces of the Mode3D mesh field | `DeviceField` | shrunk by 8 ulps (`EngineField.h`) |
| AMR block-face snapping | `EgMesh::AddCellByCoordinate` | 64 eps x nCells (`EngineMesh.h`) |
| boundary-event tolerance (default 1 m) in the box on-face test and the event-fraction window | `FindFirstEventR` | **open**: unchanged.  An exit within about 24 m of a box edge can be missed, and the next iteration then ends with INVALID_FIELD (reason 4).  Not observed on the 61 test decks.  Proposed fix: in the float build use max(tolerance, 64 ulps of the largest box bound) |
| Shue boundary tolerance | `ShueClassify*` | none needed: evaluated in double/SI |

## 7. Reproduce

```
cd srcGeospace/tools/precision_experiment
make && make -j 3 run        # about 2.5 min on this machine (HC4 dominates);
                             # writes $GEOSPACE_RUNS/results/precision_experiment/run
```
