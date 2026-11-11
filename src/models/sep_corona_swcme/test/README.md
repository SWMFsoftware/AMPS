# Corona–SWCME shared-library tests

Updated 2026-11-10. These are dependency-light C++ software/physics-reference
checks. They use real shared implementation and independent equations or
manufactured fixtures. They do not allocate PIC particles, launch a native
AMPS application, or qualify full-CME G1 propagation.

[Model APIs and equations](../README_CROSS_RAY.md) ·
[Global owner and reporting](../../../test/tests/models/cme/sep_corona_swcme_full/README.md) ·
[Complete catalog](../../../test/TEST_CATALOG.md).

## Working directory and quick builds

Run from the AMPS root, then enter the MODEL directory so relative fixture
paths such as `examples/bg3d4_piston/event.conf` resolve correctly:

```sh
cd src/models/sep_corona_swcme
make -j16 build/test_bg3d4_contact_derivative build/test_spherical_shell \
  build/test_cross_ray_transport build/test_spherical_sector \
  build/test_cross_ray_mesh build/test_cross_ray_ct \
  build/test_cross_ray_reconstruction build/test_cross_ray_fluid \
  build/test_cross_ray_contact build/test_material_transport

build/test_bg3d4_contact_derivative
build/test_spherical_shell
build/test_cross_ray_transport
build/test_spherical_sector
build/test_cross_ray_mesh
build/test_cross_ray_ct
build/test_cross_ray_reconstruction
build/test_cross_ray_fluid
build/test_cross_ray_contact
build/test_material_transport
```

The Makefile builds the shared archive and its maintained coronal dependency
through their own Makefiles. Objects stay below each model's `build` directory.
Default flags are C++17, `-O2 -Wall -Wextra -Wpedantic -Werror`.
These commands do not rebuild native AMPS, regenerate application hooks or
serve as a substitute for the mandatory native root/idle/configuration/-j16
preflight. The global adapter copies sources into a private repository-shaped
work tree before invoking these exact targets, preventing stale checkout
objects from supplying a false result.

For the repaired complete contact diagnostic and established ambient references:

```sh
make -j16 build/test_bg3d4 build/test_bg3d4_ambient_balance \
  build/test_bg3d4_spherical_piston build/test_bg3d4_append check-architecture
build/test_bg3d4
build/test_bg3d4_ambient_balance
build/test_bg3d4_spherical_piston
build/test_bg3d4_append
```

The complete CMBGU05 executable is substantial; use the derivative-control
binary for a quick independent grader check. The full diagnostic is still
required in the global owner. `make test` is fail-fast and currently stops at
the unresolved production refinement. To obtain all individual outcomes,
including controls after that failure, use the global owner described below.

## New binaries and independent oracles

Every local marker is globalized by prefix `CSWC-FULL-`. The exact individual
case descriptions are also maintained in the
[global CME adapter README](../../../test/tests/models/cme/README.md).

| Binary / source in this directory | Cases | What, how and why |
| --- | ---: | --- |
| `test_bg3d4_contact_derivative.cpp` | 1 | Shares the actual fourth-order position grader with CMBGU05; independent sphere/plane trajectories test convergence, both signs of physical slip, noncancelling absolute leakage, unit normals and malformed/out-of-support samples. |
| `test_spherical_shell.cpp` | 1 | Independent two-node Gaussian polynomial moments grade broad, thin and one-ulp shells over SI scales through 1 AU, detecting cancellation and invalid geometry without a width floor. |
| `test_cross_ray_transport.cpp` | 5 | Hand-evaluated SI ideal-MHD flux/energy/fast-wave limits; HCI rotations and ALE work; opposite extensive transfers; invalid-state/normal-B/moving-volume rejection; periodic density advection versus an analytic solution at 32/64/128 cells. |
| `test_spherical_sector.cpp` | 3 | Analytic octant volume/centroid/area vectors, proper rotation, constant-field flux closure, exact sphere totals and all 792 interfaces of a 336-cell reference shell; invalid frames/edges and thin geometry. |
| `test_cross_ray_mesh.cpp` | 3 | Unique shared topology, exact integer incidence and shell Euler characteristic; curved radial/meridional/latitude line elements including poles/seam; event/epoch/generation/cell-budget identity and rejected preparation. |
| `test_cross_ray_ct.cpp` | 3 | Independent Gaussian integration of known A and E around actual curved edges checks Stokes/Faraday; twenty common deterministic random edge updates check discrete divergence; deliberate monopoles, stale/foreign data and incomplete intervals reject. |
| `test_cross_ray_reconstruction.cpp` | 3 | Known constant Cartesian fields under rotation/scaling; independent affine divergence-free field from edge-potential integration at 4/8/16 cells per axis; geometric-rank rejection and exactly zero B. |
| `test_cross_ray_fluid.cpp` | 4 | Twelve epochs of exact constant-state integrated MHD flux; analytic opposite lateral/boundary transfers with all momentum components; conserved energy with changed CT field and explicit source work; whole-epoch rejection. The curl-E source control is deliberately manufactured and is not an ideal-MHD production law. |
| `test_cross_ray_contact.cpp` | 4 | Independent quintic startup angle/rate/acceleration, all four actual supported-ray leakage checks, material-label versus fixed-HCI ray derivatives through handoff/late time, threaded-interface Maxwell/electric/work residuals, and excluded-geometry HCI bookkeeping. |
| `test_material_transport.cpp` | 4 | Closed-form affine density/field/entropy and independently deformed signed face flux; independent temporal/Eulerian derivatives; nonlinear lateral map with explicit inverse and global mass/divergence checks; folds/unresolved rank/nonfinite output reject. |

These are 31 new model-local cases: 29 cross-ray/material/contact cases plus
the two independent numerical controls. The existing 15-row owner remains
registered, giving 46 full-volume software rows including architecture.
No API-reported conservation residual is copied into an expected answer.

## Exact marker groups

```text
CMBGU05-DERIVATIVE-CONTROLS
CMBGU05-SHELL-GEOMETRY-CONTROLS
CSWC0628-XRAY-FLUX
CSWC0628-XRAY-ALE
CSWC0628-XRAY-CONSERVATION
CSWC0628-XRAY-REJECTION
CSWC0628-XRAY-ADVECTION
CSWC0628-XRAY-GEOMETRY
CSWC0628-XRAY-MESH-CLOSURE
CSWC0628-XRAY-GEOMETRY-REJECTION
CSWC0628-XRAY-TOPOLOGY
CSWC0628-XRAY-EDGE-GEOMETRY
CSWC0628-XRAY-MESH-IDENTITY
CSWC0628-XRAY-CT-INITIALIZATION
CSWC0628-XRAY-CT-FARADAY
CSWC0628-XRAY-CT-AUTHORITY
CSWC0628-XRAY-RECONSTRUCTION-CONSTANT
CSWC0628-XRAY-RECONSTRUCTION-CONVERGENCE
CSWC0628-XRAY-RECONSTRUCTION-REJECTION
CSWC0628-XRAY-FLUID-UNIFORM
CSWC0628-XRAY-FLUID-CONSERVATION
CSWC0628-XRAY-FLUID-ENERGY
CSWC0628-XRAY-FLUID-REJECTION
CSWC0628-XRAY-CONTACT-STARTUP
CSWC0628-XRAY-CONTACT-DERIVATIVES
CSWC0628-XRAY-CONTACT-INTERFACE
CSWC0628-XRAY-CONTACT-REJECTION
CSWC0628-XRAY-MATERIAL-INVARIANTS
CSWC0628-XRAY-MATERIAL-RATES
CSWC0628-XRAY-MATERIAL-NONLINEAR
CSWC0628-XRAY-MATERIAL-REJECTION
```

Evidence-only lines such as `...-EVIDENCE` or `...-ADVECTION-LEVEL` are not
additional cases. The global parser imports only the exact registered typed
PASS/FAIL/SKIP/ERROR markers. A zero exit without all owned markers is an
ERROR; prose output or a partial grid level cannot invent a missing PASS.

## Acceptance details that must remain intact

- CMBGU05 uses h=0.1/0.05/0.025 s, the original physical leakage bound and
  convergence clause. Its reference is rho*(v-independent dX/dt).n. The
  explicitly user-approved distributed gate requires every level to satisfy
  unchanged bounds. Force/work and cross-stage model discrepancy still appear.
- Stable shell moments do not change equations, masses, viscosity, RK4/CFL or
  source samples. The independent Gaussian-centroid oracle replaces an old
  cancellation-sensitive test expression. CSWC0625 retains its 2e-13 radial
  flux bound and now measures 4.01794e-16 relative error.
- Static mesh/reference CT tests check all cells/interfaces, not a selected
  surviving subset. Arithmetic uncertainty bounds are distinct from raw
  divergence and from production physical tolerances. Invalid states are not
  cleaned, floored or silently relabeled ambient.
- The material kernel accepts a numerically resolved small positive J; its
  rejection test does not weaken the separate production minimum-J gate.
  Global map compatibility, physical force/heating/work, shock entropy and
  native field ownership are not established by a pointwise Cauchy identity.
- The refinement owner retains all five advances, all supported rays, the
  original 3600-s comparison epoch, mass/time controls, shock-branch and
  convergence thresholds, and its 180-second routine execution ceiling.

The refinement binary also has a diagnostic CLI, for example:

```sh
make -j16 build/test_bg3d4_piston_refinement
build/test_bg3d4_piston_refinement 40 0.2 3600
```

That runs ONE level and does not emit the suite PASS marker. Its zero exit is
only completion of a diagnostic. The complete original coarse run takes
471.164 s and detects two zones; a stable-geometry coarse diagnostic takes
127.902 s with the same printed shock values and two zones. Neither qualifies
the unchanged full matrix.

## Global owner, reports and machine paths

From the AMPS root:

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONA-SWCME-FULL-001 \
  --output-dir test_output/cswc-full-check
```

Use a fresh output directory. `--config` is required; another machine should
supply a profile based on `test/qualification-machine/template.json`. Compiler,
make, Python, launcher, output and temporary paths belong in that profile,
not in this test code. The complete registry is retained even when one local
binary fails or times out. `--case` narrows reported rows after exact-registry
validation; this owner still executes its complete local binary registry.
For quick isolated execution use the direct binary commands above.

The owner currently runs exclusively within the core stream and executes
its local binaries sequentially; `--jobs` does not parallelize those binaries.
The global runner's separate stream/build controls remain documented in
[`src/test/README.md`](../../../test/README.md).

## Current evidence and remaining qualification

The complete 2026-11-10 core receipt is under repository-local
`test_output/remaining-failures-repair-20261110/core-current/`:
1,067 cases = 825 PASS / 1 FAIL / 240 SKIP / 1 ERROR, in 1958.293 s.
`core-current-baseline-comparison.json` preserves every previous 792 PASS and
all 32 new global cases PASS. `CMBGU05` is the sole changed baseline result.
The remaining ERROR is `CSWC-FULL-CSWC0628-MASS-TIME`; the remaining FAIL is
`T12-COVERAGE-STRICT-001` with 26 incomplete components.

Focused optimized/debug/all-enabled sanitizer receipts and first failed or
superseded diagnostic attempts remain in the evidence root. Existing linked
archives are uninstrumented dependencies where the receipt says so. The
[comparison plot](../../../../test_output/remaining-failures-repair-20261110/bounded-comparisons.png)
uses recorded contact/advection/reconstruction errors; it does not substitute
for the production refinement or a G1 physical acceptance gate.

A reduced-front PASS, a software-reference PASS, source-only compilation,
initialization or a single diagnostic level cannot qualify spatial full-CME
regions. A production lateral map, actual ambient/thermal/source consistency,
common moving geometry, RH/contact/force/work/validity checks, BG3D-5 matching
and later native one-/four-rank zero-particle propagation remain open.

## Rule for adding or changing a test

Update the local Makefile target, local markers, global adapter `LOCAL_IDS`
and `BINARY_CASES`, central `case_registries.json`, relevant ownership and this
README/global case descriptions in the same change. Preserve exact ordered
membership; document an independent oracle, equations, units, assumptions,
validity, typed failure and retained evidence. Add new cases to the generated
catalog through registration and run `python3 src/test/generate_test_catalog.py`.
The work is incomplete until catalog/traceability checks pass. No test gate,
reference fixture or selected physical problem may be weakened to obtain PASS.
