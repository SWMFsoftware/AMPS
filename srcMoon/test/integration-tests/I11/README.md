# I11 — AMR spatial-resolution convergence

The existing I11 wrapper automatically runs `test_infrastructure.py` before
MPI launch. Five isolated tests cover both approved sequences, fixed
far-field/physical controls, exact finest-input reproduction for I13, and
rejection of unapproved custom settings. They also freeze the 500-step sample,
5000-step per-resolution cap, and disabled particle-restart build setting.
These tests generate only temporary
input files, not observations or science baselines. Infrastructure success
does not contribute a numerical-convergence PASS.

I11 uses the normal SPICE-enabled Moon executable and the reviewed spherical
M1 receipt. The approved 2026-10-14 memory campaign uses
`m1_memory200.in` and immutable 480, 240, 120, and 60 km child receipts.
The legacy FULL template retains its 40/20/10/5 km sequence when explicitly
selected. Other CUSTOM and QUICK meshes are rejected. The runner is
changing only spatial-resolution targets. Each linked run proves that
`simulation_environment_mode=frozen` holds ephemeris time, Sun/Earth states,
and rotating-frame angular velocity fixed across normal production particle
steps. Particle transport remains enabled with its positive production time
step; this is not a no-force test mode.

The vector includes all configured density profiles, six predeclared
30-degree solar-zenith-angle surface-return-flux bands and their global
integral, both Earth-view limb branches for
Na column density and both D-line brightnesses, and global inventory. Base
blocks use the production sampling interval. Batch duration doubles until
lag-one correlation is statistically insignificant; four independent means
establish three consecutive <=2% equilibrium changes, and the following five
means are scored with Student-t 95% intervals. CI half-widths must be <=1.67%
for integrated quantities and <=3.33% for pointwise bins. Successive grids
must change integrated observables by <5% and pointwise bins by <10%.

The user-approved bounded campaign uses **500 particle iterations per completed
sample** and stops each resolution after at most **5000 particle iterations**.
This produces at most ten base samples, enough to evaluate the unchanged
nine-batch statistical minimum once if the raw blocks are already independent.
It does not guarantee qualification: correlation-driven batch aggregation,
equilibration, or confidence-width failure at the cap remains `FAIL`, never
`PASS` or `SKIPPED`. The runner still executes all four resolutions and retains
their raw block/resource artifacts. This bound is an execution policy, not a
relaxation of a failed criterion.

Immediately after every completed resolution, including a capped `FAIL`, the
runner writes `stage_diagnostic.csv`, `.png`, `.eps`, and `.json` in that
resolution directory. The first-stage figure compares its raw samples with
their raw stage means. Later figures add cumulative raw-mean comparisons with
completed coarser stages. Every artifact is marked `DIAGNOSTIC ONLY - NOT
OBSERVATIONAL VALIDATION`; it is never read by formal scoring and never
substitutes for accepted `qualified_metrics.csv` or `convergence_metrics.csv`.

The four child inputs scale the surface, near-surface, angular, and polar
targets to the selected sequence (default: 480, 240, 120, and 60 km).  The far-field and outer-domain targets and
all physical extents remain fixed, because refining a remote tail region does
not estimate near-Moon discretization error and previously caused an
uncontrolled block-count increase.  All point samples, shell integrals,
surface bands, and inventories are reduced across their actual MPI owners;
only rank zero writes the common artifacts.

For the memory campaign the fixed far-field target is 480 km (not the initial
240 km allocation-only candidate): the production parser requires far >= near,
including the 480 km coarse member. Derivation checks this invariant before
MPI launch and never silently clamps an individual member. This change is
numerical, not a physical-domain or source change. Far-field discretization
error remains outside the near-surface convergence claim.

There is no wall-time or memory cap beyond the explicit 5000-iteration bound.
`process_resource_metrics.csv` covers
initialization and the mpiexec process tree, while `resource_metrics.csv`
records completed production blocks. Actual allocator exhaustion uses the
AMPS fatal exit with reduction guidance and yields `ERROR`.

Particle restart output is compiled off in `input/moon.input`. Particle-only
restart files cannot recover the in-memory I11 sampling/batch state and were
the main source of misleading disk consumption. No valid restart/resume claim
is made; a rerun begins in a new output directory from iteration zero.

`srcMoon/examples/m1_regression.in` is the bounded default linked wiring
receipt and `m1_quick.in` is the intermediate production-wiring profile.
Neither can be supplied to I11 as convergence evidence.

Run from the repository root with a clean SPICE Moon executable:

```sh
python3 srcMoon/test/integration-tests/I11/test.py \
  --amps ./amps \
  --input srcMoon/examples/m1_memory200.in \
  --output-dir test_output/srcMoon/I11
```

The output contains a directory per resolution with logs, receipt, result,
manifest, immutable input, raw block metrics, resource histories, and runtime
products. `qualified_metrics.csv` and accepted `convergence_metrics.csv` exist
only if the unchanged gates actually pass before the cap. The parent always
contains the campaign manifest and aggregate `result.json`.
