# I11 — AMR spatial-resolution convergence

The existing I11 wrapper automatically runs `test_infrastructure.py` before
MPI launch. Three isolated tests cover both approved sequences, fixed
far-field/physical controls, exact finest-input reproduction for I13, and
rejection of unapproved custom settings. These tests generate only temporary
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
must change integrated observables by <5% and pointwise bins by <10%. There is
no statistical block ceiling: while a required bin is unresolved, the run
continues rather than claiming either `PASS` or `SKIPPED`. An evaluated
successive-resolution violation is `FAIL`; interruption or resource exhaustion
before qualification is `ERROR` and leaves incomplete evidence.

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

There is no wall-time or memory cap. `process_resource_metrics.csv` covers
initialization and the mpiexec process tree, while `resource_metrics.csv`
records completed production blocks. Actual allocator exhaustion uses the
AMPS fatal exit with reduction guidance and yields `ERROR`.

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
manifest, immutable input, block/qualified metric CSVs, resource histories,
and runtime products. The parent contains the campaign manifest,
`convergence_metrics.csv`, and aggregate `result.json`.
