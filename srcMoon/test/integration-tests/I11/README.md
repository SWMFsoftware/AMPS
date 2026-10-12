# I11 — AMR spatial-resolution convergence

I11 uses the normal SPICE-enabled Moon executable and the reviewed spherical
FULL M1 receipt. The runner rejects QUICK/CUSTOM inputs before launching an
expensive model and creates immutable 40, 20, 10, and 5 km child receipts,
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
targets to 40, 20, 10, and 5 km.  The far-field and outer-domain targets and
all physical extents remain fixed, because refining a remote tail region does
not estimate near-Moon discretization error and previously caused an
uncontrolled block-count increase.  All point samples, shell integrals,
surface bands, and inventories are reduced across their actual MPI owners;
only rank zero writes the common artifacts.

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
  --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I11
```

The output contains a directory per resolution with logs, receipt, result,
manifest, immutable input, block/qualified metric CSVs, resource histories,
and runtime products. The parent contains the campaign manifest,
`convergence_metrics.csv`, and aggregate `result.json`.
