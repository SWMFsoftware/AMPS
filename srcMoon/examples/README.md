# `srcMoon` application-input examples

## 200 GB memory-budget mesh (2026-10-14)

**Subsequent approval on 2026-10-14:** The 480/240/120/60 km qualification
sequence is now approved. It is the `run_all_test.py --campaign-mode
qualification` default. I12 can return PASS for this exact memory profile;
other unapproved CUSTOM inputs remain SKIPPED. I11 accuracy/statistical
criteria and I13 independent baseline review are unchanged. The original
allocation-only result below remains historical evidence, not retroactively
relabelled. The following approval supersedes the earlier FULL-only policy
described below: select `--m1-input srcMoon/examples/m1_sphere.in` explicitly
for the legacy FULL campaign. Coarsening never implies convergence.

`m1_memory200.in` and the volume mesh in `lola_surface.in` now request
60 km near-surface cells, 30 km terminator/polar cells, 480 km far-field
cells, and 104,226 km outer target cells. These are requested sizes; the
Cartesian AMR tree realizes dyadic sizes and enforces neighbor balance.
Near/angular/polar/outer targets are six times the original FULL values;
the far-field target is twelve times FULL. The fixed 480 km far target is
required because the production parser prohibits far-field cells smaller
than the 480 km near-surface target of the coarsest campaign member. Domain extents,
angular refinement regions, physical parameters, observation sampling, the
100 km surface triangulation target, and thermal depth/time resolution are
unchanged. Coarsening the surface triangulation alone would not solve the
volume-block allocation problem.

The budget means **200,000,000,000 bytes (186.265 GiB)**, not 200 GiB.
The earlier 30 km QUICK run used 191.9 GiB; doubling its volume targets
provides additional headroom, but cubic scaling is only an estimate, not
a memory guarantee. Actual AMR allocation must be measured. This is an
operational custom profile, not a change to physics or a qualified spatial
resolution. Particle populations and different species/builds can add memory.

The existing I12 runner exercises the production parser, resolution callbacks,
tree construction, MPI distribution and block allocation on eight ranks:

```sh
python3 -B srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --output-dir test_output/srcMoon/memory200-20261014-01
```

Choose a new output directory for every subsequent run. The resource contract
is a completed linked allocation with peak `process_tree_rss_kib * 1024`
below 200,000,000,000 bytes, plus successful I12 formula/realized-cell checks.
The resource CSV samples once per second and sums RSS (including shared pages
per process); it is not an enforced allocation cap or proof of future peak use.
Logs, input snapshot, executable/source identity, SPICE hashes and science
metrics are preserved by the runner. No observational data are used as an
oracle: I12 independently evaluates the input-defined angular function and
checks actual cell sizes against its unchanged acceptance rules.

Successful wiring on the approved memory profile can return I12 PASS;
an allocation/runtime failure is ERROR and a violated mesh criterion is FAIL.
The frozen FULL `m1_sphere.in` and legacy I11 40/20/10/5 km sequence remain
available. Neither I12 nor memory success bypasses I11 convergence or I13
baseline review or authorizes I17 observational scoring. `run_all_test.py`
can select this input with `--m1-input srcMoon/examples/m1_memory200.in --only I12`.
Default qualification selects `m1_memory200.in`; explicitly selecting
`m1_sphere.in` restores the memory-intensive legacy campaign.

Measured on 2026-10-14 at source revision
`92386a37737eab45cc8a855587a881d16bcefcd1` with the existing SPICE executable
and eight MPI ranks on the initial 240 km far-field candidate: 51,136 leaves,
43,341,576 KiB peak sampled aggregate RSS
(44.382 GB / 41.334 GiB), and approximately 83 seconds. All six SZA formula
and realized-cell checks and the polar check passed; the largest formula
relative error was 2.43e-16. I12's formal result was SKIPPED because the
profile is custom (0 formal PASS, 0 FAIL, 0 ERROR, 1 SKIPPED). Evidence is in
`test_output/srcMoon/memory200-20261014-01/I12/`, including the complete
`process_resource_metrics.csv`, `result.json`, and `run_manifest.json`.
The subsequent 480 km far-field correction has not yet completed a resource
measurement; do not attribute the earlier numbers to that changed receipt.
The initial measurement covers spherical initialization and mesh diagnostics, not
a long particle simulation or the LOLA/thermal allocation. No production
source changed or recompilation was needed for these input-only changes.

`lola_surface.in` demonstrates the complete generic
`application: moon` block consumed by the lunar application. All maintained
examples use `#begin application: moon` / `#end application`; the retired
`#section` spelling is rejected. Paths are resolved relative to the file
containing each assignment, not the launch directory. The example selects
`LDEM_4`, requests a 100 km maximum
great-circle edge on the reference sphere, and writes an AMPS CEA-long mesh
plus a Tecplot `FETRIANGLE` inspection file.

The same file contains the complete adjustable M1 comparison contract: one
run-wide UTC, the three-region altitude grid, LSO SZA profiles, optional
body-fixed pole and solid-angle-weighted shell-average profiles, Earth-view
limb geometry, and the anti-lunar tail-map coordinate system, width, circular
field diameter, exposure duration, and unresolved beam-response status.
Profile altitude is in
metres above the AMPS 1737.10 km reference sphere, not above the LOLA datum or
local terrain; all angles explicitly name degrees. The production parser
materializes an 81-point common altitude grid for the committed values.

The same LOLA example is the explicit M2 receipt used by I14-I20. It enables
the moving MOON_ME_DE421 terrain boundary, terrain-horizon illumination, a
depth-resolved Hayne-et-al. thermophysical column, Na/Ar residence laws, zero
initial dynamic inventory, and versioned `<base>.Na`/`<base>.Ar` checkpoints.
The four moving-collision controls set numerical bracketing and refinement;
I14 must establish convergence. `surface_inventory_checkpoint_file=none`
disables writes. Selecting `checkpoint` as the initial state requires both
species files and fails before AMPS initialization if either is missing.

`run_all_test.py` does not require a developer to edit `input/moon.input` for
I18. When an M2 Na+Ar ID is selected without `--skip-build`, the runner creates
a uniquely named temporary Config.pl input, verifies that its only physical
compile-time changes are `SpeciesList=Na,Ar` and the `RadiogenicArgon` block,
saves the executable and build hashes inside the campaign, deletes those
temporary inputs, and restores the canonical Na-only build. With
`--skip-build`, pass that executable explicitly as `--m2-amps`; it must not
replace canonical `--amps`.

The Na sticking table and immediate non-sticking accommodation are separate
receipt choices. The example selects `temperature-table` sticking from
Tenishev et al. (2013), Fig. 3, and `kinetic-energy` accommodation with
`omega=0.2` from its Eq. 14. Changing either law or factor requires an explicit
input edit and a replacement reference string; AMPS never interprets the
coefficient silently as the historical linear-speed blend. Ar selects
constant-unity sticking and a declared full-thermal limiting law. He/Ne remain
outside M2 and are not assigned unsupported surface parameters here.

The two `diviner_eclipse_*` settings are a paired I17 evidence receipt. They
remain `none` in the general `lola_surface.in` example. The I17 runner derives
a preserved copy that names the reproducibly processed PDS4 RDR history and
the pre-holdout frozen JSON selection contract. The observation history never
initializes or forces the thermal column: production predicts kinetic
temperature, and `MoonDivinerObservation.cpp` converts its
`epsilon*sigma*T^4` emission to comparable bolometric brightness temperature.
Supplying only one of the two paths, or a missing file, is a parser error.

The same file records a continuous uniform-global radiogenic `40Ar` source at
`3.8e20 atoms s^-1`, the lower end of the Killen (2002) diffusion estimate,
but sets `argon_source_enabled=false`. This separation is intentional: the
standard executable is the frozen Na-only baseline. To run an experimental Ar
campaign, first change `SpeciesList=Na,Ar` and only the dedicated
`SourceProcessCode=RadiogenicArgon` block to `on` in `input/moon.input`, run
the normal clean Config.pl build, then copy the `.in` receipt and set
`argon_source_enabled=true`. Do not enable the historical `MySource` block.
The source launches a local-temperature Maxwell flux uniformly per lunar
reference-sphere area. `argon_source_map_file` remains `none`; the regional
D10 South Pole-Aitken product cannot be extended to a global K-weighted source
without new authoritative data and a reviewed roadmap decision.

`m1_sphere.in` is the reviewed FULL analytic-sphere template for I11-I13. It also
contains the approved M1 AMR controls explicitly: 10 km near-surface cells
through 200 km altitude, a 40 km target at 10,000 km, monotonic coarsening to
17,371 km at the 400-radius domain boundary, the six-knot SZA resolution
table, a 15 degree terminator band, and a 5 km cap poleward of 75 degrees.
It also selects the frozen-epoch environment required for spatial convergence.
The explicit `random_seed=731993` initializes the production RNG used by all
four I11 levels; the test harness does not later replace that stream.
There are no hidden angular or radial defaults. An I11 campaign must archive derived
receipts and scale every cell-size target consistently; it must not edit this
template in place. I13 derives the same finest 5 km receipt used by I11 and
requires an identical input hash before proposing the controlled spherical
baseline; it never combines finest-grid science with this coarser template's
realized mesh configuration.

`m1_quick.in` is the QUICK production-wiring profile. It preserves the FULL
domain, physical controls, epoch, sampling grid, and observation geometry, but
uses three times every FULL spatial target: 30 km near the surface, 120 km at
the edge of the 10,000 km science region, 52,113 km at the outer boundary, and
15 km minimum angular/polar refinement. The first factor-two candidate reached
1,074,788 blocks and was stopped after 28.7 minutes at 223.7 GiB RSS in an
eight-rank I12 run. The committed factor-three profile crosses the next coarser
dyadic AMR level. Its cold eight-rank run completed the wiring criteria in
1,409.9 seconds (23.50 minutes), with 191.9 GiB peak RSS, on `passion` on
2026-10-12. This timing is host/cache-dependent and not a scientific acceptance
threshold. QUICK can verify wiring but cannot satisfy I11-I13 or replace FULL
convergence evidence.

`m1_regression.in` is the default bounded developer-regression receipt used by
`srcMoon/test/run_all_test.py`.  It retains the FULL physical domain, epoch,
sources, profile grid, and observation geometry, while coarsening every
explicit spatial target by a factor of twelve. It is deliberately labelled
`mesh_profile=custom`, so a successful linked initialization/wiring check is
reported as `SKIPPED / NOT QUALIFIED`; it can never be mistaken for I11--I13
M1 evidence. QUICK remains the intermediate production-wiring profile and
FULL remains the only qualification template.

All resolution fields remain explicit in each `.in` file. For an
investigator-selected resolution, copy either receipt, change the target
values, and set `mesh_profile = custom`. Physical extents need not change when
only resolution is being studied.

All committed examples set `mesh_output_tecplot = false`. The optional full
three-dimensional AMR `mesh.dat` and `tetra_mesh*.dat` diagnostics are not
science inputs and are not required by the mover or M1 samplers; at the M1
resolution their text representation can consume hundreds of gigabytes. Set
the value to `true` only for an explicitly provisioned visualization run. The
LOLA `surface_tecplot_file` and internal MPI `mesh.msh` exchange file are
separate products and are not disabled by this setting.

Run the bounded regression wiring profile with eight MPI ranks:

```sh
mpiexec -n 8 ./amps -input srcMoon/examples/m1_regression.in
```

Run the QUICK wiring profile with:

```sh
mpiexec -n 8 ./amps -input srcMoon/examples/m1_quick.in
```

Run the FULL M1 profile with:

```sh
mpiexec -n 8 ./amps -input srcMoon/examples/m1_sphere.in
```

Build and run from the repository root:

```sh
rm -rf build
make -j test_Moon
mkdir -p test_output/srcMoon/example-run
cd test_output/srcMoon/example-run
mpiexec -n 2 /data/vtenishe/Moon-test/AMPS/run_test_Moon/amps \
  -input /data/vtenishe/Moon-test/AMPS/srcMoon/examples/lola_surface.in \
  2>&1 | tee run.log
```

The 100 km example is substantially more expensive than the 500 km U08 smoke
fixture. For a quick production-wiring check, follow
[`../test/U08_lola_geometry/README.md`](../test/U08_lola_geometry/README.md).

To retain the historical analytic collision boundary while still recording a
complete receipt, copy the example and change:

```text
surface_geometry = sphere
```

All keys remain required in sphere mode. The LOLA paths are checked but no
terrain files are generated or loaded. Runtime `spice_path` validation records
the intended SPICE tree; it does not override the build-time SPICE selection
or kernel root supplied to `Config.pl`.

Supported limb selection is `both`, `toward-projected-sun`, or
`away-from-projected-sun`. Tail coordinates can be `j2000-ra-dec` (recommended
for the native WHAM table) or `eclipj2000-lon-lat`; display can put west or
east on the right without changing the east-positive stored coordinate.
Numeric grid, exposure, angular-list, map, beam, and altitude values are
adjustable. Shell weighting remains `solid-angle`, and beam response remains
`unresolved` because the available D14 evidence does not authorize a Gaussian
or top-hat response. The parser records these values but does not make an
unimplemented observation operator valid: pole/shell output and finite-beam/
exposure convolution remain pending U26/I-series work.
The pole, shell, and two-branch near-limb products are implemented for I11.
The D14 finite-beam/exposure tail operator remains a later observational
capability and is not inferred from those numerical-convergence products.

The realistic mesh is stored in its native MOON_ME_DE421 body-fixed frame.
For an orbit-evolving run the production mover evaluates the body/LSO state
transform at each candidate collision epoch and returns both the transformed
facet normal and local moving-surface velocity. The input therefore requires a
SPICE-enabled build; the no-SPICE build fails closed instead of treating the
terrain as stationary. The analytic-sphere mode remains the M0 baseline.

This implementation removes the former first-order Nastran dispatch failure by
routing the orbit-enabled second-order mover through the production moving-
surface callback. That is implementation evidence only: I14 remains
`SKIPPED / NOT QUALIFIED` until the FULL I11--I13 predecessor chain passes and
the linked D01 criteria are executed with preserved artifacts.
