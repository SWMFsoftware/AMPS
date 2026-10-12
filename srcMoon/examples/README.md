# `srcMoon` application-input examples

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

The realistic surface is static and label-frame-native. Because the normal
production input now enables orbit evolution, this example must currently be
run from an explicitly configured no-SPICE build. The executable rejects LOLA
mode in an orbit-enabled build until the moving-boundary/frame-transform
problem has an explicit verified implementation. The analytic-sphere
production run is the current SPICE-enabled M0 baseline.

The current no-SPICE production selection uses the first-order no-force mover.
Its initial internal-boundary dispatch does not handle the registered Nastran
terrain descriptor, so a normal LOLA run presently stops at iteration 1 with
`undetermined internal boundary type`. Input parsing and terrain generation
complete before that error, but the case is not a successful linked terrain
run. Do not use it as I14 evidence until the production mover is repaired and
the linked test is repeated.
