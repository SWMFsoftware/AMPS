# `srcMoon` application-input examples

`lola_surface.in` demonstrates the complete runtime section consumed by the
lunar application. Paths are resolved relative to the input file, not the
launch directory. The example selects `LDEM_4`, requests a 100 km maximum
great-circle edge on the reference sphere, and writes an AMPS CEA-long mesh
plus a Tecplot `FETRIANGLE` inspection file.

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

The realistic surface is static and label-frame-native. Because the normal
production input now enables orbit evolution, this example must currently be
run from an explicitly configured no-SPICE build. The executable rejects LOLA
mode in an orbit-enabled build until the moving-boundary/frame-transform
problem has an explicit verified implementation. The analytic-sphere
production run is the current SPICE-enabled M0 baseline.
