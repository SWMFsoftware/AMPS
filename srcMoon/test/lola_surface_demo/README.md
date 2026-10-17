# LOLA production-surface meeting demonstration

This runner creates a presentation-quality lunar surface from the immutable
D01 `LDEM_4` product while preserving formal V&V status separately. It is not
a new permanent test ID and does not repurpose U08 or I14.

The runner first requires U08 to return `PASS`. It then calls a thin adapter
linked against the same production Moon input parser, `LolaDem` decoder,
icosphere generator, radial topography adjustment, CEA writer, and Tecplot
writer used by `srcMoon`. Python independently decodes every selected native
D01 elevation from the IMG using the detached-label layout and compares it
with every generated surface vertex. No smoothing, gap filling, fitted
correction, hand-selected site, synthetic observation, or model-generated
reference is used.

From the repository root, after the normal clean SPICE Moon build, run:

```sh
python3 -B srcMoon/test/lola_surface_demo/run.py \
  --surface-resolution-m 100000
```

An explicit new output directory may be supplied with `--output-dir`. The
runner refuses to overwrite a non-empty directory. The default creates a UTC-
stamped directory under `test_output/srcMoon/lola-surface-demo/`.

The resolution is the maximum requested pre-topography great-circle edge in
metres and is restricted to 10--1000 km. The default 100 km matches
`srcMoon/examples/lola_surface.in`. Finer values increase vertices, output
size, render time, and memory; they are not automatically stronger validation
evidence.

Artifacts include:

- `application.in`: exact derived production receipt;
- `lola_surface.cea`: AMPS-loadable one-based surface connectivity;
- `lola_surface.dat`: Tecplot `FETRIANGLE` surface;
- `surface_vertices.csv`: complete production vertex population;
- `lola_dem_comparison.csv`: every vertex and independently decoded D01 value;
- `lola_surface_demonstration.png` and `.eps`: maps, identity comparison, and
  residual distribution;
- `stdout.log`, `stderr.log`, `run_manifest.json`, and `result.json`;
- `verification/U08/result.json`: prerequisite local verification evidence.

The generated `.cea` or `.dat` surface can also be rendered from arbitrary
camera directions with `srcMoon/tools/visualize_surface.py`; usage and the
required explicit frame/reference-radius metadata are documented in
`srcMoon/tools/README.md`. That general-purpose rendering remains a separate
visualization artifact and does not alter this runner's status semantics.

`result.json` uses `execution_status=PASS` only for successful demonstration
generation. Its separate `formal_vv` object records `execution=NOT_RUN`,
`status=SKIPPED`, and `validation_status=NOT_VALIDATED` for I14. The top-level
capability remains the frozen roadmap label `EXPERIMENTAL`; the separate
artifact class is `MEETING_DEMONSTRATION`.
The figures are titled `MEETING DEMONSTRATION - FORMAL I14 NOT RUN`. I14 still
requires same-campaign I11/I13 qualification, the linked production BVH ray
queries, and second-order particle/surface collisions on the configured MPI
ranks. A good-looking surface or a small vertex residual cannot replace those
gates.
