# U08 — LOLA ingest, coordinate conversion, and geometry

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U08_lola_geometry/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The runner still executes the production-kernel subcheck when the raw local
`LDEM_4.IMG/.LBL` pair is present. It verifies native control pixels,
label-driven decoding, icosphere counts and maximum-edge resolution, outward
faces, edge uniformity, and CEA/Tecplot writers. Because the canonical D01
directory currently lacks its package README, `provenance.json`, and
`qa_report.json`, the final U08 status remains `SKIPPED / NOT VALIDATED` even
when that subcheck succeeds.

The 500 km fixture can also be used for a linked smoke run after U08 writes its
`u08.in`; that run is not U08 and does not satisfy I14. On 2026-10-10 the
production executable completed 100 steps with both one and two MPI ranks.
Those runs verify wiring and failure handling only: D01 qualification,
resolution convergence, and the complete I14 contract remain pending.

## Detailed fixture

The probe writes its own complete `moon` section beneath the selected output
directory and requests a 500,000 m maximum reference-sphere edge. The
unrelated section preceding it verifies that the production parser ignores
settings owned by another application. No test writes into the canonical
validation-data repository.

Four frozen pixel-centre controls exercise the native signed little-endian
decode and coordinate directions:

| East longitude | Planetocentric latitude | Elevation |
|---:|---:|---:|
| 0.125° | 89.875° | -119.5 m |
| 0.125° | -0.125° | -721.5 m |
| 180.125° | -0.125° | 2836.5 m |
| 359.875° | -89.875° | 91.0 m |

The topology reference for the selected refinement is 642 vertices and 1280
faces. Acceptance also requires every face to point away from the lunar
centre, maximum pre-relief great-circle edge no larger than 500,000 m, a
minimum/maximum post-relief chord ratio of at least 0.75, exact CEA record
counts, and a Tecplot `FETRIANGLE` declaration.

Run the local fixture into a new artifact directory:

```sh
python3 srcMoon/test/U08_lola_geometry/test.py \
  --output-dir test_output/srcMoon/lola-example
```

The generated linked-run input is then
`test_output/srcMoon/lola-example/U08/surface_probe/u08.in`. After a clean
`make -j test_Moon`, it can be exercised by the application from an isolated
working directory:

```sh
mkdir -p test_output/srcMoon/lola-example/linked-run
cd test_output/srcMoon/lola-example/linked-run
mpiexec -n 2 /data/vtenishe/Moon-test/AMPS/run_test_Moon/amps \
  -input /data/vtenishe/Moon-test/AMPS/test_output/srcMoon/lola-example/U08/surface_probe/u08.in \
  2>&1 | tee run.log
```

This linked command is a smoke check, not part of the local U08 score and not
a substitute for I14 resolution/convergence qualification.
