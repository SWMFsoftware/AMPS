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

The 500 km fixture can also be used to investigate linked runtime wiring after
U08 writes its `u08.in`; that run is not U08 and does not satisfy I14. A clean
no-SPICE one-rank run on 2026-10-11 parsed the receipt, loaded the generated
terrain, and printed all 81 altitude points, but stopped at iteration 1 with
`pic_mover.cpp: Error: undetermined internal boundary type`. The selected
first-order no-force mover's initial boundary-dispatch switch lacks the
Nastran-surface case that is present in the second-order path. This is a linked
production `ERROR`, not a U08 geometry-kernel failure, and supersedes the
earlier smoke-run claim until the mover path is repaired and rerun. D01
qualification, mover repair, resolution convergence, and the complete I14
contract remain pending.

## Detailed fixture

The probe writes its own complete generic `application: moon` block beneath
the selected output directory and requests a 500,000 m maximum
reference-sphere edge. It also rewrites the block delimiters to the retired
`#section` spelling and requires the production parser to reject that file;
there is no undocumented dual-grammar compatibility path. No test writes into
the canonical validation-data repository.

Because U08 already exercises the production `moon` parser, the fixture also
checks that the agreed adjustable M1 altitude receipt produces exactly 81
points with boundaries at 0, 200, 1000 and 10000 km, and retains the declared
SZA/limb/tail scalar values. Those are parser invariants only. They do not
exercise an observation operator and cannot change U26 or any I-series status.
The fixture also requires `mesh_output_tecplot = false`, verifies that the
production parser retains it, and rejects a non-Boolean value. This covers the
configuration contract for suppressing full-volume AMR visualization; the
surface Tecplot writer exercised later in U08 is a separate LOLA product.
The fixture labels its deliberately coarse mixed-resolution receipt `custom`
and rejects the retired `production` label; committed science receipts use
the explicit `quick` or `full` identities.

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

This linked command currently reproduces the first-order mover error described
above. It is not part of the local U08 score and is not a substitute for I14
resolution/convergence qualification.
