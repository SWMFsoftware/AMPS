# I12 — linked angular surface refinement

For the 200 GB resource check, pass
`--input srcMoon/examples/m1_memory200.in --mpi-ranks 8` instead of the
FULL template. The [profile contract](../../../examples/README.md#200-gb-memory-budget-mesh-2026-10-14)
defines the budget, measured peak, exact command and limits of this evidence.
This extends I12's existing configuration coverage, not its meaning: the
same production path and unchanged analytical/realized-cell tests run.
Following the 2026-10-14 approval, this exact custom memory profile can pass
I12's unchanged wiring criteria. Other unapproved CUSTOM profiles remain
SKIPPED. I12 PASS alone still does not qualify M1 or authorize I17 scoring.

I12 loads a complete Moon `.in` receipt, initializes the normal production
mesh, and queries the production boundary-refinement callback at solar zenith
angles 0, 30, 60, 90, 120, and 180 degrees. It also constructs the lunar north
pole through the SPICE `MOON_ME_DE421` to `LSO` transformation.

The executable records requested and realized boundary-adjacent cell sizes.
The independent Python checker reconstructs the input-defined piecewise-linear
selector and requires machine-precision agreement. Realized AMR cells must be
within 25 percent of the requested target. AMPS can realize only power-of-two
child sizes, so the roadmap's discrete-AMR exception accepts the nearest
adjacent dyadic level: `abs(log2(realized/requested)) <= 0.5`. AMPS also
enforces a maximum one-level difference between neighboring leaves. At a
declared knot adjoining a finer target, one additional finer level is accepted
only when the dyadic parent would satisfy the unchanged 25% criterion. Coarser
cells, non-transition samples, and more than one extra level are not excused.
`result.json` identifies the exception and parent-level error; the test never
changes the requested physical target to make the mesh pass.

Run through the M1 campaign runner. For a focused run with an already
configured executable:

```sh
python3 srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps \
  --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I12-focused
```

The runner defaults to one MPI process. Use `--mpi-ranks 8` for an eight-rank
decomposition; the exact rank count is propagated to `mpiexec`, retained in
`run_manifest.json`, and written into the production receipt. For distributed
I12, every rank participates in ownership reductions proving that each sampled
AMR leaf has exactly one allocated owner; only rank zero writes the receipt.
I11 and I13 also use reviewed distributed reductions; earlier I01--I10 linked
diagnostics retain their separately reviewed one-rank contracts.
One OpenMP thread per rank remains fixed by this runner.

The focused command uses the reviewed spherical M1 template. A campaign runner
must archive any derived convergence receipts and must not edit that template.
The template explicitly sets `mesh_output_tecplot = false`: full-volume
Tecplot serialization is not part of the I12 oracle and previously exhausted
the filesystem without changing the requested or realized AMR mesh. The
internal MPI `mesh.msh` tree exchange remains active.
The runner also preserves `process_resource_metrics.csv` at one-second cadence;
this monitoring records process-tree RSS and host-available memory but imposes
no execution or memory cap.

For a faster production-wiring check, substitute
`srcMoon/examples/m1_quick.in`. QUICK uses three times every spatial target
relative to the human-authored FULL template while retaining the same domain
and physical controls. The first factor-two candidate reached 1,074,788 blocks
and was stopped at 28.7 minutes and 223.7 GiB RSS on eight ranks. Factor three
crosses the next coarser dyadic AMR level. Its cold eight-rank run completed in
1,409.9 seconds (23.50 minutes), with 191.9 GiB peak RSS, on `passion` on
2026-10-12. That host/cache-dependent runtime is operational guidance, not an
acceptance gate. A QUICK I12 result cannot be reported as the FULL M1 I12 gate.
The runner preserves the QUICK metrics but reports `SKIPPED / NOT QUALIFIED`
when its wiring criteria pass; a criterion violation remains `FAIL`, and an
infrastructure/runtime problem remains `ERROR`.

For routine regression use `srcMoon/examples/m1_regression.in`, which uses
factor-twelve spatial targets and is the default selected by
`run_all_test.py --campaign-mode regression`.  It exercises the same linked
production initialization, callback, AMR realization, MPI ownership, and
independent Python formula checks.  Because it is explicitly `custom`, passing
wiring metrics remain `SKIPPED / NOT QUALIFIED`; a violated evaluated metric
is still `FAIL`, and a build/runtime problem is still `ERROR`.
