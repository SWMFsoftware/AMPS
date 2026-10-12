# Linked integration and numerical tests

Each `I01`–`I10` directory is a frozen test ID with an individual runner and
machine-readable contract. The runners call the actual configured Moon
executable with one MPI rank by default and one OpenMP thread. A focused run
may select another positive process count with `--mpi-ranks N`; the runner uses
that value in the `mpiexec -n N` command and records the same value in
`run_manifest.json`. The executable records raw
SI quantities after full `amps_init()`; `integration_common.py` evaluates
independent analytical, invariant, or high-resolution numerical references.
I11 applies the same pattern to a complete 40/20/10/5 km campaign. Each level
uses production transport and sampling, adaptive batch means, frozen seed and
layout controls, and archived resource histories. The independent Python
scorer evaluates the declared integrated and pointwise convergence gates. I12
applies the pattern to the
input-driven linked AMR mesh. I13 evaluates
the analytical sphere invariant and records the baseline inventory, but remains
`SKIPPED` until I11 supplies equilibrated sampled metrics and an investigator
reviews and freezes the candidate reference CSV. A wiring/preflight result
never substitutes for I11 convergence.

The aggregate runner has two explicit modes.  Its default `regression` mode
uses `srcMoon/examples/m1_regression.in`, requests eight MPI ranks, records
I11 and I13 as honest `SKIPPED` non-executions, and launches one coarse linked
I12 production initialization/wiring run.  `qualification` mode is opt-in and
runs the unbounded FULL I11--I13 campaign with the unchanged science gates:

```sh
python3 srcMoon/test/run_all_test.py --campaign-mode regression --mpi-ranks 8
python3 srcMoon/test/run_all_test.py --campaign-mode qualification --mpi-ranks 8
```

Run a correctly configured executable with:

```sh
python3 srcMoon/test/integration-tests/I02/test.py \
  --amps ./amps --output-dir test_output/srcMoon/M0
```

For example, an eight-rank focused I12 QUICK wiring run is:

```sh
python3 srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps \
  --input srcMoon/examples/m1_quick.in \
  --mpi-ranks 8 \
  --output-dir test_output/srcMoon/I12-quick-rank8
```

Changing the process count changes the numerical layout and therefore must be
treated as an explicit run control. It does not promote a QUICK profile from
`SKIPPED / NOT QUALIFIED` to a FULL M1 result.

Every normally completed run preserves `stdout.log`, `stderr.log`,
`receipt.txt`, `result.json`, `run_manifest.json`, and, when supplied,
`application.in` below `<output-dir>/IXX/`. An infrastructure ERROR retains
the subset produced before failure and identifies the missing evidence rather
than creating placeholders. `PASS`, `FAIL`, `ERROR`, and
`SKIPPED` have the meanings defined in the top-level test README. A kernel or
preflight PASS never replaces a linked-runtime PASS.

I11-I13 are registered in the top-level `run_all_test.py` campaign. The
aggregate runner passes the selected executable and reviewed M1 input through
to these same individual runners, retains all files under each per-ID artifact
directory, and consumes `result.json` with the same four-state/exit-code
protocol. It does not reproduce or reinterpret the linked science checks.

I11-I13 preserve one-second process-tree RSS and host-available memory in
`process_resource_metrics.csv`. I11 additionally records resource state at
every completed production sampling block. Monitoring imposes no wall-time,
block-count, or memory ceiling; actual C++ allocation exhaustion exits through
the detailed AMPS fatal diagnostic and is an `ERROR`.

After load balancing and before block allocation, srcMoon also computes a
host-level lower bound containing only each owned leaf's block object and
block-associated storage.  It aborts only when that unavoidable lower bound
already equals or exceeds Linux `MemAvailable`; it applies no guessed memory
cap and intentionally excludes center/corner nodes, particles, and MPI
overhead.  The diagnostic names the input-file resolutions that can be
coarsened.  The optimized debug build used by the aggregate runner passes
`OPT3=-O3 -g` to `make`, retaining GDB symbols without turning the linked
production run into an unoptimized build.

I05 requires the no-SPICE neutral-Na build. I08 requires the dedicated
no-SPICE `Na,NA_PLUS` fixture. The other dynamic tests require the DE421 SPICE
build. I04 consumes the canonical `make -j test_Moon` zero-byte diff. The
campaign runner handles these configuration boundaries and rebuilds generated
trees rather than mixing object files from different fixtures.
