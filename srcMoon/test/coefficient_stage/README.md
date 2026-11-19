# Normal environment/coefficient stage

`Moon::PrepareProductionParticleEnvironment` is the identical orbit/thermal/
cache block formerly inside `amps_time_step`. The normal step and conditional
cohort call this SAME function before particle motion. No force formula,
source injection or sampling implementation is copied. The extraction preserves
legacy evolving midpoint/end epochs, frozen environmental epochs, SI conversions,
rotation/body policy and optional thermal/inventory side effects.

Run the supporting native owner using an authentic successful normal-SPICE
producer from the same source checkout. The explicit isolated development
campaign automatically supplies this proof after ordinary restoration:

```sh
/usr/bin/python3 -B srcMoon/test/run_m3_coefficient_stage.py \
  --normal-build-evidence /absolute/normal_build_evidence.json \
  --output-dir /absolute/NEW-coefficient-stage-results \
  --wall-seconds 600
```

This wrapper does not configure the active checkout. It verifies the actual
archive/executable/generated ABI, source hashes and two successful normal build
commands, then compiles a thin adapter. The enclosing campaign keeps its own
original deadline;600s is not a new hour for this child.

The six declared jobs are frozen first/second tracing at1/8MPI and legacy evolving
at1/8MPI. Each calls the SAME normal stage three times. New tracing with evolving
external geometry remains unsupported; the evolving check preserves the existing
legacy branch. Expected named criteria:26per frozen job and25per evolving job,
154total. Direct furnished CSPICE **state** queries supply independent SI/frame/
time oracles. `spkpos` is not substituted for a bitwise state contract: dynamic
frame numerical paths can differ between those toolkit APIs.

`comparison.csv`, `comparison.png`, `comparison.eps`, named-check CSVs, compiler/
MPI logs, actual furnished-kernel hashes, source/build identities and
`result.json` remain under the selected output directory. An allocated-node
normal neutral-Na acceleration check verifies nonzero initialized coefficients;
no particle/observational qualification is inferred. Conservation and full
force/mover convergence retain their existing separate owners and source bases.

The supporting global group `coefficient-stage-controls` executes the thin input
owner's exact zero-force macro admission and actual1/8MPI status-consensus
controls. It supplies software checks only. Run it independently with:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_coefficient_stage.py \
  --output-dir /absolute/NEW-input-admission-controls
```

Missing toolkit/producer is SKIPPED; build/parser/runtime/hash/resource failures
are ERROR; evaluated numerical mismatches are FAIL. PASS always refers to these
precisely declared software/wiring criteria, never observational validation or
closure of missing physical Na states/crosswalks.
