# `srcMoon` subsystem verification tests

This directory implements the local/stand-alone `U01`-`U26` layer defined by
[`../development_and_validation.md`](../development_and_validation.md).  It does
not replace the linked `I01`-`I46` campaign.  A local-kernel pass is evidence
for that kernel only; it is not evidence that the callback is reached in a
production run, that the numerical solution converges, or that observations
agree with the model.

For milestone M0, every literal U01-U26 ID remains registered and is executed
with the four-state semantics below.  Local gates used directly by the linked
I01-I10 baseline must pass; a U-test for production physics scheduled in a
later milestone remains honestly `SKIPPED` until that capability exists.  M0
does not convert such a skip into a pass and does not pull M2-M5 physics into
Phase 0.  All linked I01-I10 tests themselves must pass before M0 is complete.

The normal M0 science configuration remains neutral Na only.  I08 will use a
separate committed `Na,NA_PLUS` configuration solely to exercise the real ion
mover in analytic uniform fields; this does not advertise Na+ as enabled in
the normal production input.

## Layout

Each frozen test ID has its own directory:

```text
UXX_descriptive_name/
  README.md
  test.py
  reference/acceptance.json
```

`acceptance.json` is the machine-readable test contract and reference
solution.  It records the purpose, production entry points, configuration,
inputs, independent oracle, units, frame/time convention, stochastic controls,
procedure, metrics, thresholds, expected result and failure modes, required
data and hashes, predecessor gates, artifacts, and status semantics.

The common C++ probe links against the generated production library
`build/libAMPS.a`.  It calls compiled `srcMoon`/exosphere functions directly;
it does not copy the production formulas into a test-only implementation.  Its
expected values come from analytical invariants or explicit control points.
Source guards are deliberately labelled as such and never reported as linked
runtime verification.

### Runner and probe design

`common/moon_testlib.py` is orchestration only. It discovers frozen IDs, reads
each `acceptance.json`, compiles the thin C++ probe when required, preserves
stdout/stderr and the compile command, assigns the four-state result, and
writes JSON artifacts. It does not contain lunar physics.

`common/production_kernel_probe.cpp` links the generated `libAMPS.a`, so calls
resolve to the same configured production functions used by `srcMoon`. The
probe contains only independent expected values and invariant calculations:
closed-form point gravity, inverse-square ratios, independently sourced table
control points, and topology/orientation checks. A probe PASS means only that
named kernel and fixture passed; it is not a full executable, convergence,
conservation, or validation result.

U05 additionally verifies the absolute one-AU Na radiation-pressure curve
against a reproducible two-resolution digitization of Combi et al. (1997),
Figure 7. The historical digitization that produced the embedded array is
unavailable and is not claimed as provenance. The runner verifies every
reference SHA-256 before comparing 14 unobscured curve points with the compiled
production kernel using a predeclared `0.5 cm s^-2` graphical uncertainty.

U08 has an additional two-part status:

1. the local production parser/decoder/geometry/writer probe runs whenever the
   exact raw IMG/LBL pair is available;
2. the enclosing U08 result checks D01 package qualification independently.

Consequently the kernel can report PASS while U08 reports
`SKIPPED / NOT VALIDATED`. This is deliberate and prevents raw directory
presence from being treated as qualified validation evidence.

## Running

### Complete campaign

`run_all_test.py` is the single entry point for the complete registered
`srcMoon` U-series campaign. By default it first removes the generated
`build/` tree, runs `make -j`, and then launches every registered test through
that test's own `test.py` runner:

```sh
python3 srcMoon/test/run_all_test.py
```

The default output directory has the form
`test_output/srcMoon/run-all-YYYYMMDDTHHMMSSZ`. An explicit `--output-dir`
must be new or empty; the runner refuses to overwrite an earlier campaign.

Use `--skip-build` only when the existing generated tree is known to match the
source and configuration. `--list` prints the explicit campaign registry,
`--only UXX` may be repeated for a focused runner check, `--stream` mirrors
individual logs to the terminal, and `--strict-skips` makes any SKIPPED test
fail the aggregate command. Without `--strict-skips`, SKIPPED remains a
recorded incomplete capability but does not change the campaign exit status.

The final terminal summary reports the number of tests executed, passed,
failed, skipped, and ended with an error. For every FAIL or ERROR it prints the
test ID, the exact individual runner, and that runner's combined stdout/stderr
log. The same information is written to
`<output-dir>/run_all_summary.json`; individual logs are under
`<output-dir>/runner_logs/UXX.log`.

The aggregate exit values are `0` when no test failed or errored, `1` for an
evaluated FAIL (or a strict-skip rejection), and `2` for ERROR. A clean-build
failure is a campaign infrastructure ERROR: no individual test is counted as
executed, and diagnostics are preserved in `<output-dir>/build.log`.

### Adding a test to the aggregate campaign

Adding a `UXX_*` directory is not sufficient. Every new permanent test must
also be added explicitly to the `TESTS` registry in `run_all_test.py` as its
stable ID and relative `test.py` path. The runner performs a preflight
comparison between that registry and all `U[0-9][0-9]_*` directories and exits
with ERROR if either side is missing. Permanent IDs must first be added to the
roadmap/test registry consistently, as required by `AGENTS.md`.

Each individual `test.py` runner must accept these arguments:

```text
--build-dir <generated-AMPS-build>
--output-dir <campaign-output-root>
```

It must write exactly `<output-dir>/<UXX>/result.json`. That JSON object must
contain `id` equal to its registered `UXX` and `status` equal to one of `PASS`,
`FAIL`, `SKIPPED`, or `ERROR`. It should also contain a descriptive `name` and,
for SKIPPED or ERROR, a precise `reason`. The process exit must agree with the
JSON status: `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. Missing, malformed,
stale, ID-mismatched, status-invalid, or exit-inconsistent output is classified
as ERROR by `run_all_test.py`; the individual runner log is retained for
diagnosis.

The aggregate removes an existing per-test `result.json` immediately before
launching that test, so a crashed runner cannot be credited with a stale
result. Tests must keep all other artifacts beneath their own
`<output-dir>/<UXX>/` directory and must never write into validation-data
`raw/` directories.

### Lower-level runners

List contracts:

```sh
python3 srcMoon/test/run_tests.py --list
```

Run all currently specified U tests:

```sh
python3 srcMoon/test/run_tests.py --all \
  --output-dir test_output/srcMoon/unit
```

Run selected tests or an individual directory runner:

```sh
python3 srcMoon/test/run_tests.py --test U03 --test U12
python3 srcMoon/test/U03_lunar_gravity/test.py
```

Use `--strict-skips` when automation must reject an incomplete campaign.  By
default, skipped capabilities are recorded but do not make the runner fail.
Exit values are `0` for PASS, `1` for FAIL, `2` for ERROR, and `77` for an
individual SKIPPED test.

The optional `--rebuild` flag first removes `build/` and then runs `make -j`, as
required by this checkout's build discipline.  It is not CMake-based.  The
canonical application regression remains:

```sh
rm -rf build
make -j test_Moon
ls -l test_Moon.diff
```

The required regression result is a zero-byte `test_Moon.diff`.  Preserve its
build log under `test_output/srcMoon/`.

To run only the currently executable local checks:

```sh
python3 srcMoon/test/run_tests.py \
  --test U02 --test U03 --test U04 --test U05 --test U07 \
  --test U08 --test U12 --test U14 --test U22 \
  --output-dir test_output/srcMoon/kernel-tests
```

U08 returns exit code 77 while D01 qualification files are absent even when
its nested kernel result is PASS. Inspect `U08/result.json` rather than
interpreting that exit as a numerical failure.

U04 additionally requires the SPICE-enabled Moon configuration and exact
frozen kernel bytes. Configure and clean-build before running it:

```sh
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
rm -rf build
make -j
python3 srcMoon/test/U04_rotating_frame/test.py \
  --output-dir test_output/srcMoon/U04
```

The runner verifies every SHA-256 declared by the U04 contract before calling
CSPICE. It writes `kernel_hashes.json`, component-wise `term_vectors.json`, and
the complete `result.json`. Missing kernels mean `SKIPPED`; bytes that disagree
with the contract mean `ERROR`.

## Current implementation status

| ID | Scope | Runner state | Meaning |
|---|---|---|---|
| U01 | LOS geometry | SKIPPED | No deterministic production ray fixture has been isolated. |
| U02 | target isolation/wiring | implemented source guard | Configuration evidence only. |
| U03 | lunar gravity | implemented linked probe | Point-mass kernel at a closed-form point. |
| U04 | rotating frame | implemented linked probe | Exact J2000/LSO/MOON_ME_DE421 names, kernel hashes, differential gravity, fictitious terms, transform derivative, and orthogonality are checked at the frozen epoch. |
| U05 | Na radiation pressure/shadow | implemented linked probe | Earth-umbra logic, inverse-square invariant, and absolute Combi-curve agreement to frozen graphical uncertainty; full mover gating remains I07. |
| U06 | Lorentz force | SKIPPED | Active species list has no ion. |
| U07 | refinement selector | implemented linked probe | Current constant callbacks only. |
| U08 | LOLA geometry | kernel subcheck implemented; campaign SKIPPED | Production parser/decoder/icosphere/writers are exercised, but D01 lacks package README/provenance/QA and I14 is pending. |
| U09 | terrain shadow | SKIPPED | Prototype is not wired to production; D01/D02 are absent. |
| U10 | Diviner temperature | SKIPPED | No Diviner-backed production model; D03 is absent. |
| U11 | thermal inertia | SKIPPED | No dynamic thermal-state production model. |
| U12 | surface interaction | implemented linked probe | Temperature/sticking control points only. |
| U13 | cold trapping | SKIPPED | No PSR/cold-trap production state. |
| U14 | Na photoionization | implemented linked probe | Legacy constant lifetime only. |
| U15 | multispecies photoionization | SKIPPED | Na is the sole active species. |
| U16 | electron impact | SKIPPED | Existing constants are not a qualified D05 implementation. |
| U17 | charge exchange | SKIPPED | No configured production reaction path. |
| U18 | plasma driver | SKIPPED | D06 reports ERROR and no qualified driver is wired. |
| U19 | helium source | SKIPPED | Source file exists but He is disabled and reservoir physics is unresolved. |
| U20 | neon source | SKIPPED | Source file exists but Ne is disabled. |
| U21 | argon source | SKIPPED | No radiogenic geography/transient model is present. |
| U22 | Na sources | implemented linked probe and configuration guard | M0 built-in impact-on/`MySource`-off invariant, impact normalization, and PSD density kernel; no night-to-day reservoir claim. |
| U23 | meteoroid driver | SKIPPED | No production forcing path. |
| U24 | H2O/OH chemistry | SKIPPED | No production species/chemistry path. |
| U25 | data provenance | SKIPPED | Requires a selected Dxx package and linked campaign context. |
| U26 | observation geometry | SKIPPED | Instrument-specific operators are not qualified. |

The status is capability-sensitive.  Downloaded files, directory presence, or
a source/preflight check never changes a missing linked capability into PASS.

## Results and reproducibility

Results go under `test_output/srcMoon/` and never under the validation-data
repository. Each test writes `result.json`; `run_tests.py` writes
`summary.json`, while the individual-runner campaign writes
`run_all_summary.json` and one log per registered runner. Probe compilation
details are retained in the result. Linked I-series runs require the larger
manifest/artifact set specified in the roadmap and are intentionally outside
these local runners.

See [`TEST_GAP_PROPOSALS.md`](TEST_GAP_PROPOSALS.md) for missing coverage.  No
new permanent IDs are assigned there.
