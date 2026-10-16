# `srcMoon` subsystem verification tests

The 200 GB allocation check reuses the registered I12 runner with
`--only I12 --m1-input srcMoon/examples/m1_memory200.in --mpi-ranks 8`.
See [the profile contract and direct command](../examples/README.md#200-gb-memory-budget-mesh-2026-10-14).
The approved memory profile can now pass I12 wiring, but M1 still requires
I11 convergence and reviewed I13 baseline evidence. Qualification defaults
to 480/240/120/60 km; explicitly select `--m1-input srcMoon/examples/m1_sphere.in`
only to run the legacy memory-intensive FULL sequence.

I17's runner now includes its launcher/plot-reader infrastructure checks
automatically and enforces the approved 2026-10-14 three-step thermal numerical
error budget before residual computation. See the
[I17 contract](integration-tests/I17/README.md) for the 10%-of-observational-sigma
criterion and the `dt_div_4/` location of the finest-step experimental figures.
Infrastructure success is not I17 validation.

This directory implements the local/stand-alone `U01`-`U26` layer and the
currently implemented linked M1 `I11`-`I13` and registered M2 `I14`-`I20`
layers defined by
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

The production input now explicitly selects
`PIC::Mover::UniformWeight_UniformTimeStep_SecondOrder`. This is required
because AMPS's global default is the no-force trajectory/boundary mover; merely
installing `Moon::TotalParticleAcceleration` does not cause that default mover
to call it. I05-I10 therefore test the same force-aware mover selected for
ordinary Moon particles, rather than a test-only integration path.

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

`run_all_test.py` is the single entry point for all 36 currently registered
`srcMoon` U-series and linked tests: U01-U26 followed by I11-I20. By default
it runs the bounded `regression` campaign with eight MPI ranks: I11 and I13
produce explicit `SKIPPED` non-execution records, while I12 runs the actual
production initialization and wiring with `srcMoon/examples/m1_regression.in`.
The default regression campaign builds only the canonical Na executable,
because I13 and hence I18-I20 are deliberately not launchable. A qualification
selection containing both I13 and any of I18-I20 first creates and preserves a
dedicated Na+Ar M2 executable and its hashed build receipt. It then removes the
generated `build/` tree again, restores the canonical Na-only Moon/SPICE
configuration, runs `make -j`, and launches every registered test
through that test's own `test.py` runner. Configuration and compilation are
both recorded in `<output-dir>/build.log` and as separate steps in
`run_all_summary.json`; the M2 build has its own log and JSON receipt:

```sh
python3 srcMoon/test/run_all_test.py
```

The default build sequence is equivalent to:

```sh
rm -rf build
./Config.pl -application=moon \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j "OPT3=-O3 -g"
```

The command-line optimization override retains debug symbols while restoring
production optimization; it does not edit the generated Makefile or
Makefile.conf. The resulting canonical `./amps` executable and the bounded
`srcMoon/examples/m1_regression.in` is passed to I11--I13 and
`srcMoon/examples/lola_surface.in` is passed to I14--I20. In default
regression mode. Request the formal, unbounded FULL campaign explicitly:

```sh
python3 srcMoon/test/run_all_test.py \
  --campaign-mode qualification --mpi-ranks 8
```

Qualification selects the reviewed `srcMoon/examples/m1_sphere.in` receipt.
Override `--amps` or `--m1-input` only when the
alternate executable/input identity is intentional and will be captured in the
linked manifests. These M1 meshes can be substantially more expensive than
the U-series kernel checks.

Relocated installations can be selected explicitly with `--spice-path` and
`--spice-kernels`. These arguments are ignored only when `--skip-build` is
used, because that option deliberately reuses the caller's existing generated
tree.

The default output directory has the form
`test_output/srcMoon/run-all-YYYYMMDDTHHMMSSZ`. An explicit `--output-dir`
must be new or empty; the runner refuses to overwrite an earlier campaign.

Use `--skip-build` only when the existing generated tree is known to match the
source and SPICE configuration. If a qualification selection includes I13 and
I18, I19, or I20 with `--skip-build`, `--m2-amps` must identify a separately
built Na+Ar executable; the canonical `--amps` remains Na-only. `--list`
prints the explicit campaign
registry, `--only UXX` or `--only IXX` may be repeated for a focused runner
check, `--stream` mirrors individual logs to the terminal, and `--strict-skips`
makes any
SKIPPED test fail the aggregate command. Without `--strict-skips`, SKIPPED
remains a recorded incomplete capability but does not change the campaign exit
status code. The machine-readable `campaign_status` remains `SKIPPED` whenever
one or more prerequisites/capabilities are skipped and no FAIL/ERROR takes
precedence; a zero developer-workflow exit never relabels incomplete evidence
as PASS.

The final terminal summary reports the number of tests executed, passed,
failed, skipped, and ended with an error. For every FAIL or ERROR it prints the
test ID, the exact individual runner, and that runner's combined stdout/stderr
log. The same information is written to
`<output-dir>/run_all_summary.json`; individual logs are under
`<output-dir>/runner_logs/UXX.log` or `IXX.log`. A normally completed linked
runner keeps `receipt.txt`, `stdout.log`, `stderr.log`, `run_manifest.json`, an
input snapshot, and `result.json` under `<output-dir>/IXX/`. An infrastructure
ERROR may retain only the subset created before failure. The aggregate JSON
records the artifact directory and manifest path without flattening or
overwriting the individual evidence. Its `artifacts` map indexes only files
that actually exist, including partial evidence from an ERROR run. A reported
PASS, FAIL, or SKIPPED missing any mandatory linked artifact is reclassified
as ERROR.

The aggregate exit values are `0` when no test failed or errored, `1` for an
evaluated FAIL (or a strict-skip rejection), and `2` for ERROR. A Moon/SPICE
configuration failure or clean-build failure is a campaign infrastructure
ERROR: no individual test is counted as executed, the failed step is named in
`run_all_summary.json`, and diagnostics are preserved in
`<output-dir>/build.log`.

### Adding a test to the aggregate campaign

Adding a `UXX_*` or linked I-series directory is not sufficient. Every new permanent
test must also be added explicitly to the `TESTS` registry in
`run_all_test.py` as its
stable ID and relative `test.py` path. The runner performs a preflight
comparison between that registry and all `U[0-9][0-9]_*` directories and exits
with ERROR if either side is missing. Permanent IDs must first be added to the
roadmap/test registry consistently, as required by `AGENTS.md`.

Each U-series `test.py` runner must accept these arguments:

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

The registered I11-I20 runners instead accept `--amps`, `--input`,
`--output-dir`, `--campaign-mode`, and `--mpi-ranks`, and write
`<output-dir>/<IXX>/result.json` using the same
status/exit-code protocol. They also preserve the linked artifacts required by
`AGENTS.md`. Successful process launch or a wiring invariant never substitutes
for science evidence. I11 now runs and scores the complete four-resolution
campaign. I13 remains `SKIPPED` until I11 passes and both its sampled-science
CSV and deterministic mesh/species/source/configuration JSON candidates are
reviewed and frozen; the runner never installs either reference automatically.

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

The reviewed four-rank reference was updated on 2026-10-11 for the production
force-aware mover after two complete runs produced identical bytes. Its
SHA-256 is
`abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
The prior no-force reference hash was
`9386ea2d353db638ce84bf1f202d55cb7276443f800938c6365807be8a8deef3`;
that physical configuration is intentionally retired rather than preserved as
a legacy option. The campaign runner only compares output and never rewrites
the reference. Any later nonzero diff therefore remains a reviewed FAIL, not
permission for an automatic baseline update.

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

### Linked M0 baseline

The linked `I01`-`I10` contracts and individual runners are under
`integration-tests/`. The campaign entry point performs separate clean builds
for the normal SPICE baseline, the no-SPICE two-body fixture, and the dedicated
no-SPICE `Na,NA_PLUS` fixture. It also runs the complete U-series registry,
reruns U06 in the ion fixture, runs the canonical `test_Moon` regression, and
restores the normal SPICE build at the end:

```sh
python3 srcMoon/test/run_tests.py --phase baseline \
  --output-dir test_output/srcMoon/M0-baseline
```

Equivalently, call `run_integration_tests.py` directly and optionally provide
relocated `--spice-path` and `--spice-kernels` roots. Every output directory
must be new or empty; scientifically meaningful earlier runs are never
overwritten. The final `baseline_summary.json` reports separate local-
prerequisite, I01-I10, and combined M0 gate states.

Individual linked runners require an executable built for their declared
configuration. For example:

```sh
python3 srcMoon/test/integration-tests/I02/test.py \
  --amps ./amps --output-dir test_output/srcMoon/I02
```

Each I-test writes `stdout.log`, `stderr.log`, `receipt.txt`, `result.json`, and
`run_manifest.json`. I01 additionally writes an identical-control repeat and a
single-control-change manifest. I09 freezes seed 731993 and classifies a finite
domain exit as physical lunar escape only when the Moon-relative inertial
specific energy is positive; LSO velocity is transformed as
`v_inertial = v_LSO + Omega × r` before applying that criterion.

U01 runs its inside-sphere guard in a separate process. Acceptance requires
the exact production fatal diagnostic and a zero-byte, incomplete receipt.
The shell return code is only recorded: the legacy AMPS exit helper invokes
`MPI_Abort` with error code zero, which OpenMPI may expose as process status
zero even though the guard terminated execution.

### M1 numerical gates

The reviewed FULL spherical input for M1 is
`srcMoon/examples/m1_sphere.in`. I11 and I13 fail preflight if given QUICK or
CUSTOM because neither may substitute for the frozen convergence campaign.
Focused linked mesh and sphere-inventory runs
use:

```sh
python3 srcMoon/test/integration-tests/I11/test.py \
  --amps ./amps --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I11
python3 srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I12
python3 srcMoon/test/integration-tests/I13/test.py \
  --amps ./amps --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I13
```

For an approximately 20-minute I12 production-wiring run, use the explicitly
coarser QUICK receipt:

```sh
python3 srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps --input srcMoon/examples/m1_quick.in \
  --output-dir test_output/srcMoon/I12-quick
```

The committed QUICK targets are three times the FULL values. The first
factor-two candidate reached 1,074,788 blocks and was stopped after 28.7
minutes at 223.7 GiB RSS in an eight-rank I12 run. Factor three crosses the
next coarser dyadic AMR level. Its cold eight-rank run completed in 1,409.9
seconds (23.50 minutes), with 191.9 GiB peak RSS; the result was correctly
`SKIPPED / NOT QUALIFIED` after every wiring criterion passed. QUICK does not
establish M1 convergence or a FULL I12 result.

They are also included in the default aggregate campaign. Run only the linked
M1 bounded-regression set while reusing a verified SPICE build with:

```sh
python3 srcMoon/test/run_all_test.py --skip-build \
  --only I11 --only I12 --only I13 \
  --campaign-mode regression --mpi-ranks 8 \
  --amps ./amps --m1-input srcMoon/examples/m1_regression.in \
  --output-dir test_output/srcMoon/M1-linked
```

This command must yield two `SKIPPED` non-executions plus the evaluated I12
status; it is not M1 completion evidence. Use `--campaign-mode qualification`
with `m1_sphere.in` for the formal gates.

I11 and I12 are complete executable gates. I13 evaluates the analytical sphere area
even before a baseline exists, but reports `SKIPPED` until a passing I11
campaign supplies equilibrated profiles, return flux, column density, and
brightness and both the science CSV and deterministic configuration JSON are
reviewed. The runner never generates or updates either reference automatically.
I11 preserves separate derived inputs,
logs, manifests, production block metrics, and resource histories for all four
levels. Unequal-epoch, unequal-seed, or unequal-layout runs are not numerical-
convergence evidence.

## Current implementation status

| ID | Scope | Runner state | Meaning |
|---|---|---|---|
| U01 | LOS geometry | implemented linked executable fixture | Production ray limits, record layout, and the inside-origin guard are evaluated after full initialization. |
| U02 | target isolation/wiring | implemented source guard | Configuration evidence only. |
| U03 | lunar gravity | implemented linked probe | Point-mass kernel at a closed-form point. |
| U04 | rotating frame | implemented linked probe | Exact J2000/LSO/MOON_ME_DE421 names, kernel hashes, differential gravity, fictitious terms, transform derivative, and orthogonality are checked at the frozen epoch. |
| U05 | Na radiation pressure/shadow | implemented linked probe | Earth-umbra logic, inverse-square invariant, and absolute Combi-curve agreement to frozen graphical uncertainty; full mover gating remains I07. |
| U06 | Lorentz force | implemented for dedicated ion fixture | Normal Na-only builds honestly SKIP; the M0 runner rebuilds with Na/Na+ and evaluates the production force callback. |
| U07 | refinement selector | implemented linked probe | Legacy no-input callbacks and the input-driven radial, SZA, terminator, and body-fixed polar production kernels are evaluated. |
| U08 | LOLA geometry | implemented linked probe; PASS | Production parser/decoder/icosphere/writers and the strict receipt are exercised; the runner parses D01 readiness metadata and independently re-hashes the immutable IMG/LBL. Linked I14 remains separate. |
| U09 | terrain shadow | implemented linked probe; PASS | Production BVH ray/facet horizon service passes flat, ridge, backface, self-hit, and blocking fixtures; I15 owns D02 comparison. |
| U10 | Diviner temperature | implemented linked probe; PASS | Strict processed-D03 schema, node reproduction, periodic axes, independent trilinear fixture, native RDR UTC/flag-mask logic, radiance-coordinate calibration inverse, Paige C3--C9 substitution, and grey bolometric operator pass; I16/I17 own observations. |
| U11 | thermal inertia | implemented linked probe; PASS | Depth-resolved nonlinear finite-volume solver passes forcing, convergence, continuity, and exact enthalpy-form energy audit; I17 owns Diviner histories. |
| U12 | surface interaction | implemented linked probe; PASS | Na/Ar sticking and residence kernels plus the selectable linear-speed/kinetic-energy accommodation identity are checked; linked stochastic I18 remains separate. |
| U13 | cold trapping | implemented linked probe; PASS | Continuous finite residence, warming release, inventory conservation, classification, and checkpoint round trip pass; I20 owns D02/D03 linkage. |
| U14 | Na photoionization | implemented linked probe | Legacy constant lifetime only. |
| U15 | multispecies photoionization | SKIPPED | Na is the sole active species. |
| U16 | electron impact | SKIPPED | Existing constants are not a qualified D05 implementation. |
| U17 | charge exchange | SKIPPED | No configured production reaction path. |
| U18 | plasma driver | SKIPPED | D06 reports ERROR and no qualified driver is wired. |
| U19 | helium source | SKIPPED | Source file exists but He is disabled and reservoir physics is unresolved. |
| U20 | neon source | SKIPPED | Source file exists but Ne is disabled. |
| U21 | argon source | local uniform kernel verified; full ID SKIPPED | The species-specific continuous source conserves its 3.8e20 atoms/s reference-sphere area partition and rejects unsupported maps. D10 is regional/QA-incomplete, Ar-linked stochastic injection and transient/K-weighted modes are absent, so this is not source validation. |
| U22 | Na sources | implemented linked probe and configuration guard | M0 built-in impact-on/`MySource`-off invariant, impact normalization, and PSD density kernel; no night-to-day reservoir claim. |
| U23 | meteoroid driver | SKIPPED | No production forcing path. |
| U24 | H2O/OH chemistry | SKIPPED | No production species/chemistry path. |
| U25 | data provenance | SKIPPED | Requires a selected Dxx package and linked campaign context. |
| U26 | observation geometry | SKIPPED | The M1 observation-grid settings are parsed strictly, but the complete two-limb, body-fixed pole/shell, finite-beam, and exposure-time operators are not yet implemented or qualified; D14 lacks the radial WHAM response and numerical full raster. |

The status is capability-sensitive.  Downloaded files, directory presence, or
a source/preflight check never changes a missing linked capability into PASS.

## Results and reproducibility

Results go under `test_output/srcMoon/` and never under the validation-data
repository. Each test writes `result.json`; `run_tests.py` writes
`summary.json`, while the individual-runner campaign writes
`run_all_summary.json` and one log per registered runner. Probe compilation
details are retained in the result. Linked I-series runs require the larger
manifest/artifact set specified in the roadmap; I11-I20 preserve that set when
launched by either their individual runners or `run_all_test.py`.

### M2 linked gates

I14-I20 are registered after I13. Each runner checks the I13 `result.json` in
the same campaign directory before launching AMPS. If that predecessor is
absent or not `PASS`, it writes `application.in`, empty `stderr.log`, an
explanatory `stdout.log`, `run_manifest.json` with `execution=not-run`, and
`result.json=SKIPPED`. This prevents the passing U09-U13 kernel layer or an old
I13 artifact from being promoted to current linked evidence.

The bounded M2 preflight with reused builds is:

```bash
python3 srcMoon/test/run_all_test.py --skip-build \
  --m2-amps /path/to/amps-moon-m2-na-ar \
  --only I14 --only I15 --only I16 --only I17 \
  --only I18 --only I19 --only I20 --mpi-ranks 8
```

Formal execution uses `--campaign-mode qualification` and includes I11-I13 in
the same command. I14 requires D01; I15 requires D01+D02; I16/I17 require the
appropriate D03 products; I20 requires D02+D03. Directory presence is not
readiness: package README, provenance, QA, hashes, native labels, and the exact
processed product are checked. I18 uses the input-selected Tenishev et al.
(2013) Na law and Kegerreis et al. (2017) Ar law through a dedicated Na+Ar
build, but remains predecessor-gated until the same FULL campaign produces
I13 PASS. Missing physics or evidence is `SKIPPED`, never `PASS`.

When an I14-I20 science criterion is actually evaluated, its runner must
preserve the row-level values as `comparison.csv` and render those exact rows
as both `comparison.png` and `comparison.eps`. The strict renderer is
`common/write_comparison_artifacts.py`; it rejects non-finite values, mixed
units within one comparison ID, negative uncertainty, and undeclared selection
labels. Multiple named comparisons are rendered as separate panels without
changing their row order. It never interpolates or changes the CSV. A
`SKIPPED`/not-run test must not
create placeholder comparison files, because a graphic without evaluated
evidence could be mistaken for validation.

See [`TEST_GAP_PROPOSALS.md`](TEST_GAP_PROPOSALS.md) for missing coverage.  No
new permanent IDs are assigned there.

### M2 reference-data inventory (audited 2026-10-13)

D02 currently contains immutable north/south PDS3 IMG/LBL pairs for
`AVGVISIB_65N/65S_240M_201608` and `LPSR_65N/65S_240M_201608`. The native
labels define 6420-by-6420 signed 16-bit little-endian grids, 240 m/pixel polar
stereographic projection, planetocentric latitude, positive-east longitude,
the DE421 mean-Earth/polar-axis body-fixed frame, and the exact DN scale and
offset. AVGVISIB covers the 2009-07-13 through 2013-07-18 calculation interval
and maps average solar-disc visibility continuously toward zero; LPSR is the
binary permanent-shadow classification. D02's immutable raw package remains
separately labelled `PARTIAL`; the committed `prepare_d02_reference.py` now
reproducibly creates campaign-owned `illumination_reference.csv` and
`d02_processed_qa.json` files on a predeclared approximately 100 km lattice
without changing `raw/`. This makes the derived reference inspectable but
does not supply the missing production extended-Sun temporal operator required
for I15.

D03 currently contains all 18 version-1.0 global cumulative-product latitude
bands (`global_cumul_avg_cyl_*_002`), both version-1.0 polar resource products
(`dlre_prp_north/south`), and all 144 ten-minute ZIP/XML RDR products for each
frozen 2011-06-15 training and 2011-12-10 hold-out date. The Diviner team
identifies these as RDR processing version 3.1; the individual PDS4 product
labels retain product version 1.0. The GCP labels provide 0.5-degree latitude/longitude
by 0.25-hour local-time bins, calibrated channel-3 through channel-9 and
bolometric brightness temperatures in kelvin, positive-east/planetocentric
coordinates, and `-9998/-9999` invalid/unknown constants. PRP supplies annual
average 2-cm temperature, annual maximum surface temperature, and ice-depth
fields on polar triangles. RDR labels provide UTC/JD, geometry, channel,
detector, radiance/brightness-temperature fields, and native quality/status
metadata. The committed GCP processor creates a hashed
`temperature_table.csv`, hold-out `benchmark_points.csv`, and a processing-QA
record; the PRP processor likewise creates a hashed polar reference. The I17
processor selects native footprints/flags, aggregates radiance before
temperature, applies the published Paige C3-C9 sensitivity/T_bol algorithm,
and derives uncertainty from calibration plus within-bin dispersion. GCP/PRP
campaign records and the canonical `diviner/processed/i17` RDR product qualify
only their deterministic transformations. They do not make D03 READY: the
canonical `diviner/` package currently has no
`README.md`, `provenance.json`, or `qa_report.json`. The linked runner checks
those three package records before reading a processed D03 table and therefore
reports I16/I17/I20 as `SKIPPED / NOT VALIDATED`. I17 now has a local
multi-channel/grey-bolometric production operator, reproducible RDR processor,
input-controlled finite-disc eclipse forcing, observation-free periodic
full-column spin-up, and a linked June-only training comparison path.
Canonical processed RDR evidence is written only by an explicit
data-maintenance command under `diviner/processed/i17`; the linked runner
verifies and stages it but never rewrites the validation-data package. Before
threshold freeze, staging exposes only the 2,039 June rows to the executable.
A completed run writes `comparison.csv`, `comparison.png`, `comparison.eps`,
`channel_diagnostics.csv`, and a review-required
`june_threshold_candidate.json`, then reports `SKIPPED / NOT VALIDATED` rather
than PASS. I17 still requires D03 package qualification, a same-campaign FULL
I13 PASS, investigator review/freeze of the June threshold, and a later
untouched December execution. Raw RDR presence or an experimental plot alone
never establishes I17.

The canonical 2026-10-13 I17 processing receipt is `READY` for the deterministic
RDR transformation itself: 77,075,901 records from 87 archives were examined,
25,848,333 passed the frozen native cuts, and 6,353 complete C3--C9 bins were
published (2,039 June training and 4,314 December hold-out).  The history,
contract, and processor SHA-256 values are respectively
`870577c82b430b098d816d210b528c714ec6cfa3a20ae9be3f6482a3639ae182`,
`88fc5196fee5863691937a505873d420990b98056db8a1e8167255cd736189b0`, and
`5a5d29d0f6dba88cd49c899e541a980eb84ce24717481b0cad471f8012459a65`.
This receipt does not override the missing D03 package qualification or linked
I17 capability gates and therefore is not observational validation evidence.
