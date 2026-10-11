# SA01 — 6–9 January 2014 RBSP/REPT validation

This directory runs the standalone AMPS T05 energetic-proton validation for the
half-open interval `2014-01-06T09:00:00Z <= t < 2014-01-10T00:00:00Z`. RBSP-A
and RBSP-B REPT are the validation targets; SOHO/ERNE supplies the independent
upstream spectrum and Qin-Denton supplies the five-minute T05 history. The
baseline boundary is isotropic and spatially uniform, and the electric field is
disabled.

The accepted data argument is either the extracted package root or its
`validation_data` directory. Despite the compatibility name `--data-archive`, it
must be a directory, not a tar file. Generated decks always contain absolute
driver, boundary-spectrum, and trajectory paths.

## Products and interpretation

The comparison covers seven differential channels (18.5–116 MeV), their seven
omnidirectional channel integrals, the total omnidirectional integral intensity
over that same finite interval, and energetic-proton density. REPT density is
the isotropic-equivalent quantity derived from the seven spectra,
`4*pi*sum(J_i*DeltaE_i/v_rel(E_i))`; it is not a direct measurement or an
independent observation dataset. Rectangular bow-tie folding is an observation
operator approximation, not a full calibrated response simulation.

Physical-zero model values remain failed predictions against positive REPT
values. Nonfinite values, duplicate times, incomplete products, missing modeled
epochs, and excessive unresolved fractions are errors/failures. Successive
static T05 snapshots with no electric field do not establish dynamic
electric-field transport. REPT points are five-minute means while AMPS evaluates
the selected instantaneous epoch.

All four products use the fixed gates in `AGENTS.md`: median absolute log10 ratio
at most 0.30, factor-two fraction at least 0.60 per spacecraft and 0.70 pooled,
at least 12 matched epochs per product/spacecraft, unresolved fraction at most
0.01, and integrated-product convergence at most 0.02 on common epochs.

## Commands

From the repository root, a full STANDARD run is:

```bash
python3 srcEarth/test/SA01/run_SA01.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --amps-root /data/vtenishe/GEO/AMPS \
  --amps /data/vtenishe/GEO/AMPS/amps \
  --output-dir /data/vtenishe/GEO/AMPS/test_output/SA01/standard \
  --profile STANDARD --mpi-ranks 4 --threads 16 --formats png eps
```

Preparation without AMPS execution:

```bash
python3 srcEarth/test/SA01/run_SA01.py --prepare-only \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --output-dir test_output/SA01/prepared-standard --profile STANDARD \
  --mpi-ranks 4 --threads 16
```

Run already-prepared branches sequentially with `--run-only`. Standalone modes
never launch or rebuild AMPS:

```bash
python3 srcEarth/test/SA01/run_SA01.py --run-only --data-root /path/to/package \
  --output-dir test_output/SA01/prepared-standard --mpi-ranks 4 --threads 16
python3 srcEarth/test/SA01/run_SA01.py --compare-only --data-root /path/to/package \
  --output-dir test_output/SA01/standard
python3 srcEarth/test/SA01/run_SA01.py --plot-only --data-root /path/to/package \
  --output-dir test_output/SA01/standard --formats png eps
```

QUICK and PUBLICATION use the same gates and physics, changing only workload:

```bash
python3 srcEarth/test/SA01/run_SA01.py --data-root /path/to/package \
  --output-dir test_output/SA01/quick --profile QUICK --mpi-ranks 4 --threads 16
python3 srcEarth/test/SA01/run_SA01.py --data-root /path/to/package \
  --output-dir test_output/SA01/publication --profile PUBLICATION --mpi-ranks 4 --threads 16
```

### Independent QUICK test and STANDARD study commands

QUICK and STANDARD are registered as independent complete workflows. They do not
wait for, read, overwrite, or supply convergence evidence for one another:

```bash
# Direct QUICK campaign; this is also the exact command in srcEarth/test/list.
python3 srcEarth/test/SA01/run_SA01.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --amps-root /data/vtenishe/GEO/AMPS --amps ./amps \
  --output-parent test_output/SA01/quick_test_runs \
  --profile QUICK --solver GRIDDED --mpi-ranks 4 --threads 6 \
  --mover HC4 --trace-dt 0.1 --formats png eps \
  --exit-code-policy WORKFLOW

# The independent command in srcEarth/studies/list.
python3 srcEarth/test/SA01/run_SA01.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --amps-root /data/vtenishe/GEO/AMPS --amps ./amps \
  --output-parent test_output/SA01/standard_study_runs \
  --profile STANDARD --solver GRIDDED --mpi-ranks 4 --threads 8 \
  --mover HC4 --trace-dt 0.1 --formats png eps \
  --exit-code-policy WORKFLOW
```

`--output-parent` is available only for the complete prepare/run/compare/plot mode.
It atomically reserves a new child named with profile, UTC start time, and PID, so a
repeated list invocation cannot reuse or erase an earlier campaign. Manual stage-only
commands continue to require the exact existing `--output-dir`.

`--exit-code-policy WORKFLOW` means exit zero only after both physical branches,
comparison evaluation, the Markdown report, and all requested plot files exist and
are hashed in `workflow_completion.json`. It does **not** turn scientific disagreement
or missing convergence into qualification PASS: `comparison_summary.json`,
`qualification_summary.json`, and the report retain the unchanged gates and actual
`PASS`, `FAIL`, or `INCOMPLETE` status. The default `QUALIFICATION` exit policy remains
appropriate when a manual caller wants scientific nonqualification to return status 2.

The QUICK master-list entry is marked `! runner: exclusive`, so `test_runner.py`
waits until its process pool is empty and launches no other test until QUICK exits.
Exclusivity is deliberately a scheduler property rather than a private command-line
prerequisite: the exact command above remains independently executable, with
standalone resource isolation under the operator's control. The STANDARD studies
command is likewise a direct invocation and is not coupled to QUICK.

`--rebuild` first validates the exact repository and non-symlink `build` target,
checks for active AMPS/build processes, removes only `build`, and runs the mandated
`make -j`, saving its log. It is incompatible with compare-only and plot-only.

Numerical refinement is evaluated independently on common epochs:

```bash
python3 srcEarth/test/SA01/compare_SA01.py --data-root /path/to/package \
  --output-dir test_output/SA01/standard \
  --convergence-with test_output/SA01/publication
```

Run implementation tests with:

```bash
python3 -m unittest discover -s srcEarth/test/SA01/tests -v
```

The permanent, data-independent regression suite remains registered in
`srcEarth/test/SA01/list` for focused development:

```bash
python3 srcEarth/test/test_runner.py -j 1 srcEarth/test/SA01/list
```

The master `srcEarth/test/list` additionally registers the full QUICK physical
workflow requested for normal code validation. Its reviewed structural digest in the
Step-12 source contract was advanced to include that command; the digest check and all
scientific gates remain active. The focused routine suite is not a replacement for
the physical campaign. It includes the production field-alias ownership regression
discovered by SA01, so the correction remains checked without the external dataset.

## Output layout

Each selected output directory contains immutable preparation/data-validation
records and separate `RBSPA/` and `RBSPB/` branches. A branch contains its
rendered input and selected trajectory, `amps.log`, `run_status.json`, the five
native AMPS products, comparison CSV files, and independent reconstruction
residuals. The root contains `comparison_summary.json`, `figures/`, environment
provenance, `SA01_validation_report.md`, `qualification_summary.json`, and, after a
complete workflow, `workflow_completion.json`. Select a new `--output-dir` for a
manual physical run or use `--output-parent` to reserve one automatically; the runner
refuses to replace prior products.

Standalone comparisons record their own analysis-source hashes. `plot_manifest.json`
binds every PNG/EPS to its source implementation and the saved comparison/context
inputs, even when rendering occurs after the physical executable run. Time-energy
maps mark physical-zero predictions with crosses; unmarked blank cells are missing.
Orbit/upstream context curves break across genuine gaps just like the model series.

## GRIDDED Mode3D branch

`--solver GRIDDED` evaluates the same selected trajectory epochs with the
mesh-backed Mode3D tracer. It is a separately labelled numerical branch, not an
interchangeable replacement for the direct GRIDLESS baseline until a mesh-resolution
convergence study satisfies the unchanged 0.02 integrated-product gate.

The runner writes `TEMPORAL_MODE SNAPSHOT_LIST` and an explicit epoch list. Mode3D
therefore builds one immutable T05 field snapshot, selects only the trajectory row
at that epoch, runs cutoff and density/flux, and then advances to the next epoch.
Using ordinary `TIME_SERIES` here would incorrectly evaluate every location at every
field epoch.

The initial QUICK diagnostic uses explicit `0.05 Re` near-Earth and `2.0 Re` outer
boundary resolution with linear coarsening. Override those values only as a recorded
mesh refinement, for example:

```bash
python3 srcEarth/test/SA01/run_SA01.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --amps-root /data/vtenishe/GEO/AMPS --amps /data/vtenishe/GEO/AMPS/amps \
  --output-dir test_output/SA01/quick_gridded_005re \
  --profile QUICK --solver GRIDDED --mpi-ranks 4 --threads 8 \
  --mode3d-mesh-res-earth-re 0.05 \
  --mode3d-mesh-res-boundary-re 2.0 \
  --mode3d-mesh-coarsening LINEAR --mode3d-mesh-exponent 1.0
```

The GRIDDED launch enables parallel field initialization and uses the OpenMP
backtracking backend. It also writes and passes
`MODE3D_TRAJECTORY_SNAPSHOT_SCHEDULER EPOCH_PARALLEL` explicitly. The application-wide
default is `COLLECTIVE`; SA01 opts into private meshes because its existing scheduling,
repeatability evidence, and one-location-per-epoch work decomposition were built for
this mode. The choice is recorded in both `prepared_run.json` and the executable
command, so a future default change cannot silently alter campaign ownership.

With `EPOCH_PARALLEL`, a world-communicator atomic counter assigns one complete
epoch/location transaction at a time to each MPI process, while every process owns a
private mesh and gridded T05 field for its claimed epoch. OpenMP then distributes that
location's energy/direction work among the requested workers. A process claims another
epoch only after both cutoff and density transactions close their uniquely suffixed
products. This outer epoch policy is separate from `MODE3D_MPI_SCHEDULER DYNAMIC`,
which divides trajectory work after a field has been prepared.

This epoch ownership is specific to timestamped trajectories. `POINTS` and locations
within one `SHELLS` snapshot share a field, so density uses finer flattened
`(location,energy)` MPI tasks there; cutoff retains its flattened scalar/directional/
rigidity task scheduler. For `SHELLS + SNAPSHOT_LIST`, field snapshots are still
advanced collectively, and fine task scheduling is repeated within each epoch.

The trajectory implementation keeps physics collectives on a singleton
`MPI_GLOBAL_COMMUNICATOR` and reserves `MPI_COMM_WORLD` for the outer epoch queue and
final status. This prevents fields from different epochs from entering the same
reduction. It increases memory to one complete mesh/compact field per MPI process and
serializes legacy mesh-cache initialization once to avoid concurrent writes; expensive
epoch field builds and trajectory calculations remain concurrent.

Before ranks beyond rank 0 allocate their private domains, Linux runs measure the
first domain's resident-memory increase and project the remaining replicas for each
shared-memory node. SA01 records the default
`MODE3D_EPOCH_PARALLEL_MEMORY_FRACTION 0.80`; if projected replicas exceed that share
of `MemAvailable`, startup fails collectively and recommends `COLLECTIVE`, fewer
ranks, or a separately justified mesh/resource change. This guard changes no SA01
physics or acceptance gate. If `/proc` is unavailable, the log says the preflight was
unavailable and makes no safety claim.

The preserved primary trajectory limit is 300 s.  To close the TIME_LIMIT population
without weakening `DS_UNRESOLVED_TOL=0.01`, SA01 enables three unresolved-only
extensions with factor two (600, 1200, then 2400 s).  Only directions still ending at
`TIME_LIMIT` or `STEP_LIMIT` are recomputed; the direction/energy grids and observation
population do not change.  SA01 also enables the C19-qualified drift-shell recurrence
detector.  A trajectory becomes physically trapped only after the configured
azimuth-resolved profile repeats for three complete drift revolutions and passes its
coverage, phase-space, and secular-radius checks.  An arbitrary timeout is never mapped
to forbidden.  The strict unresolved gate remains fatal after all bounded extensions.

If the final gate still fails, preserve that output root as a numerical failure and run
a separately labelled trace-limit or mesh-resolution refinement.  Do not disable
`DS_FAIL_ON_UNRESOLVED`, lower the tolerance, or omit the affected energy channels.

Mode3D writes one suffixed cutoff, spectrum, flux, density, termination, and manifest
file per epoch. The runner requires the count of every artifact class to equal the
prepared trajectory count; comparison mode reads all suffixed files and rejects
duplicate or missing epochs.

## Explicit numerical controls and failure diagnosis

`--mover RK4|BORIS|HC4|RK6` and `--trace-dt SECONDS` are saved with the
preparation and passed explicitly to AMPS. `--run-only` uses those saved choices.
Defaults preserve the historical RK4 / 0.1 s configuration. Boris preserves
momentum magnitude in this static magnetic-only problem; it is an explicit
numerical alternative, not a changed physical model or a relaxed trapping gate.
The 2026-10-09 investigation reproduced a TIME_LIMIT direction with 2.83%
RK4 momentum drift after 2400 s, above the unchanged 1e-4 trapping tolerance.
Matched Boris diagnostics closed both difficult epochs with zero unresolved;
see `SA01_PROGRESS.md` for subsequent refinement and full-run evidence.

For a bounded reproduction from a preserved prepared branch:

```bash
python3 -B srcEarth/test/SA01/diagnose_SA01.py \
  --source-branch test_output/SA01/quick_gridded_closure_20261009/RBSPA \
  --output-dir test_output/SA01/new-boris-diagnostic \
  --epochs 2014-01-06T10:00:00Z 2014-01-06T11:00:00Z \
  --mover BORIS --mpi-ranks 2 --threads 8
```

This keeps the full energy/direction grids and strict unresolved gate. It writes
`diagnostic_status.json` and opt-in per-direction trace records; it is not a
substitute for a complete observational branch.

Compare completed bounded experiments separately from campaign qualification:

```bash
python3 -B srcEarth/test/SA01/analyze_diagnostics.py \
  --coarse test_output/SA01/closure_diagnosis_20261009/hc4_reference \
  --fine test_output/SA01/closure_diagnosis_20261009/hc4_dt0001 \
  --output test_output/SA01/closure_diagnosis_20261009/hc4_timestep_comparison.json
```

The retained example has a **failed numerical-difference gate**, not recommended
production settings: despite zero unresolved at both timesteps, its worst channel
changes 17.5%. A subsequent identical-settings HC4 repeat changed a channel by
20.3376%; consequently the earlier timestep comparisons are confounded and cannot
be interpreted as isolated convergence evidence. The subsequent investigation
identified shared interior/ghost center-buffer writes during field initialization;
the unique-writer fix and exact physical controls are described below. The helper
records all nine integrated quantities, artifact hashes and direction-classification
changes, and labels its small population explicitly. Use the progress journal for
the latest tested settings; do not choose a mover merely because it closes traces.

`--cadence-minutes 60` can hold the QUICK population fixed while using STANDARD
or PUBLICATION energy/angular resolution. Changing cadence alone is rejected as
convergence evidence. Refinement checks all seven native channel integrals, total
flux and density, requires at least twelve common epochs per spacecraft, and
reports unshared epochs. Physical inputs must match between refinements.

The generic `DS_DIRECTION_NMU` / `DS_DIRECTION_NPHI` controls define the density
equal-solid-angle parent grid (defaults 24x48). QUICK and STANDARD retain that
parent grid; PUBLICATION explicitly uses 32x48 so all requested 1536 directions
are traced instead of being silently capped at 1152. `DS_MAX_PARTICLES` still
caps the total trajectories per observation, including every energy node.

Comparison requires input/output hashes to match the completed run, and requires
every prepared epoch in all three physical product tables. Reports distinguish
observational status from overall qualification: missing convergence cannot be
PASS. A saved `convergence_summary.csv` adds a separate refinement figure; native
spectral reconstruction is not mislabeled as resolution convergence. Physical
zeros remain in metrics and are explicitly annotated where logarithmic axes
cannot display them.

### Reproducibility before refinement

The diagnostic runner enables two read-only production probes with
`AMPS_MODE3D_FIELD_AUDIT=1`. The owner-field audit samples eight interior cell
centers per owned block and compares stored B with a fresh serial evaluation at
the exact initialization coordinates. The interpolation audit compares frozen-field
samples at 4096 deterministic interior positions across three concurrent workers
and three passes. A content fingerprint covers **all** packed B/E values, epoch,
and leaf geometry; the historical driver-derived snapshot ID stays unchanged.
Sampling probes are not proofs of correctness throughout the domain, and a
content hash is an identity check, not a magnetic-field accuracy check.

Run two fresh copies of the same bounded diagnostic, then request exact checks:

```bash
python3 -B srcEarth/test/SA01/analyze_diagnostics.py \
  --coarse test_output/SA01/audit-a --fine test_output/SA01/audit-b \
  --check-repeatability --output test_output/SA01/audit-repeatability.json
```

This mode requires complete field audits, identical field/geometry fingerprints,
zero sampled discrepancies, exact equality of all nine integrated quantities,
and matching final direction records. Log coverage is checked independently against
native spectrum sample counts at every energy; two equally truncated logs cannot
pass. The logger derives energy from rigidity, so that *identity mapping only*
allows eight machine epsilons for floating-point conversion roundoff. Product and
trajectory-repeat comparisons remain exact. Missing evidence is an error, not PASS.
The physical 2% refinement threshold is deliberately **not** used for repeats.
For ordinary bounded mover/timestep comparisons the checker also requires identical,
clean field evidence before reporting PASS. Older runs without audits are INCOMPLETE;
their raw integrated-difference gate is still recorded, but cannot establish an
isolated integration-refinement result. This helper does not perform mesh refinement.

To isolate concurrency, `diagnose_SA01.py` accepts independently:

- `--field-init SERIAL|PARALLEL`: construction of the sampled field;
- `--density-backend SERIAL|OPENMP|THREADS`: production trajectory execution;
- `--mpi-ranks N --threads N`: process/worker layout.
- `--mpi-binding DEFAULT|NONE`: retain launcher defaults, or explicitly pass
  Open MPI `--bind-to none` to avoid independent small jobs contending on its
  default bound cores. This is a recorded resource choice, not a numerical control;
  do not infer actual CPU utilization from requested thread counts alone.

For a serial control use `--field-init SERIAL --density-backend SERIAL
--mpi-ranks 1 --threads 1`. Keep source deck, epochs, mover, timestep, trace
budgets and executable fixed when changing layout. Field initialization has its
own worker mechanism, so setting just `--threads 1` is not a serial-field test.
Use new output roots for every experiment. Read the progress journal before
running costly cases or rebuilding; do not rebuild over active physical jobs.

The corrected initializer selects one writer per physical center-node buffer before
starting workers, preferring interior-owner coordinates to ghost coordinates. The
production-header regression demonstrates the old order-dependent alias behavior
and rejects duplicate interior owners. Corrected physical repeats then matched
exactly across one/two ranks, one/eight workers and serial/parallel initialization
(4800 final trajectory records and all nine integrated quantities at two epochs).
This establishes bounded reproducibility, not integration or mesh accuracy. On the
same corrected fields, HC4 dt=0.1 versus 0.001 still fails the 2% timestep gate with
18.1598% maximum channel change. The default RK4 control still exceeds the strict
unresolved limit; neither failure is hidden by selecting a trace-closing mover.

Register relevant bounded results in a campaign report using repeatable
`--diagnostic-evidence PATH` arguments (for example with `--compare-only`). This
writes `diagnostic_evidence.json`, bound to the prepared campaign hashes. Later
reporting retains those checks even when the option is omitted, and revalidates the
comparison/native artifact hashes. A registered FAIL overrides a PASS on another
refinement axis; a bounded PASS never supplies the required full-population coverage.
Registration asserts relevance explicitly: do not attach an old numerical
configuration's failure to a different configuration without explaining why it
applies. Changing settings requires a new campaign root, not deleting failed checks.

### Recorded sequential campaign continuation

For a long campaign, the full runner accepts `--wait-for-campaign QUICK_ROOT` and
repeatable `--wait-for-diagnostic DIAGNOSTIC_ROOT` arguments. It waits for each
recorded runner's atomic terminal evidence before preparing or launching new work.
PID start times guard against PID reuse; process disappearance without terminal
evidence is an error, never execution success. Existing jobs are not signalled.
A sibling `OUTPUT_ROOT.continuation.json` records the waiting process, command and
transition; exclusive creation prevents duplicate continuations targeting one root.
An existing terminal failure remains a failure but does not prohibit an independent
new campaign. The normal strict completeness gates still reject partial comparisons.

Add `--convergence-base QUICK_ROOT` to evaluate the finished STANDARD products against
the coarse campaign at all common epochs and automatically regenerate both campaigns'
PNG/EPS figures and reports. Energy/angular changes are recorded separately from
cadence; no claim of timestep or mesh convergence follows from sampling refinement.
The fine and coarse reports both retain registered bounded failures. An incomplete
coarse run produces explicit `refinement_error.json`; it does not suppress valid
fine-run observational figures or supply missing convergence credit.
`qualification_summary.json` provides the final machine-readable qualification,
separate from `physical_execution_summary.json` (execution only).

Example continuation retaining an explicitly unqualified HC4 dt=0.1 baseline:

```bash
python3 srcEarth/test/SA01/run_SA01.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --output-dir test_output/SA01/new-standard-hc4-baseline \
  --profile STANDARD --solver GRIDDED --mover HC4 --trace-dt 0.1 \
  --mpi-ranks 4 --threads 8 --formats png eps \
  --wait-for-campaign test_output/SA01/quick_unique_hc4_20261009 \
  --wait-for-diagnostic test_output/SA01/closure_diagnosis_20261009/hc4_unique_dt00005 \
  --convergence-base test_output/SA01/quick_unique_hc4_20261009 \
  --diagnostic-evidence test_output/SA01/closure_diagnosis_20261009/post_fix_timestep_comparison.json
```

This example deliberately preserves a known numerical FAIL for the baseline; it is
not a recommendation to qualify dt=0.1. Further timestep and mesh studies need new
recorded controls/populations, not a changed acceptance threshold.

The permanent workflow distinguishes regression fixtures, exact repeatability,
controlled numerical refinement, and observational agreement. A trustworthy
scientific FAIL is retained as such. No successful software test substitutes for
the full QUICK/STANDARD spacecraft populations or the fixed acceptance gates.
Reports revalidate both convergence populations and their recorded status hashes;
copying a stale convergence PASS into a newer campaign cannot qualify that run.
The full prepare/run/compare/plot command exits nonzero for FAIL or INCOMPLETE,
including missing numerical convergence. Explicit `--run-only`, `--compare-only`
and `--plot-only` modes retain their narrower execution/comparison/rendering meanings
and never claim full qualification. After adding convergence evidence, refresh the
report with `--compare-only`; its displayed qualification still enforces both ships.
