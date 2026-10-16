# SA04 — 10–13 September 2017 GLE72 frozen holdout

This directory is the independent SA04 implementation of the AMPS observational
validation workflow. It covers the half-open interval
`2017-09-10T00:00:00Z <= t < 2017-09-14T00:00:00Z`. This is the independent
high-energy-spectrum holdout for GLE72: its observations must not be used to tune
field choices, transport parameters, event selection, profile workloads, or gates.
RBSP-B REPT is the sole validation target, SOHO/ERNE supplies the upstream proton
spectrum, Qin-Denton supplies the five-minute T05 history, and OMNI supplies plotted
solar-wind and geomagnetic context. The baseline boundary is isotropic and spatially
uniform, and the electric field is disabled.

SA04 deliberately has its own runner, preparation, comparison, plotting, diagnostic,
tests, and documentation. It imports no SA01, SA02, or SA03 modules and shares none
of their output directories. Changes to this workflow therefore cannot silently
change any prior case, and SA01–SA03 remain independently executable.

RBSP-A is deliberately not a model/comparison branch. At `L >= 3.2`, the frozen
provider-quality audit contains zero flag-0 spins, 24,126 flag-2 (`bias_st off`)
spins, and two flag-3 spins; it yields no eligible trajectory, spectrum, or density
row. The raw/processed RBSP-A files remain hash-validated audit evidence. Neither a
CLI option nor a plotting choice can promote them into the SA04 validation population.

## Authoritative data contract

Pass the supplied `model_inputs/SA04` directory, the enclosing `model_inputs`
directory, the extracted package root, or its `validation_data` directory to
`--data-root`. The compatibility alias
`--data-archive` also expects one of those directories, not a tar file. The currently
registered repository commands use the location supplied for this case:

```text
/data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926/model_inputs/SA04
```

`model_inputs/SA04/RBSPB` contains a pre-rendered GRIDLESS STANDARD example with 202
epochs. It predates the strict seven-channel observation-completeness intersection
and cannot represent QUICK. The live runner therefore treats `model_inputs/SA04` only
as the package entry point, resolves authenticated sibling `validation_data`, and
regenerates the requested GRIDDED/HC4 branch with the current 51-epoch QUICK or
201-epoch STANDARD population. It never executes or silently upgrades the old deck.
`preparation_summary.json` records both the supplied entry path and resolved data root
and sets `pre_rendered_model_input_policy` to
`LOCATOR_ONLY_REGENERATE_FROM_VALIDATION_DATA`.

The workflow validates `raw/SA04/source_manifest.json` and
`processed/SA04/provenance.json` before preparing a run. It rejects changed hashes,
unexpected raw files, schema changes, a changed event window, non-monotonic time,
wrong energy bins, and recorded row-count mismatches. Generated decks use absolute
paths for the driver, boundary spectrum, and selected trajectory.

The frozen package records the following SA04 inventory:

- 1,152 consecutive five-minute T05 driver records, with a maximum internal gap of
  300 seconds.
- 784 complete ERNE boundary spectra. A real 30-hour-45-minute gap between the
  retained rows at `2017-09-10T03:45:00Z` and `2017-09-11T10:30:00Z` is not filled.
  The package reports one gap exceeding 15 minutes, no omitted incomplete spectrum,
  and only two channel values filled by its fixed, bracketed, at-most-15-minute rule.
- RBSP-B: 871 eligible trajectory samples, 6,088 spectrum rows, and 862 complete
  density rows.
- RBSP-A audit only: zero eligible trajectory samples, spectrum rows, and complete
  density rows because of the provider quality flags described above.
- The Qin-Denton provenance reports zero fixed-width `W5/W6` structural recoveries for
  this case; the runner does not repeat or extend upstream processing.

The raw CDF/text files remain the source record. The processed audit tables are
inputs to the reproducible AMPS workflow, not replacements for the raw observations.
See the package's `validation_data/README_SA04.md`, `validation_data/README.md`, and
top-level `README.md` for the upstream acquisition, units, selection rules, and case
definition.

### Comparison-epoch eligibility

The packaged spacecraft trajectory tables are deliberately broader than the valid
comparison population. A trajectory row supplies position for a five-minute bin, but
does not imply that every REPT channel survived the frozen provider-quality and
positive-sample rules. SA04 therefore determines eligibility from immutable input
structure before looking at any modeled value. An epoch is eligible only when all of
the following are present at the exact timestamp:

- the profile's fixed cadence selects the trajectory row;
- a complete SOHO/ERNE boundary spectrum exists;
- each of the seven declared REPT differential channels occurs exactly once; and
- the corresponding seven-channel spectrum-derived density row exists.

Preparation writes every rejected trajectory timestamp and one exact reason to
`prepared_run.json`; `incomplete_rept_observation` is distinct from cadence and
boundary-spectrum exclusions. It also records the complete/incomplete observational
inventory. No flux magnitude, model/reference ratio, or acceptance result participates
in this selection, so the policy cannot preferentially retain agreement.

The data-package validator checks channel identity and uniqueness per epoch and
requires the density-epoch set to equal the complete seven-channel spectrum set.
Aggregate file row counts remain checked as provenance, but are not accepted as proof
of per-epoch completeness.

Campaigns prepared by an older runner may contain additional, physically complete
modeled epochs whose REPT reference is incomplete. Comparison still authenticates
and validates all native artifacts and includes their unresolved fraction in the
unchanged physical gate. It lists those timestamps and missing channels under
`observation_ineligible_modeled_epochs`, but does not interpolate a missing reference
or insert them into an observational denominator. This permits saved expensive runs
to be analyzed without rewriting preparation evidence. Newly prepared runs omit the
same ineligible epochs before AMPS execution.

For the frozen package, QUICK selects 51 RBSP-B epochs (797 cadence exclusions and 23
missing exact boundary spectra). STANDARD selects 201 (578 cadence exclusions, 91
missing boundary spectra, and one incomplete REPT observation). These are
pre-model structural counts, not outcome-dependent tuning. The older packaged
STANDARD example has 202 epochs because it predates the seven-channel eligibility
intersection; the live workflow intentionally uses the stricter 201-epoch population.

## Profiles and unchanged scientific gates

QUICK and STANDARD execute the same physical products and the same acceptance gates.
Only their sampling workload differs. PUBLICATION remains available for separately
scheduled refinement; it is not registered in either repository list.

| Profile | Trajectory cadence | Spectrum intervals | Cutoff nodes | Directions |
|---|---:|---:|---:|---:|
| QUICK | 60 min | 24 | 24 | 96 |
| STANDARD | 15 min | 48 | 50 | 384 |
| PUBLICATION | 5 min | 96 | 120 | 1,536 |

Every profile evaluates seven REPT differential channels from 18.5 to 116 MeV,
their seven channel integrals, the total omnidirectional integral over that finite
energy interval, and energetic-proton density. REPT density is the
isotropic-equivalent value reconstructed from those seven spectra,
`4*pi*sum(J_i*DeltaE_i/v_rel(E_i))`; it is not a separate direct measurement.
Rectangular bow-tie folding is an explicit observation-operator approximation, not a
full calibrated response simulation.

The fixed gates are:

- median absolute `log10(model/reference) <= 0.30`;
- factor-two fraction at least 0.60 for eligible RBSP-B and 0.70 pooled;
- at least 12 matched epochs for every RBSP-B product;
- maximum unresolved fraction at most 0.01; and
- integrated-product numerical refinement change at most 0.02 on common epochs.

Physical-zero model predictions remain in the comparison denominator and fail
factor-two agreement against positive observations. Nonfinite values, duplicate
times, incomplete physical artifacts, missing eligible model epochs, invalid units,
excessive unresolved support, and stale or mismatched evidence fail closed. A
cadence-only change cannot establish numerical convergence.

These are static T05 snapshots with no electric field. Agreement or disagreement
therefore tests the documented quasi-static model and must not be interpreted as
validation of dynamic electric-field transport. REPT values are five-minute means;
AMPS evaluates the selected instantaneous epoch.

## Independent QUICK test and STANDARD study

The master test list owns QUICK. The studies list owns STANDARD. Each command performs
preparation, physical execution for the eligible RBSP-B branch, comparison, report
generation, and PNG/EPS plotting in a newly reserved evidence directory.

```bash
# Direct QUICK campaign; this is also the exact command in srcEarth/test/list.
python3 srcEarth/test/SA04/run_SA04.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926/model_inputs/SA04 \
  --amps-root /data/vtenishe/GEO/AMPS --amps ./amps \
  --output-parent test_output/SA04/quick_test_runs \
  --profile QUICK --solver GRIDDED --mpi-ranks 4 --threads 6 \
  --mover HC4 --trace-dt 0.1 --formats png eps \
  --exit-code-policy WORKFLOW

# Independent STANDARD campaign; this is the command in srcEarth/studies/list.
python3 srcEarth/test/SA04/run_SA04.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926/model_inputs/SA04 \
  --amps-root /data/vtenishe/GEO/AMPS --amps ./amps \
  --output-parent test_output/SA04/standard_study_runs \
  --profile STANDARD --solver GRIDDED --mpi-ranks 4 --threads 8 \
  --mover HC4 --trace-dt 0.1 --formats png eps \
  --exit-code-policy WORKFLOW
```

Both commands are directly executable. The QUICK entry's immediately preceding
`! runner: exclusive` directive tells `test_runner.py` to drain its worker pool before
starting this high-memory campaign and to start no other tests until it exits. That
serialization is runner policy, not a hidden environment-variable or command-line
precondition. A direct invocation is responsible for its own machine-level resource
isolation.

`--output-parent` is valid only for the complete workflow. The runner atomically
creates a child named with the profile, UTC start time, and process ID; it will not
reuse or erase an older campaign. Stage-only commands instead require the exact
existing `--output-dir`.

`--exit-code-policy WORKFLOW` returns success only when the physical branch,
comparison, Markdown report, requested figures, and their evidence hashes are
complete. It never converts scientific disagreement into a PASS.
`comparison_summary.json`, `qualification_summary.json`, and the report preserve the
actual scientific `PASS`, `FAIL`, or `INCOMPLETE` result under the unchanged gates.
The default `QUALIFICATION` policy remains available for manual campaigns that should
return status 2 when scientific qualification is not PASS.

QUICK and STANDARD are independent: neither waits for, reads, overwrites, nor provides
convergence evidence for the other. A convergence claim requires an explicitly
registered second campaign on a common epoch population.

## Stage-specific operation

Preparation validates the entire external data contract and renders the sole RBSP-B
branch without launching AMPS:

```bash
python3 srcEarth/test/SA04/run_SA04.py --prepare-only \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926/model_inputs/SA04 \
  --output-dir test_output/SA04/prepared-standard \
  --profile STANDARD --solver GRIDDED --mpi-ranks 4 --threads 8 \
  --mover HC4 --trace-dt 0.1
```

Continue already prepared evidence, or re-evaluate saved evidence, with:

```bash
python3 srcEarth/test/SA04/run_SA04.py --run-only \
  --data-root /path/to/package \
  --output-dir test_output/SA04/prepared-standard \
  --mpi-ranks 4 --threads 8

python3 srcEarth/test/SA04/run_SA04.py --compare-only \
  --data-root /path/to/package \
  --output-dir test_output/SA04/existing-campaign

python3 srcEarth/test/SA04/run_SA04.py --plot-only \
  --data-root /path/to/package \
  --output-dir test_output/SA04/existing-campaign --formats png eps
```

`--run-only` uses the numerical choices saved by preparation. Compare-only and
plot-only do not execute or rebuild AMPS. Numerical refinement is registered
explicitly, for example:

```bash
python3 srcEarth/test/SA04/compare_SA04.py \
  --data-root /path/to/package \
  --output-dir test_output/SA04/coarse-campaign \
  --convergence-with test_output/SA04/fine-campaign
```

`--rebuild` is deliberately narrow: it validates the exact repository and a
non-symlink `build` target, checks for active AMPS/build processes, removes only that
build directory, runs the mandated `make -j`, and records the log. It is incompatible
with comparison-only and plotting-only operation.

## GRIDDED execution policy

The registered campaigns use the mesh-backed Mode3D tracer and explicitly render
`TEMPORAL_MODE SNAPSHOT_LIST`. Each selected trajectory row is therefore evaluated
against its matching immutable T05 snapshot before the workflow advances to the next
epoch. Ordinary `TIME_SERIES` would incorrectly evaluate every location at every
field epoch.

SA04 explicitly opts into
`MODE3D_TRAJECTORY_SNAPSHOT_SCHEDULER EPOCH_PARALLEL`, matching the reviewed SA01
execution strategy requested for this workflow. The application-wide default remains
`COLLECTIVE`. With epoch-parallel scheduling, each MPI process owns one complete
epoch/location transaction and a private mesh; OpenMP distributes that epoch's
energy/direction work. The choice is recorded in both `prepared_run.json` and the AMPS
command, so a future default cannot silently change ownership.

The runner uses the documented initial GRIDDED mesh controls: 0.05 Earth radii near
Earth, 2.0 Earth radii at the outer boundary, linear coarsening, and exponent 1.0.
Those values are numerical settings, not proof of mesh convergence. Any refinement
must be separately labelled and evaluated through the unchanged 0.02 common-epoch
gate.

Epoch-parallel execution allocates a complete mesh/compact field per MPI process. On
Linux, the first rank measures the domain resident-memory increase and projects the
remaining replicas per shared-memory node. The recorded default permits at most 0.80
of `MemAvailable`; an excessive projection fails collectively and recommends fewer
ranks, `COLLECTIVE`, or a separately justified mesh/resource change. This guard does
not change physics or acceptance gates. If `/proc` is unavailable, the log reports
that the preflight was unavailable and makes no safety claim.

The primary trajectory limit remains 300 seconds. Four unresolved-only extensions
with factor two permit 600, 1,200, 2,400, and 4,800 seconds only for directions ending
at `TIME_LIMIT` or `STEP_LIMIT`. Every pass restarts only those still-unresolved seeds
with the same launch state, field, HC4 mover, timestep policy, direction grid, and trap
criteria; already-resolved trajectories are not recomputed. A trajectory also stops
before its pass ceiling as soon as it reaches a physical boundary or passes the
unchanged trap detector. The direction and energy grids and observation population
therefore do not change.

The fourth pass retains shared production behavior established by a controlled SA02
RBSP-B trajectory that remained unresolved at 2,400 seconds and physically escaped at
2,584.0747 seconds. That earlier evidence justifies the general unresolved-only upper
bound; it is not an SA04 result and does not predict which SA04 directions will need an
extension. Calling a long-lived SA04 seed trapped still requires the unchanged
positive detector, and any final unresolved fraction above 0.01 remains fatal.

The drift-shell recurrence detector still requires three complete revolutions plus its
coverage, phase-space, and secular-radius checks; a timeout is never relabelled as
forbidden. The unresolved fraction remains fatal above 0.01 after the fourth pass.

The registered mover is HC4 with a 0.1-second trace step. Both values are explicit in
the command and evidence. Changing them is a labelled numerical experiment, not an
allowed way to retune the observational gates.

### Trap and angular-convergence diagnostics

`diagnose_SA04.py` can isolate exact preserved epochs while retaining the full energy
grid, strict unresolved gate, immutable field inputs, and production trajectory
implementation. `--dt` is a timestep experiment. `--directions-per-energy` is a
separate angular-quadrature experiment: it retains the source deck's explicit
`DS_DIRECTION_NMU x DS_DIRECTION_NPHI` parent grid and changes only the deterministic,
approximately uniformly spaced subset used at every energy. Once that parent is
exhausted, `--direction-grid NMU NPHI` explicitly replaces both dimensions; it is
accepted only together with `--directions-per-energy`, and the requested population
cannot exceed the new product. The runner computes
`DS_MAX_PARTICLES = (DS_NINTERVALS + 1) * directions` and rejects a non-integral
source population or a request larger than the parent grid. The diagnostic status
records the source and requested direction counts, both parent-grid dimensions and
capacities, energy-node count, and both particle caps. This prevents an edited particle
cap or an implicit grid change from being mistaken for controlled quadrature evidence.

For example, after confirming the epoch exists in a prepared QUICK trajectory, the
following command tests 384 directions per energy without changing the 24-interval
spectrum grid or any physical gate:

```bash
python3 -B srcEarth/test/SA04/diagnose_SA04.py \
  --source-branch /path/to/prepared/SA04/RBSPB \
  --output-dir test_output/SA04/diagnostics/rbspb-angular-384 \
  --epochs 2017-09-13T21:00:00Z --amps ./amps \
  --mpi-ranks 1 --threads 6 --density-backend OPENMP \
  --field-init PARALLEL --mover HC4 --dt 0.1 \
  --directions-per-energy 384
```

This diagnostic does not weaken or replace the trap detector. A bounded trajectory
that fails the complete-turn coverage, phase-space recurrence, and secular-radius
gates remains unresolved and is eligible only for the documented time extensions.
Angular refinement measures the integrated consequence of sensitive cutoff-boundary
directions; it must not relabel those directions as trapped or dilute them merely to
make the unresolved fraction pass. Diagnostic completion is not SA04 qualification,
and a registered campaign still needs common-population convergence evidence under
the unchanged 0.02 product gate.

## Evidence and figures

Each campaign root contains immutable preparation/data-validation records, one
`RBSPB/` branch, `physical_execution_summary.json`, comparison CSVs,
`comparison_summary.json`, `SA04_validation_report.md`,
`qualification_summary.json`, `plot_manifest.json`, and—only after a complete
workflow—`workflow_completion.json`.

Each branch contains the rendered input, selected trajectory and snapshot list,
`amps.log`, `run_status.json`, native AMPS products, termination evidence, and
independent reconstruction residuals. GRIDDED trajectory mode writes uniquely
suffixed artifact sets per epoch. The runner requires every artifact class to match
the prepared trajectory count and rejects duplicate or missing epochs.

Requested PNG and EPS comparison plots are written below `<campaign>/figures/`.
They include channel time series, representative spectra, model/reference scatter,
spacecraft time-energy maps, total-flux and density histories, residuals versus L/MLT,
orbit/upstream/storm context, unresolved/refinement diagnostics, and a validation
overview. A convergence figure is added only when valid refinement evidence exists.
The plot manifest hashes every figure and binds it to the plotting implementation and
saved comparison/context inputs. Blank time-energy cells mean missing data; explicit
physical-zero model predictions are marked separately.

No prepared deck, successful executable exit, or generated plot by itself asserts an
observational PASS. Read the qualification summary and report.

## Implementation regressions

The focused suite is data-independent and does not launch the full physical campaign:

```bash
python3 -B -m unittest discover -s srcEarth/test/SA04/tests -v
python3 srcEarth/test/test_runner.py -j 1 srcEarth/test/SA04/list
```

It checks numerical operators, physical-zero accounting, strict schemas and
completeness, provenance/stale-evidence rejection, all plotting families, output
reservation, standalone parsing of the exact repository command, concurrency
contracts, and the production unique-field-writer alias regression. Passing this
focused suite does not qualify SA04; the master-list QUICK workflow supplies the
physical observational validation.

Python plotting/comparison dependencies are listed in `requirements.txt`.
