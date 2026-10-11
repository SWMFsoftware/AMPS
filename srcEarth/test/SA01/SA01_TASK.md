# SA01 implementation and execution specification

Read `AGENTS.md` and `SA01_PROGRESS.md` first. Implement and conduct SA01 observational
validation for the live AMPS `srcEarth` application. Build Python runners, compare
the physical model output with packaged observations, and generate model/reference
figures. A dry run or synthetic fixture is not completion of the physical study.

## 1. Paths and repository inspection

```text
AMPS: /data/vtenishe/GEO/AMPS
Application: /data/vtenishe/GEO/AMPS/srcEarth
Test: /data/vtenishe/GEO/AMPS/srcEarth/test/SA01
Data package: /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926
```

Read applicable repository instructions and inspect git status before changes.
Inspect existing build configuration, executable, live parser, available outputs,
and prior runners. Inspect `srcEarth/studies/rbsp_sep_validation` and
`srcEarth/test/URBSPSEPValidation` if present. Their absence must be recorded, not
silently assumed to be an implementation defect. The packaged source is a reference;
the live repository controls executable behavior.

If compilation is necessary, use the mandatory remove-`build`, then `make -j`
sequence in `AGENTS.md`. Log the build and propagate failure. Check that an existing
executable corresponds to the source/configuration before reusing it.

## 2. Data and science setup

The half-open event window is 2014-01-06 09:00 UTC through 2014-01-10 00:00 UTC.
RBSP-A/B REPT are the validation targets. SOHO/ERNE is the external proton-spectrum
input. Qin-Denton provides the five-minute T05 driver, including W1–W6. OMNI is
independent context, not a replacement for missing storm-history parameters.

Verify these files relative to the data-package root:

```text
validation_data/raw/SA01/source_manifest.json
validation_data/processed/SA01/provenance.json
validation_data/processed/SA01/ts05_driver_5min.txt
validation_data/processed/SA01/boundary_soho_erne_5min.dat
validation_data/processed/SA01/soho_erne_5min.csv
validation_data/processed/SA01/qin_denton_5min.csv
validation_data/processed/SA01/omni_hro2_5min.csv
validation_data/processed/SA01/rbspa_trajectory_geo_5min.txt
validation_data/processed/SA01/rbspb_trajectory_geo_5min.txt
validation_data/processed/SA01/rbspa_rept_energy_bins.csv
validation_data/processed/SA01/rbspb_rept_energy_bins.csv
validation_data/processed/SA01/rbspa_rept_spectrum_5min.csv
validation_data/processed/SA01/rbspb_rept_spectrum_5min.csv
validation_data/processed/SA01/rbspa_rept_density_5min.csv
validation_data/processed/SA01/rbspb_rept_density_5min.csv
```

Inspect hashes, schema, units, epochs, energy definitions, quality selections,
missing-data masks, and processing provenance. Inputs are already preprocessed;
do not redownload/regenerate them unless a specific verified problem requires it.
Inspect `model_inputs/SA01/RBSPA` and `model_inputs/SA01/RBSPB` as baseline examples.

Use the packaged domain L >= 3.2 and the documented `Data_Quality == 0` convention.
Inspect whether that L condition is implemented per sample or on a bin mean; preserve
and describe the actual frozen reference definition. Do not change it silently.
Do not claim quality filtering proves isotropy or removes every trapped proton.

## 3. Runner layout and CLI

Create these entry points, or a documented equivalent consistent with repository conventions:

```text
run_SA01.py       unified preparation/build/run/compare/plot orchestration
prepare_SA01.py   validation, trajectory selection and input generation
compare_SA01.py   independent product comparison and metrics
plot_SA01.py      plots from saved model/reference comparison artifacts
sa01_common.py   shared scientific and file-handling utilities
README.md        science, data, usage, output schema, gates and limitations
requirements.txt explicit dependencies
tests/           meaningful numerical and failure-handling tests
```

Keep the standing instructions, task specification and progress journal in this
directory too. Avoid duplicating production numerical kernels when existing modules
are correct; make their imports/dependencies explicit and testable.

The unified CLI must support:

```text
--data-root PATH                 extracted package root, or validation_data itself
--data-archive PATH              optional alias; document that this is a directory
--amps-root PATH
--amps PATH
--output-dir PATH
--profile QUICK|STANDARD|PUBLICATION
--mpi-ranks N
--threads N
--spacecraft RBSPA RBSPB
--rebuild
--prepare-only
--compare-only
--plot-only
--formats png eps
```

Make mutually exclusive modes and invalid combinations fail clearly. Validate
positive resource counts and all paths. Use absolute paths in generated decks and
subprocess commands. Run each spacecraft in its own working/output directory.
Use a new run directory or an explicit narrowly scoped overwrite option; never
delete the source/data trees or silently replace earlier physical results.

Run spacecraft sequentially by default. Print useful stage/branch progress and
stream physical-run output while saving full logs. Standalone comparison and plotting
must work without launching or rebuilding AMPS.

## 4. Input rendering

Preserve T05, the measured trajectory, time-dependent ERNE spectrum, proton species,
isotropic/uniform baseline boundary, zero electric field, combined cutoff and
density-spectrum outputs, strict temporal handling and unresolved accounting.
Verify every keyword against the live parser before execution.

Use these channel intervals after checking the reference-bin metadata:

| Channel | Minimum MeV | Maximum MeV |
| --- | ---: | ---: |
| REPT_P1 | 18.5 | 24.0 |
| REPT_P2 | 24.0 | 31.2 |
| REPT_P3 | 31.2 | 40.6 |
| REPT_P4 | 40.6 | 52.8 |
| REPT_P5 | 52.8 | 68.6 |
| REPT_P6 | 68.6 | 89.2 |
| REPT_P7 | 89.2 | 116.0 |

The 208-MeV integral channel is not a differential-spectrum bin. The study makes
no qualification claim outside the supported finite energy interval.

Prepared trajectories must follow the packaged boundary-availability policy,
retaining exact complete five-minute boundary rows. Preserve fail-on-gap/out-of-range
behavior, including the baseline maximum interpolation span of 900 seconds.
Record all excluded epochs with reasons. Verify energy interpolation/edge handling:
the highest ERNE center and the nominal instrument upper band edge are different
concepts. Confirm a valid documented boundary representation through 116 MeV;
do not silently extrapolate or invent an extra spectral node.

Use baseline workload profiles with unchanged scientific acceptance gates:

| Profile | Trajectory cadence min | Spectrum intervals | Cutoff nodes | Angular samples |
| --- | ---: | ---: | ---: | ---: |
| QUICK | 60 | 24 | 24 | 96 |
| STANDARD | 15 | 48 | 50 | 384 |
| PUBLICATION | 5 | 96 | 120 | 1536 |

Verify each option's exact meaning against current parser/implementation semantics.
Record requested/effective MPI and thread settings. Expected launch pattern, subject
to verification of current supported options:

```bash
mpiexec -n N /absolute/path/to/amps \
  -mode gridless -i input/AMPS_PARAM.in \
  -gridless-parallel THREADS -gridless-threads T \
  -gridless-mpi-scheduler DYNAMIC
```

## 5. Model/reference comparison

Inspect actual output schema. Expected physical products are:

```text
gridless_points_spectrum.dat
gridless_points_flux.dat
gridless_points_density.dat
gridless_termination_summary.dat
standalone_run_manifest.json
```

Require complete products, unique times, coherent epochs, declared units and finite
valid numerical values. A process exit of zero is execution success, not observational
PASS. Do not collapse unresolved trajectories into magnetic blocking.

### Spectrum operator

Compute the model bin average over each REPT interval:

`J_model,i = integral(J_model(E), E_min,i, E_max,i) / DeltaE_i`.

Use documented piecewise log-linear integration where appropriate. Compare the
bin average, not merely a center-energy node. Convert m^-2 intensity to cm^-2
with `1e-4`. Label rectangular energy-bin folding as an approximation, not a full
calibrated angular/energy instrument simulator.

Physical zero model intensity is meaningful. Do not remove zero/nonfinite nodes
and interpolate across the remaining support. Nonfinite output is an error; zero
against a positive observation must count as disagreement. Handle transitions and
nonpositive segments explicitly instead of applying logarithmic interpolation there.

### Flux and reconstruction

The isotropic-equivalent observed channel integral is `4*pi*J_i*DeltaE_i`.
Compare native channel outputs and the total integral over precisely 18.5–116 MeV.
Distinguish omnidirectional integral intensity from one-sided planar incident flux.
Reconstruct model integrals independently from saved spectral nodes and compare
with native AMPS products. Keep internal reconstruction residuals separate from
observational errors and identify distinct quadrature conventions.

### Density

The reference is `n_eq = 4*pi*sum(J_i*DeltaE_i/v_rel(E_i))`, in cm^-3 when speed
is cm/s and intensity is in the native centimeter units. Use relativistic proton
speed. Match the model and reference energy interval. Quantify the center-speed
bin approximation against a continuous reference using a controlled spectrum.

Label this product isotropic-equivalent density. Spectrum-derived flux and density
are correlated products of one observation set, not additional independent measurements.

### Exposure, uncertainties and selection

Observations are five-minute means; the workload cadence selects which means are
compared. Determine whether model output is instantaneous or exposure averaged and
assess the difference near rapid changes. Do not equate a 15-minute trajectory
sampling interval with a 15-minute observation accumulation window.

Retain available standard errors and angular diagnostics. Do not invent calibration
uncertainty or present a sector coefficient of variation as conclusive isotropy proof.
Discuss possible trapped/background contamination. Do not introduce target-informed
background subtraction, scaling, or exclusions into the baseline.

### Metrics and gates

Apply the fixed gates in `AGENTS.md` to all four products and both spacecraft. Report
signed log bias, median absolute log error, RMS log error, factor-two fraction and
sample counts by spacecraft/channel/product, plus pooled scores. Preserve positive-
reference/zero-model disagreements in the denominator. Handle missing reference
channels explicitly and keep product-specific coverage visible.

Where uncertainty intervals are calculated, use documented block resampling by
orbit/pass or time window; do not assume channel/time independence. Check sufficient
coverage for each product, not merely for the spectrum series. No favorable metric
can substitute for missing comparison evidence.

## 6. Visualizations

Build an independent plotting entry point using saved comparison products. Produce
PNG and EPS with readable English labels, units, legends and captions:

1. Model/reference time series for all seven channels.
2. Representative energy spectra selected by a deterministic documented rule.
3. Total 18.5–116 MeV integral-intensity time series.
4. Density time series with the reference interpretation clearly stated.
5. Model/reference scatter with unity and factor-two guides.
6. Observed/model/residual time-energy maps.
7. Residuals versus L and MLT.
8. Orbit/upstream/storm-driver context.
9. Unresolved-support and refinement/convergence diagnostics.

Distinguish spacecraft, missing data, physical zeros and numerical failure. Break
curves across genuine gaps. Do not select only favorable spectra. Provide individual
publication-suitable plots and an overview; do not invent AMPS comparison plots before
physical results exist. Avoid enormous EPS rasterizations without checking file size.

## 7. Implementation verification and physical runs

Add meaningful numerical/failure-handling checks for unit conversion, relativistic
density, analytic channel integration, physical zeros, malformed/duplicate times,
boundary availability, gate failure, pooled metrics, absolute paths, build safety,
and plot generation from a small known fixture. Run relevant existing pipeline tests.
Keep fixture PASS explicitly distinct from observational PASS.

Execute:

1. Data/provenance validation.
2. QUICK for execution diagnosis.
3. STANDARD for both spacecraft.
4. All product comparisons and figures.
5. Controlled numerical refinement on the same matched epochs.

For convergence, hold the observation population fixed while refining numerical
parameters. Compare at least two recorded resolutions; use three where estimating
order or checking nonmonotone behavior requires it. Target <=0.02 fractional change
for integrated products with a documented absolute criterion near zero. Changing
cadence alone is not an isolated numerical refinement. One PUBLICATION run is not
convergence evidence.

If a runner/parser interface fails, diagnose and correct it. A demonstrated production
defect may require a focused fix consistent with repository instructions; document
its evidence and rerun relevant checks. Scientific model/reference disagreement must
be reported, not hidden through altered observations or gates.

## 8. Evidence, documentation and reporting

Generated artifacts belong under the CLI-selected output directory, with separate
RBSPA/RBSPB branches. Preserve rendered decks, trajectories, selection reports,
raw outputs, logs, comparison CSV/JSON, plots, and a final Markdown validation report.

Capture repository revision/dirty state, executable hash, build configuration and
logs, command/working directory, dependency versions, MPI/thread layout, data/input
hashes, timestamps, profile parameters, refinement history, exclusions, counts and
actual metrics. Report process state and resume instructions after interruption.

Every meaningful stage must update `SA01_PROGRESS.md` with real evidence and the
next action. Keep the detailed logs and large data in output directories; the journal
should point to them, not embed full logs. Maintain a short current-state section and
an append-only session history. Do not commit/push unless the user requests it.

Status meanings:

- PASS: an actual comparison satisfies its frozen requirements.
- FAIL: an actual comparison violates its requirements.
- ERROR: execution, schema, units or inputs fail.
- SKIP/BLOCKED: a needed capability, dependency or asset is unavailable.

Do not silently convert SKIP/BLOCKED to coverage or claim current-revision results
from stale cached outputs.

## 9. Required usage

```bash
cd /data/vtenishe/GEO/AMPS
python3 srcEarth/test/SA01/run_SA01.py \
  --data-root /data/vtenishe/GEO/sep-in-geospace-sa01-sa04-standalone-20260926 \
  --amps-root /data/vtenishe/GEO/AMPS \
  --amps /data/vtenishe/GEO/AMPS/amps \
  --output-dir /data/vtenishe/GEO/AMPS/test_output/SA01/standard \
  --profile STANDARD --mpi-ranks 4 --threads 16 --formats png eps
```

Provide analogous prepare-only, compare-only, plot-only, QUICK, refinement, and
`--rebuild` examples in `README.md`. Examples must match implemented CLI behavior.

## 10. Completion checklist

- [ ] Applicable instructions and existing work inspected.
- [ ] Data validated without unnecessary redownload/reprocessing.
- [ ] Runners implemented and extensively commented.
- [ ] Absolute-path rendering and live parser compatibility verified.
- [ ] Required build procedure implemented and used if rebuilding was necessary.
- [ ] Relevant numerical/failure-handling tests executed with real results recorded.
- [ ] QUICK physical execution completed or its actual failure documented.
- [ ] STANDARD RBSPA/RBSPB runs completed or actual blockers recorded.
- [ ] All four observational product comparisons evaluated.
- [ ] Independent reconstruction checks reported.
- [ ] PNG/EPS model/reference figures generated and inspected.
- [ ] Numerical refinement evaluated on a common population.
- [ ] README/application documentation updated.
- [ ] Reproducible final report and machine-readable evidence provided.
- [ ] Progress journal updated with exact next action and process state.

Successful completion requires evidence for the physical comparisons and convergence;
a fully implemented runner alone is not successful scientific qualification.
