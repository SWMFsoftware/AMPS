# EV01 — calibration ensemble and parameter identifiability

## Purpose

EV01 is the campaign-level **calibration and identifiability** test defined by
the srcSEP physics-validation plan.  Its purpose is not to prove that srcSEP can
fit one selected SEP event.  Instead, it asks whether one common transport
configuration and one declared, low-dimensional calibration law can describe a
preselected ensemble of observed SEP events **without event-by-event retuning**.
It also measures whether the fitted campaign parameters are stable when the
training event set is resampled.

The current repository implementation is a nine-event pilot that exercises the
complete campaign mechanics.  The validation plan calls for a larger 15–30
event/non-event ensemble before EV01 is treated as the final statistical
campaign gate.

## Validation question

Given a fixed srcSEP transport configuration and a predeclared event split, can
one training-derived parameterization reproduce observed GOES proton event
properties across multiple events, and are the fitted calibration parameters
identifiable rather than strongly sample-dependent or degenerate?

EV01 therefore tests several things simultaneously:

1. repeated end-to-end execution of the same linked srcSEP/AMPS model over a
   multi-event observational ensemble;
2. strict separation of the training, validation, and held-out event lists;
3. use of a **single shared** calibration relation rather than independent
   normalization of every event;
4. event-level bootstrap estimates of calibration-parameter stability;
5. event classification skill, peak-flux agreement, and peak-time agreement;
6. complete traceability of each observed comparison value back to its source.

## Observational reference and exact source location

The expected values in EV01 are **real observations/observationally derived
quantities**, not synthetic, theoretical, manufactured, or assumed reference
solutions.

The committed pilot table is
`reference/ccmc_2021_goes_observations.csv`.  It is a machine-readable
transcription of the NASA CCMC **SHINE/ISWAT SEP Model Validation Challenge**
tables for GOES-13 corrected integral proton fluxes.

The exact CCMC locations are:

* master source page: `https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/input_parameters.php`;
* table **“GOES Proton Measurements for >10 MeV, 10 pfu Threshold”**;
* table **“GOES Proton Measurements for >100 MeV, 1 pfu Threshold”**;
* event-list/data-download page: `https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/sep_events.php`;
* data-set description and processing discussion: `https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/data_sets.php`;
* per-event pages: `.../campaign2020/<YYYYMMDD>.php`; the compound
  2017-09-04/06 entry is on `20170904.php`.

The CCMC event pages contain **unnumbered** GOES Proton Measurements figures.
Consequently, the code and plots do not invent a figure number.  The visual
location is the GOES Proton Measurements section containing the integral-flux
threshold panels and the adjacent threshold tables.  The numeric EV reference
values are transcribed from the tables, **not digitized from the plots**.

CCMC states that the event quantities were calculated event-by-event using
OpSEP / `operational_sep_quantities.py`; the processing code is available at
`https://github.com/ktindiana/operational-sep`.  The current EV pilot stores
those CCMC-derived event quantities rather than raw NOAA/NCEI GOES telemetry
files.  Therefore the case intentionally does **not** claim a raw NOAA file or
NCEI dataset identifier that is not present in the source bundle.

The full machine-readable source statement is retained in
`reference/provenance.json`.  Each campaign-score row also records the
spacecraft, data product, master table, per-event CCMC page, figure location,
and processing method so a detached score CSV remains auditable.

## Events and split

The immutable pilot split in `input.json` is:

* **training:** 2012-03-07, 2012-05-17, 2012-07-12, 2013-04-11,
  2014-01-07, and 2017-09-04/06;
* **validation:** 2017-07-14;
* **held out for EV02 only:** 2014-01-06 and 2017-09-10.

The split is by **physical event**, not by individual threshold/channel row.  A
>10 MeV row and a >100 MeV row from the same event are never allowed to land in
different campaign partitions.

## Observed quantities used in scoring

For each event the reference table contains two operational channel definitions:

* integral proton flux `>10 MeV`, event threshold `10 pfu`;
* integral proton flux `>100 MeV`, event threshold `1 pfu`.

The primary score inputs are threshold crossing/non-crossing, maximum observed
peak flux, and maximum observed peak time.  CCMC `N/A` entries are preserved as
real observed non-crossings; they are not replaced by artificial fluxes.

The same source table also contains flare timing and CME parameters used as
observationally constrained **model inputs**.  These inputs are not treated as
reference solutions.  In particular, the source-time history is built from the
reported flare onset/peak/end times.  When peak and end are identical (the
2012-07-12 event), the zero-duration right-hand shutdown is serialized at the
immediately next representable floating-point time only to satisfy the native
strictly increasing time contract; no finite decay time is invented.

For 2014-01-06, where the challenge table provides no flare timing, the recorded
CME 21.5-Rs epoch is used as an explicitly marked source-limited fallback.

## srcSEP configuration

Every event is run with the same transport settings from `input.json`.  The
pilot currently uses:

* injection radius: 2.5 solar radii;
* nominal solar-wind speed: 400 km/s;
* nominal shock speed: 1000 km/s;
* mean free path normalization: 0.3 AU;
* radial MFP exponent: 1;
* rigidity exponent: 1/3;
* injection momentum index: 5;
* energy grid: 10–500 MeV at the explicitly listed 11 grid energies;
* 300 particles per energy;
* 180 s time step;
* 168 h event window;
* 1 h output cadence.

These values are fixed for the case.  EV01 does not alter transport parameters
separately for individual events.

## Model-to-observation mapping

The linked executable returns a model intensity time profile on the configured
energy grid.  The campaign runner integrates model intensities above 10 and
100 MeV to form the two GOES-like integral channels.  For each event/channel it
extracts the model maximum and its time and compares them with the observed CCMC
maximum and maximum time.

A single campaign-wide amplitude relation is then fitted using **training rows
only**:

`log10(A) = a + b log10(V_CME / 1000 km/s)`.

`A` scales the uncalibrated model peak; `V_CME` is the CCMC forecaster-quality
CME speed from the standard-input table.  No per-event amplitude factor is
fitted independently.

## Parameter-identifiability test

EV01 performs 1,000 deterministic event-level bootstrap resamples by default.
Whole SEP events are resampled, keeping the >10 and >100 MeV measurements from
the same event together.  For every replicate the shared calibration intercept
and CME-speed slope are refitted.

The output reports:

* bootstrap mean and standard deviation of the intercept;
* bootstrap mean and standard deviation of the CME-speed slope;
* intercept/slope correlation;
* all individual bootstrap replicates.

Large spreads or a strong compensating correlation indicate weak parameter
identifiability even if the training residual is small.

## Metrics and interpretation

The pilot reports:

* number of observational events represented;
* log10 peak-flux RMSE;
* mean absolute peak-time error;
* hit/miss/false-alarm/correct-negative counts and True Skill Statistic.

The numerical acceptance values in `input.json` are pilot thresholds.  A good
EV01 result means the shared parameterization is reasonably stable and produces
acceptable ensemble skill.  A poor EV01 result should not automatically be
interpreted as a numerical-code defect: event-correlated residuals may instead
indicate missing source, connectivity, transport, or covariate physics.

## Comparison figures and provenance printed on the plot

`EV01_campaign_summary.png` and `.eps` contain two panels:

1. observed versus predicted peak proton flux;
2. model-minus-observed peak-time residual for each event/channel.

The figure itself now carries an observation-source footer identifying:

* NASA CCMC SHINE/ISWAT SEP Model Validation Challenge;
* GOES-13 corrected integral proton fluxes;
* the exact >10 MeV/10 pfu and >100 MeV/1 pfu table names;
* the 2021-04-29 table version;
* the per-event unnumbered GOES Proton Measurements figure/table location;
* OpSEP as the event-quantity processing method;
* `reference_provenance.json` as the location of the full URLs.

The run also writes `EV01_observation_reference.txt`, a human-readable
provenance companion that can be distributed with the figure.

## Outputs

Important products under the EV01 output directory include:

* `resolved_input.json` — exact frozen case configuration;
* `ccmc_observations.csv` — copied observational reference table;
* `reference_provenance.json` — machine-readable data provenance;
* `EV01_observation_reference.txt` — human-readable data-source citation;
* `events/<event>/...` — linked model inputs, outputs, logs, JSON/JUnit evidence;
* `EV01_campaign_scores.csv` — per-event/channel scores with source-location columns;
* `calibration.json` — training-only coefficients, bootstrap summary, split, and reference checksum;
* `EV01_calibration_bootstrap.csv` — all bootstrap replicates;
* `EV01_campaign_summary.png/.eps` — comparison figures with observation provenance printed on-figure.

## Running EV01

From the `srcSEP` directory:

```sh
test/run_tests.py --amps ../amps --validation-case EV01 --output-dir test_output/EV01
```

EV01 currently requires serial linked execution.  Do not pass `--mpi-np`.

For the lightweight registry/provenance/source-serialization checks:

```sh
test/run_ev01_ev02_tests.sh
```

## Scope limitation

This nine-event set is deliberately retained as a pilot because it is real,
auditable observational evidence already used by the CCMC challenge.  It must
not be padded with synthetic events simply to meet the target sample count.  A
formal EV01 campaign should migrate to a larger frozen CLEAR/CCMC event and
non-event cohort while retaining the same no-leakage, provenance, calibration,
and identifiability contracts.
