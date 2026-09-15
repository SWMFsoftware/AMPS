# EV01–EV02 campaign-level evidence implementation

This document describes the final two cases in the srcSEP validation plan and
the software/data contracts implemented for them.

## Relationship between EV01 and EV02

EV01 and EV02 deliberately share the same model configuration, observation
pipeline, event split, and calibration form, but they answer different
scientific questions.

* **EV01 — calibration ensemble and parameter identifiability.**  Fit one shared
  calibration law using only the predeclared training events, evaluate the
  training/validation ensemble, and quantify how stable the fitted parameters
  are under event-level bootstrap resampling.
* **EV02 — locked held-out validation.**  Reconstruct the same calibration from
  training events only, then score only the sealed holdout events.  The holdout
  observations never enter the fit.  Any post-unsealing modification creates a
  new campaign version and requires a new untouched cohort.

Neither case generates a synthetic, analytical, manufactured, or assumed
reference solution.  The expected values are observed/observationally derived
GOES quantities published by NASA CCMC.

## Observational reference

The pilot reference is
`validation/cases/EV0*/reference/ccmc_2021_goes_observations.csv`.  It contains
the nine SEP events used by the SHINE/ISWAT SEP Model Validation Challenge, with
two operational channel definitions per event:

* GOES integral proton flux `>10 MeV`, event threshold 10 pfu;
* GOES integral proton flux `>100 MeV`, event threshold 1 pfu.

The authoritative source is the NASA CCMC Model Input Parameters page:

`https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/input_parameters.php`

Specifically, the values come from the tables:

* **GOES Proton Measurements for >10 MeV, 10 pfu Threshold**;
* **GOES Proton Measurements for >100 MeV, 1 pfu Threshold**.

The CCMC Data Sets page documents that the event quantities were calculated
with OpSEP / `operational_sep_quantities.py`.  Per-event pages provide
unnumbered GOES Proton Measurements figures and adjacent numeric tables.  Since
the web figures are not numbered, the implementation records their exact
section/location rather than inventing a figure number.  Numeric values are
transcribed from the tables, not digitized from the figures.

The pilot commits CCMC-derived event quantities rather than the underlying raw
NOAA/NCEI telemetry files.  The provenance record says this explicitly and does
not claim a raw NOAA dataset/file identifier that is not actually archived with
the test.

## Model inputs versus observational targets

The same CCMC source page also supplies observationally constrained model input
quantities:

* GOES X-ray flare onset, peak, and end times;
* CME time at 21.5 solar radii;
* forecaster-quality CME 3-D speed.

These quantities are inputs/covariates, not the expected SEP solution.  The
observational target used for scoring is the GOES proton event behavior.

A triangular source history is constructed from flare timing only as an input
to the linked model.  If peak and end timestamps are identical, the right-hand
shutdown is represented at the immediately next IEEE-754 floating-point time to
meet the native strictly increasing time-file contract without inventing a
finite physical decay time.  If flare timing is absent (2014-01-06), the
recorded CME 21.5-Rs time is used as an explicitly marked source-limited
fallback.

## Immutable campaign split

The pilot split is event-based:

* training: 2012-03-07, 2012-05-17, 2012-07-12, 2013-04-11,
  2014-01-07, 2017-09-04/06;
* validation: 2017-07-14;
* holdout: 2014-01-06, 2017-09-10.

Rows belonging to the same physical event are never split across partitions.
EV01 never scores the holdout.  EV02 refits from the training partition and
scores only holdout rows.

## Shared calibration

The linked model is run with one fixed transport parameter set across all
selected events.  No event-specific transport tuning is permitted.  The one
campaign calibration relation is

`log10(A) = a + b log10(V_CME / 1000 km/s)`

where `A` multiplies the uncalibrated model peak and `V_CME` is the CCMC
forecaster-quality CME speed.  The two coefficients are estimated from training
rows only.

## Identifiability

Both cases retain 1,000 deterministic event-level bootstrap replicates by
default.  Whole events, not individual threshold rows, are resampled.  The
runner records the mean and spread of `a` and `b` and their correlation.  This
is an identifiability diagnostic; it does not change any observed reference
value.

## Metrics

The campaign score table contains, per event/channel:

* observed and predicted maximum peak flux;
* observed and modeled maximum-peak time;
* threshold crossing/non-crossing classification;
* peak-time residual;
* log10 peak-flux ratio;
* calibration factor;
* observation spacecraft, dataset, master-table location, per-event URL,
  figure/section location, and processing method.

Campaign metrics include log10 peak-flux RMSE, mean absolute peak-time error,
contingency-table quantities, and True Skill Statistic.

## Plot provenance

`EV01_campaign_summary.*` and `EV02_campaign_summary.*` now print the source of
the observations directly on the figure.  The footer identifies NASA CCMC,
GOES-13 corrected integral proton fluxes, the exact two threshold-table names,
the 2021-04-29 table version, the unnumbered per-event GOES Proton Measurements
figure/table location, and OpSEP processing.  Full URLs are retained in
`reference_provenance.json` and in the generated
`EV0*_observation_reference.txt` file.

This explicit provenance is important because a figure exported from the result
directory must remain scientifically interpretable without relying on the
person who ran the test to remember where the observations came from.

## Pilot scope

The nine-event bundle is a **pilot evidence set**.  The validation specification
calls for 15–30 events and non-events for EV01 and a larger untouched historical
and ultimately prospective cohort for EV02.  The current CCMC CLEAR benchmark
is the appropriate direction for expansion.  The implementation intentionally
keeps the pilot limitation visible instead of manufacturing extra events or
negative cases.

## Run commands

Both cases:

```sh
test/run_tests.py --amps ../amps --validation-case EV01 --validation-case EV02 --output-dir test_output/EV01-EV02
```

EV01 only:

```sh
test/run_tests.py --amps ../amps --validation-case EV01 --output-dir test_output/EV01
```

EV02 only:

```sh
test/run_tests.py --amps ../amps --validation-case EV02 --output-dir test_output/EV02
```

EV01/EV02 currently use serial linked execution; omit `--mpi-np`.

Dependency-light implementation/provenance checks:

```sh
test/run_ev01_ev02_tests.sh
```
