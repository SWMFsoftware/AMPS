# EV02 — locked held-out and forward-validation campaign

## Purpose

EV02 is the campaign-level **unseen-event validation gate**.  It is designed to
answer a different question from EV01: after the model configuration,
calibration form, event definitions, thresholds, and data-processing rules have
been frozen, how well does srcSEP perform on observations that were not used to
fit or select the calibration?

The current repository implementation is a historical holdout pilot.  It
implements the sealed-gate mechanics using two untouched CCMC challenge events.
A final release campaign requires a substantially larger sealed cohort and,
when operational inputs are available, a distinct prospective CCMC SEP
Scoreboard/CLEAR-style evaluation.

## Validation question

With the transport configuration and calibration prescription fixed before the
holdout is scored, what event-classification, peak-flux, and peak-time skill does
srcSEP achieve on previously unseen SEP events?

The defining rule is that **holdout target values never enter the fit**.  If the
holdout is inspected and then the model, split, thresholds, or calibration is
changed, that result is a new campaign version and requires a new untouched
cohort.

## Observational reference and exact source location

EV02 uses the same immutable real-observation table as EV01:
`reference/ccmc_2021_goes_observations.csv`.

The observed targets are NASA CCMC SHINE/ISWAT SEP Model Validation Challenge
GOES-13 corrected integral proton quantities from the official Model Input
Parameters page:

* `https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/input_parameters.php`;
* **“GOES Proton Measurements for >10 MeV, 10 pfu Threshold”**;
* **“GOES Proton Measurements for >100 MeV, 1 pfu Threshold”**.

Supporting locations are the CCMC SEP-event list and data-set pages and the
per-event pages under `.../campaign2020/<YYYYMMDD>.php`.  The per-event GOES
Proton Measurements figures are **not numbered** on the website.  The EV code
therefore identifies them by section/location rather than assigning a false
figure number.  Exact numeric comparison values come from the adjacent CCMC
tables and are not digitized from the figures.

CCMC documents that these event quantities were calculated with OpSEP /
`operational_sep_quantities.py` from GOES observations.  Full URLs, table names,
processing-code location, and the raw-data limitation are in
`reference/provenance.json` and in the generated `EV02_observation_reference.txt`.

## Frozen partition

The common split is defined in `input.json` before scoring:

* training: 2012-03-07, 2012-05-17, 2012-07-12, 2013-04-11,
  2014-01-07, 2017-09-04/06;
* validation: 2017-07-14;
* **EV02 holdout:** 2014-01-06 and 2017-09-10.

EV02 executes the training events because it must reconstruct the same shared
calibration from training data, but the final target list contains only the
holdout rows.  The validation member is not scored as part of EV02, and holdout
peak fluxes/times are not supplied to `_fit_amplitude()`.

## Calibration lock

The same shared law used by EV01 is fitted from training rows only:

`log10(A) = a + b log10(V_CME / 1000 km/s)`.

The event transport parameters remain those in `input.json`; no holdout-specific
normalization or coefficient changes are allowed.  The training bootstrap is
also repeated so the release artifact contains the calibration uncertainty that
was present before evaluating the holdout.

## Observed quantities and model mapping

For each holdout event EV02 scores the two operational GOES-like channels:

* `>10 MeV`, 10 pfu event threshold;
* `>100 MeV`, 1 pfu event threshold.

The linked model time profile is integrated above the corresponding energy
threshold.  The campaign then compares the model peak flux and peak time with
the CCMC observed maximum peak flux/time and checks whether model/observation
cross the operational event threshold.

The held-out rows remain real observations even when a channel did not cross a
threshold.  Missing/non-crossing entries are not filled with an assumed flux.

## Metrics

EV02 reports the same campaign skill quantities needed to compare the locked
prediction against the held-out observations:

* log10 peak-flux RMSE;
* mean absolute peak-time error;
* contingency-table counts;
* True Skill Statistic;
* observational-event-count/completeness evidence.

These metrics are computed only after the training calibration has been fixed.
The current thresholds are pilot release criteria rather than a claim that a
two-event holdout is statistically sufficient.

## Comparison figures and source attribution

`EV02_campaign_summary.png` and `.eps` show:

1. observed versus predicted holdout peak flux;
2. holdout model-minus-observed peak-time residuals.

Every figure now prints the observation provenance directly below the panels:
NASA CCMC challenge, GOES-13 corrected integral proton fluxes, exact table names,
table update date, per-event unnumbered GOES Proton Measurements plot/table
location, and OpSEP processing.  Full URLs are retained in
`reference_provenance.json` and `EV02_observation_reference.txt`.

`EV02_campaign_scores.csv` additionally contains source columns for every row:
spacecraft, dataset, master-table location, per-event URL, figure location, and
processing method.  This allows the CSV to remain traceable if it is separated
from the rest of the run directory.

## Outputs

Important products include:

* `resolved_input.json`;
* `ccmc_observations.csv`;
* `reference_provenance.json`;
* `EV02_observation_reference.txt`;
* per-event linked model inputs/outputs/logs;
* `EV02_campaign_scores.csv` with row-level observation provenance;
* `calibration.json` and `EV02_calibration_bootstrap.csv`;
* `EV02_campaign_summary.png/.eps` with provenance printed on the figure.

## Running EV02

From `srcSEP`:

```sh
test/run_tests.py --amps ../amps --validation-case EV02 --output-dir test_output/EV02
```

EV02 currently requires serial linked execution.  Do not pass `--mpi-np`.

To run EV01 and EV02 in sequence:

```sh
test/run_tests.py --amps ../amps --validation-case EV01 --validation-case EV02 --output-dir test_output/EV01-EV02
```

For dependency-light checks:

```sh
test/run_ev01_ev02_tests.sh
```

## Interpretation and campaign governance

A held-out degradation relative to EV01 is direct evidence of overfitting or
poor transfer across events.  A systematic error with CME speed, connection, or
energy can identify missing covariates/physics.  A correction made after seeing
the holdout does not invalidate the science effort, but it **does invalidate the
claim that the same holdout remains unseen**.  Such a correction must create a
new version and a new untouched evaluation cohort.

The long-term extension of EV02 is prospective scoring: a release-tagged model
is run using only information available before the prediction cutoff and is
compared later with independently processed observations.  Retrospective and
prospective scoreboards must remain distinct.
