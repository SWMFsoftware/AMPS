# EV01–EV02 campaign-level evidence implementation

This change implements the final two cases in the srcSEP validation plan.

* **EV01** runs a sealed multi-event calibration/validation ensemble. The model
  is linked srcSEP/AMPS. The reference table contains real GOES observations
  published by NASA CCMC for the SEP Model Validation Challenge. A single
  campaign-wide CME-speed amplitude law is fit on training rows; no per-event
  scaling is allowed.
* **EV02** repeats the training-only fit and evaluates only the predeclared
  holdout events. Holdout target values do not enter fitting.

The observational table contains the >10 MeV/10 pfu and >100 MeV/1 pfu event
quantities for the nine SHINE/ISWAT challenge events. N/A threshold crossings
are kept as true observational non-crossings. The table is not a theoretical
or synthetic reference. Provenance is stored beside the CSV.

The nine-event bundle is a **pilot evidence set**. The campaign specification
calls for 15–30 events and non-events; the current CCMC CLEAR benchmark now
provides the appropriate larger long-baseline source. The implementation keeps
this limitation explicit rather than inventing additional events or synthetic
negative cases.

Run both cases:

```sh
python3 test/run_tests.py --amps ../amps \
  --validation-case EV01 --validation-case EV02 \
  --output-dir test_output/EV01-EV02
```

For the dependency-light registration/provenance checks:

```sh
test/run_ev01_ev02_tests.sh
```

## Parameter stability and figures

EV01/EV02 now include a deterministic event-level bootstrap of the shared
calibration law (1,000 replicates by default).  The implementation resamples
whole SEP events rather than individual energy-threshold rows and records the
parameter standard deviations and intercept/slope correlation.  This is an
identifiability diagnostic only: all expected values remain the NASA CCMC
observations in the committed reference table.

Each run also produces PNG and EPS campaign figures from the saved score CSV,
plus a copied reference-provenance record, so plotting and provenance can be
reproduced independently of the live model process.
