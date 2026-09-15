# EV01 — calibration ensemble and parameter identifiability

EV01 is the campaign-level calibration case from the physics-validation plan. It
uses **observed NASA CCMC/GOES event quantities as the reference**; no synthetic,
theoretical, or assumed expected fluxes are generated. The committed pilot table
is a transcription of the CCMC SHINE/ISWAT challenge tables for the >10 MeV,
10 pfu and >100 MeV, 1 pfu definitions. Threshold non-crossings are retained as
non-events for the corresponding energy threshold.

The runner invokes the linked srcSEP/AMPS executable once per predeclared event.
Each event uses the same transport parameters. A triangular source-time profile
is constructed from the observed flare onset/peak/end times solely as a **model
input**. If an archive row reports identical flare peak and end timestamps (the
2012-07-12 CCMC event does), the zero-duration shutdown is represented by the
immediately next floating-point time so the native source parser remains strictly
monotone without inventing a finite decay timescale. The campaign then fits one
global amplitude relation,
`log10(A)=a+b log10(V_CME/1000 km/s)`, on the training partition only. There is
no event-by-event amplitude normalization. The validation member is scored after
the global fit; the EV02 holdout members are never used by EV01.

The repository bundle is deliberately labelled a pilot because the implementation
plan calls for 15–30 events and non-events. The nine-event CCMC 2021 challenge is
real observational evidence and is useful for testing the full campaign
machinery, but a formal campaign should replace/extend it with the current CCMC
CLEAR event/non-event package and freeze a larger split before scientific use.

Run:

```sh
python3 test/run_tests.py --amps ../amps --validation-case EV01 \
  --output-dir test_output/EV01
```

EV01 currently requires serial linked execution. Omit `--mpi-np`.

Outputs include the copied observation table, per-event linked model products,
training-only calibration coefficients, event/channel scores, JSON/JUnit
registry evidence, and provenance. `input.json` is fixed by the registry and may
not be replaced from the command line.

### Identifiability diagnostic

The campaign runner performs 1,000 deterministic **event-level** bootstrap
resamples of the training partition.  It records the intercept/slope spread and
their correlation in `calibration.json`, and writes every replicate to
`EV01_calibration_bootstrap.csv`.  Threshold rows from the same physical event
are kept together during resampling so they are not treated as independent
events.  This bootstrap quantifies stability of the shared calibration; it does
not modify or synthesize any observational reference value.

The runner also writes `EV01_campaign_summary.png` and `.eps` from the saved
campaign-score table, showing observed-versus-predicted peak flux and the peak
time residuals.
