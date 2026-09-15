# EV02 — locked held-out validation campaign

EV02 uses the same immutable NASA CCMC/GOES observational table and the same
predeclared split as EV01. It re-derives the two global calibration coefficients
from the training rows only and then scores **only** the holdout event IDs. The
holdout observations are therefore not read by the fitting routine. Any change
to the split, input, transport configuration, observation table, or metric
threshold changes the case input/checksum and constitutes a new campaign
version rather than an overwrite of the old score.

The committed pilot holdout contains 2014-01-06 and 2017-09-10. This is a
mechanical implementation of the sealed-gate workflow, not a claim that two
holdout events are sufficient for release science. The full EV02 campaign must
use a larger sealed historical cohort and, when operational inputs are
available, a prospective CCMC SEP Scoreboard cohort.

Run:

```sh
python3 test/run_tests.py --amps ../amps --validation-case EV02 \
  --output-dir test_output/EV02
```

EV02 currently requires serial linked execution. Omit `--mpi-np`.

The expected/reference values are observations only. Model source timing and
transport assumptions are explicitly separated from those observations and are
never promoted to reference data. Equal observed flare peak/end timestamps are
encoded as a zero-duration right-hand shutdown using the next representable
floating-point time; this is a serialization detail, not a synthetic decay model.

### Locked calibration and diagnostic outputs

EV02 recomputes the shared calibration from the EV01 training partition only,
including 1,000 deterministic event-level bootstrap resamples.  The sealed
hold-out targets never enter that fit.  Bootstrap replicates are retained in
`EV02_calibration_bootstrap.csv`; observed-versus-predicted peak flux and peak
time residuals are written as PNG and EPS campaign-summary figures.
