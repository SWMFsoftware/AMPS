# VP02 — Parker magnetic-vector angle and polarity

VP02 compares the SWCME Parker background with independent Helios 1/2 proton
core-fit observations between 0.29 and 1.01 AU. It is an observational
validation case, not a deterministic unit test: the solar wind contains
Alfvénic fluctuations, stream interactions, current-sheet crossings, and
transients that are intentionally absent from the stationary analytical model.

## What and why

The primary observable is the sector-folded RT-plane angle
`atan2(-Bt sign(Br), |Br|)`. Folding removes the inward/outward sector
sign but preserves whether the field winds in the Parker direction and by the
expected amount. A second metric measures the fraction of accepted samples for
which `Br*Bt < 0`, the quadrant predicted by an outward Parker spiral after
accounting for sector polarity. These tests validate field direction and
polarity without pretending that a single 1-AU magnitude predicts turbulent
instantaneous amplitudes.

## Data and reproducibility

`download_data.py` retrieves the immutable `corefit.gz` object from Zenodo
record 1009506. Size, MD5, and SHA-256 must all match `data/PROVENANCE.json`.
When VP01 already has the same verified bytes, the downloader uses a hard link
where the file system permits. Raw third-party data remain ignored; derived
CSV, JSON, PNG, and EPS products are reproducible evidence.

Only status-1 records are used. Samples must have 0.29–1.01 AU radius,
250–850 km/s radial speed, `sigma_B/|B| <= 0.35`, and `|Br|/|B| >= 0.15`.
The latter two rules avoid unstable magnetic intervals and sector boundaries.
Days with fewer than ten retained records are excluded, and each retained day
has equal statistical weight.

## How it is tested

1. Stream the archive without extracting it and form daily medians.
2. Build fourteen radial bins with observed 16th/50th/84th percentiles.
3. Run `vp02_model_driver.cpp`, which constructs the public `swcme1d::Model`
   at each bin's measured median speed and latitude and evaluates `Br,Bphi`.
4. Independently evaluate the published Parker relation in
   `reference_solution.py`; it imports no production model code.
5. Compare production/reference values and compare both with Helios.

Frozen gates require at least 2,000 daily medians and ten populated bins,
production/reference agreement within `2e-11` degree, a median daily absolute
observational discrepancy no larger than 20 degrees, and at least 65% Parker
quadrant agreement. Natural variability is reported in the figure and is not
hidden by tuning the rotation rate or source radius.

## Run and outputs

From `swcme/test`:

```sh
make vp02-data
make vp02-test
make vp02-validation
# or through the global campaign runner
make validation-case CASE=VP02 VALIDATION_ARGS="--download"
```

`run_vp02.py --download` is the direct equivalent. The runner writes
`vp02_comparison.csv`, `vp02_reference_solution.csv`, exact model input/output,
`vp02_result.json`, an artifact
manifest, and `vp02_parker_comparison.png` plus
`vp02_parker_comparison.eps`. `--no-plots` supports headless orchestration while
retaining all numeric evidence.

## Interpretation and limitations

A pass means the production Parker implementation agrees with an independent
equation and its median radial directional trend and quadrant are consistent
with the selected Helios background sample. It does not validate solar-cycle
sector timing, transient fields, turbulence, or a global magnitude forecast.
