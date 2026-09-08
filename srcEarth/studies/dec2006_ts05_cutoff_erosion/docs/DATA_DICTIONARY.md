# Generated-product data dictionary

## Universal conventions

- UTC timestamps use `YYYY-MM-DDTHH:MM:SSZ`.
- Rigidity is in GV.
- Boundary latitude is the magnitude of AACGM latitude unless a signed
  coordinate is explicitly named.
- Residual is always `model - observation`.  Positive residual therefore means
  excessive modeled shielding, or a boundary that is too poleward.
- Positive lag means that the cutoff diagnostic follows the driver.
- Unresolved access state is `2` and is never treated as forbidden.

## `paired_model_observation.csv`

| Column | Meaning |
|---|---|
| `dataset` | `PAMELA_TABLE_S1` or `NOAA_POES_METOP_SEM2` |
| `epoch_utc` | Orbit/window midpoint used for pairing |
| `rigidity_gv` | PAMELA bin center or MEPED nominal lower-threshold rigidity |
| `channel` | Empty for PAMELA; P6--P9 for MEPED |
| `hemisphere` | `ABS_MEDIAN_NS` for PAMELA or N/S for POES/MetOp |
| `mlt_hour` | Empty for PAMELA or center of the three-hour MLT sector |
| `observed_boundary_aacgm_deg` | Observation-derived T50 latitude magnitude |
| `modeled_boundary_aacgm_deg` | AMPS observation-equivalent ACCESS_T50 |
| `model_minus_observation_deg` | Signed paired residual |
| `validation_role` | Primary or diagnostic |
| `used_for_primary_metrics` | Whether the row enters the predeclared primary score |

## `morphology_boundaries.csv`

Each row is one epoch, altitude, rigidity, hemisphere, and MLT sector.  The
boundary is valid only when the resolved transmission brackets 0.5 and is at
least one degree inside the retained AACGM latitude range.

## `morphology_harmonics.csv`

`mean_latitude_deg`, `amplitude_deg`, and `phase_mlt_hour` describe the first
MLT harmonic.  `fit_rms_deg` is the RMS cell residual.  The accessible-area
fraction is the fraction of the configured 35--85 degree analysis band poleward
of the T50 boundary, averaged over valid MLT sectors.

## `cutoff_dynamics_timeseries.csv`

Adds a quiet reference, `cutoff_erosion_deg`, and centered finite-difference
boundary speed.  Negative erosion is equatorward motion.

## `lag_correlations.csv`

Contains every predeclared driver, lag, altitude, rigidity, and hemisphere.
Intervals use a moving-block bootstrap with the configured three-hour block.
The raw five-minute count is not treated as independent sample size.

## `hysteresis_pairs.csv` and `hysteresis_summary.csv`

Pairs use identical altitude, rigidity, hemisphere, and MLT.  `SYMH_ONLY`
requires the configured 10-nT match.  `STRICT` additionally requires dynamic
pressure within 20% and IMF Bz within 2 nT.  The reported contrast is recovery
minus main phase; its physical sign must be interpreted together with the event
driver and confidence interval.
