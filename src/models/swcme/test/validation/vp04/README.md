# VP04 — CME and shock-apex height–time kinematics

VP04 validates SWCME's `DATA_DRIVEN` apex kinematics against six independent
SOHO/LASCO CDAW manual leading-edge tracks. It separates exact implementation
checks from observational interpolation skill: the former asks whether the
public API implements the intended monotone PCHIP, while the latter asks
whether a curve fitted without alternate measurements predicts those withheld
measurements.

## What is compared and why

For every event, the first, last, and alternating height-time measurements are
provided as fit knots. Intervening points are kept blind. The public
`swcme1d::Model` returns apex radius and the derivative of the same curve as
speed. `reference_solution.py` independently implements the
Fritsch-Butland/Fritsch-Carlson monotone cubic Hermite construction and imports
no SWCME code. The following comparisons have distinct meanings:

- SWCME versus the independent oracle detects an implementation error in
  interval selection, tangent construction, units, radius, or derivative.
- SWCME versus supplied knots proves the documented exact-knot contract.
- SWCME versus withheld CDAW heights measures interpolation skill on real,
  irregularly sampled leading-edge tracks.

Speed is preserved in the comparison CSV for every query and the CDAW catalog
linear speed is included as context. No catalog speed is used to tune the
interpolant.

## Data and event selection

`download_data.py` retrieves six CDAW `.yht` digital tracks: 2013-04-11,
2013-05-22, 2014-01-06, 2014-01-07, 2017-09-06, and 2017-09-10. These are
quality-index-5 halo-event traces with 8–23 monotonically increasing height
measurements. Every byte count and SHA-256 digest is fixed in source and the
raw files are ignored by source control.

CDAW heights describe the fastest leading edge in the observer's plane of sky.
They are not a deprojected three-dimensional apex. VP04 therefore validates
the numerical height-time representation and within-coronagraph interpolation,
not the accuracy of a particular CME deprojection.

## Acceptance and expected result

The case requires all six events and at least 30 withheld points. Maximum
production/reference errors must be no more than `5e-12 R_sun` in radius and
`5e-9 km/s` in speed; fit-knot residuals must be below `5e-12 R_sun`.
The observational median and worst withheld errors must be no greater than
0.35 and 1.25 solar radii, respectively. The wider single-point ceiling allows
one image-scale/manual-pick excursion while the much tighter median criterion
still governs typical interpolation skill. A pass means the data-driven public API
is numerically correct and resolves these coronagraph tracks at the stated
withholding scale. It does not establish 1-AU arrival accuracy; VP05 addresses
that separate extrapolation problem.

## Running and outputs

From `swcme/test`:

```sh
make vp04-data
make vp04-test
make vp04-validation
make validation-case CASE=VP04 VALIDATION_ARGS="--download"
```

The direct runner is `python3 validation/vp04/run_vp04.py --download`. It writes
`vp04_comparison.csv`, `vp04_reference_solution.csv`, `vp04_result.json`, an
artifact manifest, and the six
event comparison panel as both `vp04_height_time_comparison.png` and
`vp04_height_time_comparison.eps`. `--no-plots` retains numeric evidence only.
