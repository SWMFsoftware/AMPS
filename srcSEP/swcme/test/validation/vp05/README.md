# VP05 — Multi-spacecraft arrival time and speed

VP05 tests whether the default SWCME drag-based kinematics can propagate an
observed inner-heliosphere ICME boundary to a near-radial outer spacecraft. It
uses HELIO4CAST LineupCAT v3.0, which identifies multipoint ICME observations
and supplies time, spacecraft radius/longitude, and catalog speed.

## What is compared and why

Each eligible pair initializes SWCME with only the inner spacecraft radius,
event time, and speed. The outer radius is a target coordinate; outer time and
speed remain held out. The default constant-wind DBM (`Vsw=400 km/s`,
`Gamma=1e-7 km^-1`) predicts transit duration and arrival speed. These are
compared with:

- an independent Python implementation of the sign-aware analytical DBM and a
  separately coded bracketed arrival solve; and
- the cataloged outer event time and speed.

The first comparison detects production formula, sign, unit, or root-solving
errors. The second measures useful propagation skill without fitting the
quantity being scored.

## Frozen selection and data

`download_data.py` retrieves the exact 50,510-byte LineupCAT v3.0 CSV and
requires its SHA-256 from `data/PROVENANCE.json`. Raw data are ignored by source
control. Pair selection is deterministic and outcome-blind:

1. both rows must be recognized in-situ spacecraft with numeric catalog speed;
2. the outer observation must be later and at least 0.08 AU farther out; and
3. HEEQ longitude separation must be at most 12 degrees.

Every qualifying pair is retained, including alternative boundary entries.
This avoids choosing the boundary that happens to agree best with SWCME.

## Acceptance and interpretation

At least 15 pairs are required. Production/reference transit and speed must
agree within `2e-8 s` and `2e-10 km/s`. Median absolute observational arrival
error must be at most 12 hours and median outer-speed error at most 25%. The
result also reports event count and the 84th-percentile arrival error.

A pass means a one-dimensional fixed-parameter DBM has useful median skill for
this narrowly aligned catalog sample. LineupCAT start time and speed are ICME
boundary/mean-speed proxies rather than a uniform re-fit of shock onset and
shock speed. Flank curvature, longitudinal evolution, interaction, and
spacecraft-specific boundary uncertainty remain outside this case.

## Running and products

From `swcme/test`:

```sh
make vp05-data
make vp05-test
make vp05-validation
make validation-case CASE=VP05 VALIDATION_ARGS="--download"
```

The direct runner is `python3 validation/vp05/run_vp05.py --download`. It writes
the frozen model input, public-API output, full pair comparison,
`vp05_reference_solution.csv`, result and artifact JSON, and
`vp05_multipoint_comparison.png` plus the equivalent EPS.
