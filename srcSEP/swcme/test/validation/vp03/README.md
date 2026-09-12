# VP03 — Background speed, pressure, and characteristic speeds

VP03 compares the default stationary SWCME background with Helios 1/2 plasma
and magnetic moments from 0.29 to 1.01 AU. It tests the primitive quantities
that set shock Mach number: bulk speed, proton pressure, sound speed, Alfvén
speed, and the scalar perpendicular fast-mode proxy.

## What is compared and why

The SWCME background has constant 400 km/s radial speed, normalized Leblanc
density, a 120,000 K proton-only pressure closure, and a Parker magnetic field.
VP03 compares those predictions with daily Helios proton-core medians. The
observed scalar proton temperature is the trace average
`(T_parallel + 2 T_perpendicular)/3`; measured density and `|B|` then determine
pressure, sound speed, Alfvén speed, and `sqrt(cs^2+vA^2)`. These quantities
directly control whether a modeled front is sub-fast or a physical shock.

The Alfvén and fast values are diagnostic scalar magnitudes. The fast proxy is
the perpendicular upper branch rather than an angle-resolved wave measurement;
that limitation is recorded in every result.

## Data, filtering, and independence

`download_data.py` obtains the immutable Zenodo record 1009506 `corefit.gz`
archive, verifies byte count, MD5, and SHA-256, and may hard-link a separately
verified VP01/VP02 cache. Raw data remain ignored by source control.

Only status-1 records in 0.29–1.01 AU are retained. Density, speed, temperature,
and field must be finite and inside broad physical guardrails. A day needs at
least ten samples, after which its median receives one vote. Fourteen radial
bins therefore describe typical background behavior rather than sampling
cadence.

The C++ driver uses the public checked `swcme1d::Model` evaluator and production
thermodynamic closure. `reference_solution.py` independently implements the
Leblanc, Parker, ideal-gas, Alfvén, and fast-speed equations without importing
production code. Production/reference relative error must not exceed `5e-12`.

## Acceptance criteria

At least 2,500 daily medians and 12 radial bins are required. Across radial-bin
medians, symmetric model/observation factors must be no worse than 1.35 for
speed, 2.0 for pressure, 1.5 for sound speed, and 1.5 for the fast proxy. These
are validation tolerances for a deliberately stationary background, not fit
parameters, and they are frozen in `run_vp03.py`.

## Running and outputs

From `swcme/test`:

```sh
make vp03-data
make vp03-test
make vp03-validation
make validation-case CASE=VP03 VALIDATION_ARGS="--download"
```

The direct runner is `python3 validation/vp03/run_vp03.py --download`. It writes
model inputs/outputs, `vp03_comparison.csv`, `vp03_reference_solution.csv`,
`vp03_result.json`, an artifact
manifest, and `vp03_background_comparison.png` plus the matching EPS file.

## Interpretation

A pass means the default closure and fields are numerically reproduced by an
independent oracle and describe typical Helios radial-bin medians within the
published factors. It does not validate electron/alpha pressure, transient
heating, stream structure, or wave propagation at a measured oblique angle.
