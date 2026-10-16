# MMS01 — MMS1/HPCA local H+ product validation

MMS01 is a calculation-kernel validation built from the frozen MMS1/HPCA
magnetosheath interval on 2015-11-10. It follows SA01's operational pattern: all
test-owned files live in a separate directory, QUICK and COMPLETE are independently
executable, every campaign gets a unique evidence directory, comparison can be
rerun without repeating calculation, and both profiles create comparison figures.

MMS01 does **not** import or modify SA01. It also does not route the data through the
AMPS trajectory solver. The current trajectory product accepts a boundary spectrum;
the downloaded material contains a local MMS distribution but no documented
boundary-to-MMS distribution mapping. Treating the local measurement as a boundary
condition would make a physically unsupported transport claim.

## What the test validates

For H+ at each archived energy centroid, MMS01 performs two calculations:

1. it converts HPCA phase-space density `f` in `cm^-6 s^3` to directional
   differential intensity `J` in `cm^-2 s^-1 sr^-1 eV^-1` using HPCA CMAD revision
   4 equation 5-13, `J = 2 E f / m^2`, and the CMAD H+ mass coefficient
   `m = 1.0453453e-12 eV s^2 cm^-2`;
2. it integrates the directional intensities over the calibrated HPCA angular
   response at each of the 63 stored energy centroids. The weights implement CMAD
   equation 5-15 using `delta_alpha=0.318 rad`, `delta_beta=0.39 rad`, the measured
   elevation centroids, and the 16 Table 5-8 anode correction factors.

The native C++ helper calls reusable production kernels in
`srcEarth/util/FluxNumerics.h`. An independent Python path integrates the archived
Level-2 H+ differential-flux values. The comparator requires complete row identity,
exact agreement of physical-zero masks, and maximum relative error no larger than
`1e-6`. That tolerance is fixed from the published eight-significant-digit HPCA mass
coefficient: equation 5-13 squares the coefficient, so its last-place representation
appears near `1e-7`. It is not tuned to an AMPS/observation residual.

This is an instrument interpretation and AMPS numerical-kernel test. The archived
flux and phase-space density share HPCA calibration and processing, so PASS is not
independent absolute-calibration validation.

## What the test intentionally does not validate

MMS01 does not claim any of the following:

- transport from the solar wind, bow shock, magnetosheath boundary, or another
  spacecraft to MMS1;
- a magnetic-field or trajectory model;
- number density or energy-integrated flux;
- a spectrum outside the 63 measured energy centroids;
- directional agreement in GSE/GSM or another external frame; or
- instrument calibration independent of HPCA Level 2.

Density and energy-integrated flux are excluded because the CDF stores energy
centroids, while the CMAD states that the required flight-model velocity-bin sizes
are supplied in separate spreadsheets that are not present in the validation bundle.
MMS01 never creates bin edges from adjacent centroids or substitutes the nominal
20-percent energy resolution. Scalar angular integration does not require attitude;
the 30-second MEC ephemeris is used only to record the location of every comparison
distribution using bracketed component-wise linear interpolation.

## Frozen data contract

The accepted bundle root is
`/data/vtenishe/GEO/observational_validation_datasets_20261013`. The runner also
accepts its `mms1_hpca_magnetosheath_20151110` child directly. Before reading values,
it checks the fixed sizes and SHA-256 hashes of:

- the MMS1 HPCA survey Level-2 ion v4.1.0 CDF;
- the MMS1 MEC survey Level-2 ephemeris v2.0.2 CDF; and
- the local HPCA calibration and measurement algorithm document.

The half-open event interval is
`2015-11-10T13:57:54Z <= t < 2015-11-10T14:04:24Z`. It must contain 624 contiguous
625-ms ion sweeps, or 39 full 16-sweep distributions. Every selected sweep must have
H+ quality 2, science mode 2, decimation factors 2/2/4, sweep table 3, and each
distribution must contain every `start_azimuth` value exactly once. Flux and PSD
must be finite, nonnegative, and have identical zero masks.

## Profiles and commands

QUICK evaluates distribution indices 0, 19, and 38: the first, middle, and last
qualified distributions, selected before calculating any residual. COMPLETE
evaluates all 39. Formulas, source gates, comparison tolerance, energy support, and
plot families are identical.

From the repository root:

```bash
python3 srcEarth/test/MMS01/run_MMS01.py \
  --data-root /data/vtenishe/GEO/observational_validation_datasets_20261013 \
  --amps-root /data/vtenishe/GEO/AMPS \
  --output-parent test_output/MMS01/quick_test_runs \
  --profile QUICK --formats png eps

python3 srcEarth/test/MMS01/run_MMS01.py \
  --data-root /data/vtenishe/GEO/observational_validation_datasets_20261013 \
  --amps-root /data/vtenishe/GEO/AMPS \
  --output-parent test_output/MMS01/complete_study_runs \
  --profile COMPLETE --formats png eps
```

NumPy and Matplotlib must be available in the interpreter used for the campaign.
The runner checks for exactly `cdflib==1.3.2`; when it is absent, it installs that
pure-Python reader into `test_output/MMS01/python_deps` and immediately continues.
The bootstrap uses `requirements-cdf-reader.txt` with pip hash enforcement, so a
same-version wheel replacement or source-distribution fallback is rejected. It never
modifies the system Python or user site-packages. To prepare all declared dependencies
manually instead:

```bash
python3 -m pip install -r srcEarth/test/MMS01/requirements.txt
```

For a network-isolated validation host, populate the local reader once and disable
automatic installation during campaigns:

```bash
python3 -m pip install --no-deps --target test_output/MMS01/python_deps \
  --require-hashes -r srcEarth/test/MMS01/requirements-cdf-reader.txt
MMS01_AUTO_INSTALL_DEPS=0 python3 srcEarth/test/MMS01/run_MMS01.py --profile QUICK ...
```

The exact reader version is recorded in `prepared_manifest.json`. A different
installed version is rejected rather than silently changing CDF interpretation.

The native helper additionally requires a C++11 compiler. It is rebuilt into each
campaign and records both the adapter-source and production-header hashes.

Individual stages are independently rerunnable against one evidence root:

```bash
python3 srcEarth/test/MMS01/run_MMS01.py --stage COMPARE \
  --output-root test_output/MMS01/quick_test_runs/<campaign> --profile QUICK
python3 srcEarth/test/MMS01/run_MMS01.py --stage PLOT \
  --output-root test_output/MMS01/quick_test_runs/<campaign> \
  --profile QUICK --formats png eps
```

`PREPARE` requires an empty output root. `CALCULATE`, `COMPARE`, and `PLOT` consume
saved evidence and fail if required provenance or hashes differ.

## Evidence and figures

Every campaign contains:

- `prepared_manifest.json`, authenticated native input, Level-2 reference tables,
  and interpolated MMS1 positions;
- `native_run_manifest.json`, compile/run logs, native directional products, and
  native angular spectra;
- `comparison/` tables plus `comparison_summary.json` with every fixed gate;
- `figures/` and `figure_manifest.json`; and
- `workflow_summary.json`.

The PNG/EPS figure pairs are:

- `mms01_representative_spectra` — the angular integral of the archived HPCA
  Level-2 flux and the production-kernel reconstruction from archived HPCA PSD for
  the predeclared first/middle/last distributions;
- `mms01_time_energy_comparison` — the same Level-2 reference, PSD reconstruction,
  and reconstruction/reference numerical residual at every selected time and stored
  energy centroid;
- `mms01_directional_scatter` — individual positive HPCA angular-energy cells versus
  their production-kernel reconstructions, with the exact zero-cell count stated
  separately; and
- `mms01_mms1_locations` — context-only GSE locations from the MEC ephemeris. It has
  no calculated particle value and is not an agreement plot.

The word **reconstruction** is essential: the Level-2 flux and PSD are two products
of the same HPCA measurement and are connected by the documented equation tested
here. Consequently close overlap demonstrates correct CDF interpretation, units,
indexing, zero handling, PSD conversion, and angular summation. It does **not** show
that an AMPS transport simulation reproduced an independent observation. No MMS01
figure or CSV column labels this reconstruction as a model or simulation result.

Plots are evidence views only. They read saved comparison tables and never replace
values, discard low-intensity cells, or choose representative records by agreement.
