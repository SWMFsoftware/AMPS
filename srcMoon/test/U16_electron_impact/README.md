# U16 — atomic Maxwellian electron-impact coefficients

U16 verifies the shared production service in
[`srcMoon/species/ElectronImpact.cpp`](../../species/ElectronImpact.cpp), using
immutable Voronov/Verner D05 coefficients. It does not enable the generic AMPS
chemistry executor or qualify an ensemble lifetime. The complete test contract
is [CONTRACT.md](CONTRACT.md), with machine-readable acceptance in
[reference/acceptance.json](reference/acceptance.json).

Run through the existing per-ID/global registration:

```bash
/usr/bin/python3 -B srcMoon/test/U16_electron_impact/test.py
```

Or execute only this standalone kernel work package, without configuring or
rebuilding AMPS:

```bash
/usr/bin/python3 -B srcMoon/test/U16_electron_impact/voronov_verification.py \
  --output-dir test_output/srcMoon/u16-local
```

`run_u16(build_dir, output_root)` is the reusable aggregate-runner hook. It
builds the production pure kernel plus a thin adapter with C++17 (`CXX`, default
`g++`), preserving compiler commands/stdout/stderr. The generated AMPS build
argument is recorded but is not used to substitute another implementation.
The normal application adapter/link and compiled neutral-to-ion mappings are
separately checked by the parent runtime controls and linked campaign.

## Physics, units and explicit limits

The selected closure is Maxwellian electrons with explicitly configured `n_e`
and `T_e`; it is not measured plasma forcing. The neutral direct-ionization
fit is `A*(1+P*sqrt(U))/(X+U)*U^K*exp(-U)`, `U=dE/Te`. The pure API takes eV
and m^-3 and returns m^3/s and s^-1. Native `A` is cm^3/s and is converted by
exact factor1e-6. Native `Tmin` is eV and `Tmax` is keV. For the archived He,
Na, Ne and Ar neutral rows the inclusive temperature interval is 1..20000 eV.
Outside-domain, nonfinite and negative inputs are refused without changing
caller outputs. Zero density yields exact zero only after validating the
species and temperature. The optional K→eV helper uses exact SI `k_B/e`.

The evaluator loads all 406 atomic stages using the native fixed-width schema,
rejects missing/duplicate/invalid rows, verifies atomic symbols/charge stages,
and transactionally publishes only the four selected neutrals. Molecular/ionic
coefficients, ion product kinematics and charge exchange are not supplied.
Its pure parser validates structure; application/data preflight independently
binds authoritative bytes and provenance. A structurally valid altered table
is never automatically an authoritative replacement.

The formula is evaluated in logarithmic form, preserving representable small
rates. It does not inherit the original Fortran's `U>80` early-zero shortcut.
Actual floating-point underflow is a valid zero limit. The real neutral domain
has `U<=24.6`, so two explicitly labeled altered-byte arithmetic fixtures test
a representable `U=100` tail and natural deep underflow. These temporary files
are not observations or altered D05, and authoritative preflight would reject
their different hashes.

## Independent checks and artifacts

The frozen reference manifest hashes `cfit.dat`, `cfit.txt`, `cfit.f`, the full
processed 406-row table and four-row neutral selection. Every run also checks
package README/provenance/QA/processing summary and native metadata hashes.
The native rows and all 406 Fortran DATA blocks must agree exactly. The selected
processed CSV preserves the native coefficients and temperature bounds.

An independent 70-digit Decimal expression evaluates a predeclared 33-point
logarithmic grid including endpoints for each neutral. Production coefficients
and rates must retain the unchanged relative **1e-12** criterion. Density 0,
1e6 and1e7 m^-3 cover zero and linear scaling; parser/input/transaction guards
are deterministic. No stochastic samples or lunar observations are created.

Where `FC` (default `gfortran`) exists, the unchanged archived `cfit.f` is also
compiled and executed. It uses single-precision REAL and rounded
`8.617385e-5 eV/K`; comparisons use its actual effective temperature and a
separately declared, conservative arithmetic-conditioning budget
`64*2^-23*(1+U)`. This supporting original-Fortran cross-check never replaces
or relaxes the 1e-12 high-precision production gate. A missing Fortran compiler
is an explicit supporting SKIPPED record; a missing required compiler/data is
U16 SKIPPED, and build/parser/provenance corruption is ERROR.

Each run creates a unique `U16/maxwellian-kernel-*/` evidence directory with
`result.json`, `rate_vs_temperature.csv`, `data_hashes.json`,
`fortran_comparison.json`, compiler logs, per-query logs and clearly named
`control-fixtures-not-data/`. Production/header/adapter/script and both roadmap
hashes are recorded. Earlier result directories are preserved. Exit codes are
`0=PASS`, `1=FAIL`, `2=ERROR`, `77=SKIPPED`.

Local U16 PASS establishes formula/rate-service correctness only. I23 still
requires actual linked dispatch, conservative particle-event/product behavior
and its separate 2% coefficient/3% lifetime gates. Formal M3 and measured plasma
claims remain unavailable until their own predecessors and capabilities pass.
