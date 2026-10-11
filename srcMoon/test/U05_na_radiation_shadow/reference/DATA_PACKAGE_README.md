# Combi (1997) Na radiation-pressure reference

Status: **READY for U05 publication-curve verification only**.

This package supports the local U05 comparison of the AMPS sodium
radiation-pressure kernel with Figure 7 of:

M. R. Combi, M. A. DiSanti, and U. Fink, “The Spatial Distribution of
Gaseous Atomic Sodium in the Comae of Comets: Evidence for Direct Nucleus and
Extended Plasma Sources,” *Icarus* 130 (1997), 336–354,
DOI `10.1006/icar.1997.5832`.

The accessible authoritative source is NASA NTRS report `19990024948`, which
reproduces the paper's Figure 7 on printed page 13/PDF page 25. Its URL and
SHA-256 are recorded in `combi_1997_figure7_pixels.json` and
`raw_inventory.json`.

## Scope and limitation

The original machine-readable digitization used to construct the production
array in `src/species/Na.cpp` is unavailable. This package does not claim to
recover that file or every 0.1 km/s production-table value. It independently
qualifies the absolute published curve at 14 annotation-free velocities to a
predeclared graphical uncertainty of `0.5 cm s^-2`.

The March/April annotation arrows and circular markers obscure the curve at
`-10`, `-5`, `+5`, `+10`, and `+15 km/s`; those locations are excluded before
comparison with AMPS. No point was selected or rejected based on an AMPS
residual.

## Reproduction

From the AMPS repository root, run:

```sh
python3 srcMoon/test/U05_na_radiation_shadow/reference/prepare_reference.py \
  --source-pdf /data/vtenishe/moon_validation_data/radiation_pressure/raw/combi_nasa_report.pdf
```

The script verifies the source PDF hash, linearly converts the recorded axis
and curve pixels from independent 150 dpi and 300 dpi rasterizations, rejects a
cross-raster disagreement above `0.5 cm s^-2`, and generates the CSV,
provenance, QA report, and checksums. It performs no interpolation, smoothing,
renormalization, or gap filling and never reads `Na.cpp` or AMPS output.

Run the production-kernel comparison with:

```sh
python3 srcMoon/test/U05_na_radiation_shadow/test.py
```

A U05 PASS is local verification at publication-figure resolution. Full
production mover application and lunar/Earth shadow gating remain I07.
