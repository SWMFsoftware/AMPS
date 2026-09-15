# XM02 — publication-informed M-FLAMPA transport reconstruction

XM02 runs a controlled srcSEP transport calculation inside the selected linked
`AMPS` executable and compares it with Zhao et al. Figure 7. No production CSV
and no `--case-input` argument are required. The registry always selects the
single source-owned `input.json`; attempts to override it are rejected.

## What came from the publication

`publication_input.json` records the values stated by Zhao et al., including
the 2013 April 11 event, 743 km/s CME, GONG/AWSoM-R/EEGGL model chain, 648
field lines seeded at 2.5 solar radii, 120 s coupling cadence, 10 keV injection,
p^-5 suprathermal tail, injection coefficient 1.25, absence of perpendicular
diffusion, and the 0.05, 0.3, and 1.0 au far-upstream mean free paths. Each
group identifies its source section, table, equation, or figure.

The reference solution is
`reference/mflampa_2013_apr11_mfp_sensitivity.csv`. It contains the three
digitized >10 MeV intensity histories from Figure 7. The exact source-PDF hash,
pixel-to-axis calibration, curve colors, sampling interval, and estimated
0.12-dex graphical uncertainty are preserved in `reference/provenance.json`.

## Controlled reconstruction required by missing publication data

Figure 7 does not provide the evolving AWSoM-R field line, CME shock history,
sample-line coordinates, plasma density/temperature, numerical grids, or the
author's tabulated solution. An exact event rerun cannot be made from the paper
alone. XM02 therefore isolates the reported transport sensitivity with the
following fixed assumptions, all stored in `input.json`:

- already-accelerated 10.1 MeV protons are released at 2.5 solar radii and
  observed at 1 au; 10.1 MeV lies unambiguously inside the >10 MeV channel;
- a causal two-stage release is the sum of exponential 0.5 h rise and 3 h
  decay clocks, standing in for the unpublished evolving shock source;
- the controlled line has zero plasma advection and zero magnetic focusing,
  so the comparison isolates field-aligned scattering and streaming;
- the outward surface flux has density `P(mu)=2*mu` on `0<=mu<=1`, the inner
  boundary reflects particles, and 1 au is a first-passage observer;
- `D_mumu=D0*(1-mu^2)` with `D0=v/(2*lambda_parallel)`, which gives each
  requested parallel mean free path through the standard diffusion integral;
- 6,000 particles per MFP use keyed seed 30202, a 120 s mover step, a 44 h
  duration, and 2 h output bins; and
- model and reference profiles are independently divided by their peak before
  scoring. Absolute pfu is not scored because its required shock/source
  normalization is absent from the publication.

The linked callback calls the production focused-transport core for every
particle and timestep. Python supplies only orchestration, the digitized
reference, scoring, plotting, and provenance. The reference is never passed to
the native model and is never used as its output.

## Running XM02

From the `srcSEP` directory:

```sh
test/run_tests.py --amps ../amps --validation-case XM02 \
  --output-dir test_output/XM02
```

The case writes the resolved input, copied publication reconstruction, native
argument manifest, linked model CSV, native JSON/JUnit, comparison metrics,
provenance, and PNG/EPS overlays. It reports PASS or FAIL rather than SKIP once
the linked executable has been rebuilt with this implementation.

Every generated comparison image is self-attributing: its title includes
`Reference: Zhao et al. (2024), arXiv:2309.16903; extracted from Figure 7`.
This short label is stored as `reference.plot_citation` plus
`reference.figure` in `input.json`; full title, URL, PDF hash, source page,
digitization calibration, and uncertainty remain in `publication_input.json`
and `reference/provenance.json`. Thus a detached PNG/EPS still identifies its
reference, while the machine-readable artifacts retain complete provenance.

The provisional controlled-reconstruction gates are complete coverage of all
three MFP series, unit-peak log-intensity RMSE no greater than 0.6 dex, and a
maximum peak-time difference no greater than 6 h. These deliberately include
model-form uncertainty from the missing shock history; they must not be cited
as validation of the full 2013 April 11 CME simulation.

## Interpretation

A PASS means the srcSEP production streaming/scattering kernel reproduces the
published ordering and approximate temporal response to the three mean free
paths under the declared controlled reconstruction. It does not validate the
unavailable AWSoM-R background, EEGGL CME, shock acceleration, absolute pfu
normalization, or the exact unidentified M-FLAMPA sample line.
