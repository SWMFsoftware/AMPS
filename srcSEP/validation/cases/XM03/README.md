# XM03 — published M-FLAMPA 2013 April 11 event reproduction

XM03 compares linked srcSEP products with Liu et al. (2025), DOI
`10.3847/1538-4357/adc4e3`. It includes a structured reconstruction of the
published model input (`publication_input.json`) as well as the reference
solution. The reference contains Earth 2.0–2.5 and 20–25 MeV
time-intensity curves for lambda0=0.3 au from Figure 15(a,b), the Earth-connected
radial mean-free-path profile from Figure 14(b), and the reported Earth fluence
spectral index -1.78 from Figure 15(c). The event interval begins at 2013-04-11
06:00 UTC. Differential intensity units are `(cm2 s sr MeV)^-1`, radius is in
solar radii, mean free path is in au, and the fitted index is dimensionless.

Export the production comparison products in long form:

```text
observable,coordinate,value
```

The required observable names are exactly those in
`reference/mflampa_2013_apr11_event.csv`. Time-series coordinates are elapsed
hours, mean-free-path coordinates are solar radii, and the spectral-index
coordinate is zero. The linked callback validates finite coordinates/values,
permits a negative value only for the signed spectral index, and writes its
normalized artifact transactionally.

## Extracted publication input

The article is unusually detailed, and `publication_input.json` records the
reported setup at section/table/equation granularity. Major items include:

- the 2013-04-11 06:04 UTC GONG map, weak-field transform from Equation (1),
  PFSS source surface at 2.5 solar radii, harmonic order 180, and 5-by-5-pixel
  EEGGL smoothing;
- all four AWSoM-R free parameters from Table 1 and the approximate SC/IH grid
  extents, block sizes, resolutions, cell counts, and CME-path refinement;
- the GL pole locations, radius, stretching, height, field strength, and the
  675 km/s CME speed used by EEGGL;
- 648 M-FLAMPA lines seeded at 2.5 solar radii over 360 degrees longitude and
  +/-85 degrees latitude, 120 s coupling, the -120 km/s shock criterion, the
  Parker/Poisson-bracket formulation, and the absence of perpendicular
  diffusion;
- the upstream mean-free-path law and 0.1/0.3/1.0 au sensitivity values, the
  downstream turbulence prescription and diffusion floor, and the 10 keV,
  p^-5 injection with unit coefficient followed by a 1.2 flux scaling; and
- the Earth/STA/STB coordinates from Table 2 plus the energy channels,
  mean-free-path snapshot, three-day fluence interval, and 1–50 MeV fit range.

Reported values are not the same as a complete executable input. The exact map
bytes, PARAM files, SWMF revisions, restart, AMR state, particle grids, solver
controls, all 648 seed coordinates, observer interpolation, and author-produced
Figure 14/15 tables are unavailable in the paper. The manifest therefore sets
`reproduction_status` to `partial` and enumerates these blockers in
`missing_required_inputs`. A future author data release can fill those fields
without changing the scored reference schema.

Save the production result at the single fixed location
`model/srcsep_output.csv`. It is a model output, not an alternate input. The
registered `input.json` automatically selects `publication_input.json`, the
reference, acceptance gates, and that result path. Missing production output
yields SKIP while retaining the extracted input, reference, and figures.

```sh
python3 test/run_tests.py --amps /path/to/linked/srcSEP/AMPS \
  --validation-case XM03 \
  --output-dir test_output/XM03
```

Passing `--case-input` is an error for XM03. This keeps every archived XM03
command tied to the same paper-derived physical setup.

Provisional gates are a factor-three (0.477 dex) intensity RMSE, 3 h peak-time
error, 0.3 spectral-index error, 0.2 dex mean-free-path RMSE, and complete
series coverage. The plan's 1 h onset metric requires an agreed background and
threshold definition and is recorded in input but is not scored prematurely.
Model-to-M-FLAMPA and model-to-spacecraft results must be reported separately.
Digitization calibration and uncertainty are in `reference/provenance.json`.
Every run also copies the source-reviewed reconstruction into the output as
`XM03_publication_input.json`; completed comparisons include its SHA-256 in
`provenance.json` so the exact literature interpretation is recoverable.
