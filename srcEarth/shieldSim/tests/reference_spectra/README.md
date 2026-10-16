# Published spectra and observation fixtures

These files are immutable, offline test fixtures. They test provenance,
schema, units, uncertainty retention, and the dedicated Stage-2 adapters. They
are not accepted `shieldSim --spectrum` inputs: the legacy three-column parser
cannot represent bin edges, energy-per-nucleon, uncertainties, direction,
epoch, or general species.

Run the dependency-free integrity tests with:

```bash
python3 tests/test_reference_spectra.py
```

Every file is pinned in `manifest.sha256`. A changed upstream product must be
reviewed and committed as a new version; it must not silently replace a
fixture.

## GCR: AMS-02 Bartels-rotation average

Files:

- `gcr/ams02_20131118_hapi_info.json`
- `gcr/ams02_20131118_hapi_data.csv`
- `gcr/ams02_20131118_interpretation.v1.json`
- `gcr/raw/SM-p-He-vs-BR-REV.pdf`

These are an unmodified NASA CDAWeb HAPI metadata response and one data row for
the interval beginning 2013-11-18. The requested fields are `BR`, `ProFlux`,
`ProFlux_TotErr`, `HeFlux`, and `HeFlux_TotErr`. The metadata gives 45 proton
and 40 helium values in `1/(m^2 sr s GeV/n)`. The source product is calibrated,
27-day AMS-02 data and is described by the AMS Collaboration publication.

- Dataset: [NASA SPASE record](https://spase-metadata.org/NASA/NumericalData/ISS/AMS-02/P27D)
- HAPI endpoint: `https://cdaweb.gsfc.nasa.gov/hapi`
- HAPI dataset ID: `ISS_27DAY-AVERAGES_AMS-02`
- Publication: [Observation of Fine Time Structures in the Cosmic Proton and
  Helium Fluxes](https://doi.org/10.1103/PhysRevLett.121.051101)
- Rights: CC0-1.0 according to the SPASE record; acknowledge V. Bindi, the
  principal investigators/data providers, and CDAWeb.

The raw fixtures deliberately preserve the HAPI representation. The metadata
calls the coordinate arrays `centers` while naming them `BinProLow` and
`BinHeLow`. The interpretation record resolves that ambiguity against Table SM
XXXV in the pinned primary supplement: they are lower edges transformed from
the published rigidity bins. `AMS02HapiAdapter` applies the relativistic
rigidity-to-kinetic-energy-per-nucleon transform and converts
`1/(m2 sr s GeV/n)` to `1/(cm2 sr s (MeV/n))`. Analytic tests independently
reproduce the first proton and helium Jacobians.

This single near-Earth spectrum remains ingestion-only, not a complete deep-
space GCR environment: it includes only protons and inclusive mixed-isotope
helium, one Bartels rotation, a limited energy range, and no mission boundary
definition. Inclusive helium is intentionally not assigned a Geant4 alpha
particle until an isotope treatment is approved.

The exact requests were:

```text
https://cdaweb.gsfc.nasa.gov/hapi/info?id=ISS_27DAY-AVERAGES_AMS-02
https://cdaweb.gsfc.nasa.gov/hapi/data?id=ISS_27DAY-AVERAGES_AMS-02&time.min=2013-11-18T00:00:00Z&time.max=2013-12-15T00:00:00Z&parameters=BR,ProFlux,ProFlux_TotErr,HeFlux,HeFlux_TotErr&format=csv
```

## SEP: GOES-16 SGPS September 10-11, 2017 GLE

Files:

- `sep/raw/se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc`
- `sep/goes16_sgps_20170911T0730Z.csv`
- `sep/goes16_sgps_gle72_timeseries.csv`
- `sep/goes16_gle72_event_definition.v1.json`

The single-record fixture is a loss-minimized extraction of the NOAA GOES-16 SGPS
Level-2 five-minute record starting 2017-09-11 07:30 UTC, the time used for the
September 2017 event spectrum comparison in the SGPS publication. The event
fixture contains 1,181 contiguous sample times from 2017-09-10 16:25 UTC
through 2017-09-14 18:50 UTC, while retaining the single record as an immutable
spot check. The data retain both sensors, all 13 differential proton channels, lower and upper
energy bounds, total flux uncertainty, valid sample counts, all three
quality-flag sums, and the yaw-flip flag. The original product units were
`protons/(cm^2 sr keV s)`; the committed file records the exact keV-to-MeV
conversion.

- Product directory: [NOAA/NCEI SGPS September 2017 event data](https://www.ngdc.noaa.gov/stp/space-weather/satellite-data/satellite-systems/goesr/solar_proton_events/sgps_sep2017_event_data/)
- Instrument/data paper: [GOES-R Series Solar and Galactic Proton Sensor
  calibration](https://doi.org/10.1029/2021SW002750)
- Source netCDF SHA-256:
  `f5da29fe3a7dc3092a311a191ee015d6cb97773c6425e64214db2d0dab5fafb2`
- Source netCDF license: its global metadata permits use and redistribution
  without restriction, with NOAA/NCEI attribution.

Regenerate the CSV after downloading the source netCDF with:

```bash
python3 scripts/extract_goes16_sgps_fixture.py \
  tests/reference_spectra/sep/raw/se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc \
  tests/reference_spectra/sep/goes16_sgps_20170911T0730Z.csv \
  --time=2017-09-11T07:30:00Z

python3 scripts/extract_goes16_sgps_fixture.py \
  tests/reference_spectra/sep/raw/se_sgps-l2-avg5m_g16_s20172440000000_e20172732355000_v2_0_0.nc \
  tests/reference_spectra/sep/goes16_sgps_gle72_timeseries.csv \
  --start=2017-09-10T16:25:00Z --end=2017-09-14T18:50:00Z
```

Regeneration alone requires Python `netCDF4`; tests of the committed fixture do
not. The two sensor directions must not be averaged or labeled isotropic
without a documented angular model. For the Stage-2 event fluence they remain
separate; missing or rejected samples are not interpolated, only valid
five-minute sample intervals are integrated, coverage is reported for every
sensor and channel, and incomplete fluences are labeled incomplete. The event
definition records each event onset as an inclusive UTC start and each event
end as an exclusive UTC stop, together with the source and threshold/event
definition. The primary `>100 MeV` window is 2017-09-10 16:25 through
2017-09-12 22:40 UTC; the sensitivity `>10 MeV` window is 2017-09-10 16:45
through 2017-09-14 18:50 UTC. The NMDB 16:08 UTC GLE72 onset is retained as an
independent marker, not substituted for either integration boundary.

NOAA defines the raw SGPS timestamp as the averaging-interval start. The CSV
retains it and derives `midpoint = start + 150 s` and `end = start + 300 s`.
The canonical half-open support is therefore both `[start,start+300 s)` and
`[midpoint-150 s,midpoint+150 s)`. A boundary sample is weighted only by its
exact duration overlap. Treating the source coordinate itself as the midpoint
would shift every record and is rejected by the adapter. The resulting output
is a directional geosynchronous partial event fluence with explicit coverage,
not an unmodified interplanetary or isotropic SEP spectrum.

## Observation: CRaTER GCR LET spectra

File: `observations/crater_gcr_let_spectra_v1_3.csv`

This is the unmodified CSV payload of Mark D. Looper's CRaTER GCR LET spectra,
version 1.3. It contains 34 quiet-time, approximately half-year spectra from
2009 through 2026. The published record describes measurements through the
CRaTER telescope behind 9.9 g/cm2 of mostly tissue-equivalent plastic.

- Version DOI: [10.5281/zenodo.22883077](https://doi.org/10.5281/zenodo.22883077)
- Concept DOI: [10.5281/zenodo.16880787](https://doi.org/10.5281/zenodo.16880787)
- Method paper: [Long-Term Observations of Galactic Cosmic Ray LET Spectra in
  Lunar Orbit by LRO/CRaTER](https://doi.org/10.1029/2020SW002543)
- License: CC-BY-4.0.

This dataset is reserved for an instrument-response validation stage. CRaTER's
reported variable is derived from energy deposited along accepted telescope
tracks and is commonly presented as LET; it is not interchangeable with a
generic slab's unrestricted stopping power. A comparison is valid only after
implementing CRaTER detector materials and dimensions, coincidence/track cuts,
acceptance, shielding, lunar obstruction/albedo, time matching, and the same
LET/lineal-energy estimator.

## What the tests establish

`tests/test_reference_spectra.py`, `stage2_source_adapters`, and
`test_stage2_source_tool.py` establish only that:

- committed bytes match reviewed checksums;
- published metadata, units, shapes, timestamps, and uncertainties survive;
- energy coordinates are positive and ordered where their meaning is known;
- invalid/fill/non-finite values and silent direction loss are detected;
- the approved event windows, marker, sample timing, and quality policy remain
  explicit; and
- the implemented unit/Jacobian and event-integration transforms agree with
  analytic and independently calculated spot checks.

It does not establish transport correctness or observational agreement. Those
gates are specified in `development_and_validatoin.md`.

The CTest `validation_figures` target visualizes the AMS and GOES adapter
comparisons and the CRaTER reference product in PNG/SVG form under the build
directory. CRaTER remains explicitly observation-only: no model curve is
drawn until the instrument-faithful response and preregistered comparison gate
are implemented.
