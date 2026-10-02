# Downloaded references for CME3D02: 12 July 2012 CME

## Install and run

Extract both the updated source overlay and
`SWCME_20120712_reference_data_20261001.tar.gz` from the AMPS repository root.
The data archive installs this fixed, repository-owned location:

```
srcSEP3D/validation/reference_data/CME3D02/manifest.json
```

No data-path flag is needed. From the AMPS root, choose a fresh output directory:

```bash
python3 srcSEP3D/validation/run_swcme_coupling_validation.py --output-dir test_output/swcme-reference-check
python3 srcSEP3D/validation/run_validation.py --case CME3D02 --output-dir test_output/swcme-phase-v
python3 srcSEP3D/test/run_tests.py --all --output-dir srcSEP3D/test_output/all
```

The standard `--all` runner discovers CME3D02 automatically. Existing
`--validation-data /path/to/evidence` (standard runner) and
`--evidence-root /path/to/evidence` (Phase-V runner) override the installed data;
each root must contain `CME3D02/manifest.json`. The dedicated runner's
`--bundle` accepts the case directory itself. The earlier
`test/run_coupled_sep_corona.py` aggregates the shared coronal-model tests and
native initialization callbacks; it does not run this observational campaign.

With references alone the expected result is **SKIP, `reference_ready=true`**,
with the message that native shock history and `native-manifest.json` are still
required. Every declared reference checksum and required CSV contract is
verified before that SKIP. Missing/corrupt declared files are ERROR. This makes
data installation verifiable without manufacturing a model run. The dedicated
`--require-evidence` option returns exit 2 for the incomplete campaign; its
ordinary reference check and the general runners allow explicit SKIP.

## Frozen data snapshot

The observation interval is **2012-07-12 00:00 UTC to 2012-07-17 00:00 UTC**,
with the exclusive upper boundary. Download URLs, selected HAPI parameters,
timestamps, HTTP metadata where available, byte counts and SHA256 hashes are in
`raw/download_receipt.json`. Cached downloads retain their original receipt;
an initial cache without a receipt uses its file modification time, explicitly
labeled as such. The manifest hashes the raw downloads, selected profiles,
normalized products, this description and preview figures.

| Product | Dataset or source | Records in this snapshot | Purpose |
|---|---|---:|---|
| `wind_protons.csv` | NASA `WI_H1_SWE` | 4,193; 3,503 accepted | Proton bulk speed, density, scalar temperature; upstream context |
| `wind_magnetic_1min.csv` | NASA `WI_H0_MFI@0` | 7,200 | Definitive MFI magnitude and GSE components |
| `wind_magnetic_3sec_shock.csv` | NASA `WI_H0_MFI@1`, July 14 15:00–20:00 UTC | 6,000 | Resolve the local shock signature |
| `wind_pesa_ion_moments.csv` | NASA `WI_PLSP_3DP` | 16,982 | Separate ground-computed PESA Low ion moments; context only |
| `wind_orbit.csv` | NASA `WI_OR_PRE` | 720 | Predicted orbit, heliocentric Cartesian position and radius |
| `stereo_a_positions.csv`, `stereo_b_positions.csv` | NASA `STA/STB_COHO1HR_MERGED_MAG_PLASMA` | 120 each | Hourly radius and declared heliographic longitude/latitude |
| `helcats_time_elongation.csv` | HELCATS HIGeoCAT v06 raw profiles | 109 | Repeated brightness-feature measurements in angular coordinates |
| `arrival.json`, `ipshocks_wind_event.json` | IPShocks Zenodo v1 and CfA event 00525 | One selected Wind shock | Withheld spacecraft arrival marker |

All normalized dimensional columns use SI units except columns ending in
`_deg`, which use degrees. Metadata and raw values are retained for auditing.
No resampling, smoothing, interpolation through instrument gaps or replacement
of rejected measurements is performed. Wind MFI vectors retain **GSE**; STEREO
longitude/latitude retain the product's declared heliographic coordinates.
No silent conversion into the AMPS event coordinate system is made.

SWE fit flags must exceed 1, and speed, density and thermal width must be finite
and positive. Fill values are masked before unit conversion. Scalar proton
temperature is `m_p * W_nonlin^2 / (2*k_B)`; it is not total plasma pressure.
Rejected records keep their timestamp and flag with empty physical columns.
SWE has serious fit/coverage loss around this event: **the declared two-hour
downstream context contains no accepted SWE samples**. Its statistics are JSON
null, not invented medians. `WI_PLSP_3DP` ion moments are provided as a separately
identified product with its own `VALID==1` and fill checks; they are not spliced
into SWE or promoted to plasma acceptance evidence. MFI invalid records are
also retained and flagged.

## Arrival and feature identity

The frozen Wind marker is **2012-07-14T17:39:09Z**, selected from the immutable
IPShocks v1 CSV, DOI <https://doi.org/10.5281/zenodo.19730292>.
The downloader verifies its publisher MD5
`7432ac3cd68094b32fad2a0b01e33c28`; the manifest additionally uses SHA256.
CfA event 00525 independently reports **17:39:07.5 UTC, +/-20 s**. We adopt
20 s as an explicitly declared analysis timing uncertainty; it is not an
uncertainty field published with the IPShocks timestamp. Neither timestamp is
bow-shock shifted OMNI data.

Wind's radius at the shock is **1.00773304 AU**, obtained by linearly
interpolating the Cartesian `HEC_POS` from `WI_OR_PRE` and taking its norm.
This is the **Predicted Orbit** product. A usable Definitive Orbit HAPI response
was not available for this preparation. The CfA page's spacecraft-position
numbers disagree with the NASA orbit and modern catalog, so those numbers are
not used. A publication campaign should review the orbit accuracy and mapping
of the spherical model front to the Wind measurement.

Selected HELCATS events are `HCME_A__20120712_02` at PA75 and
`HCME_B__20120712_01` at PA260. The A event ending `_01` precedes the afternoon
eruption and is excluded. Raw repeat-trace IDs and position angles are retained.
These tracks measure **CME brightness-feature elongation**, not an independently
identified radial shock front. HIGeoCAT's geometric fits are retained as context.
They must not be used to create purported held-out shock radii from a fit to the
entire same track. No radial shock track or 20-Rsun launch fit is fabricated.

## Native evidence still required

The supplied manifest uses `srcsep3d-swcme-reference-v1` and an **arrival-only**
comparison scope. Keep the downloaded references intact. A later campaign adds
`native-manifest.json` and the native files it describes under a new subdirectory
such as `native-run/`. Copy the existing
`validation/templates/swcme_coupling_evidence.template.json` to that native
manifest; populate the schema, native producer, input/log/history descriptors,
executable identity, actual launch epoch, clock, MPI ranks, source-disabled
state, frozen thresholds and independent parameter provenance. Its observation
track descriptor is unused in this arrival-only mode. Input/log/history file
names are relative to this case directory. Native file hashes are checked
separately; rerunning reference preparation does not checksum native products as
observations.

The runner takes the arrival only from the checksum-verified reference marker.
It does not derive the launch clock, speed or drag from Wind arrival. A native
arrival comparison produces the 600-dpi PNG and vector EPS, JSON metrics and
caption. A successful arrival diagnostic does not qualify a radial shock track,
CME plasma/IMF evolution or a production release. For the full radial evolution
test, independently identify the shock near 20 Rsun, freeze its launch fit,
and prepare a separate reviewed full-track evidence bundle under the original
native evidence schema. See `validation/SWCME_COUPLING_VALIDATION.md`.

Native propagation-only/source-off input support and installed-provider history
export still need implementation. The current model supplies shock kinematics
on an analytic Parker background; evolving CME plasma/IMF requires a further
background adapter. This data package does not change those capabilities.

## Refresh and reproduce

The supplied archive works offline. Python 3.9+ is sufficient for download,
normalization and checksum checks; Matplotlib/NumPy are needed only to regenerate
figures or plot a native comparison. To repeat acquisition at the installed path:

```bash
python3 srcSEP3D/validation/cases/2012-07-12/prepare_reference_data.py
```

Existing raw files are reused only with matching checksums and source URLs.
`--refresh` retrieves new service responses and produces a new snapshot;
preserve the bundled snapshot for a reproducible publication. `--no-figures`
allows preparation without Matplotlib. `--output-dir` can create a separate case
directory; `--raw-dir` imports a downloaded raw cache.

`wind-reference.png/.eps` show only measured Wind SWE/MFI data with a shock
marker; invalid fits and telemetry gaps remain breaks in the curves.
`helcats-elongation-reference.png/.eps` show the raw repeated angular tracks.
These are reference previews with **no model curves**, not completed validation
figures. Publication captions must retain the stated instrument/feature roles.

## Sources and acknowledgements

* NASA/GSFC Space Physics Data Facility, Wind MFI/SWE/3DP teams and STEREO teams:
  <https://cdaweb.gsfc.nasa.gov/hapi/>,
  <https://cdaweb.gsfc.nasa.gov/misc/NotesW.html>.
* IPShocks v1 snapshot: <https://doi.org/10.5281/zenodo.19730292>.
  Use the requested acknowledgement: This paper uses data from the Heliospheric
  Shock Database, generated and maintained at the University of Helsinki.
* CfA Interplanetary Shock Database, Michael L. Stevens and Justin C. Kasper:
  <https://lweb.cfa.harvard.edu/shocks/wi_data/00525/wi_00525.html>.
* HELCATS consortium and STEREO/SECCHI teams:
  <https://www.helcats-fp7.eu/catalogues/wp3_cat.html>;
  Barnes et al. (2019), <https://doi.org/10.1007/s11207-019-1444-4>.
* Hu et al. (2016), event interpretation:
  <https://doi.org/10.3847/0004-637X/829/2/97>.

Retain source acknowledgements and consult each source's citation/usage terms
when publishing the comparison; this package does not replace them with a new
license.
