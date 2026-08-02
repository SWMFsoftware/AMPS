# C10 observational data and reference-solution construction

## 1. Purpose of this directory

This directory contains the **external observational inputs** used by the C10
validation test and the provenance records needed to reproduce the numerical
reference solution.

C10 compares AMPS storm-time proton access boundaries with boundaries extracted
from the NOAA/NCEI historical POES/MetOp SEM-2 Level-2 archive. The publications
listed in [`../REFERENCES.md`](../REFERENCES.md) define the instrument, channel
response, data format, quality issues, and boundary-extraction method, but they
do not contain the complete pass-level numerical reference. The numerical C10
reference must therefore be generated from the archived measurements.

The normal directory layout is:

```text
srcEarth/test/C10/
├── data/
│   ├── README.md                         this document
│   ├── reference_source/                 downloaded NCEI daily files
│   │   └── download_manifest.json        URL/size/SHA-256 inventory
│   └── ts05_driving.txt                  optional local copy of the C9 driver
├── C10_poes_boundary_crossings.csv       generated pass-level observations
├── reference_C10_poes_meped_boundary.csv generated runner reference
├── C10_reference_manifest.json           generated scientific provenance
└── C10_reference_summary.json            generated extraction summary
```

Raw NCEI daily files should remain below `data/reference_source/`. Generated
crossing and reference products are written in the C10 directory by default so
that `run_C10.py` can find them without an additional path argument.

## 2. What is—and is not—the C10 reference solution

The authoritative observational product is:

```text
C10_poes_boundary_crossings.csv
```

Each row represents one measured inbound or outbound satellite-pass crossing
of the 50%-of-polar-cap-flux boundary for one MEPED omnidirectional proton
channel.

The file consumed directly by the C10 runner is:

```text
reference_C10_poes_meped_boundary.csv
```

It is derived from the crossing table by aggregating measurements into
overlapping time windows, hemispheres, channels, and MLT bins. It is not a
synthetic curve and it must not be replaced by digitized points from a paper or
by a published empirical fit.

The checked-in five-line `reference_C10_poes_meped_boundary.csv` is only a
placeholder. `run_C10.py --validate-references` deliberately rejects it until
`build_poes_reference.py` successfully writes archive-derived rows.

## 3. Measurement source

Use the NOAA National Centers for Environmental Information historical
POES/MetOp Space Environment Monitor archive:

```text
https://www.ncei.noaa.gov/data/poes-metop-space-environment-monitor/
```

For the December 2006 C10 event, the intended processed product is the
historical Level-2, version `v01r00`, 16-second product:

```text
access/l2/v01r00/txt/2006/
access/l2/v01r00/cdf/2006/
```

The ASCII product is recommended because it is human-readable and does not
require a CDF library. The CDF reader is available for users who maintain a CDF
mirror.

The standard C10 download includes these historical spacecraft tokens:

| Token | Spacecraft |
|---|---|
| `n15` | NOAA-15 |
| `n16` | NOAA-16 |
| `n17` | NOAA-17 |
| `n18` | NOAA-18 |
| `m02` | MetOp-02 / MetOp-A |

The preferred source interval is 5–16 December 2006. This provides event
context around the main 14–16 December validation interval and allows the
source measurements to be inspected independently of the final extraction
window.

## 4. Python setup without administrator privileges

### 4.1 Data download requires no extra packages

`download_poes_sem2.py` uses only the Python standard library. The raw archive
files can be downloaded before creating an environment or installing any
packages.

### 4.2 Recommended `/bin/tcsh` environment

The following procedure installs everything under the user's home directory and
does not require `sudo`, `apt`, or write access to `/usr`.

Run from `srcEarth/test/C10`:

```tcsh
python3 -m pip install --user virtualenv
setenv PATH "${HOME}/.local/bin:${PATH}"
python3 -m virtualenv .venv
source .venv/bin/activate.csh
python -m pip install "pip<25.1"
python -m pip install -r requirements.txt
```

The `pip<25.1` constraint is useful on systems that still provide Python 3.8.
After activation, verify the scientific dependencies:

```tcsh
python -c "import numpy, matplotlib, aacgmv2; print('C10 dependencies are available')"
```

For ASCII input, `cdflib` is not used during reference construction even though
it is listed as an optional dependency. For CDF input, verify it as well:

```tcsh
python -c "import cdflib; print('CDF support is available')"
```

To leave the environment:

```tcsh
deactivate
```

### 4.3 Standalone `virtualenv` fallback

If `python3 -m pip` cannot install `virtualenv`, obtain the standalone
application in the user's home directory:

```tcsh
curl -L https://bootstrap.pypa.io/virtualenv.pyz -o "$HOME/virtualenv.pyz"
python3 "$HOME/virtualenv.pyz" .venv
source .venv/bin/activate.csh
python -m pip install "pip<25.1"
python -m pip install -r requirements.txt
```

Do not run `apt install python3-venv` on a system where you do not have root
access. The C10 reference can be built entirely in a user-local environment.

## 5. Downloading the NCEI source data

All commands in the remainder of this document are run from:

```text
srcEarth/test/C10
```

### 5.1 Automatic archive discovery and download

The recommended command is:

```tcsh
python download_poes_sem2.py --start 2006-12-05 --end 2006-12-16 --satellites n15,n16,n17,n18,m02 --format txt --output-dir data/reference_source
```

The downloader:

1. starts at the official NCEI Level-2 2006 directory;
2. recursively reads ordinary directory-index pages;
3. selects files matching the requested date and spacecraft tokens;
4. downloads the matching daily products without changing their contents; and
5. writes `data/reference_source/download_manifest.json`.

Each manifest entry records:

```text
original URL
local path
file name
byte count
SHA-256 digest
download/reuse status
```

The manifest is part of the scientific provenance. Do not delete or manually
edit it after a successful download.

### 5.2 Re-download existing files

Use `--overwrite` to replace existing files:

```tcsh
python download_poes_sem2.py --start 2006-12-05 --end 2006-12-16 --satellites n15,n16,n17,n18,m02 --format txt --output-dir data/reference_source --overwrite
```

This is appropriate after an interrupted transfer or when the file inventory
contains zero-byte files.

### 5.3 Manual URL-list fallback

Some institutional networks allow direct file downloads but block recursive
indexing. In that situation:

1. open the official NCEI `txt/2006/` directory in a browser;
2. identify the exact daily-file URLs for the requested spacecraft and dates;
3. save one exact HTTPS URL per line in `poes_urls.txt`; and
4. run:

```tcsh
python download_poes_sem2.py --url-list poes_urls.txt --format txt --output-dir data/reference_source
```

Blank lines and lines beginning with `#` are ignored. This route produces the
same checksum manifest as automatic discovery.

### 5.4 Expected number of source files

A request for 12 days and five spacecraft can produce as many as 60 daily files,
but a smaller number is not automatically an error. A spacecraft may not have a
published product for every requested day. Scientific usability is determined
from the actual observation and crossing counts, not from an assumed exact file
count.

## 6. Mandatory source-file checks

Perform these checks before building the reference.

### 6.1 List the inventory

```tcsh
find data/reference_source -type f | sort
```

### 6.2 Inspect the download manifest

```tcsh
python -m json.tool data/reference_source/download_manifest.json | less
```

### 6.3 Detect zero-byte files

```tcsh
find data/reference_source -type f -size 0 -print
```

A zero-byte file contains no measurements. Depending on the installed C10
revision, the builder either skips it with a warning or stops with a diagnostic
such as:

```text
File preview: <empty file>
```

The preferred correction is to remove and redownload the file:

```tcsh
find data/reference_source -type f -size 0 -delete
python download_poes_sem2.py --start 2006-12-05 --end 2006-12-16 --satellites n15,n16,n17,n18,m02 --format txt --output-dir data/reference_source
```

When the same satellite-day product remains unavailable, remove the empty file
and continue with the rest of the constellation. Record the missing product in
the run notes and evaluate its effect on crossing and MLT coverage.

### 6.4 Detect HTML or archive error pages saved as data

```tcsh
foreach f (`find data/reference_source -type f -name '*.txt'`)
    grep -qi -E '<html|access denied|not found|forbidden' "$f"
    if ($status == 0) echo "possible non-data response: $f"
end
```

Any identified response page must be deleted and downloaded again from the
correct file URL.

### 6.5 Inspect one ASCII product

```tcsh
set f = `find data/reference_source -type f -name '*.txt' -size +0c | head -1`
echo "$f"
head -5 "$f"
```

Historical products may contain either:

- an explicit column-name header; or
- a documented fixed-order numerical layout without a header.

Recent C10 parsers recognize both forms. The explicit-header form must contain
fields equivalent to:

```text
year mo dy hr mi second sslat sslon mepomp6 mepomp7 mepomp8 mepomp9
```

For a headerless product, inspect the first nonblank record:

```tcsh
awk 'NF && $1 !~ /^#/ {print "fields in first data row:", NF; exit}' "$f"
```

The documented headerless 16-second layout contains at least 41 numerical
fields. A nonempty file that matches neither layout is treated as malformed and
must not be silently included in the reference.

### 6.6 Verify local files against the manifest

The builder writes a second scientific manifest, but it is useful to verify the
raw download inventory independently. A simple spot check is:

```tcsh
sha256sum "$f"
```

Compare the value with the corresponding `sha256` field in
`download_manifest.json`.

## 7. Building the archive-derived reference

### 7.1 Standard reference-build command

```tcsh
python build_poes_reference.py --input-dir data/reference_source --event-start 2006-12-14T00:00:00Z --event-end 2006-12-16T12:00:00Z --crossings-output C10_poes_boundary_crossings.csv --reference-output reference_C10_poes_meped_boundary.csv --manifest-output C10_reference_manifest.json --summary-output C10_reference_summary.json
```

A successful run ends with messages naming all four outputs and returns status
zero. In `tcsh`, check the status immediately:

```tcsh
echo "build exit status = $status"
```

### 7.2 Capture a complete build log

To retain both standard output and errors:

```tcsh
python build_poes_reference.py --input-dir data/reference_source --event-start 2006-12-14T00:00:00Z --event-end 2006-12-16T12:00:00Z --crossings-output C10_poes_boundary_crossings.csv --reference-output reference_C10_poes_meped_boundary.csv --manifest-output C10_reference_manifest.json --summary-output C10_reference_summary.json |& tee C10_build_reference.log
```

With a pipeline, the shell status may reflect `tee` rather than the Python
process. Therefore, verify the output files rather than assuming the build
succeeded merely because the prompt returned.

### 7.3 Default extraction settings

| Option | Default | Role |
|---|---:|---|
| `--default-altitude-km` | 850 | fallback when altitude is absent |
| `--minimum-pass-abs-lat-deg` | 45 | minimum absolute AACGM latitude retained in a polar pass |
| `--polar-plateau-abs-lat-deg` | 75 | latitude threshold for estimating polar-cap flux |
| `--minimum-polar-samples` | 4 | minimum valid samples in the plateau estimate |
| `--minimum-plateau-to-low-ratio` | 2 | rejects passes without adequate polar enhancement |
| `--minimum-leg-samples` | 4 | minimum samples on an inbound/outbound leg |
| `--window-hours` | 2 | width of the aggregation window |
| `--step-hours` | 1 | separation of adjacent window centers |
| `--mlt-bin-hours` | 3 | magnetic-local-time bin width |
| `--minimum-crossings-per-cell` | 1 | observations required for a populated reference cell |

Display the exact options supported by the installed code with:

```tcsh
python build_poes_reference.py --help
```

Do not relax quality thresholds merely to force a nonempty reference. First
inspect individual passes, channel values, source formats, and event coverage.

## 8. Scientific extraction performed by the builder

The builder follows this sequence.

### 8.1 Read and normalize measurements

For every supported daily file, the reader:

1. infers the spacecraft identity from the file name;
2. reads the 16-second epoch and geographic sub-satellite position;
3. reads MEPED omnidirectional proton channels P6–P9;
4. rejects negative fill values, nonfinite values, and unphysical large archive
   values;
5. uses the file altitude when available or the explicit fallback altitude;
6. recomputes AACGM latitude and MLT with `aacgmv2`; and
7. retains source-file identity and checksum information.

### 8.2 Map channel thresholds to nominal proton rigidities

The nominal channels are:

| Channel | Nominal lower threshold | Assigned rigidity |
|---|---:|---:|
| P6 | 16 MeV | 0.174013525 GV |
| P7 | 36 MeV | 0.262395866 GV |
| P8 | 70 MeV | 0.369131538 GV |
| P9 | 140 MeV | 0.531334344 GV |

The rigidity is calculated from the proton kinetic-energy threshold. This is a
transparent nominal mapping, not a full convolution of the integral detector
response with an event spectrum.

### 8.3 Separate polar passes and pass legs

Measurements are grouped into single-spacecraft, single-hemisphere polar
passes. Each pass is split near its highest absolute AACGM latitude into:

```text
inbound/equator-to-pole leg
outbound/pole-to-equator leg
```

The two legs are analyzed independently because they occur at different times
and often at different MLT.

### 8.4 Estimate the polar-cap plateau

For one spacecraft pass, channel, and hemisphere, valid measurements at:

```text
|AACGM latitude| >= polar-plateau threshold
```

are used to estimate the polar-cap channel level. The implemented reference
level is based on the median of the retained polar samples, making it less
sensitive to isolated spikes than a simple mean.

### 8.5 Locate the half-flux boundary

The boundary level is:

```text
0.5 × polar-cap plateau
```

For each leg, the algorithm locates adjacent records bracketing that level and
linearly interpolates the crossing time, geographic position, AACGM latitude,
MLT, and altitude. Each accepted crossing becomes one row in
`C10_poes_boundary_crossings.csv`.

### 8.6 Aggregate crossings for the AMPS comparison

Pass-level crossings are then grouped by:

```text
overlapping time window
channel/nominal rigidity
hemisphere
MLT bin
```

Cells lacking enough observations are retained explicitly with `missing=TRUE`.
The runner excludes missing cells from numerical metrics; it does not fill or
interpolate them.

## 9. Generated products

### 9.1 `C10_poes_boundary_crossings.csv`

This is the primary measurement-level product. It includes fields such as:

```text
satellite
channel
energy_threshold_mev
assigned_rigidity_gv
hemisphere
pass_id
leg
crossing_time_utc
geographic_lat_deg
geographic_lon_deg
altitude_km
aacgm_lat_deg
mlt_hour
polar_plateau_flux
half_plateau_flux
boundary_uncertainty_deg
quality_flags
source_file
source_sha256
```

Use this file to inspect individual satellite passes, spacecraft subsets,
hemispheric differences, and MLT coverage.

### 9.2 `reference_C10_poes_meped_boundary.csv`

This is the numerical table consumed by `run_C10.py`. Important columns include:

```text
interval_midpoint_utc
interval_start_utc
interval_end_utc
rigidity_gv
energy_threshold_mev
channel
hemisphere
mlt_hour
boundary_aacgm_lat_deg
sigma_deg
altitude_km
n_crossings
satellites
missing
source
source_ref
notes
```

A valid file must contain numerical rows after the header. A five-line file
beginning with `C10 reference placeholder` has not been built.

### 9.3 `C10_reference_manifest.json`

This scientific manifest records:

```text
extraction settings
source-file paths, sizes, and SHA-256 values
channel threshold and rigidity mapping
number of pass-level crossings
number of total and populated reference cells
generation time
```

Archive this manifest with every C10 result. It is distinct from the download
manifest: the download manifest records network acquisition, while the
reference manifest records the exact sources and settings used to create the
numerical reference.

### 9.4 `C10_reference_summary.json`

The summary provides a quick integrity check, including:

```text
number of source files
number of observations in the event window
number of extracted crossings
number of total reference cells
number of nonmissing reference cells
spacecraft represented
channels represented
first and last observation times
output paths and checksums
```

## 10. Verifying a successful build

Run all of the following:

```tcsh
ls -l C10_poes_boundary_crossings.csv reference_C10_poes_meped_boundary.csv C10_reference_manifest.json C10_reference_summary.json
wc -l C10_poes_boundary_crossings.csv reference_C10_poes_meped_boundary.csv
head -5 reference_C10_poes_meped_boundary.csv
python -m json.tool C10_reference_summary.json
```

The reference CSV must contain more than the placeholder/header lines. The
summary must report nonzero `n_observations`, `n_crossings`, and
`n_nonmissing_reference_cells`.

Validate the generated table and the shared December-2006 TS05 driver:

```tcsh
python run_C10.py --validate-references --reference reference_C10_poes_meped_boundary.csv
python run_C10.py --validate-driver
```

Or as one command:

```tcsh
python run_C10.py --validate-references --reference reference_C10_poes_meped_boundary.csv && python run_C10.py --validate-driver
```

## 11. Troubleshooting

### 11.1 `no reference rows parsed`

Example:

```text
C10 reference could not be loaded: no reference rows parsed
```

Cause: the runner is still reading the five-line placeholder or a build that
failed before writing rows.

Checks:

```tcsh
wc -l reference_C10_poes_meped_boundary.csv
head -5 reference_C10_poes_meped_boundary.csv
ls -l C10_reference_summary.json
```

Correction: successfully rerun `build_poes_reference.py` before validation.

### 11.2 `File preview: <empty file>`

Cause: a downloaded daily product has zero bytes.

Locate it:

```tcsh
find data/reference_source -type f -size 0 -print
```

Delete and redownload it, or remove it and continue with the remaining
constellation if the archive has no usable product for that spacecraft-day.
Never create a fabricated replacement.

### 11.3 Header cannot be recognized

Example:

```text
could not locate an explicit 16-second Level-2 header and could not recognize a headerless NCEI data row
```

Check whether the file is empty, HTML, truncated, compressed under an unexpected
suffix, or from a different data product. For a nonempty ASCII file:

```tcsh
head -10 path/to/file.txt
awk 'NF && $1 !~ /^#/ {print NF; exit}' path/to/file.txt
```

Use the documented 16-second Level-2 product rather than an unrelated summary
or telemetry file.

### 11.4 `No valid observations fall inside the requested event interval`

Check:

- the source-file dates;
- `--event-start` and `--event-end`;
- whether the downloaded products contain the intended spacecraft and year;
- whether time columns were parsed correctly; and
- whether all relevant channel records were rejected as fill values.

### 11.5 `No boundary crossings were extracted`

Do not immediately loosen thresholds. First inspect:

- whether an SEP enhancement is present in P6–P9;
- whether each pass reaches the configured polar-plateau latitude;
- the number of valid polar samples;
- the plateau-to-low-latitude ratio;
- P8/P9 alternating telemetry values; and
- individual flux-versus-AACGM-latitude profiles.

Only after those checks should sensitivity runs change plateau latitude, sample
counts, or contrast thresholds.

### 11.6 `aacgmv2` import failure

Activate the environment and verify installation:

```tcsh
source .venv/bin/activate.csh
python -m pip install -r requirements.txt
python -c "import aacgmv2; print(aacgmv2.__version__)"
```

### 11.7 CDF input failure

CDF is optional. For CDF source files:

```tcsh
python -m pip install cdflib numpy
python -c "import cdflib, numpy; print('CDF dependencies available')"
```

When possible, use the ASCII products first because that is the primary tested
path for C10.

## 12. TS05 driver used by C10

C10 also requires the checksum-verified December-2006 TS05 driver used by C9.
The runner searches first for:

```text
data/ts05_driving.txt
```

and then for the sibling C9 copy:

```text
../C9/data/ts05_driving.txt
```

A separate local copy is optional when the sibling C9 file exists. The runner
verifies the driver's SHA-256 before a scientifically eligible run.

## 13. Running C10 after the reference is built

A routine two-branch run from the C10 directory is:

```tcsh
python run_C10.py --solver BOTH --profile ROUTINE --reference reference_C10_poes_meped_boundary.csv --boundary-cutoff effective --interval-samples 1 --output-root test_output/C10_real --amps /home/vtenishe/T11/AMPS/amps --cutoff-scan-n 120 --shell-lon-res-deg 15 --shell-lat-res-deg 2 --access-abs-lat-min-deg 45 --access-abs-lat-max-deg 85 --dynamic-chunk 32 -np 4 -nt 16
```

The runner produces C9-style time-series and residual plots, an
observed-versus-modeled scatter plot, an MLT comparison plot, comparison CSVs,
and branch result JSON files. See [`../README.md`](../README.md) for complete
runner behavior and acceptance criteria.

## 14. Reproducibility and archiving

A scientific C10 result should retain:

```text
all original NCEI daily files
data/reference_source/download_manifest.json
C10_poes_boundary_crossings.csv
reference_C10_poes_meped_boundary.csv
C10_reference_manifest.json
C10_reference_summary.json
C10_build_reference.log
the TS05 driver and checksum
all generated AMPS input files and logs
raw AMPS shell and penumbra products
comparison CSVs, plots, and result JSON files
AMPS git commit and local diff
Python version and installed-package inventory
```

Record the Python environment with:

```tcsh
python --version
python -m pip freeze > C10_python_environment.txt
```

Do not overwrite a frozen baseline without retaining the previous source files,
manifests, extraction settings, and checksums.

## 15. Related documentation

- [`../README.md`](../README.md): scientific definition, runner, metrics, and
  complete C10 workflow.
- [`../INSTALL.md`](../INSTALL.md): compact installation and run commands.
- [`../REFERENCES.md`](../REFERENCES.md): dataset, instrument, calibration, and
  cutoff-method references.
- [`../VALIDATION.md`](../VALIDATION.md): completed software checks and required
  archive validation.
- [`../CHANGES.md`](../CHANGES.md): implementation history.
