# Archived observations of the 29 November 2020 CME shock

The legacy `SEP--CME-shock-evolution` test reads this
[observations.json](observations.json) through the shared Python shock fitter.
It contains 59 manual marks on actual STEREO-A/SECCHI COR2 and SOHO/LASCO C2
and C3 white-light exposures. These are sight lines toward the candidate
outer compression edge. An image measures line-of-sight brightness. It
cannot supply the three-dimensional depth of a mark or identify whether
each part of the front is locally superfast. Kollhoff et al. identify the wave
in C2 and COR2 in Fig. 3. Kouloumvakos et al. combine these instruments
with other viewpoints in their ellipsoidal reconstruction in Figs. 2, 4
and 5. The annotations here were made on the archived images. They are
independent of the papers' fitted coordinates.
[Kollhoff et al. 2021](https://doi.org/10.1051/0004-6361/202140937),
[Kouloumvakos et al. 2022](https://doi.org/10.1051/0004-6361/202142515).

## Archive retrieval and feature review

I used the 2020-11-29 SECCHI/COR2-A daily `d4c2A` sequence and the SOHO/LASCO
level 0.5 daily C2 and C3 inventories. The manifest records every input and
difference exposure by its download URL, local filename, SHA-256, primary
FITS astrometry and time cards. The local aliases `AHHMMSS.fts` correspond
to the SECCHI archive names `20201129_HHMMSS_d4c2A.fts`. The LASCO aliases
are the archive names. The raw files can be kept outside Git.
[SECCHI daily archive](https://stereo-ssc.nascom.nasa.gov/pub/ins_data/secchi/L0/a/img/cor2/20201129/),
[LASCO C2 daily inventory](https://umbra.nascom.nasa.gov/pub/lasco_level05/201129/c2/img_hdr.txt),
[LASCO C3 daily inventory](https://umbra.nascom.nasa.gov/pub/lasco_level05/201129/c3/img_hdr.txt).

| Instrument | FITS alias | Mid-exposure UTC | Subtracted exposure | Marks | Use |
| --- | --- | --- | --- | ---: | --- |
| SOHO/LASCO C2 | `22793508.fts` | 13:26:02.179050 | `22793507.fts` | 9 | Fit |
| SOHO/LASCO C2 | `22793509.fts` | 13:36:38.368400 | `22793508.fts` | 8 | Reserved |
| STEREO-A/COR2 | `A133900.fts` | 13:39:13.006000 | `A132400.fts` | 9 | Reserved |
| SOHO/LASCO C2 | `22793510.fts` | 13:48:18.092100 | `22793509.fts` | 5 | Reserved |
| STEREO-A/COR2 | `A135400.fts` | 13:54:13.005000 | `A133900.fts` | 10 | Fit |
| SOHO/LASCO C3 | `32643517.fts` | 13:54:14.906650 | `32643516.fts` | 8 | Fit |
| STEREO-A/COR2 | `A142400.fts` | 14:24:13.005000 | `A135400.fts` | 10 | Fit |

The four fit exposures cover the early C2 profile, near-simultaneous COR2
and C3 views separated by about 58° in HCI longitude, and a later COR2 profile.
The reserved exposures retain intermediate evolution and the C2 upper arc
near the image boundary. All marks from a single image belong to one role.
The current model was rejected during fitting, so the reserved marks have
no validation score.

For feature review, I applied the primary FITS `BSCALE` and `BZERO`, divided
by exposure seconds and formed the fractional difference
`(event_count_rate - preceding_count_rate) / max(preceding_count_rate, 1)`.
The floor is one image-count unit per second and keeps the display defined
at zeros in the occulted area. The per-instrument median window and
grayscale stretch in [sources.json](sources.json) suppress isolated points
and show the surrounding brightening. The review rasters retain the
original detector grids. These quicklooks guide position marking and
supply no absolute brightness, density or compression measurement.
Compare the moving outer edge with the preceding view and published
outlines. Leave gaps where the occulter, a streamer, detector artifacts,
or the bright ejecta obscure the feature. The lower arc in the late COR2
view is more diffuse, so its marks have a wider pixel uncertainty.

[picks.csv](picks.csv) stores exactly four fields per manual mark.

| Column | Definition |
| --- | --- |
| `frame` | Exact local FITS alias in `sources.json`. |
| `x_px` | Zero-based column centre in the unchanged FITS raster, increasing rightward. |
| `y_px` | Zero-based row centre in the unchanged FITS raster, increasing upward in the review plot. |
| `sigma_px` | Estimated transverse uncertainty of identifying the edge, in pixels. The 2020 marks use 14–18 pixels for C2, 4–6 for C3 and 20–28 for COR2. |

## Astrometry, time and input units

- SECCHI/COR2 supplies an already rectified image with primary
  `HPLN-TAN/HPLT-TAN` WCS. Use `DATE-AVG` for the summed exposure and the
  header `HCIX_OBS/HCIY_OBS/HCIZ_OBS` in metres. Use `DATE-OBS` and
  `DATE-END` to identify the first and last exposures.
  [SECCHI FITS keyword definitions](https://sohowww.nascom.nasa.gov/solarsoft/stereo/secchi/doc/FITS_keywords.pdf).
- The selected rectified LASCO C2/C3 level 0.5 files have `SOLAR-X/SOLAR-Y`
  with equal scales in arcseconds per pixel. Their legacy `CROTA2` roll
  provides the matrix orientation. In these files `CRVAL1/2` are zero.
  The mid-exposure UTC is `DATE-OBS` plus `TIME-OBS` plus half of
  `EXPTIME` seconds. This uses the archived nominal level 0.5 pointing.
  Qualification for more precise science would require calibrated
  pointing and its uncertainty.
  [NRL LASCO FITS keyword definitions](https://lasco-www.nrl.navy.mil/keywords.html).
- Both adapters add one to the CSV pixel coordinates before subtracting
  FITS `CRPIX`. The pixel displacement is mapped by
  `diag(CDELT) × PC` in arcseconds, converted to radians, and transformed
  by the inverse tangent-plane projection about `CRVAL`. The final unit
  vector is rotated into HCI using local solar west, projected solar north,
  and the direction from the spacecraft toward the Sun. This preserves
  the finite observer distance and the correct sign of east/west.
  [Thompson 2006 solar image WCS](https://doi.org/10.1051/0004-6361:20054262).
- I requested Horizons targets `-21` (SOHO) and `-234` (STEREO-A), with
  centre `500@10` (Sun body centre), geometric vectors, `VEC_TABLE=2`,
  `REF_SYSTEM=ICRF`, `REF_PLANE=ECLIPTIC`, `KM-S`, UTC and one-minute
  spacing. The manifest stores the exact URLs, original response hashes and
  the list of `[UTC, x_km, y_km, z_km]` knots. For SOHO the packer linearly
  interpolates the bracketing Cartesian positions at exposure time. The
  J2000 ecliptic-to-HCI rotation uses the IAU76 obliquity (84381.448
  arcsec) and the solar pole (ICRF right ascension 286.13°, declination
  63.87°). HCI +X lies along the ascending node of the solar equator on
  that ecliptic. HCI +Z points along the solar pole. The transformed
  STEREO-A positions independently agree with all three COR2 headers
  within 345 km. Horizons is used as a cross-check for COR2.
  [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html),
  [HCI frame definition](https://docs.sunpy.org/en/stable/generated/api/sunpy.coordinates.frames.HeliocentricInertial.html).
- `observed_utc` is the instrument's mean receive time. The physical front
  time is approximated by subtracting `norm(observer_hci_m) / c` from
  the receive time. `elapsed_s` is that corrected time minus
  `2020-11-29T13:12:00Z`. This retains the approximately 8.8 s difference
  between SOHO and STEREO-A light travel times. Tangent depth is unresolved
  from a single image ray. The manifest budgets 40 s of depth-related
  timing error, equivalent to a line-of-sight scale of 17.2 solar radii,
  for this short inner-corona sequence. The packer combines the marking
  angle and `(40 s × 3000 km/s) / norm(observer_hci_m)` in quadrature
  as the effective `sigma_rad`. The 3000 km/s envelope is above the
  published apex speed for this interval. These weights represent a
  conservative marking and timing budget with correlated errors within
  an exposure.
  [Event kinematics, Kouloumvakos et al. 2022, Fig. 5](https://doi.org/10.1051/0004-6361/202142515).
- `available_utc=2026-11-24T12:00:00Z` is the retrospective release cutoff
  for this prepared annotation set. For an operational replay, use the
  actual publication, processing or annotation timestamp for each
  input asset. A fit seeds `initial_fit` with the geometry and speed
  scale of the observed region. This case uses one axis, a fixed width
  and a constant speed over the interval. The seed initializes the solver,
  which fits the speed alongside the surface geometry.

The [event schema](../EVENT_INPUT.md) defines every replay key. The packer
adds `frame_file` to keep each ray linked to its FITS exposure and
`holdout_tracks` to store the reserved measurements outside the fitted
`tracks` array.

## Reproduce the checked-in input

From the AMPS root, download the manifest files to a separate directory,
using the manifest alias as the filename. This fetches both difference
exposures and the two Horizons answers. A regenerated Horizons response can
have new descriptive text or a new generation date. The verification
compares the numeric knots against the frozen answer.

```bash
case_dir=srcSEP/test/cme_shock_evolution
obs="$case_dir/observed_2020_11_29"
raw=/tmp/amps-cme-2020-11-29
mkdir -p "$raw"
/usr/bin/python3 - "$obs/sources.json" <<'PY' > "$raw/download.tsv"
import json
import sys
from pathlib import Path
m = json.loads(Path(sys.argv[1]).read_text())
for item in [*m["frames"], *m["backgrounds"], *m["ephemerides"].values()]:
    print(item["file"], item["url"], sep="\t")
PY
while IFS=$'\t' read -r file url; do
    curl -fL "$url" -o "$raw/$file" || exit 1
done < "$raw/download.tsv"

/usr/bin/python3 "$case_dir/pack_observations.py" \
  --sources "$obs/sources.json" --picks "$obs/picks.csv" \
  --verify-archive-dir "$raw" --output "$raw/repacked-observations.json"
cmp "$obs/observations.json" "$raw/repacked-observations.json"

/usr/bin/python3 "$case_dir/preview_observations.py" \
  --sources "$obs/sources.json" --archive-dir "$raw" \
  --picks "$obs/picks.csv" --output-dir "$raw/preview"
```

The image preview requires matplotlib and SciPy. Packing requires NumPy.
SHA-256 for the frozen [observations.json](observations.json) is
`eeb8f6acfaa714066af1316a09a44f0f0d7f051db2a906c242fcc04052359303`.

## Offline script regression and references

The preview test uses a pair of 13-column by 9-row cutouts in each of C2,
C3 and COR2. Each pair uses the same detector window in the event
exposure and the preceding exposure listed above. Coordinates in the table
are zero-based columns and rows of the original FITS images.

| Event file | Preceding file | Window origin | Archived mark |
| --- | --- | --- | --- |
| `22793508.fts` | `22793507.fts` | (97, 452) | (103, 455) |
| `32643517.fts` | `32643516.fts` | (360, 540) | (365, 544) |
| `A135400.fts` | `A133900.fts` | (601, 957) | (605, 962) |

[preview-cutouts.json](preview-cutouts.json) stores the signed 16-bit
samples, original exposure, scale and zero-point cards, and upstream file
hashes. The test materializes small FITS rasters from these samples and
translates each selected mark to its cutout's origin. The preview uses the
cutout pixel axes. The packer still reads the full-image astrometry and all
marks from `sources.json` and `picks.csv`.

With `--save-panel-data`, the preview writes NPZ files containing the plotted
image and marks beside its PNGs. The test compares every pixel against
[preview-reference.json](preview-reference.json) at an absolute tolerance
of `2e-15` in fractional count rate. The reference uses exact rational
exposure normalization and a separate scalar median calculation with
half-sample reflection at the cutout boundary. The checker compares the
origin, extent, display limits and mark coordinates exactly. It also
decodes the PNGs to verify the mark rendered. This comparison is independent
of font rasterization and PNG encoding.

Run the preview check directly from the AMPS root with a fresh directory.

```bash
case_dir=srcSEP/test/cme_shock_evolution
run_dir=$(mktemp -d)
/usr/bin/python3 "$case_dir/preview_fixture.py" prepare "$run_dir/preview-input"
/usr/bin/python3 "$case_dir/preview_observations.py" \
  --sources "$run_dir/preview-input/sources.json" \
  --archive-dir "$run_dir/preview-input" --picks "$run_dir/preview-input/picks.csv" \
  --output-dir "$run_dir/preview" --save-panel-data
/usr/bin/python3 "$case_dir/preview_fixture.py" check "$run_dir/preview"
```

To reproduce both preview fixture and reference from the verified downloads,
generate candidates outside the checkout and compare them.

```bash
/usr/bin/python3 "$case_dir/preview_fixture.py" record "$raw" \
  --fixture "$raw/preview-cutouts.json" --reference "$raw/preview-reference.json"
cmp "$obs/preview-cutouts.json" "$raw/preview-cutouts.json"
cmp "$obs/preview-reference.json" "$raw/preview-reference.json"
```

The cutout SHA-256 is
`4beb425fecf1edf975f442328b4e81364923c3bcb78b11160454efb67d80248d`.
The preview reference SHA-256 is
`4e6d459f091bf3131db43aa13cdc51ffbcb70cc4a4e5b6c7e880c7d08f04cfca`.
The [full-input reproduction commands](#reproduce-the-checked-in-input)
also run the packer directly. The legacy target runs that CLI from the
checked-in astrometry and marks, compares the resulting JSON byte for byte
against its own reference and passes the generated JSON into the fitter.

## Use the same workflow for a different event

1. Locate the event in a published shock study or joint instrument
   inventory. Obtain science FITS on at least two well-separated lines
   of sight and cover three or more front epochs. Record a reproducible
   reason for feature choices, which frames are reserved, and any gaps.
2. Copy `sources.json` into the new event directory. Set `replay` fields,
   issue UTC and uncertainty budget. Replace the `frames`,
   `backgrounds` and `ephemerides` with the new archive products and
   Horizons queries. Each frame needs a local `file`, `instrument`
   (one of `STEREO-A/COR2`, `SOHO/LASCO-C2`, `SOHO/LASCO-C3`),
   download `url`, `role` (`fit` or `holdout`), `available_utc`,
   `difference_frame`, `preview` (`median_px`, `vmin`, `vmax`),
   and an `ephemeris` key for LASCO. Each background needs its local
   `file`, `instrument` and `url`. Each ephemeris entry needs its
   local `file` and query `url`. The recording command fills in
   `sha256`, `header` and `knots`. If using a different product or
   WCS, verify its astrometry and generate the rays with the instrument's
   calibrated transform.
3. Download these assets to a new raw directory. Run
   `preview_observations.py` without `--picks` for unannotated panels.
   Save the reviewed marks and their uncertainty to a new `picks.csv`
   and rerun the preview with the overlay. Keep all flags and missed
   frames in the event's provenance.
4. Run the packer on the new files. Use `--record-archive-dir` once to
   record the FITS headers, hashes and numeric ephemeris samples in the
   new manifest. Use `--verify-archive-dir` for subsequent reproductions.

   ```bash
   /usr/bin/python3 srcSEP/test/cme_shock_evolution/pack_observations.py \
     --sources /path/to/event/sources.json --picks /path/to/event/picks.csv \
     --record-archive-dir /path/to/event/raw --output /path/to/event/observations.json
   ```

5. Pass the new JSON path through `--data` and the appropriate UTC cutoff
   through `--issue` to `regression.py prepare`. For a real native
   five-hour SEP shock trace, you also need admissible geometry, a
   qualified dated upstream realization and data coverage or a physically
   justified continuation for the full interval. The
   [event-input guide](../EVENT_INPUT.md) gives the receipt and native
   CLI invocation.

## Regression result and limits

The legacy Table check runs offline on the checked-in package. It repacks
the observations from the source metadata and picks, checks the STEREO-A
frame against Horizons, and runs the ordinary fitter. The optimizer reaches
the supported upper half-width bound of `π/2 - 0.0001` rad with a trial RMS
image residual of about 489 arcseconds. The expected result is
`REJECTED` in `observed-2020-11-29.json` and a passing software check that
the invalid fit cannot become a frozen C++ input. This result constrains
this data selection and reduced model. The 90° value is the optimizer's
limit. The published reconstruction uses a more flexible
ellipsoid with evolving geometry. Any future positive fit for this event
needs independent review using the reserved exposures.
[Kouloumvakos et al. 2022, Sect. 2.2 and Fig. 4](https://doi.org/10.1051/0004-6361/202142515).

The existing five-hour native trace and final comparison exercise the
prescribed numerical benchmark in the same Table test. The archival case
has an empty `ambient_members` array and no outer propagation law. A
five-hour real-event comparison of normal speed, Mach number or compression
requires a qualified dated upstream and propagation model as well as
admissible front geometry.
