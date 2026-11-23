# 2020 December 14: eclipse and K-Cor low-corona validation package

`case.json` registers authentic reference availability and missing inputs for
this proposed comparison. It is not an observational PASS. No density, shock
surface, uncertainty, model forcing or scoring tolerance has been invented.
The central case selector must retain `NOT_RUN`/incomplete scientific status
until the listed forward operator and observations are supplied.

The small official numerical reference is
`reference/20201214.kcor.radial.intensity.txt`: 53 azimuthally averaged polarized
brightness samples at heliocentric radii 1.07–2.11 solar radii. Its two brightness
columns are in B/Bsun. The second is variation over solar azimuth, **not** the
measurement uncertainty of the mean. The original header and bytes are
preserved, and the file has no exact exposure epoch. It cannot replace a local
streamer/coronal-hole profile or be treated directly as electron density.

The [official date selector](https://mlso.hao.ucar.edu/mlso_data_summary.php?date=2020-12-14&inst=kcor)
links the [numerical pB scan](https://mlso.hao.ucar.edu/hao/acos/engineering/2020/12/14/20201214.kcor.radial.intensity.txt)
and its [instrument-team description](https://mlso.hao.ucar.edu/files/web_page_pB_radial_plot_info.pdf).
`reference/provenance.json` records their meaning and the actual acquisition.
`reference/observing_windows.csv` and `reported_pos_speed_ranges.csv` retain
published metadata only; the approximate projected speeds are ineligible for
quantitative density or three-dimensional shock scoring.

The two eclipse views discussed by [Boe et al. (2021)](https://arxiv.org/html/2106.04027v1)
precede the later K-Cor interval. The [HAO event description](https://www2.hao.ucar.edu/mlso/event-gallery/2020-12-14-total-solar-eclipse-day-cmes)
reports several CMEs and weather gaps. No component association or simultaneous
cross-calibration is assumed. The two-site images show multiple evolving
components and nonradial motion; one spherical front is not thereby qualified.

`reference/electron_density_metadata.json` separately records the CIP inversion
study by [Edwards et al. (2023)](https://doi.org/10.1007/s11207-023-02231-5):
its observing site, directional regions, limited radial range, later K-Cor
cross-calibration and data-availability statement. The actual density arrays
and their calibration/error information have not been obtained. The paper's
spherical inversion assumption and sensitivity to polarizer-angle selection
must remain explicit. Neither is a numerical acceptance criterion chosen here.

Reproduce only the checksum-pinned small product into a fresh directory:

```sh
python3 -B reference/fetch_reference.py --output-dir /absolute/fresh/evidence
```

Restricted execution environments may block DNS inside Python children even
when a direct top-level `curl` works. In that case acquire the exact URL above
with `curl -fsSL --max-filesize 2000000 -o /cache/20201214.kcor.radial.intensity.txt`,
then run `fetch_reference.py --source-cache /cache --output-dir /fresh/evidence`.
The receipt labels cache replay separately, verifies size/SHA256, and never
rewrites the committed reference. Acquisition success qualifies provenance
only. Scientific comparison still needs exposure metadata, calibrated local
images/density profiles with errors, geometry/line-of-sight treatment and the
reviewed forward operator identified in `case.json`.

K-Cor dataset DOI: `10.5065/D69G5JV8`. Credit the Mauna Loa Solar Observatory,
High Altitude Observatory and NSF NCAR when using these data. The
[MLSO acknowledgement instructions](https://www2.hao.ucar.edu/mlso/data-use-acknowledgements)
and [UCAR data terms](https://www.ucar.edu/terms-of-use/data) describe attribution
and the applicable CC BY4.0 terms. The original numerical table is unchanged;
manual transcriptions and incomplete metadata are labeled separately.

A single per-exposure calibrated L2 product was subsequently verified at
[the official exposure URL](https://mlso.hao.ucar.edu/hao/acos/2020/12/14/20201214_175641_kcor_l2_pb.fts.gz).
`reference/kcor_level2_exposure.json` and the original primary-header cards
record17:56:41–17:56:56UTC, calibrated pB in B/Bsun, helioprojective WCS,
solar-north orientation, bandpass and tangential-polarization pipeline.
The listing begins at17:56:41UTC; no eclipse-time exposure is claimed.
The raw2.56MB compressed image stays outside the source package. Inspect a
fresh independently acquired copy with:

```sh
curl -fsSL --max-time 60 --max-filesize 6000000 \
  'https://mlso.hao.ucar.edu/hao/acos/2020/12/14/20201214_175641_kcor_l2_pb.fts.gz' \
  -o /absolute/evidence/20201214_175641_kcor_l2_pb.fts.gz
python3 -B reference/inspect_level2.py \
  --input /absolute/evidence/20201214_175641_kcor_l2_pb.fts.gz \
  --output-dir /absolute/fresh/header-check
```

The current archive product was reprocessed in2025 with software2.2.6;
identity with the2023 density study's calibration is not established. Its
`QUALITY=ok` keyword does not establish a quiet baseline, and no error image
extension is present. Header integrity cannot supply the missing directional
density inversion or its uncertainty and cannot earn observational PASS.
