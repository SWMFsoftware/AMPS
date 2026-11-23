# 1997 May 27 SOHO/CDS and Ulysses connectivity screen

This standalone analysis evaluates whether the two CDS prominence rasters
sample low coronal material that could feed the wind measured by Ulysses.
It uses the historical observations and the existing AMPS magnetic field
solver. Stage 13 qualification requires a calibrated spectrum, plausible
emission attribution, and an independent in situ connection.

## Sources and epoch

The [retrieval manifest](retrieval.json) pins archive URLs and SHA-256 values.
Keep the source files outside the repository under their manifest names.

- SOHO/CDS `s7997r00.fits` and `s7997r01.fits` are Level 0 legacy NIS
  distribution copies. Preserve their telemetry null values, detector
  background windows, wavelength windows and exposure times. The raster
  starts on its rightmost column, so observation times decrease as the
  native column index increases. The supplied preliminary WCS is for geometry
  reconnaissance. [SOHO archive](https://soho.nascom.nasa.gov/data/catalogues/) · [SOHO product descriptions](https://umbra.nascom.nasa.gov/soho/sr08/soho_map.pdf)
- The NSO/Kitt Peak `m1923f.fits` map provides 360 Carrington longitudes
  increasing from zero in degrees and 180 equal-area latitude rows from
  south to north. The `f` product includes the provider's line-of-sight to
  radial correction and polar fill. Stored integers have scale 0.1 gauss.
  The pipeline treats the rows as equal-area bin centers in sine latitude,
  converts the field to tesla, and reports the signed monopole removed
  before the PFSS solve. This is a rotation-averaged magnetic boundary.
  Some contributing longitudes were observed after the raster date. The
  NSO instantaneous series for 1997 begins on November 7. [Map format and acknowledgement](https://nispdata.nso.edu/ftp/kpvt/synoptic/README) · [Instantaneous map index](https://nispdata.nso.edu/ftp/kpvt/synoptic/mag.instantaneous/1997/)
- SOHO's hourly ephemeris gives distance, heliographic latitude and
  rotating Carrington and inertial heliographic longitudes at the source.
  Ulysses' hourly merged measurements give distance, inertial coordinates,
  signed RTN magnetic field and proton flow speed and direction in the
  delayed June 14–29 arrival window. Use explicit native fills and preserve
  the sign of \(B_r\). The ±12-hour source-time envelope around the raster
  is a coverage choice. It is not a test that two parcels match. [SOHO ephemeris metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOHO_HELIO1HR_POSITION) · [Ulysses metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_COHO1HR_MERGED_MAG_PLASMA)
- Suess et al. place May 27 within I2, the variable, fanned-out streamer
  interval. Their LASCO and UVCS analysis supplies morphology and a
  \(2.5\,R_\odot\) source-surface reference. The independent Kitt Peak He I
  coronal-hole atlas is a rotation-average photospheric context map.
  Neither supplies the line-of-sight emitting location of these CDS
  rasters. [Suess et al., Sections 2 and 3](https://ntrs.nasa.gov/citations/20000032969) · [Kitt Peak He I atlas](https://nispdata.nso.edu/ftp/kpvt/synoptic/choles/ch_rot1923.pdf)

## Reproduce the geometry screen

From the AMPS root, supply a directory containing the eight named files
from `retrieval.json`. NumPy and a C++17/OpenMP compiler are required.
The geometry program builds the small adapter against the canonical AMPS
PFSS evaluator in a temporary directory.

```sh
OMP_NUM_THREADS=8 /usr/bin/python3 \
  srcSEP/validation/cases/stage13/soho-ulysses-1997/reproduce.py \
  --input-dir "$CDS1997_INPUT" --output-dir "$CDS1997_OUTPUT"
PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=2 /usr/bin/python3 \
  srcSEP/validation/cases/stage13/soho-ulysses-1997/test_mapping.py -v
```

Each off-limb ray that intersects \(1.05\,R_\odot\) has a near and far
intersection. Their Carrington Cartesian positions use the measured
distance and heliographic latitude of SOHO at the individual exposure
midpoint. The tracer integrates both field directions with fourth-order
Runge–Kutta and admits an open tube only when one end reaches the
photosphere and the other reaches \(2.5\,R_\odot\). Its spherical step,
length cap and boundary-distance checks are diagnostic numerical
controls. Two step sizes at \(\ell_{\max}=12\) and a second truncation
\(\ell_{\max}=8\) expose changes in topology and source locations. These
checks sample one magnetogram and one source-surface radius.

The outer map uses the measured radial flow component, travel time from
the source surface and the Carrington coordinate phase at departure.
This constant-speed ballistic reconstruction omits the unknown travel
time below the source surface and any outer stream interaction. The output
lists all candidate arrival rows and both signs. For open rays it records
the nearest source-surface candidate of matching sign, without interpreting
a small angular distance as a confirmed connection. The actual columns and
units appear in the CSV headers. `result.json` records the input hashes,
geometry runs, coverage and qualification status.

## Radiometric and pointing calibration

`calibrate_raster.pro` wraps SolarSoft's `READCDSFITS`, `GT_SOLAR_XY`,
`CDS_NEW_SPIKE` and `NIS_CALIB` routines. Run it in a CDS-enabled
SolarSoft IDL session with versioned pointing, wavelength, detector and
radiometry databases. Its `/ADJUST` setting applies SolarSoft's
helium-limb to photospheric-limb pointing correction. The exported
corrected coordinates are specific to the Mg IX window. SolarSoft
also retains the original detector array and window descriptors in the
saved product. Radiometry is in erg cm\(^{-2}\) s\(^{-1}\) sr\(^{-1}\)
per spectral pixel. A subsequent science fit needs compatible coordinates
for each line, detector and wavelength uncertainties, a foreground
treatment and the covariance of fitted intensities. No temperature,
pressure or composition enters the Stage 13 validation cases from this
screen. [CDS instrument calibration note](https://soho.nascom.nasa.gov/softops/cds/lrg_data/info/swnote/vds_calibration.pdf) · [CDS analysis guide](https://solar.bnsc.rl.ac.uk/software/uguide/uguide.shtml)

```idl
.compile calibrate_raster.pro
calibrate_raster, '/archive/cds1997/s7997r00.fits', '/calibrated/cds1997/s7997r00'
calibrate_raster, '/archive/cds1997/s7997r01.fits', '/calibrated/cds1997/s7997r01'
```

For a corrected geometry replay, pass `--calibration-dir` pointing at the
directory with both `.sav` and both `.wcs.csv` files. Generate the receipt
with `/usr/bin/python3` and SciPy. The inspector opens both saved arrays,
checks the detector-to-radiance units and validates the exported Mg IX
coordinates against the full wavelength-dependent WCS. It records the
input and output hashes and the installed database inventory. Retain the
SolarSoft logs with the receipt.

```sh
PYTHONDONTWRITEBYTECODE=1 /usr/bin/python3 \
  srcSEP/validation/cases/stage13/soho-ulysses-1997/inspect_calibration.py \
  --input-dir "$CDS1997_INPUT" --calibration-dir "$CDS1997_CALIBRATED" \
  --cds-data-dir "$SOLARSOFT_ROOT/soho/cds/data" --software "$CDS_SOFTWARE_VERSION"
OMP_NUM_THREADS=8 /usr/bin/python3 \
  srcSEP/validation/cases/stage13/soho-ulysses-1997/reproduce.py \
  --input-dir "$CDS1997_INPUT" --calibration-dir "$CDS1997_CALIBRATED" \
  --output-dir "$CDS1997_CORRECTED_OUTPUT"
```

## Results on the two archived rasters

Both rasters passed the SolarSoft preprocessing and archive-array checks
using the LMSAL CDS, SOHO and general software bundles with GNU Data
Language 0.9.9. The inspector verified all 8,580 corrected spatial positions
and radiance units in each science window. The fitted-line and foreground
requirements still apply to any temperature or pressure estimate.

| Position input | Off-limb pixels, raster 00 | Off-limb pixels, raster 01 | Open / closed / unresolved shell intersections |
| --- | ---: | ---: | ---: |
| Preliminary archive WCS | 1,707 | 1,705 | 0 / 6,824 / 0 |
| SolarSoft Mg IX WCS | 1,687 | 1,685 | 0 / 6,744 / 0 |

The topology counts agree across all three PFSS runs for each position
input. The Ulysses candidates reaching the spacecraft June 18–20 UTC
include both radial-field polarities. These rasters supply no open-stream
sample at \(1.05\,R_\odot\) under the dated synoptic reconstruction and
the chosen source surface. A ray can contain emission from several radii
and structures, so the result does not locate the actual emitter. The
Stage 13 wind qualification remains SKIPPED. The raw calibration logs,
database inventory and geometry receipts are retained in
`/data/vtenishe/TEST/stage13-cds1997.D2PjlV`.
[SolarSoft CDS bundle](https://www.lmsal.com/solarsoft/ssw/offline/swmaint/tar/ssw_soho_cds.tar.Z) · [SolarSoft SOHO bundle](https://www.lmsal.com/solarsoft/ssw/offline/swmaint/tar/ssw_soho_gen.tar.Z) · [SolarSoft general bundle](https://www.lmsal.com/solarsoft/ssw/offline/swmaint/tar/ssw_ssw_gen.tar.Z)
