# Stage 13 data sources and provenance

This inventory covers the seven maintained spherical SWCME3D `SHOCK_ONLY`
campaign inputs, the independent wind constraints used to design the thermal
solver, and candidate observations for qualifying the wind. The physical
particle domain is 1.05 solar radii to 5 AU. A wind observation or circulation
control beyond 5 AU may constrain provider support separately from the 5 AU
particle escape surface. The build, case-selection commands and validation
rules are in the [Stage 13 campaign guide](../../STAGE13_SPHERICAL_VALIDATION.md).

## Pinned Stage 13 case inputs

Each case `case.json` records the source IDs, exact product URLs, artifact
paths, byte counts, checksums, and case readiness. Its README describes data
processing, units, quality information and coverage in more detail.
`retrieval.json` and the case downloader or reproducer record exact archive
queries, time subsetting and retained source bytes where applicable. Keep
observed data, derived observables, model-informed forcing, and provenance
metadata distinct when using the manifests.

### 2020 December 14 eclipse and K-Cor

- The Mauna Loa Solar Observatory daily radial profile contains 53
  azimuthally averaged polarized-brightness entries over approximately
  1.07–2.11 solar radii. Use its brightness units and processing description
  when synthesizing line-of-sight Thomson scattering. Azimuthal spread
  describes structure around the limb and must be interpreted with the data
  provider's definition before use in a measurement-error model. [MLSO date-specific index](https://mlso.hao.ucar.edu/mlso_data_summary.php?date=2020-12-14&inst=kcor) · [radial profile](https://mlso.hao.ucar.edu/hao/acos/engineering/2020/12/14/20201214.kcor.radial.intensity.txt) · [radial-profile documentation](https://mlso.hao.ucar.edu/files/web_page_pB_radial_plot_info.pdf)
- One calibrated K-Cor Level 2 FITS exposure spans 17:56:41–17:56:56 UTC.
  The case retains the header, WCS, processing and reprocessing provenance,
  units and checksum. Compare at the recorded exposure, separately from
  the earlier eclipse images. K-Cor is cited as dataset DOI `10.5065/D69G5JV8`.
  [Exact Level 2 exposure](https://mlso.hao.ucar.edu/hao/acos/2020/12/14/20201214_175641_kcor_l2_pb.fts.gz) · [MLSO acknowledgements](https://www2.hao.ucar.edu/mlso/data-use-acknowledgements) · [UCAR terms](https://www.ucar.edu/terms-of-use/data)
- Eclipse event timings and reported plane-of-sky speeds come from the
  original papers and instrument event report. The Edwards et al. density
  estimate comes from polarization inversion. The line of sight can contain
  both open and closed material. The retained products require a forward
  brightness operator and a separate temperature or pressure constraint.
  [Boe et al.](https://arxiv.org/html/2106.04027v1) · [Edwards et al.](https://link.springer.com/article/10.1007/s11207-023-02231-5) · [Abramson et al.](https://www.papersinphysics.org/papersinphysics/article/view/718) · [MLSO event report](https://www2.hao.ucar.edu/mlso/event-gallery/2020-12-14-total-solar-eclipse-day-cmes) · [case and manifest](eclipse-2020-12-14/README.md)

### 2012 May 17 ground-level enhancement

- The official CCMC GOES-13 JSON contains processed integral proton
  event quantities for five channels. The two retained official comparison
  rows preserve their (>10) MeV/10 pfu and (>100) MeV/1 pfu definitions.
  The JSON names the underlying time series separately. The Oulu GLE and
  NMDB databases are routes for station counts, which also require response,
  station corrections and geometry. PAMELA is an independent energetic
  particle reference. [CCMC event](https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/20120517.php) · [GOES-13 JSON](https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/docs/json/GOES-13_integral.2012-05-17T000000Z.json) · [Oulu](https://gle.oulu.fi/) · [NMDB](https://www.nmdb.eu/nest/) · [PAMELA analysis](https://ntrs.nasa.gov/api/citations/20180005611/downloads/20180005611.pdf)
- The 1.38 and 2.32 solar-radii constraints are published heliocentric
  proxy heights associated with the corrected metric type-II onset and
  inferred release, respectively. Source and electromagnetic comparison
  clocks and the COR1 exposure time remain separate. The alternative
  multi-species timing analysis carries its own path-length model and
  uncertainty. [Gopalswamy et al.](https://cdaw.gsfc.nasa.gov/publications/gopal/gopal2013ApJL.pdf) · [Rouillard et al.](https://arxiv.org/pdf/1605.05208v2) · [case and manifest](gle-2012-05-17/README.md)

### 2013 April 11 near-Earth continuity case

- OV01 uses 80 points digitized from the three spectra in Figure 12 of
  Liu et al. The instruments are ACE/EPAM, GOES-13/EPEAD and SOHO/ERNE.
  The observing clock starts at 07:24 UTC. The model-informed Earth-connected
  thermal-energy source proxy starts at 06:00 UTC. Keep the original XM03
  common source and global normalization for the continuity comparison.
  [Liu et al.](https://doi.org/10.3847/1538-4357/adc4e3) · [publisher provenance and license](https://authors.library.caltech.edu/records/qwsnd-jmx38) · [CCMC event](https://ccmc.gsfc.nasa.gov/assessment/topics/SEP/campaign2020/20130411.php) · [case and manifest](sep-2013-04-11/README.md)

### 2001 September 24 returning particles

- The file is the original CDAWeb HAPI `SOHO_ERNE-HED_L2-1MIN` CSV subset,
  from 10:00 through before 16:00 UTC, with `Time`, `est`, `PH` and `PHC`.
  Preserve energy bounds, fill values and channel definitions from HAPI.
  The one-minute intensity and count product contributes event context.
  The directional return comparison also needs sector counts or rates,
  response matrices and the field or pitch-angle mapping used by the paper.
  [Exact HAPI request](https://cdaweb.gsfc.nasa.gov/hapi/data?id=SOHO_ERNE-HED_L2-1MIN&parameters=Time,est,PH,PHC&time.min=2001-09-24T10:00:00Z&time.max=2001-09-24T16:00:00Z&format=csv) · [HAPI metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOHO_ERNE-HED_L2-1MIN) · [Tan et al.](https://epact2.gsfc.nasa.gov/don/09TanApJ.pdf) · [SOHO event](https://soho.nascom.nasa.gov/hotshots/2001_09_24/) · [case and manifest](sep-2001-09-24/README.md)

### 1978 September 23 and 1979 March 1 multi-spacecraft SEP events

- Reames, Ng and Tylka provide the spacecraft arrangement and published
  event constraints. The retained Voyager 2 CRS daily and six-hour excerpts
  audit real data availability. The selected event intervals have gaps in
  these products. A channel-matched independent comparison requires the
  calibrated flux products and reported response for the actual observer.
  [Reames, Ng and Tylka](https://epact2.gsfc.nasa.gov/don/12Helios.pdf) · [Voyager daily product guide](https://voyager.gsfc.nasa.gov/yearly/voyager-2/flux/P1D) · [1978 daily](https://voyager.gsfc.nasa.gov/yearly/voyager-2/flux/P1D/1978P1Y.asc) · [1978 six-hour](https://voyager.gsfc.nasa.gov/yearly/voyager-2/flux/PT6H/1978P1Y.asc) · [1979 daily](https://voyager.gsfc.nasa.gov/yearly/voyager-2/flux/P1D/1979P1Y.asc) · [1979 six-hour](https://voyager.gsfc.nasa.gov/yearly/voyager-2/flux/PT6H/1979P1Y.asc) · [1978 case](sep-1978-09-23/README.md) · [1979 case](sep-1979-03-01/README.md)

### Ulysses October–November 2003

- NASA SPDF COSPIN/HET daily PHA flux and uncertainty contribute 61
  original daily records for the larger-domain control. McKibben et al.
  describe this multiple-event period. The 8 AU comparison control in the
  deck and its observational coverage are distinct from the 5 AU physical
  boundary. [PHA record](https://spdf.gsfc.nasa.gov/pub/data/ulysses/particle/cospin/het/het_pha_flux_h_he_1day/cospin_het_phaflux_h_he_1day_1990_20090630.asc) · [PHA format](https://spdf.gsfc.nasa.gov/pub/data/ulysses/particle/cospin/het/het_pha_flux_h_he_1day/cospin_het_phaflux_format.txt) · [NASA dataset metadata](https://data.nasa.gov/dataset/ulysses-cosmic-ray-and-solar-particle-investigation-cospin-high-energy-telescope-het-proto) · [McKibben et al.](https://doi.org/10.1029/2005JA011049)
- CDAWeb `UY_COHO1HR_MERGED_MAG_PLASMA` supplies 1,464 hourly records
  from 2003-10-01 through before 2003-12-01 UTC. It supplies position, RTN
  field, flow, proton and alpha density, and the two legacy ion-temperature
  alternatives. `UY_PROTON-MOMENTS_SWOOPS` supplies 401 records from
  2003-10-28 through before 2003-10-30 UTC, with RTN velocity and field,
  proton number density and parallel/perpendicular proton moments. Values
  retain native nT, km/s, cm(^{-3}), K and AU. No interpolation, fill
  replacement, scalar-temperature assumption or coordinate conversion was
  applied to these files. [COHO HAPI metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_COHO1HR_MERGED_MAG_PLASMA) · [SWOOPS proton metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_PROTON-MOMENTS_SWOOPS) · [columns, cadence and source URLs](ulysses-2003-oct-nov/reference/thermal_context_provenance.json)
- The SPDF `UY_M0_BAE` electron-moments 2003 directory ends in May for
  the retained product. An October query produced an out-of-window
  May sample, which was rejected. The electron guide also describes
  dependence on interpolated ion densities in the core-density calculation.
  These contextual data leave dated mapping, stream/ICME masks, electron
  pressure, and the wind mass-loading calibration open. [SPDF electron archive](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swoops_cdaweb/bae_m0/2003/) · [coverage audit](ulysses-2003-oct-nov/reference/electron_coverage_audit.json) · [case and manifest](ulysses-2003-oct-nov/README.md)

## Candidate wind data from the archive search

The two SOHO–Ulysses quadratures are direct candidates for tying coronal
observations to wind observations near and beyond 5 AU. Keep the documented
coronal exposure dates separate from the delayed spacecraft arrival times.
Both campaigns require the actual pointing and field mapping, calibration
and independent lower-corona thermodynamics before they can supply a
qualified Stage 13 input.

- **1997 May–June SOHO–Ulysses low-latitude streamer wind.** The May 21 to
  June 1 coronal coverage includes UVCS O VI and H I measurements around
  3.5 and 4.5 solar radii and LASCO imagery. Figure 6 covers Ulysses DOY
  165–180 near 5.1 AU, about 10 degrees north of the solar equator. The
  wind speed over the campaign was approximately 375 ± 25 km/s. It is
  a promising quiet slow-wind control after excluding the published
  transient windows. The contemporaneous SUMER catalogue lists a May 24
  off-limb Lyα program, May 25 and 27 off-limb Lyβ/O VI programs and a
  May 29 Mg IX program. The pointing check below places those particular
  scans away from the east-limb stream.
  [Suess et al. 2000](https://doi.org/10.1029/2000JA000044) · [NASA full text](https://ntrs.nasa.gov/citations/20000032969) · [SUMER chronological catalogue](https://www2.mps.mpg.de/projects/soho/sumer/text/chronolist.html)
- **1998 November–December SOHO–Ulysses quadrature.** The SOHO observing
  program is #5490, with UVCS H I Lyα and O VI at 3.5 and 4.5 solar radii
  from November 29 through December 7 and LASCO C2 polarized brightness
  for density context. Ulysses was at approximately 5.2 AU, 17.4 degrees
  south, with SWOOPS proton and alpha measurements, VHM/FGM magnetic field,
  and SWICS charge states. The paper distinguishes a rapidly expanding
  low-latitude coronal-hole stream and a slower streamer-associated stream
  around a transient and sector transition. Its coronal proton speeds
  use the Ulysses mass-flux mapping, so those speeds and the same mass
  flux cannot be counted as independent constraints. Its reported oxygen
  line widths, charge-state freeze-in temperatures and assumed/set electron
  temperatures cannot supply a measured proton and electron pressure at
  1.05 solar radii. Preserve separate selection and mapping for the two
  streams. [Poletto et al. 2002, sections 3–6 and Tables 1–3](https://doi.org/10.1029/2001JA000275) · [SOHO schedule](https://soho.nascom.nasa.gov/soc/calendar1998_old.html) · [SOHO science archive](https://soho.nascom.nasa.gov/data/archive/) · [CDAWeb Ulysses product catalogue](https://cdaweb.gsfc.nasa.gov/misc/NotesU.html)
- **Confirmed outer archive products.** The in-period HAPI queries use
  `UY_M0_BAE` for total/core/halo electron density and scalar temperature,
  `UY_PROTON-MOMENTS_SWOOPS` for proton RTN velocity, density and temperature
  tensor, `UY_COHO1HR_MERGED_MAG_PLASMA` for hourly position and proton/alpha
  density context, and `UY_1MIN_VHM` for one-minute RTN magnetic vectors.
  The 1997 June 14–29 UTC window contains 579 electron, 1,955 proton,
  384 hourly and 22,808 magnetic records. The 1998 December 19–27 UTC
  window contains 337 electron, 1,076 proton, 216 hourly and 12,794 magnetic
  records. These are archived rows before scoring or a quiet-stream mask.
  The largest gaps between electron timestamps are about 43.5 and 8.3 hours,
  respectively. The 1997 paper flags a small shock around June 17 and
  possible CME material. The 1998 paper flags a CIR at the beginning of
  the comparison interval. The archive electron-core density uses ion
  density with an assumption of neutrality, so it is correlated with the
  ion-based mass flux. The SPDF SWICS 3-hour ratio products also supply
  128 rows in 1997 and 72 rows in 1998, with archived alpha speed,
  C(^{6+})/C(^{5+}), O(^{7+})/O(^{6+}), Fe/O, mean Fe charge and
  accumulated-cycle count. One and two O(^{7+})/O(^{6+}) entries are
  missing in the two respective intervals. Preserve
  the original zero-missing convention and the documented charge-state
  ratio uncertainty. Retain each raw time base and missing flag.
  [Electron metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_M0_BAE) · [Proton metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_PROTON-MOMENTS_SWOOPS) · [Hourly merged metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_COHO1HR_MERGED_MAG_PLASMA) · [One-minute magnetometer metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_1MIN_VHM) · [1997 SWICS ratios](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swics/3hr_ratios/uswichst97.asc) · [1998 SWICS ratios](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swics/3hr_ratios/uswichst98.asc)
- **Historical magnetic boundary and spectral calibration.** The 1998
  coronal observing interval falls in Carrington rotation 1943. Use a
  dated, versioned photospheric map from the NSO historical archive or an
  equivalent supported synoptic product. Preserve the photospheric field
  component, harmonic normalization, map window and source-surface
  geometry. The VSO lists SOHO UVCS data. The SOHO archive distinguishes
  current calibrated products from historical files with obsolete FITS
  headers. Keep native times, spacecraft attitude, slit orientation,
  radiometry and measured line intensities when building the comparison.
  Catalogue entries establish a retrieval route. Usable pointed spectra
  require verification of the exposure metadata and spatial coverage.
  [Carrington date catalogue](https://data.darts.isas.jaxa.jp/pub/yohkoh/YLA4ISAS/sxt_synop_maps.html) · [NSO historical maps](https://nso.edu/data/historical-archive/) · [GONG synoptic map catalogue](https://magmap.nso.edu/archive.html) · [VSO UVCS catalogue](https://vso.nascom.nasa.gov/cgi/show_details?instrument=UVCS) · [SOHO calibrated product descriptions](https://umbra.nascom.nasa.gov/soho/sr08/soho_map.pdf)

The archive audit verified all eight HAPI time subsets and 19 artifact
hashes. Exact queries, retained raw data, fill counts and SWICS subsets are
recorded in `audit.json` and `README.md` under
`/data/vtenishe/TEST/stage13-quadrature-observation-audit-20261123`.

## Independent inner-corona and inner-heliosphere controls

- **Hinode low corona.** Hahn et al. used EIS in a 2007 solar-minimum
  polar coronal hole and report (n_e) and (T_e) from 1.03 to 1.20 solar
  radii. Their line widths include ion and nonthermal contributions. These
  profiles constrain that polar wind population. The 1997 and 1998 cases
  have different epochs and coronal topology. Cranmer's UVCS archive
  supplies independent high-latitude H I and thermal estimates farther
  out and has uncertainty and machine-readable data. [Hahn et al.](https://arxiv.org/abs/1010.4277) · [Cranmer 2020](https://arxiv.org/abs/2007.13180) · [UVCS numerical data](https://zenodo.org/records/3908519)
- **2021 January Metis and Parker Solar Probe.** Metis observed on
  January 17 at 16:30 UTC over 3.5–6.3 solar radii using polarized
  brightness and H I Lyα. Parker sampled the southern coronal-hole
  extension on January 18, 17:49–20:12 UTC at about 0.11 AU. January 23
  is an opposite-polarity source at about 0.26 AU. The inferred coronal
  outflow uses assumptions about (T_e), helium, and inward-scaled PSP
  proton properties. The preliminary CDAWeb `PSP_FLD_L3_SQTN_RFS_V1V2`
  check found 2,432 rows in the January 18 window. Its electron-core
  temperature and density need quality and population reconciliation
  with the corresponding SWEAP measurements before combining them in
  a pressure EOS. [Telloni et al.](https://doi.org/10.3847/1538-4357/ac8103) · [ESA Solar Orbiter Archive](https://soar.esac.esa.int/soar/) · [NASA PSP FIELDS/SWEAP archive](https://cdaweb.gsfc.nasa.gov/pub/data/psp/00readme.html) · [PSP assessment](ulysses-2003-oct-nov/THERMODYNAMIC_QUALIFICATION.md)
- **2022 March Solar Orbiter Slow Wind SOOP.** The March 3–6 and
  17–22 observing windows coordinate Solar Orbiter, Hinode/EIS and IRIS.
  The March 18–21 active region 12967 case provides tracked
  near-source spectroscopy and solar wind at about 0.35 AU. SWA-PAS/HIS
  and EAS products can constrain composition and separate thermodynamics.
  The March 2–7 electron temperature paper uses a different observing
  window. Treat changing source connectivity and ICME intervals explicitly.
  [Yardley et al. observing plan](https://doi.org/10.3847/1538-4365/acd24b) · [March 18–21 S-web study](https://arxiv.org/abs/2303.12192) · [March 2–7 electron study](https://doi.org/10.3847/1538-4357/ae3c7b) · [ESA SWA datasets](https://esdcdoi.esac.esa.int/doi/html/data/heliophysics/SolarOrbiter/SWA.html)
- **2023 March near-sonic slow wind.** PSP observed March 16 11:54:58 to
  March 17 05:57:20 UTC from 22 to 15.6 solar radii. Wind and MMS provide
  a candidate linked interval near 1 AU. This is a direct sonic-region
  constraint for an active-region outflow. Its photospheric connectivity
  and the wind observed near 1 AU have mapping and evolution uncertainty.
  [Ervin et al.](https://doi.org/10.3847/1538-4357/ad57c4) · [Wind measurements](https://wind.nasa.gov/data_sources.php)
- **Population comparisons.** Dakeyo et al. derive separate electron and
  proton iso-poly trends from PSP and Helios between roughly 0.1 and 1 AU.
  The inferred basal temperatures are extrapolations, and their supported
  thermal-only regimes differ by wind population. The Helios–Ulysses
  thermal study extending to 5.4 AU analyzes fast wind. Its slopes cannot
  establish the outer slow-wind law without a separate population match.
  [Dakeyo et al.](https://arxiv.org/abs/2207.03898) · [Cranmer et al. 2009](https://arxiv.org/abs/0907.2650)

## Follow-up archive audit, 2026-11-23

The best historical connection reaches Ulysses at 5.1 AU. The best modern
low-corona energy-budget candidate follows a different stream and reaches
Solar Orbiter at 0.335 AU. Their plasma parameters belong to different
epochs and sources. Each comparison keeps its own field mapping and
measurement uncertainties.

### SOHO–Ulysses quadratures

- The MSSL CDS catalogue contains 28,892 entries for 1997 and 13,779 for
  1998. A preliminary footprint screen retained 41 of 439 entries in
  1997 May 21–June 1 UTC and 34 of 787 in 1998 November 29–December 7 UTC,
  inclusive. The screen intersects unrotated footprints with helioprojective
  boxes in arcseconds. In 1997 the box is `x = [-1250, -800]`,
  `y = [-150, 500]`. In 1998 it is `x = [800, 1250]`, `y = [-650, 50]`.
  The 1998 hits are all from campaign 5485, with 27 catalogued for the
  gravitational focusing cone and seven for the full-Sun mosaic. Those
  windows still need pointing, line and plasma-connectivity analysis before
  use as a thermal constraint for the separate UVCS quadrature campaign 5490.
  [MSSL CDS catalogue](https://www.mssl.ucl.ac.uk/www_solar/gis/data_archive.html) · [1997 catalogue bytes](https://www.mssl.ucl.ac.uk/www_solar/gis/cat97.save.gz) · [1998 catalogue bytes](https://www.mssl.ucl.ac.uk/www_solar/gis/cat98.save.gz)
- Two *retrieved* 1997 May 27 NIS rasters are the closest spectral leads.
  Their FITS headers confirm the preliminary center `(-980.6, +113.0)`
  arcsec and size `121.9 × 240.2` arcsec. `s7997r00.fits` spans
  `08:08:39.132–09:02:40.096 UTC` and `s7997r01.fits` spans
  `09:02:50.785–09:56:52.810 UTC`. Both are CDS study 44, campaign 1500,
  `POBS1_2/v18`, a prominence program. The line windows include Mg IX
  368.07 Å, Mg X 624.94 Å, Si IX 345.13 Å and He I/II. The legacy window
  named `SI_12_320_66` has a FITS wavelength range near 320 Å. There is
  no recorded window at the 520.66 Å Si XII wavelength used in commonly
  published CDS line lists. The complete published Si IX 350/342 Å
  density diagnostic is also absent. The retained lines' value depends on
  finding calibrated off-limb pixels that trace to the open stream and
  quantifying prominence and foreground contamination. Retain calibrated
  line intensities and their covariance for an emissivity-based electron
  temperature fit. An ion line
  width needs its own thermal and nonthermal inversion, and an observed
  helium intensity needs an abundance and ionization model before it can
  constrain alpha density.
  [First FITS](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedCDS/cds/cdsfits11/s7997r00.fits) · [Second FITS](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedCDS/cds/cdsfits11/s7997r01.fits) · [CDS calibration guide](https://solar.bnsc.rl.ac.uk/software/uguide/uguide.shtml) · [CDS line list with Si XII](https://discovery.ucl.ac.uk/9749/1/9749.pdf) · [Gibson et al. diagnostic analysis](https://doi.org/10.1029/98JA02681)
- The specific May 24, 25, 27 and 29 SUMER scans cited above start at
  `(+639, +840)`, `(+639, +999)`, `(+639, +900)` and `(+600, +1020)`
  arcsec in the respective day catalogues, on the northern/western side.
  This excludes those pointings from the targeted low-latitude east-limb
  thermal comparison. The instrument chronology records closing on 1998
  November 12 and reopening on 1999 February 23, with no December observing
  entries. The catalogue does not establish an alternate SUMER spectrum
  for either target close to 1.05 solar radii.
  [Coordinate conventions](https://www2.mps.mpg.de/projects/soho/sumer/text/sff_cat1.html) · [May 24](https://www2.mps.mpg.de/projects/soho/sumer/text/an24May97.lst) · [May 25](https://www2.mps.mpg.de/projects/soho/sumer/text/an25May97.lst) · [May 27](https://www2.mps.mpg.de/projects/soho/sumer/text/an27May97.lst) · [May 29](https://www2.mps.mpg.de/projects/soho/sumer/text/an29May97.lst) · [November 1998 chronology](https://www2.mps.mpg.de/projects/soho/sumer/text/sfNov98.html) · [December 1998 chronology](https://www2.mps.mpg.de/projects/soho/sumer/text/sfDec98.html) · [February 1999 chronology](https://www2.mps.mpg.de/projects/soho/sumer/text/sfFeb99.html)

### 2024 March–April eclipse campaign

Rivera et al. track a small, negative-polarity coronal hole with remote
sensing from March 29 to April 9. Current-sheet crossings, electron strahl
and PFSS footpoints identify its wind at Parker Solar Probe on March 31
at 11.60 solar radii (0.054 AU) and at Solar Orbiter at 72.06 solar radii
(0.335 AU) in the April 7–14 context. Table 2 uses spacecraft Carrington
longitude ranges 85–121° and 70–105°, respectively. The event can test the
acceleration region with observed densities, plasma moments and magnetic
fields. The connection and expansion carry uncertainty from ballistic
mapping and the coronal field model. The archive copy is the manuscript
identified by the authors as accepted in *The Astrophysical Journal*.
[Rivera et al., Sections 2–4, Tables 1–3 and Appendix A](https://arxiv.org/abs/2610.04669)

- **Coronal density and temperature.** The reported electron-density
  scale near 1.05 solar radii is `(2.29 ± 0.71) × 10^8 cm^-3`. This is
  a remote electron diagnostic, even where the preprint uses it as an
  estimate of proton density. On-disk AIA/SPICE differential emission
  measure estimates give about `1.10–1.23 MK` for electrons, with a
  reported basal `T_e = 1.23 ± 0.04 MK`. Table 2 repeats that value for
  low-corona protons and alphas as a proxy. Independent basal ion
  temperatures and composition are required to obtain the ion contribution
  to pressure. The April 9 NSO Cryo-NIRSP query returned 25 public
  metadata records. Three sequential pairs covering the study region in
  the southwestern limb are `TFPUCV/NDTDBO` (17:39/17:42 UTC),
  `KYLOTD/HPRKOY` (17:54/17:57) and `SMTRNK/WRAEVA` (18:09/18:12),
  ordered as Fe XIII 1079.8/1074.7 nm. Their paired archive bounding
  boxes agree to roughly 0.1 arcsec at each pointing. Use the L1 spectra
  and exact WCS to build exposure-matched line ratios and carry along the
  line-of-sight contribution of nearby structures. Only the NSO metadata
  were retrieved in this audit. MLSO UCoMP line ratios and K-Cor
  polarized brightness provide separate coronal observables.
  [Exact NSO metadata request](https://api.dkistdc.nso.edu/datasets/v1?startTimeMin=2024-04-09T00%3A00%3A00Z&startTimeMax=2024-04-10T00%3A00%3A00Z&instrumentNames=CRYO-NIRSP&pageSize=100) · [MLSO archive](https://www2.hao.ucar.edu/mlso/mlso-data-and-movies) · [Rivera et al., Appendix A.2–A.3](https://arxiv.org/abs/2610.04669)
- **Wind speed, composition and mass loading.** The EIS Fe XII Doppler
  speed reported near the base is `8.09 ± 3.36 km/s`. This is a projected
  ion-line diagnostic. The March 29 `05:15:01 UTC` Metis polarized
  brightness and H I Lyα provide density and Doppler-dimming speed
  constraints at 2–4 solar radii. The ESA catalogue lists the corresponding
  Level 2 neighborhood as `solo_L2_metis-vl-pb_20240329T051456` and
  `solo_L2_metis-uv-image_20240329T051503`. The catalogue start times are
  `05:14:56.291` and `05:15:03.289 UTC`, respectively. The speed depends
  on the assumed neutral-hydrogen kinetic temperature (`1.2 MK`) and
  electron profile.
  The LASCO speed inversion uses the PSP flux and modeled expansion, so
  it shares the information used for the mass-loading anchor. At PSP,
  proton and alpha radial speeds are about `266 ± 28` and `378 ± 49 km/s`.
  The mass-loading observation should use their separate radial fluxes,
  `(m_p n_p U_p,r + m_alpha n_alpha U_alpha,r) / |B_r|`, with time, frame,
  fill, instrument quality and signed sector preserved. The preprint's
  4% helium assumption at the base and at Solar Orbiter does not establish
  composition at either location.
  [ESA Solar Orbiter Archive](https://soar.esac.esa.int/soar/) · [Rivera et al., Table 2 and Appendix A.7](https://arxiv.org/abs/2610.04669)
- **Transverse momentum and work.** UCoMP phase speeds and
  emissivity-weighted plane-of-sky magnetic field estimates constrain
  expansion and curvature with the density and projection uncertainties
  carried through the comparison. In Table 3 the inferred basal
  Alfvénic contribution is about `12 W/m²` and the alternative slow or
  fast compressive contributions are about `7` or `23 W/m²`, all scaled to
  the common reference area used in the paper. The fast-mode interpretation
  gives the closest agreement between the coronal and in-situ energy
  budgets. These fluxes bound the work ledger conditionally on mode and
  expansion. A transverse force requires the divergence of a stress
  tensor or a self-consistent field and pressure balance across neighboring
  tubes. The light-of-sight diagnostics and the scalar energy fluxes leave
  that force underdetermined. The dated Downs et al. MAS realization can
  supply a distinct model comparison for geometry and stress, with the
  observable synthesized using the instrument's viewing direction and
  emissivity. Treat its shared inputs and modeled expansion as correlated
  with the sections of the preprint that use MAS.
  [Rivera et al., Section 3 and Appendix A.4](https://arxiv.org/abs/2610.04669) · [Downs et al. model data and boundary inputs](https://zenodo.org/records/14889337) · [MAS download catalogue](https://www.predsci.com/corona/apr2024eclipse/download_data.php)

The CDAWeb HAPI coverage probes below preserve the instrument units and
timestamps. UTC ranges are inclusive at the start and exclusive at the
end. Counts precede quality filtering, synchronization, strahl and sector
selection, and application of the study's stream-mapping masks.

| Product and metadata | UTC date range | Rows returned | Fields requested |
| --- | --- | ---: | --- |
| [PSP SWEAP proton moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=PSP_SWP_SPI_SF00_L3_MOM) | March 29–April 2 | 125,178 | `QUALITY_FLAG,DENS,VEL_RTN_SUN,TEMP` |
| [PSP SWEAP alpha moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=PSP_SWP_SPI_SF0A_L3_MOM) | March 29–April 2 | 125,174 | `QUALITY_FLAG,DENS,VEL_RTN_SUN,TEMP` |
| [PSP FIELDS quasi-thermal noise](https://cdaweb.gsfc.nasa.gov/hapi/info?id=PSP_FLD_L3_SQTN_RFS_V1V2) | March 29–April 2 | 89,472 | `electron_density,density_quality_flag,electron_core_temperature` |
| [PSP FIELDS RTN](https://cdaweb.gsfc.nasa.gov/hapi/info?id=PSP_FLD_L2_MAG_RTN_1MIN) | March 29–April 2 | 5,760 | `psp_fld_l2_mag_RTN_1min` |
| [Solar Orbiter SWA-PAS moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOLO_L2_SWA-PAS-GRND-MOM) | April 8–15 | 265,643 | `quality_factor,N,V_RTN,T` |
| [Solar Orbiter MAG RTN](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOLO_L2_MAG-RTN-NORMAL-1-MINUTE) | April 8–15 | 10,080 | `B_RTN` |
| [Solar Orbiter SWA-EAS pitch-angle product](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOLO_L3_SWA-EAS-NMPAD-PSD) | April 7–15 | 63,616 | `QUALITY_FLAG,QUALITY_BITMASK` |

The CDAWeb SWA-HIS Level 3 composition probe returned no data in the
April 7–15 window. The ESA science catalogue query for overlapping SWA
products also returned zero `swa-his` entries. The CDAWeb RPW 10-second
density probe returned a fill-only record outside the window, which was
rejected. ESA catalogue hits include EAS velocity distributions and RPW
spacecraft-potential files in the same window. The HAPI EAS probe assessed
flags only. A derived electron pressure at Solar Orbiter would require a
reviewed electron-moment inversion using the distributions, spacecraft
potential and separation of core, halo, strahl and spacecraft populations.
[HIS catalogue coverage](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOLO_L3_SWA-HIS-COMP-10MIN) · [RPW metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOLO_L3_RPW-BIA-DENSITY-10-SECONDS) · [ESA Solar Orbiter Archive](https://soar.esac.esa.int/soar/) · [SWA product descriptions](https://esdcdoi.esac.esa.int/doi/html/data/heliophysics/SolarOrbiter/SWA.html) · [EAS correction method](https://doi.org/10.1051/0004-6361/202658881)

### Published constraints for the remaining gaps

- **Low-corona proton and electron pressure.** Landi and Testa's
  streamer survey and Gibson et al.'s multi-instrument streamer analysis
  supply electron diagnostic methods and population comparisons. The
  May 8, 1997 SUMER spectrum in Landi and Testa is at 1.04 solar radii
  with a different pointing and epoch from the May 27 east-limb rasters.
  Cranmer's UVCS polar-hole hydrogen analysis demonstrates how neutral
  hydrogen line profiles can constrain proton temperature and its
  anisotropy, with uncertainties and numerical tables. It describes polar
  fast wind at larger radii. Use these sources to check the diagnostic and
  physical plausibility, with plasma population and observation geometry
  explicit in any prior. The dated basal ion pressure for either candidate
  stream is still missing. [Landi and Testa](https://doi.org/10.1088/0004-637X/787/1/33) · [Gibson et al.](https://doi.org/10.1029/98JA02681) · [Cranmer 2020](https://arxiv.org/abs/2007.13180) · [Cranmer numerical tables](https://zenodo.org/records/3908519)
- **Basal composition and helium.** The 2009 September 14 HERSCHEL
  sounding rocket measured paired He II 30.4 nm and H I 121.6 nm coronal
  emission. Moses et al. derive an abundance map out to 2.2 solar radii
  and find significant variation across magnetic structures. This is a
  method and population comparison for composition. For the 2024 eclipse,
  Molnar et al. detected He I 1083 nm on the occulting lunar disk as well
  as above its limb, implicating terrestrial scattered light. Those line
  intensities cannot establish basal alpha abundance in the tracked 2024
  wind. [Moses et al. HERSCHEL analysis](https://www.nature.com/articles/s41550-020-1156-6) · [Full text](https://flore.unifi.it/retrieve/7d895e14-eff8-4ac3-a708-abf171bc2bd9/2020_Nature_Moses_Herschel.pdf) · [Molnar et al. eclipse spectroscopy](https://arxiv.org/abs/2501.01009)
- **Outer wind near 5 AU.** The Ulysses quadratures remain the dated
  observations near the upper boundary. Ebert et al. provide slow-wind
  radial proton trends from Ulysses for a population-level check. Wilson
  et al. demonstrate near-5-AU proton and alpha moments from Juno during
  the 2016 Jupiter approach. The public PDS JADE calibrated archive
  includes the 2024 epoch in its mission coverage. An upstream plasma
  interval, polarity, travel time and magnetic connection to the
  2024 eclipse stream have yet to be established.
  [Ebert et al. Ulysses analysis](https://doi.org/10.1029/2008JA013631) · [Wilson et al. Juno approach](https://doi.org/10.1002/2017JA024860) · [NASA PDS JADE catalogue](https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=JNO-J%2FSW-JAD-3-CALIBRATED-V1.0)

The raw legacy FITS, archive metadata, screened catalogue rows and HAPI
responses are retained beside their query and checksum receipt in
`/data/vtenishe/TEST/stage13-archive-followup-20261123.sP1RzB/audit.json`.
The current Stage 13 case decks contain the previously pinned datasets.

The 1997 follow-up now has a reproducible geometry screen. It uses the
actual CR 1923 Kitt Peak magnetic map, the CDS exposure times and the
archived SOHO and Ulysses ephemerides and plasma moments. The initial
archive-pointing run classified all 6,824 near- and far-side intersections
at 1.05 solar radii as closed for all three trace settings. SolarSoft
processed both rasters and the corrected Mg IX pointing replay classified
all 6,744 intersections as closed. Suess et al. also identify May 27 as
part of the variable streamer interval I2. [Calibration, geometry and fit requirements](soho-ulysses-1997/README.md)

## Use in wind qualification

Freeze the source region, population, epoch, field and flow frames, sector and
transient masks, uncertainties, calibration and qualification roles before
fitting. With co-located ion densities, normal bulk flow and a resolved signed
magnetic sector, the radial mass flux divided by |B_r| can anchor the mass
loading of a steady matched open tube. If proton and alpha drift are measured,
use their separate radial fluxes in the observation operator. A local
electron-core or ion freeze-in temperature is a different observable from
the scalar thermal pressure used in the model. Keep independent forward
observables for qualification, including signed vector field and thermal
channels. A missing required channel or an unresolved mapping stays visible
in the manifest. The calibrated 1997 case has no open target pixels in
these rasters at the inner particle boundary in the tested reconstruction.
A foreground and emission-location analysis would still be needed to use
the measured intensities. For the 2024
candidate, the recovered archive products support a narrower coronal and
inner-heliosphere comparison with one declared mass-loading anchor, an
independent set of thermal and field observables, and a comparison of
transverse force budgets. Separate basal ion pressure, a constrained
transverse force closure, and independent same-stream support to 5 AU
remain open for the live wind qualification.
