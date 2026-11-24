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
  includes the 2024 epoch in its mission coverage. The ephemeris screen
  below evaluates the arrival geometry for the April 2024 stream.
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

### Additional archive search for basal thermodynamics and outer-wind coverage

#### SOHO and Ulysses in 2003 and 2004

The search used the SOHO campaign calendar, the SUMER planning and FITS
catalogues, the NASA SOHO FITS archive, the MSSL annual CDS catalogues
and NASA CDAWeb Ulysses plasma and field products. The 2003 CDS catalogue
has 130 records tagged with campaign 6942 and 716 tagged with campaign
7007. The 2004 catalogue returned zero records under campaign IDs 7054
and 7101. These are campaign-tag counts. Individual records still require
pointing, line and exposure checks. The SUMER archive supplied the
original FITS files for the 2004 pointings below.
[SOHO 2003 calendar](https://soho.nascom.nasa.gov/soc/calendar2003_old.html) · [SOHO 2004 calendar](https://soho.nascom.nasa.gov/soc/calendar2004_old.html) · [SUMER campaign catalogue](https://www.sidc.be/users/dammasch/mssl/campaigns.html) · [MSSL CDS catalogue](https://www.mssl.ucl.ac.uk/www_solar/gis/) · [NASA SUMER archive](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/)

| Campaign and recovered SUMER data | What the archive establishes |
| --- | --- |
| 7054, May 28–31 and June 1, 2004. [Study 2468](https://www.sidc.be/users/dammasch/mssl/studies.html#2468) and [study 2536](https://www.sidc.be/users/dammasch/mssl/studies.html#2536) | The [May 28 raw spectrum](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/may/sum_040528_000604.fts) has a 50-pixel window centered at 972.640 Å with nominal coverage of 971.53–973.71 Å, which includes H I Lyγ. The nominal slit footprint projects to 1.0945–1.1822 solar radii. The [May 30 Lyα setting](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/may/sum_040530_223226.fts) is at 1215.66 Å with a nominal footprint of 1.1129–1.1482 solar radii. The checked [May 31 setting](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/may/sum_040531_000104.fts) at 749.57 Å has a nominal footprint of 1.0835–1.1467 solar radii. |
| 7101, November 27–December 3, 2004. [Study 2606](https://www.sidc.be/users/dammasch/mssl/studies.html#2606) | The [November 27 raw file](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/november/sum_041127_073904.fts) identifies the Mg IX temperature and Ne V density observing program. This slit intersects the projected 1.05-solar-radius shell in the nominal geometry (range 1.0223–1.0996). The [December 3 example](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/december/sum_041203_101020.fts) projects to 1.0094–1.0282 solar radii. Each exposure still requires the actual spectral windows and pointing, plus a foreground and magnetic-connectivity analysis. |

The projected ranges use the raw primary-header slit center, width,
height and roll, and the apparent solar radius from the nearest hourly
SOHO heliocentric distance, using 695,700 km as the solar radius. They
are a nominal footprint calculation before pixel-level optical pointing
correction and radiometric calibration. A line of sight whose projected
radius exceeds 1.05 solar radii cannot intersect that shell. The hydrogen
settings therefore offer a proton diagnostic above the basal shell. An
inversion would need instrument and scattered-light corrections, separation
of radiative and collisional emission, neutral–proton coupling and
separation of thermal and nonthermal broadening. The published streamer
H I Lyβ/Lyγ analysis describes the emission separation and illustrates
how rapidly usable line coverage can end with height. The SUMER helium
abundance measurements and off-limb line-width studies also document
population and stray-light sensitivity. Their values are for their own
observed regions and dates.
[SOHO ephemeris metadata](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOHO_HELIO1HR_POSITION) · [SUMER archive and calibration](https://arxiv.org/abs/1309.1314) · [Landi and Testa H I diagnostic](https://doi.org/10.1088/0004-637X/800/2/110) · [Laming and Feldman helium diagnostic](https://arxiv.org/abs/astro-ph/0008427) · [Dolla and Solomon line-width corrections](https://doi.org/10.5194/angeo-27-3551-2009)

Ulysses records provide a timing screen before any magnetic or
compositional pairing. For each proton row with positive density,
radial speed and temperatures and a present radial field, the script
back-propagates from its archived distance at its measured radial
speed to 2.5 solar radii. The source interval is the scheduled SUMER
observing dates in the table. Arrival envelopes below encompass the
earliest and latest proton candidates lying in the downloaded windows.
Electron counts within the envelope require positive reported density
and temperature in each of the total, core and halo channels. The
method assumes constant radial speed for timing. Longitude, latitude,
sector continuity, intervening stream interactions and calibration
uncertainties are still required to associate a flux tube.
[SWOOPS proton moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_PROTON-MOMENTS_SWOOPS) · [SWOOPS electron moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_M0_BAE) · [SWOOPS proton and alpha densities](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_M0_BAI) · [Hourly Ulysses fields and plasma](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_COHO1HR_MERGED_MAG_PLASMA)

| SUMER source dates (UTC) | Ulysses time-only candidate envelope (UTC, rounded to day) | Distance (AU) | Proton candidates | Electron rows inside envelope |
| --- | --- | ---: | ---: | ---: |
| 2003 May 19–23 (6942) | June 5–8 | 4.94 | 287 | 0 |
| 2003 November 25–December 4 (7007) | December 13–24 | 5.28–5.30 | 1,111 | 0 |
| 2004 May 28–June 1 (7054) | June 19–22 | 5.41 | 256 | 150 |
| 2004 November 27–December 3 (7101) | December 20–25 | 5.31–5.32 | 677 | 0 |

The November 2003 electron query included 421 valid January 2004 rows
after the candidate envelope. The December 2004 query returned only
electron fill records. The full June 8–July 12, 2004 query returned
1,343 electron rows, 3,552 usable proton rows under the numeric screen,
8,245 SWOOPS rows with positive proton and alpha densities and 811
hourly rows with plasma and vector field values present. The archive
contains both radial-field signs during the 256 June timing candidates.
June 2004 is the strongest of these new near-5-AU leads for a multi-species
comparison. Dated mapping and a quality-controlled matched interval are
still required before computing a mass-loading anchor. The separate
SWOOPS alpha-moments and SWICS kinetic-temperature products provide
additional June coverage, documented below.
[2004 SWOOPS electron archive](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swoops_cdaweb/bae_m0/2004/) · [Ulysses instrument and fill notes](https://cdaweb.gsfc.nasa.gov/misc/NotesU.html#UY_M0_BAE)

#### June 2004 alpha temperatures and coronal connectivity

Two additional archive products supply outer alpha temperatures at 5.41
AU. The SWOOPS alpha distribution moments provide density, radial velocity
and field-aligned parallel and perpendicular kinetic temperatures. The
SWICS main-charge-state moments provide a separate alpha kinetic
temperature and speed for each multi-hour accumulation. SWICS normalizes
its alpha density channel to O VI. Use a measured proton and alpha
density pair for He/H. Both products were retrieved for June 18 through
24 UTC. In the exact proton arrival envelope, June 19 19:56:35 to June
22 16:48:28 UTC, 43 SWOOPS rows contain positive alpha density, radial
speed and both temperatures, and 22 SWICS rows contain positive alpha
speed and kinetic temperature. These are archive coverage counts before
instrument quality or single-sector selection.
[SWOOPS alpha moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_ALPHA-MOMENTS_SWOOPS) · [SWICS kinetic moments](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_M1_SWI) · [2004 native SWICS ZIP](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swics/3p5h_moments/zip/swics04.zip) · [Native archive format](https://spdf.gsfc.nasa.gov/pub/data/ulysses/plasma/swics/new_readme.txt)

The follow-up traces both sides of each nominal SOHO line of sight on
the 1.05, 1.10, 1.15 and 1.20 solar-radius shells for the three recovered
campaign 7054 FITS pointings. The central slit column was sampled every
10 arcseconds at the midpoint of the primary-header observing interval.
None of the rays intersects the basal shell. The 316 intersections above
it were traced with the existing potential-field source-surface (PFSS)
evaluator and a 2.5-solar-radius source surface. The SOLIS/VSM maps are
synoptic pseudo-radial fields in gauss, assembled over a rotation. The
complete 2017 and 2018 maps were tested with harmonic degrees 8 and 12
and integration steps of 0.01 and 0.02 solar radii. The CR 2016 map
contains 5,400 nonfinite cells and was omitted from the trace. The
later CR 2018 map tests sensitivity to magnetic reconstruction and
epoch. Neither reconstruction resolves source variability at the
exposure time.
[SOLIS map processing](https://solis.nso.edu/0/vsm/aboutmaps.html) · [SOLIS rotations 2006–2025](https://solis.nso.edu/0/vsm/crmaps/6302l_cr2006-2025.html)

All 34 intersections at 1.10 solar radii were closed in every run. The
May 28 Lyγ sightline at 1.20 solar radii had open negative-polarity
intersections in both complete maps. In the constant-radial-speed
back-map of the Ulysses proton rows, every run identified the same
nearest arrival at June 20 05:59:58 UTC, with a computed departure at
May 28 01:29:22 UTC. For that sightline at 1.20 solar radii, the closest
same-sign directions on the source surface were 6.72–7.61° apart using
CR 2017 and 1.32–4.02° using CR 2018. The direction search admitted
departures within 24 hours of the exposure midpoint. The angle ranges
describe map and truncation sensitivity and give a candidate interval.
The screen leaves line-of-sight emission attribution, unsteady transport,
and the magnetic connection uncertain.

Native data near this candidate include a SWOOPS proton sample at
05:59:58 UTC, an electron sample at 05:56:01 UTC, an alpha sample at
05:43:55 UTC, and a SWICS multihour alpha moment timestamped 06:05:43 UTC.
The proton sample has radial speed 403.1 km/s and radial magnetic field
−0.0752 nT. The SWOOPS alpha sample has parallel and perpendicular
temperatures of 74,449 K and 99,894 K. The SWICS alpha moment lists 24,980 K as
its kinetic temperature. The samples have different times, integration
intervals and temperature definitions. Keep them separate in the
comparison and carry the difference into any moment and mass-flux analysis.
The retrieval and calculation receipts are in
`/data/vtenishe/TEST/stage13-wind-match.VH7zds/report.json`,
with individual rays in `ray-results.csv` and the replay command
`/usr/bin/python3 check.py` in that directory.

#### Juno arrival geometry for the April 2024 eclipse stream

The archive returned 101 hourly Solar Orbiter positions inside the
published April 7–14 context and 70–105° Carrington range, before
applying the sector, strahl or instrument masks. They span April 9
01:00 to April 13 05:00 UTC at 0.3138–0.3567 AU. The screen uses 840
hourly Juno positions from April 15 through May 19 UTC at 4.9985–5.0137 AU.
For each Solar Orbiter row, the ballistic travel time is the difference
in heliocentric distance divided by speed, with the Juno radius evaluated
at the resulting arrival epoch on the hourly grid. The angular test
compares the measured source and arrival directions in the same
heliographic inertial (HGI) frame. All dates are in UTC.
[Solar Orbiter positions](https://cdaweb.gsfc.nasa.gov/hapi/info?id=SOLO_HELIO1HR_POSITION) · [Juno positions](https://cdaweb.gsfc.nasa.gov/hapi/info?id=JUNO_HELIO1HR_POSITION) · [Rivera et al., Table 2](https://arxiv.org/abs/2610.04669)

| Radial speed (km/s) | Arrival date range | HGI angular separation |
| ---: | --- | ---: |
| 363.54 | May 1–5 | 90.04–114.30° |
| 443.03 | April 27–May 1 | 89.65–113.93° |
| 522.52 | April 24–28 | 89.40–113.69° |

Even allowing every pair of observation times in the full ephemeris
windows, the minimum direction separation is 88.29°. This excludes a
direct ballistic match of the tracked parcels to Juno. A differently
dated recurrent stream would require its own mapping and evidence. The
published Juno upstream and bow-shock measurements on December 10–11,
2024 establish instrument coverage in a separate interval.
[Juno JADE instrument archive](https://pds.nasa.gov/ds-view/pds/viewProfile.jsp?dsid=JNO-J%2FSW-JAD-3-CALIBRATED-V1.0) · [Juno December bow-shock study](https://www.nature.com/articles/s41467-026-76223-x)

#### April 2024 transverse magnetic constraint

The six matched April 9 DKIST metadata records list Stokes I alone.
The MLSO UCoMP archive supplied two April 9 pairs of Level 2 Fe XIII
FITS. The early 1079.8 and 1074.7 nm start times were 18:00:09.92 and
18:07:48.30 UTC, close to the phase-speed and density observations
listed at 18:00 UTC by Rivera et al. The later start times were
21:03:22.87 and 21:07:52.08 UTC. All four files contain weighted
Stokes I, Q and U and a noise mask. The early pair also contains two
azimuth products for each line. The real-data inversion requires
co-registration, time alignment, the occulter and noise masks and
polarization calibration. The phase-speed and density maps must be
compared across the source hole and neighboring structures using the
same geometry and observing window.
[NSO metadata request](https://api.dkistdc.nso.edu/datasets/v1?startTimeMin=2024-04-09T00%3A00%3A00Z&startTimeMax=2024-04-10T00%3A00%3A00Z&instrumentNames=CRYO-NIRSP&pageSize=100) · [MLSO April 9 archive](https://mlso.hao.ucar.edu/mlso_data_summary.php?date=2024-04-09&inst=ucomp) · [Early 1079.8 nm FITS](https://mlso.hao.ucar.edu/hao/acos/2024/04/09/20240409.180009.ucomp.1079.l2.fts) · [Early 1074.7 nm FITS](https://mlso.hao.ucar.edu/hao/acos/2024/04/09/20240409.180748.ucomp.1074.l2.fts) · [Late 1079.8 nm FITS](https://mlso.hao.ucar.edu/hao/acos/2024/04/09/20240409.210322.ucomp.1079.l2.fts) · [Late 1074.7 nm FITS](https://mlso.hao.ucar.edu/hao/acos/2024/04/09/20240409.210752.ucomp.1074.l2.fts)

Paraschiv's published method combines the two lines' linear polarimetry
with a seismology-derived plane-of-sky field. Its reported performance
is based on synthetic observations. Check signal, ambiguity, viewing
geometry and the single-dominant-emitter approximation in the observed
source before using an inversion to constrain the magnetic part of the
transverse force. The full momentum and work comparison still requires
cross-field gradients, thermal and wave stresses, and their joint
uncertainty on neighboring tubes. The dated MAS calculation can provide
a model comparison for those quantities.
[Paraschiv inversion, Sections II–III](https://arxiv.org/html/2607.12846v2) · [Rivera et al.](https://arxiv.org/abs/2610.04669) · [Downs et al. MAS reference data](https://zenodo.org/records/14889337)

The retrieved files, exact URLs, SHA-256 checksums, footprint calculations,
UTC envelopes, and FITS header summaries are retained in
`/data/vtenishe/TEST/stage13-physical-input.Se1RHf/audit.json`. The
companion `audit.py` can be rerun with `/usr/bin/python3 audit.py` in
that directory. No candidate from this search yet supplies basal
proton and alpha pressure and helium abundance on a confirmed open
tube, a constrained transverse force, and independent outer comparison
for the same stream.

### Derived pressures and additional hydrogen spectra

#### June 2004 outer wind

The June 20 05:43:55 UTC proton and alpha SWOOPS records have the same
timestamp and heliocentric distance, 5.41 AU. This ion pair is 16 minutes
and 3 seconds before the closest ballistic candidate. The density and
radial speed pairs are 0.12853 and 0.00451 cm⁻³, and 407.7 and 412.2 km/s,
respectively. They give a local number-density ratio of **3.51%** and
a radial alpha-to-proton number-flux ratio of **3.55%**. In an H–He
plasma with ion charges +1 and +2, the corresponding electron density
would be 0.13755 cm⁻³. For each ion species, the two measured temperature
components give pressures `n k_B T_parallel` and `n k_B T_perp`.
The scalar pressure reported here is one third of the trace of the
gyrotropic pressure tensor.
[Matched proton moments source](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_PROTON-MOMENTS_SWOOPS) · [Matched alpha moments source](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_ALPHA-MOMENTS_SWOOPS)

| Species | Parallel pressure (Pa) | Perpendicular pressure (Pa) | Trace / 3 (Pa) |
| --- | ---: | ---: | ---: |
| Proton | 7.30 × 10⁻¹⁴ | 3.18 × 10⁻¹⁴ | 4.56 × 10⁻¹⁴ |
| Alpha | 4.64 × 10⁻¹⁵ | 6.22 × 10⁻¹⁵ | 5.69 × 10⁻¹⁵ |

The measured radial ion mass flux is `1.00e-16 kg/m²/s`, including
each ion's own radial speed. Four independently retrieved 1-minute VHM
samples within two minutes of the ion timestamp give radial field
values from −0.132 to −0.054 nT. Dividing the ion flux by their
absolute field values gives `0.76e-6–1.85e-6 kg/(s Wb)`. This is a
local sensitivity calculation with the same ion flux in every quotient.
The archive intervals and unresolved stream linkage preclude using
it as a qualified mass-loading anchor.
[Ulysses VHM one-minute magnetic field archive](https://cdaweb.gsfc.nasa.gov/hapi/info?id=UY_1MIN_VHM)

Steady separate conservation of the proton and alpha fluxes would carry
the *number-flux* ratio along a proven tube if both species remain in
their stated ionization states and there are no sources or losses.
The corresponding basal number-density ratio would still be
`(F_alpha / F_p) (U_p / U_alpha)`, using the local normal speeds.
It cannot be obtained from the 5.41-AU density ratio alone. The
temperatures also require species energy equations with heat exchange and
transport. The outer pressure estimates therefore support a local
comparison only.

Two independent SOHO/UVCS O VI-channel raw FITS files confirm additional
May 28 campaign 7054 wavelength coverage. Observation 21242 at
00:28:04–00:48:40 UTC has nominal slit center 1.75 solar radii at PA 95°
and commanded windows that cover H I Lyγ and Lyα. Observation 21244
at 11:36:53–18:13:37 UTC has nominal slit center 1.70 solar radii at
PA 95° and a window covering H I Lyβ. The FITS primary headers establish
dates, pointing and commanded wavelengths. Calibrating and modeling an
actual line profile would give a constraint higher up the tube.
[UVCS 2004 archive](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedUVCS/2004/) · [UVCS instrument and calibration metadata](https://esdcdoi.esac.esa.int/doi/html/data/heliophysics/soho/UVCS.html) · [Dolei et al. hydrogen-width analysis](https://doi.org/10.1051/0004-6361/201628753)

#### The March–April 2024 source

The published basal electron density and temperature give
`p_e = n_e k_B T_e = 3.89e-3 Pa` near 1.05 solar radii. The range
`2.60e-3–5.26e-3 Pa` uses the corners of the published marginal
density and temperature ranges. It is a sensitivity interval without
an assumed joint error distribution. The same paper's in-situ Parker
Solar Probe summary means give a radial ion mass flux per absolute
radial field of `2.17e-6 kg/(s Wb)` when the measured alpha drift
is included. The quotient formed from the tabulated means is a scale
estimate. Fitting an anchor requires products of simultaneous
quality-controlled, sector-selected samples. The coronal sound
speed and enthalpy in that analysis use the adopted ion temperatures.
Reversing those calculations would recover the adopted temperatures.
Even an independent total pressure would constrain only a sum of
species pressures.
[Rivera et al., Sections 3.2 and A.12 and Tables 2–3](https://arxiv.org/abs/2610.04669)

The ESA Solar Orbiter Archive returned 3,675 Level 2 SPICE footprint
records, including multiple releases of some observations, during the
March 28–April 9 search. Two March 29 Probe Quadrature rasters start
at 07:12:33 and 12:45:33 UTC, after the 05:15:01 Metis observation of
the tracked region. The public MEDOC Data Release 6 files both contain
a 31-channel H I Lyβ spectral image with calibrated spectral radiance,
and both record a separate pointing-correction array for each sky axis.
The audit retrieved the actual Lyβ data and the pointing extensions.
Using the nominal WCS in an exploratory region spanning projected
1.04–1.06 solar radii and position angles 96–107°, 3,796 and 1,735 pixels
have finite radiance in the seven wavelength bins around Lyβ. This is a
coverage check around the published Metis outflow direction. The
maximum recorded pointing variations are 24.6 and 6.9 arcseconds, so
the exact aperture and connection require correction and remapping.
[ESA SOAR TAP](https://www.cosmos.esa.int/web/soar/how-to-use-tap) · [MEDOC March 29 FITS directory](https://spice.osups.universite-paris-saclay.fr/spice-data/release-6.0/level2/2024/03/29/) · [SPICE DR6 notes](https://spice.osups.universite-paris-saclay.fr/spice-data/release-6.0/release-notes.html) · [Rivera et al., Figure 14](https://arxiv.org/abs/2610.04669)

The Lyβ window supports a candidate proton diagnostic through a
hydrogen line-formation model. It requires instrument line-spread
calibration, removal of scattered disk emission, a neutral–proton
coupling model, and independent constraints on unresolved motions
and the emitting structure. He II at 1025.27 Å is only 0.0452 nm
from H I Lyβ, about 2.35 bins in these FITS files and smaller than the
published representative SPICE long-wave line width. The archive
check yields no separate coronal helium line measurement in the two
rasters. A helium blend fit would require independent calibration and
detectability evidence before using its intensity or width.
[SPICE in-flight spectra](https://doi.org/10.1051/0004-6361/202141221) · [Instrument PSF study](https://doi.org/10.1051/0004-6361/202245582) · [NIST wavelengths](https://physics.nist.gov/PhysRefData/Handbook/Tables/findinglist1_a.htm)

The Downs et al. archive identifies 3D MAS snapshot `002396` at
19:00 UTC on April 9, about one hour after the published coronal
seismology interval. It includes density, one temperature and all
three components of velocity, field and current. These fields support
an independent model comparison of the Lorentz and thermal-gradient
terms. The archived April 9 field list does not include wave-stress
tensors or separate ion temperatures. An instantaneous transverse
momentum and work comparison also needs the time derivative of the
flow or a demonstrated steady approximation. Neither a scalar wave
energy flux nor the Lorentz term alone fixes the cross-field force.
[Downs et al. data and run index](https://zenodo.org/records/14889337) · [Stage 13 force and work convention](../../../../src/models/sep_corona_swcme/docs/STAGE13_SPHERICAL_SHOCK_ONLY.md#kinematic-material-api-and-physical-ownership)

Run `/usr/bin/python3 derive.py` in
`/data/vtenishe/TEST/stage13-derived-input.Ny8HG2` to replay the
calculations. The resulting `derived.json` records source URLs, byte
ranges, checksums, native units, and the local comparison and
sensitivity values. The observations provide separate basal electron
pressure and outer ion pressure constraints for their own epochs.
The same-stream basal ion pressures, helium abundance and transverse
stress remain unqualified.

### Targeted search for helium at the inner boundary

The NSO inventory for March 28 through April 9, 2024 returned 66
Cryo-NIRSP records. Three April 9 He I 1083 nm records have a
wavelength axis and Stokes I. Their sequential pointings correspond
to the three Fe XIII density-ratio pairs for the tracked southwestern
coronal hole. The largest difference between corresponding archive
bounding-box coordinates is 0.11 arcseconds. The public ASDF
metadata contain the 181 exposure headers for each raster. The
footprint below is computed from the per-exposure FITS WCS,
the gnomonic sky projection and the archived angular solar radius.
The radius is a line-of-sight impact parameter.
[NSO inventory query](https://api.dkistdc.nso.edu/datasets/v1?startTimeMin=2024-03-28T00%3A00%3A00Z&startTimeMax=2024-04-10T00%3A00%3A00Z&instrumentNames=CRYO-NIRSP&pageSize=100) · [Rivera et al., Figures 8 and 12](https://arxiv.org/abs/2610.04669)

| He I raster | April 9 start UTC | Paired Fe XIII rasters | Nominal projected radius (solar radii) | Spatial samples in 1.04–1.06 solar radii |
| --- | --- | --- | ---: | ---: |
| [ICODSO](https://dkist.data.nso.edu/datasetview/ICODSO) | 17:46:01.67 | TFPUCV, NDTDBO | 0.965–1.262 | 27,189 |
| [JTOBIO](https://dkist.data.nso.edu/datasetview/JTOBIO) | 18:00:51.81 | KYLOTD, HPRKOY | 1.148–1.445 | 0 |
| [VLGIHK](https://dkist.data.nso.edu/datasetview/VLGIHK) | 18:15:26.24 | SMTRNK, WRAEVA | 1.351–1.645 | 0 |

ICODSO is the direct observing candidate for this boundary. Its
headers describe the 1080.68–1084.48 nm interval at 0.004214 nm per
spectral sample. The sample count is geometric and precedes
occulter, data quality and stray-light screening. The ASDF metadata
were retrieved and checked. The calibrated Level 1 FITS spectra
require an authenticated Globus transfer and were unavailable in
this session. The He I triplet traces neutral helium. Converting
any detected solar signal to coronal He/H requires a
collisional-radiative model with disk illumination, photoionization,
and a constrained line-of-sight plasma distribution. The neutral
line width also requires a velocity-formation model to constrain the
He²⁺ kinetic temperature. The independent April 8 eclipse spectroscopy
found terrestrial scattered light in its own He I observations,
including over the occulting lunar disk. The DKIST spectra would
require a separate scattering analysis.
[DKIST Level 1 access](https://docs.dkist.nso.edu/projects/python-tools/en/stable/topic_guides/net.html) · [He I and He II emissivity calculations](https://doi.org/10.3847/1538-4357/ab9d84) · [April 8 eclipse scattering analysis](https://arxiv.org/abs/2501.01009)

In the Solar Orbiter view on March 29, the public SOAR archive
provides Level 2 FSI 304 Å context around the 05:15:01 UTC Metis
observation. The retrieved V02 images started at 05:10:20.216 and
05:20:20.217 UTC, with 10-second integrations in the
`FSI_SYNOPTIC_B` mode. The FITS files store one 304 Å passband image
each in DN/s. They provide a near-simultaneous spatial comparison.
That passband combines He II with coronal emission including
Si XI, and the images contain no resolved line width.
[EUI image before Metis](https://soar.esac.esa.int/soar-sl-tap/data?product_type=SCIENCE&RETRIEVAL_TYPE=PRODUCT&data_item_id=solo_L2_eui-fsi304-image_20240329T051020216) · [EUI image after Metis](https://soar.esac.esa.int/soar-sl-tap/data?product_type=SCIENCE&RETRIEVAL_TYPE=PRODUCT&data_item_id=solo_L2_eui-fsi304-image_20240329T052020217) · [FSI 304 Å passband analysis](https://doi.org/10.1051/0004-6361/202244020)

The 401 indexed SUMER records in the two campaign 7054 studies
contain no listed wavelength range bracketing the dedicated
1084.94 Å He II diagnostic. Four native FITS files do carry other
He II wavelength coverage. The May 28 Lyγ window includes 972.11 Å,
and the May 28 Lyα setting includes 1215.09 and 1215.17 Å. The
May 31 and June 1 Lyβ settings include 1025.27 Å at 0.04436 Å per
spectral pixel. These are commanded spectral windows with native
WCS. Calibrated line detection and separation from the adjacent
hydrogen and scattered-light contributions remain unverified.
The four nominal slit footprints start at 1.085–1.112 solar radii,
outside the basal shell. Any confirmed helium signal could
constrain the outer part of the coronal profile for this epoch.
[SUMER campaign 7054](https://www.sidc.be/users/dammasch/mssl/campaigns.html#7054) · [May 28 Lyγ](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/may/sum_040528_000604.fts) · [May 28 Lyα](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/may/sum_040528_083641.fts) · [May 31 Lyβ](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/may/sum_040531_070309.fts) · [June 1 Lyβ](https://soho.nascom.nasa.gov/data/Deprecated/deprecatedSUMER/2004/june/sum_040601_070255.fts) · [NIST He II wavelengths](https://physics.nist.gov/PhysRefData/Handbook/Tables/heliumtable2_a.htm)

The replay command is `/usr/bin/python3 audit.py` from
`/data/vtenishe/TEST/stage13-helium-search.I9L6fE`. Its `report.json`
contains the native header geometry and the exact queries and file
checksums. The search located a co-pointed helium spectrum at the
inner boundary. It did not establish a basal helium abundance or
alpha pressure for a stream with verified outer coverage.

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
