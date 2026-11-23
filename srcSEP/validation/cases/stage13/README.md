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
  wind speed over the campaign was approximately (375pm25) km/s. It is
  a promising quiet slow-wind control after excluding the published
  transient windows. The contemporaneous SUMER catalogue lists a May 24
  off-limb Lyα program, May 25 and 27 off-limb Lyβ/O VI programs and a
  May 29 Mg IX program. Their slit footprints, line diagnostics and
  applicability to the east-limb streamer have yet to be established.
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
in the manifest. The available data leave the low-corona proton/electron
thermal profile for either historical quadrature and the spherical material
transport validation open.
