# `srcMoon` clarification questions with code evidence

This file expands the six clarification questions raised by the `srcMoon`
audit.  Each question includes the relevant source excerpts, the behavior that
can be established from those excerpts, and the information needed before the
corresponding production code or validation test can be changed.

The initial audit excerpts describe source revision
`17e86acbafa264a43f2f06b3d2113b29466b5a20`; later resolved-decision and M0
sections describe the working tree based on revision
`fe6153ed5e238359e01f14e13a5588a0176b89bf`.  Excerpts from `build/` are shown
only as evidence of parser output.  `build/` is generated and must not be
patched as the primary fix.  Indentation and configuration line-continuation
markers are normalized for readability; identifiers and active values are
unchanged.  Where long comments or tables are shortened, the omission is
explicitly marked.

## 1. Is the duplicate impact-vaporization source intentional?

### Configuration evidence

`moon.input:174-186` enables a user-defined source whose rate and particle
generator are exactly the impact-vaporization functions.  The built-in impact
source is enabled separately:

```text
Source:UserDefined=on
  SourceProcessCode=MySource
  SourceRate=Exosphere::SourceProcesses::ImpactVaporization::GetTotalProductionRate
  GenerateParticleProperties=Exosphere::SourceProcesses::ImpactVaporization::GenerateParticleProperties
  ModifySurfaceSpeciesAbundance=false

Source:ImpactVaporization=on
  HeliocentricDistance=1.0*_AU_
  SourceRatePowerIndex=0.0
  SourceRate(NA)=1.69E22, SourceTemperature(NA)=6000.0
```

The generated `build/models/exosphere/Exosphere.cpp` consequently adds both
entries in `totalProductionRate()`:

```cpp
#if _EXOSPHERE_SOURCE__IMPACT_VAPORIZATION_ == _EXOSPHERE_SOURCE__ON_
  res += Exosphere::SourceProcesses::ImpactVaporization::
      GetTotalProductionRate(spec, BoundaryElementType,
          (cInternalSphericalData *)(BoundaryElement));
#endif

#if _EXOSPHERE__USER_DEFINED_SOURCE_MODEL__MODE_ == _EXOSPHERE_SOURCE__ON_
  res += Exosphere::SourceProcesses::ImpactVaporization::
      GetTotalProductionRate(spec, BoundaryElementType, BoundaryElement);
#endif
```

It also assigns both source IDs the same rate:

```cpp
FluxSourceProcess[_EXOSPHERE_SOURCE__ID__IMPACT_VAPORIZATION_] =
    Exosphere::SourceProcesses::ImpactVaporization::
        GetTotalProductionRate(spec, BoundaryElementType,
            (cInternalSphericalData *)(BoundaryElement));

FluxSourceProcess[_EXOSPHERE_SOURCE__ID__USER_DEFINED__0_MySource_] =
    Exosphere::SourceProcesses::ImpactVaporization::
        GetTotalProductionRate(spec, BoundaryElementType, BoundaryElement);
```

### Established behavior

For Na at the configured reference distance, each call returns
`1.69e22 s^-1`.  Therefore the aggregate rate contains two equal impact terms,
and particles can be labelled with either the built-in impact ID or `MySource`.
`ModifySurfaceSpeciesAbundance=false` changes surface depletion bookkeeping for
the user source; it does not prevent the second rate from being added.

### Resolved M0 decision (2026-10-10)

The built-in `Source:ImpactVaporization` is the sole active M0 impact source,
and `1.69e22 s^-1` is its total Na rate.  `MySource` is disabled in the input
but retained for possible later development.  Its reported intent is to track
Na adsorbed on the nightside and release it on the dayside; the current
callbacks shown above do not implement that reservoir or transport law, so the
source cannot be re-enabled until those physics and I19 conservation are
defined and verified.

The alternatives below are retained as the audit trail that led to the
decision:

The considered alternatives were:

1. One physical impact source: disable/remove `MySource` and keep the built-in
   impact source.
2. One custom impact source: disable the built-in source and keep `MySource`,
   explaining why custom dispatch is required.
3. Two physically distinct impact populations: retain both, give each a
   distinct name, rate, generator, provenance, and depletion policy.  The sum
   must then be the declared physical total rather than an accidental doubling.

The answer should also state whether `1.69e22 s^-1` is meant to be the total
impact rate or the rate of each of two components.

## 2. Resolved: inertial frame and SPICE kernel set for lunar rotation

### Confirmed decision

The confirmed contract is `J2000` inertial, Moon-centred `LSO` rotating,
`MOON_ME_DE421` body-fixed, geometric aberration correction `NONE`, and fixed
epoch `2009-01-24T00:00:00 UTC`. The exact kernel set is `de421.bsp`,
`moon_pa_de421_1900-2050.bpc`, `moon_080317.tf`, `moon_assoc_me.tf`,
`naif0012.tls`, `pck00011.tpc`, and `OTHER/Moon.LSO.tf`.

The production constants now centralize that choice in `srcMoon/Moon.h`:

```cpp
namespace Frames {
  static const char Inertial[]="J2000";
  static const char SolarOrbital[]="LSO";
  static const char BodyFixed[]="MOON_ME_DE421";
  static const char ForceAberrationCorrection[]="NONE";
}
```

The production time-step path no longer contains `MSGR_HCI`. It queries the
beginning, midpoint, and end lunar states in J2000, freezes external force
states at the midpoint, and gets the instantaneous rotating-frame vector from
the SPICE state-transform derivative:

```cpp
sxform_c(Moon::Frames::Inertial,Moon::Frames::SolarOrbital,etMiddle,
    J2000ToLSO);
Moon::OrbitalDynamics::GetSolarOrbitalAngularVelocityLSO(etMiddle,
    Moon::RotationVector_SO_FROZEN);
```

`GetSolarOrbitalAngularVelocityLSO()` uses `xf2rav` rather than the former
finite-angle `acos`/`sin(angle)` construction. U04 independently reconstructs
the angular velocity with `W=(dR/dt)R^T`, checks every differential-gravity and
fictitious-force component, verifies rotation orthogonality, and verifies the
exact SHA-256 of all seven kernels before CSPICE is called.

### Current status

This question no longer blocks U04: the SPICE-enabled linked-kernel probe is
PASS-capable at the frozen epoch. It does not by itself complete I06 or I10;
those still require actual AMPS mover trajectories and convergence evidence.

## 3. Is surface abundance prescribed or dynamically conserved?

### Configuration evidence

`moon.input:131-132` explicitly selects user-defined surface content:

```text
define _EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(spec,el) (3.0E16)
define _EXOSPHERE__SURFACE_CONTENT_ _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_
```

Later, `moon.input:220-222` declares a different constant for Na:

```text
SurfaceVolatileDensity=const
  function=?????
  const(Na)=2.3E16
```

The parser-generated definition uses the later species table, so the effective
Na value in the clean build is `2.3e16 m^-2`, not `3.0e16 m^-2`.

At a sticking impact, the production callback in
`src/models/exosphere/Exosphere.cpp` records downward statistical weight and
deletes the particle:

```cpp
if (rnd() < StickingProbability(spec, ReemissionParticleFraction,
        SurfaceTemp)) {
  Sphere->GetSurfaceElementProjectionIndex(x_LOCAL_IAU_OBJECT,
      nZenithElement, nAzimuthalElement);
  el = Sphere->GetLocalSurfaceElementNumber(nZenithElement,
      nAzimuthalElement);

  Sphere->SurfaceElementAdsorptionFluxDOWN[spec][el] +=
      ReemissionParticleFraction*ParticleWeight;
  Exosphere::Sampling::PlanetSurfaceStickingRate[spec] += ParticleWeight;
  return _PARTICLE_DELETED_ON_THE_FACE_;
}
```

`src/models/exosphere/Exosphere_Parallel.cpp:24-31` constructs the MPI-reduced
net flux:

```cpp
memcpy(TotalFlux_LOCAL, Sphere->SurfaceElementAdsorptionFluxDOWN[spec],
    TotalSurfaceElementNumber*sizeof(double));
memcpy(TotalFlux_GLOBAL, Sphere->SurfaceElementDesorptionFluxUP[spec],
    TotalSurfaceElementNumber*sizeof(double));
for (el=0; el<TotalSurfaceElementNumber; el++)
  TotalFlux_LOCAL[el] -= TotalFlux_GLOBAL[el];

MPI_Allreduce(TotalFlux_LOCAL, TotalFlux_GLOBAL,
    TotalSurfaceElementNumber, MPI_DOUBLE, MPI_SUM,
    MPI_GLOBAL_COMMUNICATOR);
```

It adds that net flux to the population, but user-defined mode immediately
overwrites the result (`Exosphere_Parallel.cpp:101-110`):

```cpp
t = TotalFlux_GLOBAL[el];
Sphere->SurfaceElementPopulation[spec][el] += t;

#if _EXOSPHERE__SURFACE_CONTENT_ == \
    _EXOSPHERE__SURFACE_CONTENT__USER_DEFINED_
  Sphere->SurfaceElementPopulation[spec][el] =
      SurfaceElementArea[el] *
      (_EXOSPHERE__SURFACE_CONTENT_DENSITY__USER_DEFINED__FUNCTION_(spec,el));
#endif
```

Thermal desorption currently uses an Arrhenius-like production rate from the
instantaneous surface population (`src/models/exosphere/Exosphere.h`):

```cpp
res = ThermalDesorption_VibrationalFrequency[spec] *
    SodiumSurfaceElementPopulation *
    exp(-ThermalDesorption_uThermal[spec]/(Kbol*temp));
```

### Established behavior

Surface counters exist and net adsorption/desorption is calculated, but the
effective user-defined branch resets Na to a prescribed density on every
exchange.  It is therefore an externally replenished fixed boundary condition,
not a closed accumulating reservoir.  No separate residence-time clock,
surface diffusion state, PSR trap state, or delayed-release queue is visible in
this path.

### Approved decision (2026-10-13)

Please specify whether the intended baseline is:

- prescribed surface density, in which case identify the authoritative value,
  spatial map, units, and replenishment interpretation; or
- a conserved reservoir, in which case define initial inventory, capacity,
  adsorption fraction, residence-time law, thermal/PSD/sputtering depletion,
  lateral migration, cold trapping, release, and boundary conditions.

For a dynamic reservoir, also specify which surface variables and random state
must survive checkpoint/restart and the conservation tolerance expected by
I19/I20.

## 4. What is the intended production species table?

### Active-list evidence

The executable currently requests only Na (`moon.input:4`):

```text
SpeciesList=Na  !, NA_Plus
```

`species.input` contains definitions that could be selected, including:

```text
#component=Na
  mass=3.82E-26
  charge=0

#component=NA_PLUS
  mass=3.82E-26
  charge=1

#component=Ar
  mass=40*_AMU_
  charge=0

#component=AR_PLUS
  mass=40*_AMU_
  charge=1.0

#component=H2O
  mass=18*_AMU_
  charge=0

#component=OH
  mass=17*_AMU_
  charge=0

#component=He
  mass=4.0*_AMU_
  charge=0.0

#component=He_PLUS
  mas =4.0*_AMU_
  charge=1.0
```

The `He_PLUS` definition contains `mas` rather than `mass`; it must be checked
against parser behavior before He+ is enabled.

The Na photochemical function in `srcMoon/Moon.h:411-449` accepts Na only and
can create Na+ only if that species is present:

```cpp
if (spec != _NA_SPEC_) {
  PhotolyticReactionAllowedFlag = false;
  return -1.0;
}

if ((spec == _NA_SPEC_) && (_NA_PLUS_SPEC_ >= 0)) {
  PIC::ParticleBuffer::SetI(_NA_PLUS_SPEC_, ParticleData);
  return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;
}

return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
```

The surface dispatch supports a broader but fixed set
(`srcMoon/Moon.cpp:195-212`):

```cpp
switch (spec) {
case _NA_SPEC_: case _NA_PLUS_SPEC_:
  res = SodiumStickingProbability(ReemissionParticleFraction, Temp);
  break;
case _AR_SPEC_: case _AR_PLUS_SPEC_:
  res = ArgonStickingProbability(ReemissionParticleFraction, Temp);
  break;
case _HE_SPEC_: case _HE_PLUS_SPEC_:
case _NE_SPEC_: case _NE_PLUS_SPEC_:
  res = 0.0;
  break;
default:
  exit(__LINE__, __FILE__, "the option is not implemented");
}
```

### Established behavior

Definitions in `species.input` do not activate a species.  The active
`SpeciesList` controls the compiled indices and reachable reaction/source
paths.  Adding a name without completing every mover, chemistry, surface,
weight, source, and diagnostic branch can cause a build failure, particle
removal instead of conversion, or a runtime `exit()`.

### Decision needed

Please provide a species matrix with, for every intended production species:

- canonical configuration name and chemical identity/isotope;
- mass in kg and charge in coulombs or documented AMPS charge units;
- neutral/ion parent and daughter reactions with branching ratios;
- source processes and surface sticking/accommodation/residence behavior;
- particle-weight and timestep ratio relative to the reference species;
- required diagnostics and observation operator;
- baseline-enabled, optional, or deliberately unsupported status.

At minimum, confirm whether the target list is Na alone, Na+Na+, or the roadmap
set Na/Na+/He/He+/Ne/Ne+/Ar/Ar+/H2O/OH.

## 5. Which LOLA product and triangulated-surface convention should be used?

### Active production geometry

`srcMoon/main_lib.cpp:1009-1014` creates a spherical boundary:

```cpp
cInternalSphericalData::SetGeneralSurfaceMeshParameters(60, 100);
PIC::BC::InternalBoundary::Sphere::Init(ReserveSamplingSpace, NULL);
SphereDescriptor = PIC::BC::InternalBoundary::Sphere::RegisterInternalSphere();
Sphere = (cInternalSphericalData *)SphereDescriptor.BoundaryElement;
Sphere->SetSphereGeometricalParameters(sx0, rSphere);
```

Thus the two programs under `srcMoon/Topography/` are not part of the active
application.

### LOLA prototype evidence

`srcMoon/Topography/lola.cpp:31-37` hard-codes filenames:

```cpp
string fNameLON = "ldem_4_float_LON.txt";
string fNameLAT = "ldem_4_float_LAT.txt";
string fNameALT = "ldem_4_float.txt";
char fNameSphere[100] = "sphere-fine.mail";
```

Altitude lookup selects one grid cell by integer truncation, with endpoint
clamping (`lola.cpp:127-143`):

```cpp
GetLonLat(Lon, Lat, x);

iLon = (Lon-LonDataStart)/dLon;
if (iLon < 0) iLon = 0;
else if (iLon >= vLON.size()) iLon = vLON.size()-1;

iLat = (Lat-LatDataStart)/dLat;
if (iLat < 0) iLat = 0;
else if (iLat >= vLAT.size()) iLat = vLAT.size()-1;

return vALT[iLat][iLon];
```

The file reader multiplies altitude values by 1000, and the mesh is displaced
radially (`lola.cpp:106-114,145-156`):

```cpp
vRes[i].push_back(1.0E3*number);

l = Vector3D::Length(CutCell::BoundaryTriangleNodes[iVertex].x);
if (l > 1.0) {
  alt = GetAltitude(CutCell::BoundaryTriangleNodes[iVertex].x);
  Vector3D::MultiplyScalar((l+alt)/l,
      CutCell::BoundaryTriangleNodes[iVertex].x);
}
```

The result is written in a CEA long format:

```cpp
PIC::Mesh::IrregularSurface::SaveCEASurfaceMeshLongFormat(
    "surface-mesh.mail");
```

The separate shadow prototype has machine-specific kernels and paths
(`srcMoon/Topography/shadow_calc.cpp:32-40`):

```cpp
char fNameMesh[100] = "surface-mesh.mail";
char SpicePath[100] = "/nobackup/vtenishe/SPICE/Kernels";

const int nKernels = 4;
char Kernels[100][100] = {
  "NAIF/naif0010.tls", "NAIF/de430.bsp",
  "NAIF/pck00010.tpc", "NAIF/sat375.bsp"
};
```

### Established behavior

The prototype assumes text longitude/latitude/altitude arrays, kilometres for
altitude, radial displacement of an existing sphere mesh, nearest-cell
sampling, and CEA long-format output.  It does not record a LOLA product ID,
datum/radius, longitude direction, pixel registration, polar rule, missing
value, uncertainty, or interpolation method.  The shadow program is not
portable or wired into the production illumination callback.

### Decision needed

Please identify:

- exact D01/LOLA product name, archive identifier/version, file hash, and
  whether radius, elevation, or topography is stored;
- datum/reference radius and all input/output units;
- body-fixed frame, longitude direction/range, latitude type, pixel-centre or
  pixel-edge registration, pole and seam handling;
- fill-value and data-gap policy—real gaps must remain gaps;
- allowed interpolation/resampling method and target mesh resolution;
- required triangle winding, normal direction, watertightness, and CEA format
  version;
- whether production should replace the sphere globally or use terrain only
  for illumination/temperature while keeping a spherical particle boundary;
- authoritative SPICE kernels and solar-ray convention for terrain shadowing.

### Resolved surface decision (2026-10-10)

The requested production choice is now encoded rather than assumed. The
application file names the product, raw detached IMG/LBL files, maximum mesh
edge in metres, and both outputs (`srcMoon/examples/lola_surface.in`):

```text
#section begin: moon
surface_geometry = lola
surface_mesh_resolution_m = 100000
lola_product_id = LDEM_4
lola_image_file = /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.IMG
lola_label_file = /data/vtenishe/moon_validation_data/lola/raw/LDEM_4.LBL
surface_cea_file = ../../test_output/srcMoon/surface/lola_surface.cea
surface_tecplot_file = ../../test_output/srcMoon/surface/lola_surface.dat
#section end
```

`LunarSurface.cpp` rejects an incompatible label rather than guessing its
interpretation:

```cpp
if (candidate.projection != "SIMPLE CYLINDRICAL" ||
    candidate.coordinateSystemType != "BODY-FIXED ROTATING" ||
    candidate.positiveLongitudeDirection != "EAST" ||
    candidate.sampleType != "LSB_INTEGER" || candidate.sampleBits != 16 ||
    candidate.unit != "METER") {
  *error = "unsupported LOLA label ...";
  return false;
}
```

The native label reports `MEAN EARTH/POLAR AXIS OF DE421`; the corresponding
NAIF frame name is `MOON_ME_DE421`. The implementation uses the label-native
east longitude and planetocentric latitude directly and does not substitute
`IAU_MOON`. It builds an icosphere to the requested maximum great-circle edge,
bilinearly samples centre-registered pixels with longitude wrapping and polar
row clamping, adjusts every vertex radius, writes CEA and Tecplot, and reloads
CEA through the normal AMPS loader:

```cpp
generated = dem.Load(configuration, &localError) &&
    BuildIcosphere(dem.metadata().referenceRadiusM,
        configuration.surfaceMeshResolutionM, &mesh, &localError) &&
    ApplyLolaTopography(dem, &mesh, &localError) &&
    WriteCeaSurface(mesh, configuration.surfaceCeaFile, &localError) &&
    WriteTecplotSurface(mesh, configuration.surfaceTecplotFile, &localError);

PIC::Mesh::IrregularSurface::ReadCEASurfaceMeshLongFormat(
    configuration.surfaceCeaFile.c_str(), 1.0);
```

The loaded triangulation replaces the analytic sphere as the geometric
particle boundary in LOLA mode. The established 60×100 spherical grid remains
only for source and inventory indexing. Orbit-on moving terrain and terrain
shadowing remain unresolved; the code therefore rejects LOLA mode when orbit
calculation is compiled on. D01 also remains unqualified because its local
directory lacks package README/provenance/QA files, so this implementation is
EXPERIMENTAL and not observational validation.

## 6. Which authoritative references control the Na parameters?

### Source-parameter evidence

The active values in `moon.input:183-196` have abbreviated and uncertain
comments:

```text
Source:ImpactVaporization=on
  HeliocentricDistance=1.0*_AU_
  SourceRatePowerIndex=0.0
  SourceRate(NA)=1.69E22, SourceTemperature(NA)=6000.0
  ! Killen-2012-JGR; alternatives remain in the comment

Source:ThermalDesorption=on
  uThermal(NA)=1.85*eV2J, VibrationalFrequency(NA)=1.0E13

Source:PhotonStimulatedDesorption=on
  PhotonFlux_1AU=2.0E14*1.0E4
  ! Killen-2012-JGR, Yakshinskii+Madey-1999-?
  CrossSection(NA)=3.0E-21*1.0E-4,
  InjectionVelocityRange(NA)=10,10.0E3
  ! "Satantos-2010-?" and alternative values remain in the comment
```

The Na sticking table in `srcMoon/Moon.cpp:152-166` is labelled only as a
digitization:

```cpp
static const double tmin = 100.0;
static const double tmax = 500.0;
static const double dt = 5.0;
static const int nPoints = 80;

struct cDataStruct { double t, s; };

static const cDataStruct data[nPoints] = {
  // digitized from Yakshinskiy-2005-SS
  {100.00000, 0.99830}, {105.00000, 0.96768},
  {110.00002, 0.94512}, /* ... */ {494.99997, 0.18605}
};
```

The radiation-pressure function in `src/species/Na.cpp:73-87` identifies a
paper and units but not the table extraction provenance:

```cpp
// the solar radiation pressure as a function of heliocentric distance and velocity
// taken from Combi-1997-icarus.pdf
// input: heliocentric velocity (m/s), heliocentric distance (m)
// return: radiation pressure acceleration (m/s^2)

double SodiumRadiationPressureAcceleration__Combi_1997_icarus(
    double HeliocentricVelocity, double HeliocentricDistance) {
  static const int
      sodiumRadiationPressure_TableLength__Combi_1997_icarus = 1201;
  static const double
      sodiumRadiationPressure_HeliocenticVelocity_Min__Combi_1997_icarus =
          -60.0E3;
  static const double
      sodiumRadiationPressure_HeliocenticVelocity_Max__Combi_1997_icarus =
          60.0E3;
  static const double
      sodiumRadiationPressure_dHeliocenticVelocity__Combi_1997_icarus =
          0.1E3;
  // embedded values follow
}
```

The implementation clamps values outside the tabulated velocity interval,
linearly interpolates inside it, and scales by inverse-square heliocentric
distance.  U05 verifies that scaling invariant but not the absolute table.

### Established behavior

The compiled values and formulas can be regression-tested, but the comments do
not provide enough bibliographic and transformation detail to reproduce the
science inputs independently.  In particular, alternate values remain in
comments, at least one citation name is misspelled/uncertain, and digitization
methods and errors are not recorded.

### Decision needed

For impact vaporization, PSD, thermal desorption, radiation pressure/g factors,
and Na sticking, please provide:

- complete citation: authors, title, journal, year, volume/pages, DOI or archive
  identifier;
- exact source table/equation/figure and version;
- original units and the explicit conversion to AMPS units;
- digitization file/script, calibration points, and uncertainty where a figure
  was digitized;
- valid temperature, velocity, heliocentric-distance, illumination, and surface
  composition ranges;
- extrapolation/clamping rule;
- whether each number is a fixed literature input, calibration-eligible
  parameter, or held-out validation quantity.

## 12. M2: which frame contract controls an evolving LOLA boundary?

### Production-code evidence

The production LOLA mesh is stored in the label-native
`MOON_ME_DE421` body-fixed axes and loaded once as a static AMPS cut-cell
surface.  `srcMoon/LunarSurface.cpp` therefore rejects the combination with
the currently compiled orbit/evolving-frame path:

```cpp
#if _EXOSPHERE__ORBIT_CALCUALTION__MODE_ == _PIC_MODE_ON_
  // The cut-cell triangulation is static in the AMR frame.  Enabling lunar
  // rotation without moving/rebuilding it would mix the label body-fixed frame
  // with LSO and produce physically incorrect impacts, so fail explicitly.
  *error = "LOLA surface mode currently requires orbit calculation off: the "
      "static AMR triangulation is in the LDEM body-fixed frame";
  return false;
#endif
```

The normal Moon mover, force, source, and observation paths instead use the
Moon-centred solar-orbital `LSO` frame.  For example, `srcMoon/Moon.h` defines:

```cpp
namespace Frames {
  static const char Inertial[]="J2000";
  static const char SolarOrbital[]="LSO";
  static const char BodyFixed[]="MOON_ME_DE421";
}
```

This is not a coordinate-label-only issue.  A static body-fixed triangle has a
time-dependent location and surface velocity in LSO, so a physically correct
evolving collision must account for both rotation and the relative incident
velocity.  Reinterpreting the same numeric vertex coordinates as LSO would be
wrong.

### Decision needed

Please select the production contract for M2 evolving runs:

1. **Recommended:** retain LSO for AMPS particles and the downstream sodium
   tail, but implement a time-dependent body-fixed LOLA collision/intersection
   service.  Each particle segment is transformed at the applicable ephemeris
   times, the crossing with the body-fixed triangulation is solved to a declared
   tolerance, and accommodation uses velocity relative to the rotating facet.
   The static cut-cell surface can remain an initialization/AMR aid but cannot
   be treated as the evolving physical boundary.
2. Change the entire AMPS mesh/mover frame to `MOON_ME_DE421` and generalize
   gravity, rotating-frame forces, sources, diagnostics, and the long sodium
   tail domain accordingly.
3. Limit LOLA to frozen-epoch runs.  This can support part of I14--I16 but
   cannot complete I17 or M2 and therefore is not recommended as the M2
   contract.

For option 1, please also confirm that collision geometry is evaluated at the
particle-segment time, not frozen at the beginning of a global AMPS step.

## 13. M2: may the missing D02 average-visibility products be acquired?

### Data evidence

The canonical D02 raw directory currently contains only the two permanent-
shadow products:

```text
LPSR_65N_240M_201608.IMG/.LBL
LPSR_65S_240M_201608.IMG/.LBL
```

The label states that these are binary permanent-shadow classifications:

```text
PRODUCT_ID = "LPSR_65N_240M_201608"
NAME = PERMANENT_SHADOW
DESCRIPTION = "... permanently shadowed regions ... binary ..."
MAP_SCALE = 240 <m/pix>
COORDINATE_SYSTEM_NAME = "MEAN EARTH/POLAR AXIS OF DE421"
```

I15 additionally requires independent visibility/classification and
sunrise/sunset evidence.  The roadmap names the matching average-visibility
products `AVGVISIB_65N_240M_201608` and `AVGVISIB_65S_240M_201608`; they are not
present.  A PSR mask cannot be converted into time-dependent illumination
evidence without inventing information.

### Decision needed

May the two matching `AVGVISIB` IMG/LBL pairs be downloaded from the
authoritative LOLA archive into the immutable D02 raw package?  The recommended
answer is yes.  If not, I15 must remain `SKIPPED / NOT VALIDATED`; the existing
LPSR products can still support the independent PSR portion of I20.

## 14. M2: which Diviner thermal model and frozen event split control I17?

### Data and test-contract evidence

The canonical D03 package contains all 18 GCP latitude bands, both north and
south PRP tables, and complete ten-minute RDR packages for the lunar-eclipse
dates 2011-06-15 and 2011-12-10.  However, U11 deliberately has no invented
solver tolerance or parameter set:

```json
"independent_reference": "Equilibrium solution, linear periodic response, and independently integrated one-dimensional heat equation.",
"acceptance_criteria": "Thresholds must be frozen after a solver and independent reference are selected; no threshold is inferred here.",
"skip_reason": "No dynamic thermal solver or D03 eclipse history is present."
```

The last part of that skip reason is now stale because the raw RDR packages
have been downloaded, but the solver/parameter decision remains unresolved.
The roadmap permits either a one-dimensional conduction model or an equivalent
reviewed reduced model; those choices are not numerically or physically
interchangeable.

### Decision needed

Please confirm the recommended contract:

- implement a depth-resolved, energy-conserving one-dimensional conduction
  solver with explicit SI albedo, emissivity, density, heat-capacity, and
  conductivity/depth parameters rather than an empirical relaxation curve;
- use 2011-06-15 only for parameter training/initialization and freeze
  2011-12-10 as the untouched I17 hold-out event;
- use GCP/PRP for `DIVINER_MAP`, not as a synthetic eclipse history;
- require every material parameter and its citation/range in a versioned input
  table rather than hiding lunar-regolith constants in C++.

The remaining required information is the authoritative thermophysical model
or publication controlling those material parameters and which parameters, if
any, may be calibrated on the June training event.  Without that source and
calibration list, I can implement and analytically verify a general solver but
cannot claim the lunar thermal model or I17 is qualified.

## 15. M2: which species and adsorption/residence parameters control I18--I20?

### Production-code evidence

The current build activates Na only.  Its generic surface callback has
temperature-dependent sticking hooks for Na and Ar, forces He/Ne sticking to
zero, and has no residence-time state:

```cpp
switch (spec) {
case _NA_SPEC_: case _NA_PLUS_SPEC_:
  res=SodiumStickingProbability(ReemissionParticleFraction,Temp);
  break;
case _AR_SPEC_: case _AR_PLUS_SPEC_:
  res=ArgonStickingProbability(ReemissionParticleFraction,Temp);
  break;
case _HE_SPEC_: case _HE_PLUS_SPEC_: case _NE_SPEC_: case _NE_PLUS_SPEC_:
  res=0.0;
  break;
default:
  exit(__LINE__,__FILE__,"the option is not implemented");
}
```

On a sticking event the current implementation increments a downward counter
and deletes the particle immediately:

```cpp
Sphere->SurfaceElementAdsorptionFluxDOWN[spec][el]+=
    ReemissionParticleFraction*ParticleWeight;
return _PARTICLE_DELETED_ON_THE_FACE_;
```

The generic MPI exchange can update `SurfaceElementPopulation`, but the active
user-defined surface-content mode overwrites it with a prescribed density.
Therefore the current model is not a conserved reservoir and cannot satisfy
I18--I20.

### Decision needed

I19 requires Na and at least one noble gas.  Please confirm the recommended M2
scope of **Na plus Ar**, with ions using their own explicitly supplied surface
law rather than silently inheriting the neutral law.  For each enabled species
or charge state, provide or approve a literature review that freezes:

- sticking probability law and valid temperature/composition range;
- binding-energy value or distribution in joules/eV;
- attempt frequency in s^-1;
- thermal accommodation and re-emission angular/speed kernel;
- adsorption-site capacity, if finite;
- residence-time extrapolation policy outside the measured range;
- whether surface diffusion exists in addition to ballistic re-emission.

The current abbreviated `Yakshinskiy-2005-SS` and piecewise Ar comments are not
sufficient provenance for new M2 physics.  I will not promote those regression
values into a qualified residence model without an authoritative source.

## 16. M2: how is finite residence distinguished from a cold trap?

### Roadmap and state evidence

The roadmap explicitly rejects equating geometric PSR with permanent trapping:

```text
Combine terrain/illumination, thermal state and species binding kinetics;
distinguish temporary adsorption from cold-trapped inventory; allow re-release
if a site warms.
```

No production cold-trap state, release queue, restart representation, or
classification timescale currently exists.  The physically computed
Arrhenius residence time can be retained without a binary threshold, but I20
also requires reporting a species-specific thermal-retention classification.

### Decision needed

Please define or approve the recommended configuration-driven policy:

- the reservoir always uses the continuous residence law and releases material
  when its sampled desorption time is reached; `cold_trapped` is a diagnostic,
  never an irreversible sink;
- classification requires both the production terrain-illumination history and
  residence time greater than a declared species-specific reference timescale;
- the reference timescale is explicit in the input/parameter table rather than
  hard-coded;
- M2 starts from zero dynamic adsorbed inventory unless an authoritative
  initial-inventory map is supplied;
- inventory, scheduled release state, random state, and cumulative budget
  counters are included in checkpoint/restart.

Please provide the physical reference timescale for Na and Ar, or the cited
rule from which it is derived.  Choosing the simulation duration merely because
it is convenient would make the classification run-dependent and is not an
acceptable unstated assumption.

## 17. M2 predecessor-gate status

The formal M2 milestone cannot be declared complete until FULL-profile
I11--I13 pass.  The current eight-rank regression evidence is correctly
`SKIPPED / NOT QUALIFIED` because it uses `mesh_profile=custom`.  M2 production
and test infrastructure can be developed meanwhile, but I14--I20 results must
remain phase-gated and cannot substitute for the missing M1 qualification.

The authoritative material should be stored as immutable raw input with hashes,
and every processed table should be regenerated by a committed script.  AMPS
residuals must not be used to select among the alternative values already
present in comments.

## Requested response format

A concise answer can use this template:

```text
Q1 duplicate impact source:
  intended option:
  physical total rate:
  depletion policy:

Q2 frames/SPICE:
  inertial frame:
  LSO definition:
  kernels + SHA-256:
  epoch range/reference cases:

Q3 surface inventory:
  prescribed or conserved:
  authoritative equations/parameters:
  restart and conservation requirements:

Q4 species:
  enabled baseline species:
  optional species:
  reactions, weights, and timesteps:

Q5 LOLA/terrain:
  product/version/hash:
  datum/coordinates/units:
  interpolation and mesh conventions:
  production integration choice:

Q6 Na provenance:
  complete references:
  source equations/tables/figures:
  conversions, validity ranges, and uncertainty:
  fixed versus calibration-eligible parameters:

Q12 evolving LOLA frame:
  option 1, 2, or 3:
  segment-time geometry confirmed:

Q13 D02 acquisition:
  authorize AVGVISIB north/south download:

Q14 thermal model:
  approve depth-resolved solver:
  approve June-training/December-holdout split:
  controlling publication/parameter table:
  calibration-eligible parameters:

Q15 surface species:
  approve Na+Ar scope:
  ion policy:
  authoritative parameter sources or authorize literature review:

Q16 cold trapping:
  approve continuous-release diagnostic policy:
  Na reference timescale + source:
  Ar reference timescale + source:
  zero initial dynamic inventory:
  restart-state requirements:
```

## 7. Which stand-alone gates are prerequisites for declaring M0 complete?

### Roadmap evidence

The milestone table defines M0 in terms of the linked baseline
(`srcMoon/development_and_validation.md:50-53`):

```text
| M0 | Linked baseline qualified | I01–I10 run in the actual linked AMPS Moon
executable; no unresolved compile/API mismatches; provenance frozen. |
```

The R0.1 acceptance text is broader (`development_and_validation.md:1764`):

```text
Run U01–U26 and I01 preflight; repeat identical manifest generation twice;
deliberately change one control.

Acceptance / exit criteria

U01–U26 remain PASS; identical-control manifests are equivalent; a changed
control is recorded; no unknown/unhashed external file is allowed in a release
run.
```

However, the same roadmap schedules capabilities exercised by many U-tests in
later milestones.  For example, M2 requires the LOLA boundary, terrain
illumination, Diviner/dynamic thermal state, inventory, and cold trapping; M3
requires multi-species chemistry; and M4 requires measured drivers and
multi-species source modules.  Those production capabilities do not currently
exist, so their corresponding U-tests correctly report `SKIPPED` rather than
inventing evidence.

The implemented runner also explicitly distinguishes registered tests from
available capabilities (`srcMoon/test/run_tests.py`):

```python
if counts["ERROR"]:
    return 2
if counts["FAIL"]:
    return 1
if args.strict_skips and counts["SKIPPED"]:
    return 1
```

### Established behavior

All literal U01–U26 IDs are registered and executed, but the latest campaign
contains 8 PASS and 18 SKIPPED.  No I01–I10 linked runners exist yet.  Requiring
every U-test to PASS before M0 would move substantial M2–M5 production physics
into Phase 0 and make the documented sequential milestone order impossible.

### Resolved M0 decision (2026-10-10)

The phase-appropriate baseline gate (option 1 below) controls M0.  All literal
U01-U26 IDs remain registered and execute with honest status; local
prerequisites used by I01-I10 must PASS, later-capability U-tests may remain
explicitly SKIPPED, and I01-I10 must all PASS.  The alternatives are retained
below as the decision record:

1. **Phase-appropriate baseline gate (recommended):** all U01–U26 remain
   registered and execute with honest status; the local prerequisites directly
   used by I01–I10 must PASS; later-capability U-tests may remain explicitly
   `SKIPPED`; I01–I10 must all PASS with frozen provenance.
2. **Strict all-U gate:** every U01–U26 must PASS before M0.  This makes the M0
   work package include later LOLA/thermal/chemistry/driver/species/operator
   implementation and changes the stated milestone ordering.

If option 1 is selected, the initial local M0 prerequisite set will be recorded
explicitly in the roadmap/test registry instead of being inferred during a
campaign.

## 8. May I use a dedicated Na/Na+ integration-test configuration for I08?

### Code evidence

The production input enables neutral sodium only (`input/moon.input:4`):

```text
SpeciesList=Na  !, NA_Plus
```

The production force callback applies Lorentz acceleration to any configured
charged species (`srcMoon/Moon.h`):

```cpp
if (PIC::MolecularData::ElectricChargeTable[spec] != 0.0) {
  // Obtain E and B through the configured production field callbacks, then
  // add q(E + v cross B)/m to the particle acceleration.
}
```

Na photoionization already contains a production Na-to-Na+ dispatch, but it
removes the particle when Na+ is not compiled:

```cpp
if ((spec == _NA_SPEC_) && (_NA_PLUS_SPEC_ >= 0)) {
  PIC::ParticleBuffer::SetI(_NA_PLUS_SPEC_, ParticleData);
  return _PHOTOLYTIC_REACTIONS_PARTICLE_SPECIE_CHANGED_;
}

return _PHOTOLYTIC_REACTIONS_PARTICLE_REMOVED_;
```

### Resolved M0 decision (2026-10-10)

The confirmed M0 choice is to keep the normal baseline `SpeciesList=Na` and
create a committed, dedicated I08 configuration with `SpeciesList=Na,NA_PLUS`.
That test configuration would use analytical uniform E/B fields and the actual
production mover; it would not advertise Na+ as enabled in the normal science
configuration.

## 9. Which impact-vaporization registration controls the M0 baseline?

Section 1 shows that the current input registers the same impact-vaporization
rate and generator twice: once as the named built-in process and once as
`MySource`.  I02 must audit actual dispatch and I03 must verify particle weights
and rates, so they cannot define one physical baseline while this duplication
is unresolved.

### Resolved M0 decision (2026-10-10)

The confirmed M0 choice is to disable `MySource`, retain the named
`Source:ImpactVaporization=on` path, and treat `1.69e22 s^-1` as the one total
Na impact-vaporization rate.  This preserves the production implementation and
unambiguous source ID while eliminating accidental double counting.  The
definition remains in the input because its intended future role is
night-side-adsorbed Na release on the dayside.  That role remains unverified:
the present callbacks alias impact vaporization and contain no conserved
surface transport/release model.

## 10. What physical-time contract should control I11 equilibration and comparison?

### Code evidence

The normal production step advances the SPICE epoch by the mesh-dependent
global particle time step before moving particles (`srcMoon/main_lib.cpp`):

```cpp
const double SimulationTimeStep=
    PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
const SpiceDouble etBegin=Moon::OrbitalMotion::et;
const SpiceDouble etMiddle=etBegin+0.5*SimulationTimeStep;
const SpiceDouble etEnd=etBegin+SimulationTimeStep;
// ... update Sun, Earth, frame, and force state from CSPICE ...
Moon::OrbitalMotion::et=etEnd;
PIC::TimeStep();
```

The application also changes its sampling length in **iteration counts**, not
physical seconds (`srcMoon/main.cpp`):

```cpp
if ((PIC::DataOutputFileNumber!=0) &&
    (PIC::DataOutputFileNumber!=LastDataOutputFileNumber)) {
  PIC::RequiredSampleLength*=2;
  if (PIC::RequiredSampleLength>40000)
    PIC::RequiredSampleLength=40000;
}
```

I11 deliberately changes cell size. The production local time-step callback
therefore changes `GlobalTimeStep[0]`, so equal iteration counts at 40, 20, 10,
and 5 km end at different TDB epochs and sample different illumination,
gravity, frame rotation, and source conditions. A burn-in rule alone cannot
make those outputs a spatial-discretization comparison.

### Resolved M1 decision (2026-10-11)

The approved contract is option 1, frozen-environment steady state. The Moon
receipt now requires `simulation_environment_mode = frozen | evolving`.
`frozen` keeps all SPICE-derived state at `simulation_start_utc` while the
normal particle mover, sources, surface interactions, and particle clocks
continue to advance. I11 uses `frozen`; later transient campaigns may select
`evolving` only with matched physical-time windows.

The alternatives considered were:

Choose one physical contract before I11 science sampling is implemented:

1. **Frozen-environment steady state (recommended for M1):** keep the full
   M0 force/source model, but evaluate SPICE, illumination, and the rotating
   frame at the one input epoch while particle transport advances. Compare
   statistically stationary outputs after the approved 2%-for-three-blocks
   burn-in. This requires an explicit runtime input mode and a linked test that
   proves the environment remains fixed; it is not the legacy no-force mode.
2. **Matched evolving transient:** retain normal epoch advancement, declare an
   exact physical burn-in duration and science-window duration in seconds, and
   make every resolution end/sample at identical TDB times. The required
   durations cannot be derived from the currently qualified observational
   packages and must be supplied as a physical/model decision.

Until one contract is selected, I11 and the sampled portion of I13 remain
`SKIPPED / NOT VERIFIED`; initialization-only or unequal-epoch output will not
be frozen as a reference solution.

## 11. How should the approved 40 km far-field target be limited in the 400-radius tail domain?

### Code evidence

The production domain is not a 10,000 km near-Moon box. It retains a compiled
400-lunar-radius tail extent and 30-radius transverse minimum
(`srcMoon/main_lib.cpp`):

```cpp
const double xMaxDomain=400.0; //modeling of the tail
const double yMaxDomain=30;
```

The generic AMR builder evaluates the volume callback at every block centre
and probe and refines whenever the requested cell size is smaller than the
current dyadic cell size (`src/meshAMR/meshAMRgeneric.h`):

```cpp
requredResolution=localResolution(blockMiddlePoint);
// ... probe the block and retain the smallest requested value ...
c=localResolution(xProbe);
if (c<requredResolution) requredResolution=c;
```

The currently approved M1 receipt declares 40,000 m as the far value. The
implemented radial callback approaches that value at the outer domain, so it
requests approximately 40 km cells throughout an enormous tail volume. The
linked I12 run reached only AMR level 6 (422 km cells) before the cost began
growing rapidly; it still needed several more three-dimensional dyadic levels.
This is not a viable way to test the intended near-Moon 0--10,000 km profile
grid, and calling an unfinished mesh a PASS would be invalid.

### Resolved M1 decision (2026-10-11)

The recommended two-stage profile is approved. The receipt explicitly sets
`mesh_far_field_extent_m = 10000000` and
`mesh_outer_cell_size_m = 17371000`. The production selector reaches 40 km at
10,000 km altitude, then coarsens monotonically to 17,371 km at the existing
400-radius outer boundary. M1 convergence claims are restricted to the common
0--10,000 km grid; tail validation requires its own convergence evidence.

The alternatives considered were:

Recommended: retain the 400-radius domain, but add two explicit input controls:

```text
mesh_far_field_extent_m = 10000000
mesh_outer_cell_size_m = 17371000
```

Interpret the existing 40 km value as the science-region target through
10,000 km altitude, then transition monotonically to the explicit 17,371 km
outer-tail target at the domain boundary. The second value preserves the old
configured outer scale (`10 * lunar radius`) but makes it visible and
adjustable. M1 convergence claims would be restricted to the declared
0--10,000 km common grid; later tail validation would require its own spatial
convergence study and may choose a finer outer value.

Alternatives are to reduce the M1 domain explicitly (which no longer tests the
normal tail-domain wiring) or retain 40 km through all 400 radii (which is not
computationally feasible with this Cartesian AMR topology). No value will be
introduced silently.

## 12. Resolved: what law controls immediate re-emission of non-sticking neutral Na?

### Code evidence

The approved M2 Na table deliberately permits a non-sticking branch above its
coldest node (`srcMoon/examples/lola_surface.in`):

```text
sodium_sticking_model = temperature-table
sodium_accommodation_speed_model = kinetic-energy
sodium_accommodation_coefficient = 0.2
```

The shared production interaction first tests that probability, then combines
a thermal draw with the incident speed using an accommodation coefficient
(`src/models/exosphere/Exosphere.cpp`):

```cpp
if (rnd()<StickingProbability(spec,ReemissionParticleFraction,SurfaceTemp))
  return _PARTICLE_DELETED_ON_THE_FACE_;

vf=sqrt(accommodation*vt*vt+(1.0-accommodation)*vi*vi);
```

Historically the coefficient came from a one-element array
`AccomodationCoefficient[]={0.0}`. Indexing it by any species other than zero
was unsafe, and using zero for M2 Na would preserve incident speed without an
identified physical source. The user identified the local primary source
`/home/vtenishe/1-s2.0-S0019103513003692-main.pdf`: Tenishev et al. (2013),
doi:10.1016/j.icarus.2013.08.021. Its Eq. 14 supplies the kinetic-energy
interpolation, the text evaluates `omega=0.2`, and the scattered direction is
isotropic. Its Fig. 3 supplies the adopted temperature-dependent sticking
curve and identifies Yakshinskiy and Madey (2005) as the underlying source.

### Implemented decision (2026-10-13)

The shared Exosphere path now obtains both the factor and its mathematical
interpretation from the application callback. M2 Na selects the cited
kinetic-energy law and `omega=0.2`; M2 Ar selects its separately declared
full-thermal limit. The generic no-callback path retains the old linear-speed
behavior for other AMPS applications. I18 still must evaluate the linked
stochastic CDF/moments and conservation budget before this capability can be
called verified.

## 16. I17: which eclipse irradiance and thermal spin-up contracts control the comparison?

### Code evidence

The production thermal solver currently reduces Earth eclipse to a binary
infinite-cylinder shadow (`srcMoon/MoonSurfaceTemperature.cpp`):

```cpp
const double axisLength=Norm(earthToSun);
double parallel=Dot(earthToPoint,earthToSun)/axisLength;
bool earthEclipse=false;
if (parallel<=0.0) {
  // ... distance from the Earth--Sun axis ...
  earthEclipse=perpendicular2<
      std::pow(_RADIUS_(_EARTH_)/1000.0,2);
}
const bool visible=terrainVisible && !earthEclipse;
gProduction.columns[element].Advance(timeStep,cosine,visible,
    heliocentricDistanceAu,error);
```

This has no penumbra, finite solar disc, solar limb-darkening prescription, or
Earth-atmosphere refracted/transmitted irradiance. A longer run cannot repair
that forcing law. NASA's LRO eclipse material explicitly reports that sunlight
filtered through Earth's atmosphere gives the eclipsed Moon a red hue, so
treating that contribution as exactly zero is a physical model choice, not an
inference that can be made from the present source.

The I17 runner is deliberately blocked rather than pretending the missing
capabilities are present (`srcMoon/test/integration_common.py`):

```python
"I17": (
    "the calibrated C3-C9/bolometric production operator and "
    "reproducible RDR processor are implemented, but a reviewed "
    "observation-free periodic thermal spin-up and June-frozen "
    "numerical acceptance-threshold artifact are not yet available; "
    "December hold-out residuals must remain unread"),
```

There is consequently no linked I17 production dispatch yet
(`srcMoon/MoonIntegrationTests.cpp`):

```cpp
if (testId=="I14") return RunI14(&receipt,outputFile);
if (testId=="I16") return RunI16(&receipt,outputFile);
if (testId=="I18") return RunI18(&receipt,outputFile);
```

The current thermal columns are initialized at one configured temperature and
are advanced only over the application time interval. No code repeats a
pre-eclipse forcing cycle or proves that the complete subsurface temperature
profile has forgotten that arbitrary initial condition.

### Decision needed

Please approve or replace both parts of this recommended first-comparison
contract:

1. Implement exact apparent-disc overlap using SPICE Sun/Earth distances and
   radii for the direct solar beam, initially treating the Sun as a uniform
   disc and the solid Earth as opaque. Keep Earth-atmosphere
   refraction/transmission explicitly disabled by configuration, record that
   omission in every receipt, and label the resulting plot
   `EXPERIMENTAL / NOT VALIDATED`. Formal I17 validation would remain blocked
   until an authoritative atmospheric-transmission and, if required,
   limb-darkening prescription is selected. Alternatively, provide or select
   those prescriptions now and require them in the first plot.
2. Spin up each production thermal column by replaying an explicitly declared
   pre-eclipse synodic forcing interval until the full column profile converges
   at the configured absolute-temperature and energy tolerances. Put the
   interval and tolerances in the `.in` file, impose no arbitrary cycle-count
   limit, and retain the already approved resource monitor/fatal diagnostic if
   resources are exhausted. The June training event would define and freeze
   the numerical threshold before the December hold-out residuals are read.

Without these choices I can implement file plumbing, but I cannot claim that
the forcing or initialization used by the comparison is physically specified.

Both recommended parts were approved. The implementation now uses analytical
circular apparent-disc overlap, explicit radii, uniform solar radiance, opaque
solid Earth, and an explicitly disabled atmosphere contribution. It replays
the declared mean-synodic pre-training interval until the complete production
column set meets the input temperature and enthalpy tolerances. The June-only
runner rejects non-training rows and creates an experimental comparison plus a
review-required threshold candidate; it cannot promote that artifact to PASS
or inspect the December hold-out before the threshold is separately frozen.

## 2026-11-12 — corrected I04 reference review

Decision requested: approve the corrected Na-chemistry nightly candidate, or
retain the current reference and its evaluated I04 FAIL. Candidate/provenance:
`test_output/srcMoon/m1-budget-audit-20261112T064009Z/reference-review/`.
Two four-rank runs agree byte-for-byte at SHA-256
`4bb811b3e05d72fa708359f1bea0344efd1a9f7e89609d7e8a3ca7eade808cba`.
The existing pinned reference is
`abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.

Relevant former production path (`Exosphere_Chemistry.cpp`):

```cpp
ParentLifeTime=TotalLifeTime(NULL,spec,ptr,PhotolyticReactionAllowedFlag,node);
```

That private table uses `PhotolyticReactions::ModelAvailable`, whose species
switch contains H2O/O2/H2/H/OH/O and excludes Na. The corrected canonical input
now selects the existing Moon lifetime/shadow law through:

```text
ReactionProcessor=Moon::Chemistry::ProcessParticle
LifeTime=Moon::ExospherePhotoionizationLifeTime
```

Linked dispatch with dt=tau*ln(2), tau=130500 s, leaves 2062/4096 sunlit Na
and 4096/4096 shadow Na on eight ranks, within the predeclared statistical
criterion. Its actual-event model-count residual is zero. Mean nightly speed
changes from 2224.937 to 2237.848 m/s. Approval changes reference bytes/hash
and provenance only; numerical tolerances and the I13 freeze remain intact.

The user approved this exact candidate on November 12. It is now installed
with pinned reviewed provenance; the old reference is archived, and the strict
checker produces a zero-byte diff. No further approval for this same reference
installation is needed.

## 2026-11-12 — finite inventory release (approved)

The production desorption rate in `srcMoon/MoonSurfaceTemperature.cpp` is:

```cpp
const double population=sphere->SurfaceElementPopulation[species][surfaceElement];
const double residence=SurfacePhysics::ResidenceTimeS(
    *parameters,gProduction.temperaturesK[surfaceElement]);
return population/residence;
```

The generic injector previously sampled that fixed rate over a whole particle
step and rejected candidates after depletion. With Ar residence time much
shorter than the particle step, this did not integrate the finite reservoir.
The requested choice was continuous finite-inventory depletion with integrated
competing thermal/PSD hazards and explicitly weighted partial packets.
The user approved that choice. The shared production kernel now uses:

```text
released_j = N0 * (-expm1(-sum(k_j)*dt)) * k_j/sum(k_j)
remaining = N0 * exp(-sum(k_j)*dt)
```

`src/models/exosphere/FiniteSurfaceRelease.h` and
`Exosphere_FiniteInventory.cpp` implement it; U12 and the supporting I19 ledger
exercise the production path. Rates are held constant during each split
interval. Fractional packets represent a continuous stock; microscopic
atom-count shot noise is not claimed. Varying-temperature/deposition time-step
convergence is still required. No further approval is needed for this choice.

## 2026-11-12 — D02 illumination calculation epochs (initial question; superseded below)

What authoritative source gives the actual start and end epochs of the
18.6-year hourly illumination calculation used for the 201608 AVGVISIB/LPSR
products? A source document or provider metadata is needed, rather than a
choice of dates for convenience.

Reference processor before the follow-up research,
`srcMoon/test/data_processing/prepare_d02_reference.py`:

```python
"illumination_operator": {
    "reference_cycle_years": 18.6,
    "reference_time_step_hours": 1.0,
    "horizon_azimuth_spacing_deg": 0.5,
    "solar_disc": "extended source; any visible fraction defines visibility",
    "source": "Mazarico et al. 2011; NASA PGDA product 69",
},
```

The native label records:

```text
START_TIME = 2009-07-13T17:33:17.246
STOP_TIME  = 2013-07-18T14:12:13
```

Those dates do not span 18.6 years and cannot silently define the illumination
integration interval. The processor retains the published duration/cadence,
but supplies no invented epoch. Until the actual interval and matching
production extended-disc operator are established, D02 stays PARTIAL and
formal I15/I20 remain SKIPPED / NOT VALIDATED. Independent M2 work can continue.

## 2026-11-12 — D02 research result and I19/I20 proposal

The authorized primary-source research did not establish the complete
generation contract for the exact 201608 files. General PGDA hourly sampling
and older six-hour calculations do not supply version-specific epochs. The
old processor excerpt above is historical: its unsupported timing claims are
now removed from QA. All 366 spatial CSV science rows remain byte-identical.

The remaining factual question is: what authoritative provider record defines
the exact averaging epochs, time system, cadence/endpoints, solar-disc geometry,
observer height and Earth-occultation policy, and what independent timestamped
reference can evaluate I15? Averaged maps cannot reconstruct epoch labels.
Evidence and primary links are in
`srcMoon/test/D02_ILLUMINATION_RESEARCH.md`.

Current code deliberately preserves the gap:

```python
# srcMoon/test/data_processing/prepare_d02_reference.py
"reference_cycle_years": None,
"reference_time_step_hours": None,
"simulation_start_epoch": None,
"simulation_stop_epoch": None,
```

The read-only native audit records the label fields without treating them as
simulation bounds:

```python
# srcMoon/test/data_processing/audit_d02_temporal_contract.py
"native_start_time": label_value(text, "START_TIME"),
"native_stop_time": label_value(text, "STOP_TIME"),
"label_times_used_as_simulation_bounds": False,
```

For the requested concrete implementation decision, review
`srcMoon/test/M2_I19_I20_PROTOCOL_PROPOSAL.md` and its exact JSON. They propose
3x3 independent timestep/weight refinement, 16 fixed independent seeds,
simultaneous confidence intervals and explicit reversible-retention budgets.
No runner consumes these values. The relevant existing production diagnostic
writer, `srcMoon/MoonM2Budget.cpp`, currently writes step/phase/species fields;
formal synchronization additionally needs physical epoch and per-species
elapsed time. `Exosphere_FiniteInventory.cpp` uses the production species
global timestep to integrate each surface depletion interval. Thus a common
physical clock and independently pinned weights must be implemented and
verified before the proposed convergence comparisons are meaningful.

The protocol is a review proposal, not an approval already granted by the
research request. Existing accepted conservation criteria, M1 iteration limits,
formal predecessors and frozen observations remain unchanged. No long formal
I19/I20 run or December residual inspection was launched.


## 2026-11-12 — detailed I19/I20 numerical decisions for review

The requested detailed recommendation is now version 2 of
`srcMoon/test/M2_I19_I20_PROTOCOL_PROPOSAL.md` and its adjacent JSON. It is
unapproved and unconsumed. No active gate has been changed. Review the common
clock, independent dt/weight grid, two physical cases, proposed numerical
budgets, statistical approximation, cohort metrics and staged resource policy
as a single explicitly documented decision package.

Relevant production evidence for the clock decision:

```cpp
// srcMoon/main_lib.cpp:2048; target Ar, source Na, multiplier 0.1
PIC::ParticleWeightTimeStep::copyLocalTimeStepDistribution(
    _AR_SPEC_,_NA_SPEC_,0.1);
```

The retained I19 synchronized-step receipt is 0.6855865301040126 s, rather than
its earlier initialization print of about 6.855865 s. The proposed reporting
horizon therefore means millions of steps, not a small extension of the
20-step fixture. The proposed clock alignment is H/(8*ceil(H/(8*dt_bound))),
with integer checkpoint indices and both species on that clock.

Relevant production evidence for the seed/replicate decision:

```cpp
// srcMoon/MoonIntegrationTests.cpp, current supporting fixture
rnd_seed(731993+PIC::ThisThread);
```

A new formal path must actually accept and verify the selected seed; editing
only a Python manifest does not create independent replicates. The configured
RNG and stream assumptions also need explicit qualification.

Relevant production evidence for the retention interpretation:

```cpp
// srcMoon/MoonSurfacePhysics.cpp
state.thermalRetention = residence >= retentionTimescaleS;
```

This is a timescale flag. For constant first-order thermal loss, residence=H
means exp(-1), about 36.8%, remains after H. The proposed first-release survival,
current surface occupancy after returns, and restricted mean residence must
therefore be reported separately. No irreversible sink is introduced.

Recommendations remain 1% of a fixed atom scale for resolution differences,
one percentage point for fractions, 1% of H for restricted residence, 0.1%
maximum numerical-domain loss in new long cases, and simultaneous confidence
half-widths no greater than half the corresponding margin. Existing exact
count, roundoff and 0.5% conservation limits remain binding. These numerical
budgets are proposed engineering allocations, not observed error bars.

Recommend approving implementation and bounded resource/precision planning
first. The complete two-case 3x3 grid with sixteen seeds means 288 long runs;
its launch remains conditional on measured cost, statistical support and all
formal predecessors. No data gap, physics coefficient or acceptance change is
silently resolved by preparing this review package.


## 2026-11-12 — development hour / deferred cluster execution (user-directed)

The complete development test command must fit within one hour, including
required builds, plots and cleanup. Long tests may be deferred to a cluster
for a few hundred CPU-core hours. The v3 proposal documents this operational
requirement; scientific thresholds remain unapproved. The proposed concrete
cluster planning target is 300 aggregate core-hours, adjustable before launch.

Current code shows why a new aggregate guard is still needed:

```python
# srcMoon/test/run_m2_development.py, inside the I14/I18/I19 loop
codes=execute_levels({test_id:{"directory":directory,"command":cmd}},
    monitor,1,a.test_runtime_seconds,render_previews=False)
```

Each call starts a fresh scheduler deadline. Configuration/builds, local tests,
QA, plots and restoration are outside that child guard. The next implementation
must use one monotonic deadline across all phases, propagate remaining time,
reap MPI/build groups and reserve reporting/cleanup time. No shell-timeout
wrapper or per-child one-hour cap is claimed to satisfy the requirement.

The candidate short refinement workload is two physical cases x nine grid
cells x two seeds, with 24/48/96 steps over the same short physical horizon.
These are development diagnostics, not substitutes for long kinetic/thermal
convergence or the proposed sixteen-seed qualification. Deferred cluster tests
remain SKIPPED / NOT QUALIFIED. Source execution policy is amended; no model
code, active acceptance criterion or external observation was changed.

## November 16 M2 completion decisions and code evidence

These decisions apply to the working tree based on
`17d49e00d1d2ccf33b06d11a3e6e5ffef34ac810`. Earlier excerpts above remain
historical evidence, not a description of the revised implementation.

### D02 producer metadata — no additional information available

Question: are authoritative calculation epochs/time system, cadence/endpoints,
observer height, solar/ephemeris/Earth-occultation conventions and an independent
timestamped illumination benchmark available for the frozen 201608 products?
The user answered **no additional metadata available**. These values cannot
be inferred from a directory or acquisition date. The native operator exposes
separate quantities (`MoonSurfaceEnvironment.h`, `SolarDiscSample`):

```cpp
double terrainVisibleAreaFraction=0.0;
double earthVisibleAreaFraction=0.0;
double jointlyVisibleAreaFraction=0.0;
```

The formal comparison must bind those meanings and actual calculation epochs
to the independent reference contract. Geometry-only production checks are
permitted; formal D02 comparison remains SKIPPED until that evidence exists.

### December I17 initialization — independent explicit spin-up approved

Question: should December carry the June state forward, or initialize by the
same observation-free periodic method with a separately explicit interval?
The user approved **independent December spin-up**. The runtime inputs include:

```text
diviner_eclipse_selection = holdout
thermal_spinup_cycle_start_utc = <explicit investigator-selected UTC>
thermal_spinup_cycle_end_utc = <explicit investigator-selected UTC>
```

The placeholders above are explanatory syntax, not runnable epochs or an
approved interval. `RunI17` and `i17_holdout.py` require selected-only history
staging, a reviewed hash-bound freeze and native authorization before reading
the hold-out file. Approval of the initialization method does not supply the
missing December bounds or observational acceptance criteria.

### Diviner brightness versus kinetic temperature — grey facet approved

The audit found this pre-correction mismatch in
`MoonSurfaceTemperature.cpp`:

```cpp
// ProductionTemperatureK: the table's brightness value was returned directly.
return temperature;
// ProductionThermalDesorptionRate: the unchanged initialization vector was used.
const double residence=SurfacePhysics::ResidenceTimeS(
    *parameters,gProduction.temperaturesK[surfaceElement]);
```

Question: explicitly interpret the map as an isothermal grey facet, or keep it
only as brightness diagnostics and use predictive thermal physics for kinetics?
The user approved **explicit grey-facet approximation**. The new input is:

```text
surface_temperature_mode = diviner-map
diviner_map_kinetic_interpretation = isothermal-grey
```

The configured `thermal_infrared_emissivity` supplies epsilon. Production now
uses `GreyKineticTemperatureK(T_bol, epsilon)` and refreshes logical source
temperatures at their explicit epoch. I16 converts kinetic temperature through
the matching forward brightness operator before comparison. Neither archive
values nor uncertainties are relabelled as measured kinetic temperature.
Anisothermality, emissivity uncertainty and unresolved footprint structure
remain scientific approximation/qualification limits. Missing map support
remains unavailable; it is never replaced by an initialization temperature.

## M3 next-package review — 2026-11-17

The user approved early M3 development, explicit Maxwellian/Voronov atomic
rates, and literature proposals BEFORE unresolved charge-exchange/molecular
physics choices are implemented. The atomic rate-service work has its own
verification; the following proposed next scope is not automatically approved
and cannot qualify M3. Full research and primary references are in
`srcMoon/test/M3_CHARGE_EXCHANGE_PROPOSAL.md` and
`srcMoon/test/M3_MOLECULAR_KINETICS_PROPOSAL.md`.

### Q-M3-CX: initial charge-exchange scope

**Recommendation for review:** implement independently verified total-capture
rate kernels for ground-state H+ + He and H+ + Na(3s), strictly within the
reviewed published domains, while particle reactions/daughter transport remain
disabled. Use original ORNL A-32 coefficients for He; acquire/QA the NIST Na
table and freeze its interpolation rule before code. Do not extrapolate Na
below 1 keV/u. Preserve source physical uncertainties separately from numerical
implementation accuracy. Ne/Ar/alpha channels remain explicitly unsupported.
A reduced neutral-loss bath is NOT automatically selected by this proposal.

The current generic library is not a ready lunar model:

```cpp
// src/models/charge_exchange/ChargeExchange_LindsayStebbings.cpp
if (spec != _H_SPEC_ && spec != _H_ENA_V1_SPEC_ &&
    spec != _H_ENA_V2_SPEC_ && spec != _H_ENA_V3_SPEC_)
  return 1E+100;
```

The actual lunar U17 contract confirms the gap:

```json
"implementation": "not-implemented",
"skip_reason": "No Moon charge-exchange channels or authoritative channel table are configured."
```

**Decision:** approve this bounded rate-kernel next package, or wait until
complete joint product-state/angular evidence is available before any new
charge-exchange implementation. Neither choice authorizes invented cross
sections or promotes a total-rate test to I24 product/conservation PASS.

### Q-M3-MOL: molecular rates versus physical daughter transport

**Recommendation for review:** implement the complete frozen H2O/OH photon
rate/channel/state inventory and H/O-nucleus/charge bookkeeping first,
retaining both quiet/active endpoints explicitly. Keep daughter transport
DISABLED / NOT QUALIFIED until a reviewed joint distribution/state/decay
contract exists. Do not drop unsupported branches or renormalize their rates;
parent-survival/channel counting is a separately labelled rate verification,
not a transport simulation or observational result.

Relevant current product mapping:

```cpp
// src/models/photolytic_reactions/PhotolyticReactions_Huebner1992ASS.cpp
int PhotolyticReactions::H2O::Huebner1992ASS::ReactionProducts
    [nReactionChannels][nMaxReactionProducts]={
    {_OH_SPEC_,_H_SPEC_,-1},
    {_O_SPEC_,_H2_SPEC_,-1},
    {_O_SPEC_,_H_SPEC_,_H_SPEC_},
    // ... ionization branches ...
};
```

This maps excited and ground-state oxygen to the same `_O_SPEC_`, while frozen
D04 distinguishes O(3P), O(1D), O(1S). Current generic kinematics also substitutes
mean channel energy and isotropic directions for distributions, and products
and parent loss are sampled separately. The proposal includes exact code
locations and authoritative distribution/state references. Experimental-OH
sub-threshold wavelength tails need source clarification; raw values must
remain unchanged, not clipped into apparent energy consistency.

**Decision:** approve the bounded rate/topology/bookkeeping next package, or
wait for the full daughter-distribution/state evidence before molecular code.
This decision supplies no absent energy, angular PDF, lifetime or threshold.
U24/I22/I31 transport and the formal M3 gate remain unqualified either way.



## M3 completion clarifications and approvals — 2026-11-18

The following choices were asked before new qualification execution. All three
were approved in the session; approvals supply no missing source values.

1. **Numerical/statistical protocol.** Original code/data contract:
   `srcMoon/test/M3_QUALIFICATION_PROTOCOL.json`:

   ```json
   "additional_statistical_design": "PROPOSED_REQUIRES_PREEXECUTION_REVIEW",
   "numerical_error_budget": null
   ```

   Approved: numerical bias budgets0.002 forI21/I22 and0.003 forI23/I24,
   endpoint-survival primary estimators, separate branch/partial-rate diagnostics,
   simultaneous95% exact branch intervals with10% relative half-width for every
   positive source branch, original-parent keyed rank-independent streams, and
   cluster deferral without dropped branches when development resources do not
   suffice. The immutable original is retained; the explicit approved addendum
   is `srcMoon/test/M3_QUALIFICATION_ADDENDUM_20261118.json`.

2. **Moving-parent photon composition.** Existing
   `srcMoon/species/ReducedKinematics.cpp` supplies `PhotonTwoBody`, while
   `srcMoon/MoonReactionController.cpp` passes the actual event parent:

   ```cpp
   context.parent.velocityMS[k]=velocity[k];
   channel.kinematics(context,&policy,&error);
   ```

   Approved: extend the source-mean/isotropic experimental closure with a
   consistently transformed incoming photon and complete electron/radiation
   ledgers. An explicit source-threshold/state crosswalk and independently
   declared model-error bounds remain required. A simple velocity boost plus
   a fitted residual bath is not an admissible substitute.

3. **Missing physical inputs.** Frozen approved decisionP06 explicitly says:

   ```text
   retain selected scalar Maxwellian/Voronov rates and keep impact daughter
   kinetics disabled until source-consistent joint data/model exist
   ```

   Approved: research authoritative sources while completing the software;
   preserve strict data gates and disable unsupported channels. D04 currently
   provides no molecular threshold table; its OH_exp label is not a resolved
   rovibrational population. Neither this approval nor a successful API test
   creates those absent inputs or a universal ion-regolith material law.
