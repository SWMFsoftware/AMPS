#ifndef SRCMOON_MOON_INPUT_H
#define SRCMOON_MOON_INPUT_H

#include <string>
#include <vector>

#include "MoonThermalModel.h"

namespace Moon {
namespace Runtime {

// Selects the geometric boundary used by the particle mover.  Sphere keeps
// the historical cInternalSphericalData boundary; Lola replaces only the
// collision geometry with a generated triangulation.  A logical spherical
// grid is retained in both modes for the existing source and inventory code.
enum class SurfaceGeometry {
  Sphere,
  Lola
};

// All lengths are SI.  Paths are absolute after parsing; relative paths are
// resolved against the directory containing the .in file, never against the
// process working directory.
struct Configuration {
  SurfaceGeometry surfaceGeometry = SurfaceGeometry::Sphere;
  std::string inputFile;             // Absolute path of the parsed receipt.
  std::string spiceRoot;             // Root containing cspice/ and Kernels/.
  std::string spiceToolkitDirectory; // <spiceRoot>/cspice.
  std::string spiceKernelDirectory;  // <spiceRoot>/Kernels.
  std::string lolaProductId;         // Must equal PRODUCT_ID in the PDS label.
  std::string lolaImageFile;         // Detached binary PDS image.
  std::string lolaLabelFile;         // Authoritative detached PDS label.
  std::string surfaceCeaFile;        // Generated AMPS CEA-long mesh.
  std::string surfaceTecplotFile;    // Generated Tecplot FETRIANGLE mesh.
  double surfaceMeshResolutionM = 0.0; // Maximum pre-topography arc length.

  // M1 volume-AMR controls. ``meshProfile`` is a provenance label, not a
  // source of compiled defaults: every dimensional value below remains
  // mandatory in QUICK, FULL, and CUSTOM receipts. QUICK and FULL are the
  // committed example campaigns; CUSTOM makes an investigator-chosen set of
  // explicit values machine-readable without pretending it is either one.
  // The callback holds the requested near-surface size through
  // meshNearSurfaceExtentM, then interpolates monotonically in log(radius)
  // toward meshFarFieldCellSizeM at meshFarFieldExtentM altitude.  Beyond
  // that science region it coarsens again toward meshOuterCellSizeM at the
  // compiled AMPS outer radial boundary.  Keeping the two extents distinct
  // prevents the M1 40 km science target from being applied throughout the
  // 400-lunar-radius legacy tail box, which is computationally intractable.
  std::string meshProfile;
  double meshNearSurfaceCellSizeM = 0.0;
  double meshNearSurfaceExtentM = 0.0;
  double meshFarFieldCellSizeM = 0.0;
  double meshFarFieldExtentM = 0.0;
  double meshOuterCellSizeM = 0.0;

  // Optional full-volume visualization output.  This controls the diagnostic
  // Tecplot exports ``mesh.dat`` and ``tetra_mesh*.dat`` only; it does not
  // disable the compact ``mesh.msh`` AMR-tree exchange file required by the
  // existing multi-rank initialization path.  A fine global/tail AMR tree can
  // expand to hundreds of gigabytes when serialized as text, so configured
  // runs must request this output explicitly instead of producing it as an
  // undocumented side effect of mesh construction.
  bool meshOutputTecplot = false;

  // Boundary-adjacent refinement is prescribed by an explicit piecewise-
  // linear table in solar zenith angle (degrees, metres).  A separate
  // terminator band applies the table's value at 90 degrees, while the polar
  // cap is evaluated in the label-native MOON_ME_DE421 body frame.  Taking
  // the minimum of applicable targets preserves the finest requested scale
  // where physical regions overlap.
  std::vector<double> meshAngularSzaDeg;
  std::vector<double> meshAngularResolutionM;
  double meshTerminatorHalfWidthDeg = 0.0;
  double meshPolarCapStartLatitudeDeg = 0.0;
  double meshPolarResolutionM = 0.0;

  // The M1 comparison grid is expressed as physical altitude above the AMPS
  // lunar reference sphere (_RADIUS_(_TARGET_)), in metres.  It is not height
  // above the 1737.4 km LOLA datum or above local terrain.  The first regions
  // are linear and share their boundary point.  The far region is logarithmic
  // and begins strictly above profileAltitudeMiddleMaximumM.  The parser
  // materializes profileAltitudeGridM once so every production diagnostic
  // consumes exactly the same grid rather than independently reconstructing
  // it with possibly different endpoint/round-off rules.
  double profileAltitudeNearMaximumM = 0.0;
  double profileAltitudeNearStepM = 0.0;
  double profileAltitudeMiddleMaximumM = 0.0;
  double profileAltitudeMiddleStepM = 0.0;
  double profileAltitudeFarMaximumM = 0.0;
  int profileAltitudeFarPointsPerDecade = 0;
  std::vector<double> profileAltitudeGridM;

  // Solar-zenith-angle profiles are defined in the LSO equatorial plane.
  // Body-fixed pole profiles use +/-Z in MOON_ME_DE421; longitude is
  // undefined at a pole and is therefore deliberately not an input value.
  std::vector<double> profileSolarZenithAnglesDeg;
  bool profileIncludeBodyFixedPoles = false;
  bool profileIncludeShellAverage = false;
  // A shell average is an integral over direction, so the only scientifically
  // approved weighting is physical solid angle.  This explicit receipt field
  // prevents a future implementation from silently averaging AMR cells or
  // surface facets, whose nonuniform areas would bias the result.
  std::string profileShellAngularWeighting;

  // This UTC is authoritative for the complete orbit-enabled run, including
  // dynamics and diagnostic geometry.  A separate observation-only epoch
  // would make model state and line-of-sight geometry physically inconsistent.
  // CSPICE performs the authoritative UTC-to-ET interpretation after the
  // configured leap-second and ephemeris kernels have been furnished.
  std::string simulationStartUtc;

  // ``frozen`` holds all ephemeris-driven environment quantities at
  // simulationStartUtc while PIC particle time continues to advance.  This
  // supplies the stationary physical problem required for spatial AMR
  // convergence. ``evolving`` advances the SPICE epoch once per PIC step and
  // remains available for later transient/observation campaigns.
  std::string simulationEnvironmentMode;

  // Random-number control for reproducible production particle histories.
  // This is an input receipt rather than a test-runner override: sources,
  // surface interactions, and the mover all consume AMPS's single production
  // RNG stream after amps_init() installs this seed.
  int randomSeed = 0;

  // Near-Moon limb geometry.  The currently approved M1 contract uses an
  // Earth observer and the Sun-Moon-Earth plane, with the receipt requesting
  // both apparent-limb intersections or either named branch.  String selectors
  // ensure unsupported future choices fail explicitly instead of silently
  // falling back to one geometry.
  std::string limbObserver;
  std::string limbPlane;
  std::string limbSides;
  // The two branches are named by the projection of the Moon-to-Sun direction
  // into the observer's sky plane.  This avoids the ambiguous legacy words
  // "east/west" and remains well defined except at the explicitly detectable
  // zero-projection degeneracy.
  std::string limbBranchReference;
  double limbTangentMaximumAltitudeM = 0.0;

  // Extended-tail map geometry.  ``anti-lunar`` means the sky direction from
  // the Earth opposite the instantaneous Earth-to-Moon vector; it is not an
  // alias for anti-solar.  The beam diameter describes a circular aperture,
  // not map-pixel spacing.  Until aperture convolution is implemented, the
  // value is configuration/provenance only and cannot establish U26/I33 PASS.
  std::string tailCenterDirection;
  // D14 tabulates native WHAM pointings in J2000 RA/Dec and also provides
  // independently transformed ECLIPJ2000 longitude/latitude.  Either receipt
  // can therefore be selected without manufacturing observation geometry.
  // Both use positive-east longitude/RA and positive-north latitude/Dec in
  // model storage; display orientation is a separate visualization choice.
  std::string tailMapCoordinateSystem;
  std::string tailMapDisplayOrientation;
  double tailMapHalfWidthDeg = 0.0;
  double tailBeamDiameterDeg = 0.0;
  double tailExposureDurationS = 0.0;
  // The publication establishes a circular 1-degree field but not a radial
  // response curve.  Keep that absence machine-readable instead of assuming
  // a Gaussian or uniform disk.  The production operator must remain skipped
  // while this is "unresolved".
  std::string tailBeamResponse;
  // An observational comparison must average the evolving model over the
  // finite exposure, rather than evaluating only one instantaneous epoch.
  // This field records that required operator even though the full D14
  // observation operator is not implemented yet.
  std::string tailTemporalOperator;

  // M2 surface/environment selection.  ``disabled`` preserves the reviewed
  // M0/M1 problem exactly.  ``enabled`` requires every field in the M2 group;
  // there are no compiled fallback constants for a science run.
  bool m2EnvironmentEnabled = false;
  std::string surfaceTemperatureMode;
  std::string divinerTemperatureTableFile;
  // I17 time-resolved Diviner evidence is separate from the climatological
  // I16 interpolation table.  The history contains binned PDS4 RDR radiances,
  // native brightness temperatures, actual footprint centroids/times, flags,
  // and propagated uncertainties.  The JSON contract freezes selection and
  // binning before hold-out scoring.  ``none`` keeps non-I17 production runs
  // independent of observational files; the linked I17 runner replaces both
  // values in its preserved application receipt.
  std::string divinerEclipseHistoryFile;
  std::string divinerEclipseSelectionContractFile;
  double illuminationRayOriginOffsetM = 0.0;
  // Moving LOLA collision controls.  The segment length controls conservative
  // bracketing (metres in LSO), while time/distance tolerances and the
  // iteration limit control root refinement.  I14 varies these explicit
  // values to establish numerical convergence.
  double surfaceCollisionMaximumSegmentLengthM = 0.0;
  double surfaceCollisionTimeToleranceS = 0.0;
  double surfaceCollisionDistanceToleranceM = 0.0;
  int surfaceCollisionMaximumIterations = 0;
  Thermal::ColumnConfiguration thermal;
  double thermalTimeStepS = 0.0;
  double thermalInitialTemperatureK = 0.0;

  // I17 illumination and initialization controls.  These fields are part of
  // the science receipt because a binary point-Sun shadow and an arbitrarily
  // initialized regolith column cannot reproduce an eclipse cooling curve.
  // The approved first-comparison model integrates the unocculted area of a
  // uniformly bright apparent solar disc behind an opaque solid Earth.  It
  // deliberately contains neither solar limb darkening nor refracted or
  // transmitted terrestrial-atmosphere light; both omissions are explicit
  // model limitations and keep I17 EXPERIMENTAL / NOT VALIDATED.
  std::string thermalEclipseDirectIrradianceModel;
  std::string thermalSolarDiscRadianceModel;
  std::string thermalEarthAtmosphereTransmissionModel;
  double thermalSolarDiscRadiusM = 0.0;
  double thermalOpaqueEarthRadiusM = 0.0;

  // A converged periodic initial condition is obtained without consulting a
  // Diviner temperature.  The declared UTC interval is replayed until both
  // the largest node-temperature change [K] and largest whole-column
  // enthalpy change [J m^-2] between successive cycle boundaries meet these
  // tolerances.  There is intentionally no maximum cycle count: allocation
  // failure and external job-resource termination remain ERROR conditions.
  std::string thermalSpinupCycleStartUtc;
  std::string thermalSpinupCycleEndUtc;
  double thermalSpinupTemperatureToleranceK = 0.0;
  double thermalSpinupEnergyToleranceJm2 = 0.0;

  // Species surface physics.  Energy units are stated in field names because
  // Na laboratory work reports per-particle eV whereas the approved Ar model
  // reports molar kJ/mol.  Input parsing converts Na eV to joules only at the
  // production adapter boundary.
  double sodiumActivationEnergyEv = 0.0;
  double sodiumAttemptFrequencyHz = 0.0;
  // Sticking and immediate non-sticking accommodation are independent
  // choices.  A temperature-table sticking model controls the Bernoulli
  // adsorption event, while the accommodation model controls only the speed
  // of an atom that was not adsorbed.  Keeping both selections in the run
  // receipt prevents a literature-specific Na law from becoming a hidden
  // compiled default for every simulated species.
  std::string sodiumStickingModel;
  std::vector<double> sodiumStickingTemperatureK;
  std::vector<double> sodiumStickingProbability;
  std::string sodiumAccommodationSpeedModel;
  double sodiumAccommodationCoefficient = 0.0;
  std::string sodiumSurfaceInteractionReference;
  double argonDesorptionEnergyJmol = 0.0;
  std::string argonStickingModel;
  std::vector<double> argonStickingTemperatureK;
  std::vector<double> argonStickingProbability;
  std::string argonAccommodationSpeedModel;
  double argonAccommodationCoefficient = 0.0;
  std::string argonSurfaceInteractionReference;

  // Radiogenic 40Ar source receipt.  The total rate is atoms s^-1 over the
  // complete Moon, not a per-area flux.  ``uniform-global`` is the only
  // presently qualified spatial kernel.  A potassium-weighted option is not
  // accepted until D10 supplies a global, QA-clean map; the current regional
  // SPA product must never be extrapolated or gap-filled silently.
  bool argonSourceEnabled = false;
  std::string argonSourceSpatialModel;
  double argonSourceRateAtomsS = 0.0;
  std::string argonSourceTimeModel;
  std::string argonSourceMapFile;
  std::string argonSourceReference;
  std::string ionSurfaceInteractionMode;

  // A trap remains a finite-residence inventory.  This declared comparison
  // timescale controls only the thermal-retention diagnostic and never deletes
  // particles irreversibly.  A zero initial inventory is used unless a
  // versioned checkpoint is explicitly selected.
  double coldTrapRetentionTimescaleS = 0.0;
  std::string surfaceInventoryInitialState;
  std::string surfaceInventoryCheckpointFile;
  double surfaceInventoryCheckpointIntervalS = 0.0;
};

// Parse the unique, flat ``#begin application: moon`` block through AMPS's
// shared runtime-input grammar. Unknown root blocks, nested Moon blocks, and
// legacy ``#section`` directives are fatal so no input is silently ignored.
// The operation is transactional: result is unchanged when false is returned.
bool ParseApplicationInput(const std::string& path, Configuration* result,
                           std::string* error);

// Install the immutable process configuration before amps_init().  Installation
// is intentionally one-shot: silently replacing paths or geometry after MPI
// and mesh initialization would make the run receipt disagree with the model.
bool InstallConfiguration(const Configuration& configuration,
                          std::string* error);
bool HasConfiguration();
// When no application input was provided, return a read-only default object.
// Callers must test HasConfiguration() before interpreting its empty paths.
const Configuration& GetConfiguration();

const char* SurfaceGeometryName(SurfaceGeometry geometry);

// Production M1 mesh kernels.  The angular evaluator consumes an already
// transformed absolute MOON_ME_DE421 latitude so frame transformation remains
// the caller's responsibility and can be checked independently by I12.
double EvaluateAngularMeshResolutionM(const Configuration& configuration,
                                      double solarZenithAngleDeg,
                                      double absoluteBodyLatitudeDeg);
double EvaluateVolumeMeshResolutionM(const Configuration& configuration,
                                     double radiusM,
                                     double referenceRadiusM,
                                     double outerRadiusM);

}  // namespace Runtime
}  // namespace Moon

#endif  // SRCMOON_MOON_INPUT_H
