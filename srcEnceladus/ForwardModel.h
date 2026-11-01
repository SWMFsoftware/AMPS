#ifndef SRCENCELADUS_FORWARD_MODEL_H
#define SRCENCELADUS_FORWARD_MODEL_H

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace Enceladus {
namespace ForwardModel {

// All public quantities use SI units.  The core in this file has deliberately
// no PIC or MPI dependency so parser, source-normalization, reaction-rate, and
// numerical-planning tests exercise independent application logic.
struct Status {
  bool success = true;
  std::string message;

  static Status Ok() { return Status{}; }
  static Status Error(const std::string& text) {
    Status result;
    result.success = false;
    result.message = text;
    return result;
  }
};

enum class SourceMode { Shape, Absolute };
enum class SourceGeometry { PointJet, PolylineStripe };
enum class StripeSurfaceMapping { LocalVerticalToReferenceSphere };
enum class SolarActivity { Quiet, Active };
enum class PopulationControl { Off, VelocityShift, Scatter };
enum class ReactionKind { Photo, ChargeExchange };
// Schema 2 implements parent removal only.  Product tracking is deliberately
// not an accepted runtime mode until daughter transport is re-authorized and
// independently qualified.
enum class ReactionProducts { Discard };
enum class PlasmaBackgroundMode { None, CacheOnly, CacheOrIces };

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct Source {
  std::string id;
  SourceGeometry geometry = SourceGeometry::PointJet;
  double temperatureK = 0.0;
  double bulkSpeedMPerS = 0.0;

  // Point-jet geometry.  The conventions are signed planetocentric latitude,
  // west-positive longitude, tilt from the outward radial normal, and azimuth
  // clockwise from local north in IAU_ENCELADUS.
  double latitudeDeg = 0.0;
  double westLongitudeDeg = 0.0;
  double tiltDeg = 0.0;
  double azimuthDeg = 0.0;

  // Polyline sources are the published IAU_ENCELADUS Cartesian nodes [m].
  // The tenishev2014_hybrid_8 law represents the line by the same
  // trapezoidal arc-length node weights used by the standalone fit.
  std::vector<Vec3> verticesM;

  // Shape mode: dimensionless fractions. Absolute mode: molecules s^-1.
  std::map<std::string, double> strength;
};

struct SourceDistribution {
  int schemaVersion = 0;
  std::string name;
  std::string sourceCase;
  std::string frame;
  std::string sourceLaw;
  StripeSurfaceMapping stripeSurfaceMapping =
      StripeSurfaceMapping::LocalVerticalToReferenceSphere;
  double stripeReferenceRadiusM = 0.0;
  SourceMode mode = SourceMode::Shape;
  std::vector<Source> sources;
  std::string inputFile;
};

struct Configuration {
  int schemaVersion = 0;
  std::string inputFile;
  std::string sourceDistributionFile;
  std::vector<std::string> simulatedSpecies;
  std::vector<std::string> injectedSpecies;
  std::vector<std::string> chemistrySpecies;
  std::vector<std::string> selectedReactions;
  std::map<std::string, double> totalSourceRatePerS;
  std::map<std::string, double> characteristicSpeedMPerS;
  std::map<std::string, double> courantNumber;
  std::map<std::string, std::uint64_t> primaryMacroparticlesPerStep;
  std::map<std::string, double> plasmaIonFraction;
  std::uint64_t randomSeed = 0;
  long sampleLengthSteps = 0;
  long totalIterations = 0;
  // Optional observation operator evaluated from AMPS' completed sampling
  // buffer at every normal output.  The trajectory is a frozen SPICE CSV in
  // IAU_ENCELADUS; its path is resolved relative to this run deck.
  bool flybyDensityOutput = false;
  std::string flybyId;
  std::string flybyTrajectoryFile;
  PopulationControl populationControl = PopulationControl::Off;
  ReactionProducts reactionProducts = ReactionProducts::Discard;
  // Generic PIC SplitWithVelocityShift perturbs position by this fraction of
  // the local cell diagonal and velocity by this fraction of the cell's
  // measured thermal speed. These are dimensionless API parameters, not SI
  // distances or speeds; naming them with physical units would make a run
  // configuration quantitatively wrong as the mesh or population changes.
  double splitPositionShiftCellDiagonalFraction = 0.0;
  double splitVelocityShiftThermalSpeedFraction = 0.0;
  int minParticlesPerCell = 0;
  int maxParticlesPerCell = 0;
  bool nonuniformSplitPositionShift = false;
  SolarActivity solarActivity = SolarActivity::Quiet;
  PlasmaBackgroundMode plasmaBackgroundMode = PlasmaBackgroundMode::None;
  std::string plasmaCacheDirectory;
  std::string icesInstallationDirectory;
  std::string icesCase;

  // These two values are used only to plan fixed daughter statistical
  // weights before transport. Runtime charge-exchange frequencies always use
  // the local coupled density and relative velocity instead.
  double chemistryWeightReferenceNumberDensityM3 = 0.0;
  double chemistryWeightReferenceRelativeSpeedMPerS = 0.0;
};

struct Product {
  std::string species;
  int multiplicity = 1;
};

struct Reaction {
  std::string id;
  ReactionKind kind = ReactionKind::Photo;
  std::string parent;
  std::string plasmaReactant;
  std::vector<Product> products;
  double quietRatePerS = 0.0;
  double activeRatePerS = 0.0;
  double excessEnergyQuietJ = 0.0;
  double excessEnergyActiveJ = 0.0;
  double crossSectionAt15KmS_M2 = 0.0;
  double crossSectionAt26KmS_M2 = 0.0;
};

struct NumericalPlan {
  double globalTimeStepS = 0.0;
  std::string limitingSpecies;
  std::map<std::string, double> physicalSourceRatePerS;
  std::map<std::string, double> expectedMacroparticlesPerStep;
  std::map<std::string, double> statisticalWeight;
  std::map<std::string, double> reactionRatePerS;
};

Status ParseConfiguration(const std::string& file, Configuration* output);
Status ParseSourceDistribution(const std::string& file,
                               SourceDistribution* output);
Status ValidateConfiguration(const Configuration& configuration,
                             const SourceDistribution& distribution);

const std::vector<Reaction>& Burger2007Reactions();
const Reaction* FindReaction(const std::string& id);

// Photo rates are tabulated by Burger et al. (2007) at 9.5 AU. Ion-neutral
// rates use nu=n_i*sigma(v_rel)*v_rel. Burger tabulates sigma only at 15 and
// 26 km/s; values outside that measured bracket are rejected, never silently
// extrapolated.
Status EvaluateReactionRate(const Reaction& reaction,
                            SolarActivity activity,
                            double totalPlasmaNumberDensityM3,
                            double plasmaReactantFraction,
                            double relativeSpeedMPerS,
                            double* ratePerS);

Status BuildNumericalPlan(const Configuration& configuration,
                          const SourceDistribution& distribution,
                          double minimumAllocatedCellSizeM,
                          NumericalPlan* output);

Vec3 PointJetPositionM(const Source& source, double radiusM);
Vec3 PointJetNormal(const Source& source);
Vec3 PointJetAxis(const Source& source);
Status StripeNodeWeights(const Source& source, std::vector<double>* weights);

// Dimensionless flux moment and conditional speed quantile for the
// tenishev2010_point source law.  With t=v/sqrt(2 k_B T/m), s=u/v_th, and
// c=v_hat.axis, emitted molecules have joint event density
//
//   p(t,v_hat) proportional to t^3 exp(-t^2+2 s c t-s^2)
//
// on the outward surface hemisphere.  The moment is its integral over t.
// Keeping these routines PIC-independent lets the production sampler be
// checked against independent quadrature rather than against itself.
Status Tenishev2010FluxMoment3(double directionCosine, double speedRatio,
                              double* moment);
Status Tenishev2010SpeedQuantile(double directionCosine, double speedRatio,
                                double probability, double* speedOverThermal);

const char* Name(SourceMode mode);
const char* Name(SourceGeometry geometry);
const char* Name(StripeSurfaceMapping mapping);
const char* Name(SolarActivity activity);
const char* Name(PopulationControl control);
const char* Name(ReactionProducts products);
const char* Name(PlasmaBackgroundMode mode);

}  // namespace ForwardModel
}  // namespace Enceladus

#endif
