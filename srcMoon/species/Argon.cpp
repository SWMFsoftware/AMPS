#include "pic.h"

#include "Argon.h"
#include "Moon.h"
#include "MoonInput.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>

namespace Moon {
namespace Species {
namespace Argon {
namespace {

// The AMPS user-defined source dispatcher calls these production functions
// after Config.pl assigns the source a normal Exosphere source ID.  Keeping
// the distribution here prevents generic Exosphere code from acquiring lunar
// 40Ar assumptions.  It is configured once after the logical surface exists.
bool gConfigured = false;
bool gEnabled = false;
double gTotalRateAtomsS = 0.0;
cInternalSphericalData* gLogicalSurface = nullptr;
std::vector<double> gElementRatesAtomsS;
cSingleVariableDiscreteDistribution<int> gSurfaceDistribution;

bool PositiveFinite(double value) {
  return std::isfinite(value) && value > 0.0;
}

double ElementRate(int element, int*) {
  if (element < 0 || static_cast<std::size_t>(element) >=
          gElementRatesAtomsS.size()) {
    exit(__LINE__, __FILE__,
         "radiogenic 40Ar source requested an invalid surface element");
  }
  return gElementRatesAtomsS[static_cast<std::size_t>(element)];
}

}  // namespace

bool ValidateSourceParameters(const SourceParameters& parameters,
                              std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  if (parameters.spatialModel != "uniform-global") {
    *error = "argon source spatial model must be uniform-global; "
             "the regional, QA-incomplete D10 product cannot define a "
             "global potassium-weighted production source";
    return false;
  }
  if (parameters.timeModel != "continuous") {
    *error = "enabled radiogenic 40Ar source currently requires the "
             "literature-supported continuous time model";
    return false;
  }
  if (!PositiveFinite(parameters.totalRateAtomsS)) {
    *error = "enabled radiogenic 40Ar source rate must be finite and positive";
    return false;
  }
  if (parameters.reference.empty()) {
    *error = "enabled radiogenic 40Ar source requires a non-empty reference";
    return false;
  }
  return true;
}

bool BuildUniformElementRates(const std::vector<double>& elementAreasM2,
                              double totalRateAtomsS,
                              std::vector<double>* elementRatesAtomsS,
                              std::string* error) {
  if (elementRatesAtomsS == nullptr || error == nullptr) return false;
  error->clear();
  if (!PositiveFinite(totalRateAtomsS) || elementAreasM2.empty()) {
    *error = "uniform 40Ar source requires a positive rate and nonempty area "
             "table";
    return false;
  }
  long double totalAreaM2 = 0.0L;
  for (double areaM2 : elementAreasM2) {
    if (!PositiveFinite(areaM2)) {
      *error = "uniform 40Ar source encountered a non-positive or non-finite "
               "surface-element area";
      return false;
    }
    totalAreaM2 += static_cast<long double>(areaM2);
  }
  if (!(totalAreaM2 > 0.0L) || !std::isfinite(totalAreaM2)) {
    *error = "uniform 40Ar source total surface area is invalid";
    return false;
  }
  elementRatesAtomsS->resize(elementAreasM2.size());
  for (std::size_t element = 0; element < elementAreasM2.size(); ++element) {
    (*elementRatesAtomsS)[element] = totalRateAtomsS *
        static_cast<double>(static_cast<long double>(elementAreasM2[element]) /
                            totalAreaM2);
  }
  return true;
}

bool ConfigureProduction(const Runtime::Configuration& configuration,
                         cInternalSphericalData* logicalSurface,
                         std::string* error) {
  if (error == nullptr) return false;
  SourceParameters parameters;
  parameters.enabled = configuration.argonSourceEnabled;
  parameters.spatialModel = configuration.argonSourceSpatialModel;
  parameters.timeModel = configuration.argonSourceTimeModel;
  parameters.totalRateAtomsS = configuration.argonSourceRateAtomsS;
  parameters.reference = configuration.argonSourceReference;
  if (!ValidateSourceParameters(parameters, error)) return false;

  gConfigured = false;
  gEnabled = parameters.enabled;
  gTotalRateAtomsS = 0.0;
  gLogicalSurface = nullptr;
  gElementRatesAtomsS.clear();
  if (!gEnabled) {
    gConfigured = true;
    return true;
  }

#if _AR_SPEC_ < 0
  *error = "radiogenic 40Ar source is enabled by the .in receipt, but Ar is "
           "absent from the compile-time SpeciesList; configure an Ar-enabled "
           "srcMoon build before running this receipt";
  return false;
#else
  if (logicalSurface == nullptr) {
    *error = "radiogenic 40Ar source requires the initialized logical lunar "
             "surface";
    return false;
  }
  const int count = logicalSurface->GetTotalSurfaceElementsNumber();
  if (count <= 0) {
    *error = "radiogenic 40Ar source found no logical surface elements";
    return false;
  }
  std::vector<double> areasM2(static_cast<std::size_t>(count));
  for (int element = 0; element < count; ++element)
    areasM2[static_cast<std::size_t>(element)] =
        logicalSurface->GetSurfaceElementArea(element);
  if (!BuildUniformElementRates(areasM2, parameters.totalRateAtomsS,
                                &gElementRatesAtomsS, error)) return false;
  gLogicalSurface = logicalSurface;
  gTotalRateAtomsS = parameters.totalRateAtomsS;
  gSurfaceDistribution.SetLimits(0, count - 1, ElementRate);
  gConfigured = true;
  return true;
#endif
}

double GetTotalProductionRate(int species, int boundaryElementType,
                              void* boundaryElement) {
#if _AR_SPEC_ < 0
  (void)species;
  (void)boundaryElementType;
  (void)boundaryElement;
  return 0.0;
#else
  if (!gConfigured || !gEnabled || species != _AR_SPEC_) return 0.0;
  if (boundaryElementType != _INTERNAL_BOUNDARY_TYPE_SPHERE_ ||
      boundaryElement != gLogicalSurface)
    exit(__LINE__, __FILE__,
         "radiogenic 40Ar source received an unexpected boundary object");
  return gTotalRateAtomsS;
#endif
}

void InitializeSurfaceSourceDistribution() {
  if (gConfigured && gEnabled) {
    // The generic discrete-distribution API passes an optional species token
    // into its weight callback.  ElementRate does not need it because this
    // module owns only 40Ar, but supplying the real compiled species keeps the
    // adapter conformant with the same API used by built-in source modules.
    int species = _AR_SPEC_;
    gSurfaceDistribution.Init(&species);
  }
}

bool GenerateParticleProperties(
    int species, PIC::ParticleBuffer::byte* temporaryParticleData,
    double* positionSoM, double* positionBodyM, double* velocitySoMPerS,
    double* velocityBodyMPerS, double* sphereCenterM, double sphereRadiusM,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>*& startNode,
    int boundaryElementType, void* boundaryElement) {
  (void)temporaryParticleData;
  (void)sphereCenterM;  // srcMoon's Moon-centred logical sphere is at zero.
#if _AR_SPEC_ < 0
  (void)species;
  (void)positionSoM;
  (void)positionBodyM;
  (void)velocitySoMPerS;
  (void)velocityBodyMPerS;
  (void)sphereRadiusM;
  (void)startNode;
  (void)boundaryElementType;
  (void)boundaryElement;
  return false;
#else
  if (!gConfigured || !gEnabled || species != _AR_SPEC_) return false;
  if (boundaryElementType != _INTERNAL_BOUNDARY_TYPE_SPHERE_ ||
      boundaryElement != gLogicalSurface)
    exit(__LINE__, __FILE__,
         "radiogenic 40Ar generator received an unexpected boundary object");

  // Draw a surface element with probability proportional to its physical
  // area, then draw a uniform position inside that logical element.  In LOLA
  // mode the existing production postprocessor moves this radial point to the
  // resolved terrain and rotates the velocity to the facet normal.
  const int element = gSurfaceDistribution.DistributeVariable();
  int zenith = 0, azimuth = 0;
  double outwardNormal[3];
  gLogicalSurface->GetSurfaceElementIndex(zenith, azimuth, element);
  gLogicalSurface->GetSurfaceElementRandomDirection(outwardNormal, zenith,
                                                     azimuth);
  for (int dimension = 0; dimension < 3; ++dimension)
    positionBodyM[dimension] = sphereRadiusM * outwardNormal[dimension];

  const double (*transform)[6] =
      Exosphere::OrbitalMotion::IAU_to_SO_TransformationMartix;
  for (int row = 0; row < 3; ++row) {
    positionSoM[row] = 0.0;
    for (int column = 0; column < 3; ++column)
      positionSoM[row] += transform[row][column] * positionBodyM[column];
  }
  startNode = PIC::Mesh::mesh->findTreeNode(positionSoM, startNode);
  if (startNode == nullptr || startNode->Thread != PIC::Mesh::mesh->ThisThread)
    return false;

  const double cosineSubsolar =
      Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[0]*outwardNormal[0] +
      Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[1]*outwardNormal[1] +
      Exosphere::OrbitalMotion::SunDirection_IAU_OBJECT[2]*outwardNormal[2];
  const double temperatureK =
      Exosphere::GetSurfaceTemperature(cosineSubsolar, positionSoM);
  if (!PositiveFinite(temperatureK))
    exit(__LINE__, __FILE__,
         "radiogenic 40Ar source received an invalid surface temperature");

  // InjectMaxwellianDistribution emits opposite its supplied external normal.
  // Supply the inward normal so the resulting distribution is the outward
  // Maxwell-Boltzmann flux distribution used by Kegerreis et al. (2017), not
  // a volume Maxwellian and not an isotropic direction distribution.
  double inwardNormal[3] = {-outwardNormal[0], -outwardNormal[1],
                            -outwardNormal[2]};
  const double zeroBulkMPerS[3] = {0.0, 0.0, 0.0};
  PIC::Distribution::InjectMaxwellianDistribution(
      velocityBodyMPerS, zeroBulkMPerS, temperatureK, inwardNormal, species);

  // Apply the complete 6x6 body-to-LSO state transform.  Its derivative block
  // adds the velocity of the rotating surface, so the sampled Maxwell flux is
  // relative to the lunar regolith rather than an inertial zero-velocity wall.
  for (int row = 0; row < 3; ++row) {
    velocitySoMPerS[row] = 0.0;
    for (int column = 0; column < 3; ++column) {
      velocitySoMPerS[row] +=
          transform[row+3][column] * positionBodyM[column] +
          transform[row+3][column+3] * velocityBodyMPerS[column];
    }
  }
  return true;
#endif
}

}  // namespace Argon
}  // namespace Species
}  // namespace Moon
