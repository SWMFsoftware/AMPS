#ifndef SHIELDSIM_CRATER_BENCHMARK_HH
#define SHIELDSIM_CRATER_BENCHMARK_HH

/*
 * Foreground-only axial CRaTER benchmark: Geant4-independent definitions.
 *
 * This interface deliberately separates published facts (the axial stack,
 * detector aperture, D1-D6 geometry factor, isotope identities, and event
 * cuts) from the Geant4 transport harness.  The separation lets unit tests
 * protect the scientific boundary without initializing a physics list and
 * makes it difficult for an engineering-mass-model assumption to leak into
 * this explicitly limited benchmark.
 */

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace shield_sim {

inline constexpr const char* kCraterBenchmarkScope =
    "foreground_only_axial";
inline constexpr double kCraterDetectorRadiusMm = 17.5;
inline constexpr double kCraterD1D6GeometryFactorCm2Sr = 0.569;
inline constexpr double kCraterFullAngleFieldOfViewDeg = 31.4;

struct CraterAxialElement {
  std::string id;
  std::string material;
  double thicknessMm = 0.0;
  bool sensitive = false;
  int detectorNumber = 0;
};

struct CraterPrimaryIon {
  int z = 0;
  int a = 0;
  std::string symbol;
};

struct CraterDeposits {
  // Normalized energy deposits in keV/um, in detector order D1 ... D6.
  std::array<double, 6> kevPerUm{};
};

struct CraterBenchmarkCase {
  std::string caseId;
  std::string sourcePath;
  std::string sourceSha256;
  std::string observationPath;
  std::string observationSha256;
  std::string observationColumn;
  std::string outputDirectory;
  std::string physicsList;
  std::string emPhysics;
  std::string sourceInterpolation;
  std::int64_t histories = 0;
  std::int64_t randomSeed = 0;
  double productionCutUm = 0.0;
  double targetRelativeStatisticalUncertainty = 0.10;
};

struct CraterBonSpectrum {
  std::string modelVersion;
  std::string runDates;
  double modulationPotentialMv = 0.0;
  std::vector<double> energyMeVPerNucleon;
  // fluxByEnergy[energy_index][z-1], in the units stated by fluxUnits.
  std::vector<std::vector<double>> fluxByEnergy;
  std::string fluxUnits;
};

struct CraterObservation {
  std::string column;
  std::vector<double> centersKevPerUm;
  std::vector<double> differentialIntensity;
};

// Exact values transcribed from Spence et al. (2010), Tables 2--5.
const std::vector<CraterAxialElement>& PublishedCraterAxialStack();

// Exact primary-ion Z,A mapping recovered from the GLACE v1.0 immutable ZIP
// member names ion_ZZ_AA.json.  BON2020 supplies one flux column per Z and
// excludes the separately tabulated He-3 kernel.
const std::vector<CraterPrimaryIon>& PublishedCraterPrimaryIons();

CraterBenchmarkCase LoadCraterBenchmarkCase(const std::string& path);
CraterBonSpectrum LoadCraterBonSpectrum(const std::string& path,
                                        const std::string& expectedSha256);
CraterObservation LoadCraterObservation(const std::string& path,
                                        const std::string& expectedSha256,
                                        const std::string& column);

// Log-log interpolation is explicit in the benchmark case.  Values outside
// the supplied energy grid are rejected, never extrapolated or clamped.
double InterpolateCraterFlux(const CraterBonSpectrum& spectrum, int z,
                             double energyMeVPerNucleon);

// The published v1.3 table provides centers but not edges.  This constructs
// explicit geometric-midpoint edges separately on the D6 and D5 sides of the
// exact 20 keV/um splice, and extrapolates only the two outer half bins.  The
// reconstruction is recorded as a benchmark limitation rather than described
// as publisher-supplied bin metadata.
std::vector<double> CraterObservationBinEdges(
    const std::vector<double>& centersKevPerUm);

bool CraterCoincidence(const CraterDeposits& deposits);
bool CraterD5PublishedCut(const CraterDeposits& deposits);
bool CraterD6PublishedCut(const CraterDeposits& deposits);

// D6 supplies the spectrum at and below 20 keV/um; D5 supplies it above 20.
// Returns a negative value when the published coincidence is not satisfied.
double CraterCombinedObservable(const CraterDeposits& deposits);
bool CraterCombinedPublishedCut(const CraterDeposits& deposits);

}  // namespace shield_sim

#endif
