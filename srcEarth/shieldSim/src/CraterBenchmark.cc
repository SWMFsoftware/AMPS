#include "CraterBenchmark.hh"

#include "ContentHash.hh"
#include "StrictJson.hh"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace shield_sim {
namespace {

std::string ReadAll(const std::string& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("cannot open CRaTER input '" + path + "'");
  std::ostringstream bytes;
  bytes << input.rdbuf();
  if (!input.good() && !input.eof())
    throw std::runtime_error("failed while reading CRaTER input '" + path + "'");
  return bytes.str();
}

std::string NormalizeSha(std::string value) {
  if (value.rfind("sha256:", 0) != 0) value = "sha256:" + value;
  return value;
}

void VerifySha(const std::string& bytes, const std::string& expected,
               const std::string& path) {
  const std::string actual = Sha256(bytes);
  if (actual != NormalizeSha(expected))
    throw std::runtime_error("SHA-256 mismatch for '" + path + "': expected " +
                             NormalizeSha(expected) + ", received " + actual);
}

std::string Resolve(const std::filesystem::path& parent, const std::string& value) {
  std::filesystem::path path(value);
  if (path.is_relative()) path = parent / path;
  return std::filesystem::absolute(path).lexically_normal().string();
}

std::int64_t ExactInteger(const Json& value, const std::string& context) {
  const double number = value.AsNumber(context);
  if (!std::isfinite(number) || number < 1.0 || std::floor(number) != number ||
      number > static_cast<double>(std::numeric_limits<std::int64_t>::max()))
    throw JsonError(context + " must be a positive exact integer");
  return static_cast<std::int64_t>(number);
}

std::vector<std::string> SplitCsv(const std::string& line) {
  // The Zenodo artifact is an unquoted numeric CSV.  Rejecting quotes here is
  // safer than silently implementing only a subset of RFC 4180.
  if (line.find('"') != std::string::npos)
    throw std::runtime_error("quoted fields are unsupported in the pinned CRaTER table");
  std::vector<std::string> result;
  std::size_t begin = 0;
  for (;;) {
    const std::size_t comma = line.find(',', begin);
    result.push_back(line.substr(begin, comma == std::string::npos
                                           ? comma : comma - begin));
    if (comma == std::string::npos) return result;
    begin = comma + 1;
  }
}

double ParseFinite(const std::string& text, const std::string& context) {
  std::size_t consumed = 0;
  double value = 0.0;
  try { value = std::stod(text, &consumed); }
  catch (const std::exception&) { throw std::runtime_error(context + " is not numeric"); }
  if (consumed != text.size() || !std::isfinite(value))
    throw std::runtime_error(context + " is not a finite number");
  return value;
}

}  // namespace

const std::vector<CraterAxialElement>& PublishedCraterAxialStack() {
  // Spence et al. (2010), Table 4.  Voids are retained explicitly so their
  // published separations are tested and the detector centers cannot drift
  // when the geometry code changes.
  static const std::vector<CraterAxialElement> elements = {
      {"S1", "Al6061_T6", 0.8128, false, 0},
      {"gap_S1_D1", "vacuum", 2.363, false, 0},
      {"D1", "silicon", 0.148, true, 1},
      {"gap_D1_D2", "vacuum", 3.360, false, 0},
      {"D2", "silicon", 1.000, true, 2},
      {"gap_D2_A1", "vacuum", 3.943, false, 0},
      {"A1", "A150_TEP", 53.9670, false, 0},
      {"gap_A1_D3", "vacuum", 3.554, false, 0},
      {"D3", "silicon", 0.149, true, 3},
      {"gap_D3_D4", "vacuum", 3.568, false, 0},
      {"D4", "silicon", 1.000, true, 4},
      {"gap_D4_A2", "vacuum", 11.900, false, 0},
      {"A2", "A150_TEP", 26.9800, false, 0},
      {"gap_A2_D5", "vacuum", 11.900, false, 0},
      {"D5", "silicon", 0.149, true, 5},
      {"gap_D5_D6", "vacuum", 3.060, false, 0},
      {"D6", "silicon", 1.000, true, 6},
      {"gap_D6_S2", "vacuum", 2.362, false, 0},
      {"S2", "Al6061_T6", 0.81026, false, 0},
  };
  return elements;
}

const std::vector<CraterPrimaryIon>& PublishedCraterPrimaryIons() {
  // These are not guessed atomic weights.  They are the ion_ZZ_AA filenames
  // present in dry_10cm_protons.zip from GLACE v1.0 (Zenodo 8343472,
  // md5:9bf407522d4a6529b2f0e86f3900a9be).  He-3 is a kernel in that archive,
  // but the GLACE README explicitly states that BON2020 does not include it.
  static const std::vector<CraterPrimaryIon> ions = {
      {1,1,"H"},   {2,4,"He"}, {3,7,"Li"},  {4,9,"Be"},  {5,11,"B"},
      {6,12,"C"},  {7,14,"N"}, {8,16,"O"},  {9,19,"F"},  {10,20,"Ne"},
      {11,23,"Na"},{12,24,"Mg"},{13,27,"Al"},{14,28,"Si"},{15,31,"P"},
      {16,32,"S"}, {17,35,"Cl"},{18,40,"Ar"},{19,39,"K"}, {20,40,"Ca"},
      {21,45,"Sc"},{22,48,"Ti"},{23,51,"V"}, {24,52,"Cr"},{25,55,"Mn"},
      {26,56,"Fe"},{27,59,"Co"},{28,58,"Ni"},
  };
  return ions;
}

CraterBenchmarkCase LoadCraterBenchmarkCase(const std::string& path) {
  const Json root = Json::ParseFile(path);
  const auto& object = root.AsObject("CRaTER benchmark case");
  RejectUnknownFields(object, {"schema_version", "benchmark_scope", "case_id",
      "histories", "random_seed", "accuracy", "transport", "source",
      "observation", "output"}, "CRaTER benchmark case");
  if (root.At("schema_version").AsString() !=
      "shieldSim.crater_foreground_axial_case.v1")
    throw JsonError("unsupported CRaTER benchmark case schema_version");
  if (root.At("benchmark_scope").AsString() != kCraterBenchmarkScope)
    throw JsonError("benchmark_scope must be 'foreground_only_axial'");

  const std::filesystem::path parent =
      std::filesystem::absolute(path).parent_path();
  CraterBenchmarkCase result;
  result.caseId = root.At("case_id").AsString();
  if (result.caseId.empty()) throw JsonError("case_id must not be empty");
  result.histories = ExactInteger(root.At("histories"), "histories");
  result.randomSeed = ExactInteger(root.At("random_seed"), "random_seed");

  const auto& accuracy = root.At("accuracy").AsObject("accuracy");
  RejectUnknownFields(accuracy, {"target_relative_statistical_uncertainty"},
                      "accuracy");
  result.targetRelativeStatisticalUncertainty =
      root.At("accuracy").At("target_relative_statistical_uncertainty").AsNumber();
  if (!(result.targetRelativeStatisticalUncertainty > 0.0 &&
        result.targetRelativeStatisticalUncertainty < 1.0))
    throw JsonError("accuracy target must be strictly between zero and one");

  const auto& transport = root.At("transport").AsObject("transport");
  RejectUnknownFields(transport, {"physics_list", "em_physics", "production_cut_um"},
                      "transport");
  result.physicsList = root.At("transport").At("physics_list").AsString();
  result.emPhysics = root.At("transport").At("em_physics").AsString();
  result.productionCutUm =
      root.At("transport").At("production_cut_um").AsNumber();
  if (result.physicsList != "Shielding" ||
      result.emPhysics != "G4EmStandardPhysics_option4")
    throw JsonError("CRaTER v1 benchmark requires Shielding plus "
                    "G4EmStandardPhysics_option4 (the current Shielding_EMZ equivalent)");
  if (!(result.productionCutUm > 0.0))
    throw JsonError("production_cut_um must be positive");

  const auto& source = root.At("source").AsObject("source");
  RejectUnknownFields(source, {"path", "sha256", "interpolation"}, "source");
  result.sourcePath = Resolve(parent, root.At("source").At("path").AsString());
  result.sourceSha256 = root.At("source").At("sha256").AsString();
  result.sourceInterpolation =
      root.At("source").At("interpolation").AsString();
  if (result.sourceInterpolation != "log_log")
    throw JsonError("only explicit source interpolation 'log_log' is supported");

  const auto& observation = root.At("observation").AsObject("observation");
  RejectUnknownFields(observation, {"path", "sha256", "column"}, "observation");
  result.observationPath =
      Resolve(parent, root.At("observation").At("path").AsString());
  result.observationSha256 = root.At("observation").At("sha256").AsString();
  result.observationColumn = root.At("observation").At("column").AsString();

  const auto& output = root.At("output").AsObject("output");
  RejectUnknownFields(output, {"directory"}, "output");
  result.outputDirectory =
      Resolve(parent, root.At("output").At("directory").AsString());
  return result;
}

CraterBonSpectrum LoadCraterBonSpectrum(const std::string& path,
                                        const std::string& expectedSha256) {
  const std::string bytes = ReadAll(path);
  VerifySha(bytes, expectedSha256, path);
  const Json root = Json::Parse(bytes, path);
  const auto& object = root.AsObject("BON2020 source");
  RejectUnknownFields(object, {"BON_CODE_VERSION", "RUN_DATES", "PHI_UNITS",
      "ENERGY_UNITS", "FLUX_UNITS", "RUN_PHI", "GCR_Z", "GCR_ENERGIES",
      "GCR_FLUXES"}, "BON2020 source");
  if (root.At("PHI_UNITS").AsString() != "MV" ||
      root.At("ENERGY_UNITS").AsString() != "MeV/nucleon" ||
      root.At("FLUX_UNITS").AsString() != "ions / (MeV/nuc m^2 s sr)")
    throw JsonError("unexpected BON2020 units");

  CraterBonSpectrum result;
  result.modelVersion = root.At("BON_CODE_VERSION").AsString();
  result.runDates = root.At("RUN_DATES").AsString();
  result.fluxUnits = root.At("FLUX_UNITS").AsString();
  result.modulationPotentialMv = root.At("RUN_PHI").AsNumber();
  const auto& z = root.At("GCR_Z").AsArray();
  if (z.size() != 28) throw JsonError("BON2020 GCR_Z must contain 28 species");
  for (std::size_t i = 0; i < z.size(); ++i)
    if (z[i].AsNumber() != static_cast<double>(i + 1))
      throw JsonError("BON2020 GCR_Z must be the ordered sequence 1..28");
  const auto& energies = root.At("GCR_ENERGIES").AsArray();
  const auto& fluxes = root.At("GCR_FLUXES").AsArray();
  if (energies.size() < 2 || fluxes.size() != energies.size())
    throw JsonError("BON2020 energy/flux dimensions are inconsistent");
  for (std::size_t i = 0; i < energies.size(); ++i) {
    const double energy = energies[i].AsNumber();
    if (!(energy > 0.0) || (i && !(energy > result.energyMeVPerNucleon.back())))
      throw JsonError("BON2020 energies must be positive and strictly increasing");
    result.energyMeVPerNucleon.push_back(energy);
    const auto& row = fluxes[i].AsArray();
    if (row.size() != 28) throw JsonError("each BON2020 flux row must contain 28 species");
    std::vector<double> values;
    values.reserve(28);
    for (const auto& item : row) {
      const double value = item.AsNumber();
      if (!std::isfinite(value) || value < 0.0)
        throw JsonError("BON2020 flux must be finite and nonnegative");
      values.push_back(value);
    }
    result.fluxByEnergy.push_back(std::move(values));
  }
  return result;
}

CraterObservation LoadCraterObservation(const std::string& path,
                                        const std::string& expectedSha256,
                                        const std::string& column) {
  const std::string bytes = ReadAll(path);
  VerifySha(bytes, expectedSha256, path);
  std::istringstream input(bytes);
  std::string line;
  if (!std::getline(input, line)) throw std::runtime_error("empty CRaTER observation");
  if (!line.empty() && line.back() == '\r') line.pop_back();
  const auto header = SplitCsv(line);
  const auto match = std::find(header.begin(), header.end(), column);
  if (match == header.end())
    throw std::runtime_error("observation column '" + column + "' is absent");
  const std::size_t selected = static_cast<std::size_t>(match - header.begin());
  CraterObservation result;
  result.column = column;
  std::size_t rowNumber = 1;
  while (std::getline(input, line)) {
    ++rowNumber;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    const auto fields = SplitCsv(line);
    if (fields.size() != header.size())
      throw std::runtime_error("wrong field count in CRaTER observation row " +
                               std::to_string(rowNumber));
    const double center = ParseFinite(fields[0], "CRaTER LET center");
    const double value = ParseFinite(fields[selected], "CRaTER LET intensity");
    if (!(center > 0.0) || value < 0.0 ||
        (!result.centersKevPerUm.empty() && center <= result.centersKevPerUm.back()))
      throw std::runtime_error("invalid or unordered CRaTER observation row " +
                               std::to_string(rowNumber));
    result.centersKevPerUm.push_back(center);
    result.differentialIntensity.push_back(value);
  }
  if (result.centersKevPerUm.empty()) throw std::runtime_error("empty CRaTER observation table");
  (void)CraterObservationBinEdges(result.centersKevPerUm);
  return result;
}

double InterpolateCraterFlux(const CraterBonSpectrum& spectrum, int z,
                             double energyMeVPerNucleon) {
  if (z < 1 || z > 28) throw std::out_of_range("CRaTER source Z is outside 1..28");
  const auto& energies = spectrum.energyMeVPerNucleon;
  if (energyMeVPerNucleon < energies.front() || energyMeVPerNucleon > energies.back())
    throw std::out_of_range("CRaTER source interpolation would require extrapolation");
  auto upper = std::lower_bound(energies.begin(), energies.end(), energyMeVPerNucleon);
  if (upper == energies.begin()) return spectrum.fluxByEnergy.front()[z - 1];
  if (upper == energies.end()) return spectrum.fluxByEnergy.back()[z - 1];
  const std::size_t hi = static_cast<std::size_t>(upper - energies.begin());
  if (*upper == energyMeVPerNucleon) return spectrum.fluxByEnergy[hi][z - 1];
  const std::size_t lo = hi - 1;
  const double f0 = spectrum.fluxByEnergy[lo][z - 1];
  const double f1 = spectrum.fluxByEnergy[hi][z - 1];
  if (f0 == 0.0 || f1 == 0.0) {
    // Log interpolation is undefined at zero.  Linear interpolation between
    // exact tabulated zeros is not silently substituted; a zero adjacent to a
    // positive value indicates that the source is outside a usable domain.
    if (f0 == 0.0 && f1 == 0.0) return 0.0;
    throw std::runtime_error("log-log BON2020 interpolation encountered an isolated zero");
  }
  const double fraction = std::log(energyMeVPerNucleon / energies[lo]) /
                          std::log(energies[hi] / energies[lo]);
  return std::exp(std::log(f0) + fraction * (std::log(f1) - std::log(f0)));
}

std::vector<double> CraterObservationBinEdges(
    const std::vector<double>& centersKevPerUm) {
  if (centersKevPerUm.empty()) throw std::runtime_error("CRaTER bin centers are empty");
  for (std::size_t i = 1; i < centersKevPerUm.size(); ++i) {
    const double observed = centersKevPerUm[i] / centersKevPerUm[i - 1];
    const bool d6ToD5Splice = centersKevPerUm[i - 1] < 20.0 &&
                              centersKevPerUm[i] > 20.0;
    // The paper describes 80 plotted bins/decade, while the released rounded
    // centers have a ratio near 10^(1/81).  Since no edge array accompanies
    // the table, do not silently choose either convention.  Require only a
    // fine logarithmic grid and reconstruct explicit midpoint edges below.
    if (!d6ToD5Splice && !(observed > 1.0 && observed < 1.05))
      throw std::runtime_error("CRaTER observation centers are not a fine logarithmic grid");
  }
  std::vector<double> edges;
  edges.reserve(centersKevPerUm.size() + 1);
  if (centersKevPerUm.size() == 1)
    throw std::runtime_error("at least two CRaTER observation centers are required");
  edges.push_back(centersKevPerUm.front() /
                  std::sqrt(centersKevPerUm[1]/centersKevPerUm[0]));
  for (std::size_t i = 1; i < centersKevPerUm.size(); ++i) {
    const bool d6ToD5Splice = centersKevPerUm[i - 1] < 20.0 &&
                              centersKevPerUm[i] > 20.0;
    edges.push_back(d6ToD5Splice ? 20.0
                                 : std::sqrt(centersKevPerUm[i - 1] * centersKevPerUm[i]));
  }
  edges.push_back(centersKevPerUm.back() * std::sqrt(
      centersKevPerUm.back()/centersKevPerUm[centersKevPerUm.size()-2]));
  return edges;
}

bool CraterCoincidence(const CraterDeposits& deposits) {
  // Spence et al. report approximate singles thresholds of 35 keV for thick
  // and 110 keV for thin detectors.  Expressing those in normalized-deposit
  // units uses the as-delivered D2/D4/D6=1000 um and D5=149 um thicknesses.
  const double thick = 35.0 / 1000.0;
  const double d5 = 110.0 / 149.0;
  return deposits.kevPerUm[1] >= thick && deposits.kevPerUm[3] >= thick &&
         (deposits.kevPerUm[5] >= thick || deposits.kevPerUm[4] >= d5);
}

bool CraterD5PublishedCut(const CraterDeposits& d) {
  const double d2=d.kevPerUm[1], d4=d.kevPerUm[3], d5=d.kevPerUm[4];
  return (((d5 < 10.0 || d4 > 1.7) && (d5 < 25.0 || d2 > 4.0) &&
           (d5 < 100.0 || d2 > 15.0)) ||
          (d5 > 10.0 && d4 < 1.7 && d2 < 0.6)) &&
         (d2 > 10.0 || 3.0*d2 + 4.0 > 2.0*d4);
}

bool CraterD6PublishedCut(const CraterDeposits& d) {
  const double d2=d.kevPerUm[1], d4=d.kevPerUm[3], d6=d.kevPerUm[5];
  return (((d2 > 4.0 || d6 < 10.0*d2 - 27.5 || 2.0*d6 < 5.0*d2 - 4.0 ||
            d6 < 4.0 || d6 + 10.0*d2 < 22.5 || d2 < 1.6) &&
           d2 > 0.45 && (d2 > 4.0 || d6 < 13.0) &&
           (d2 > 0.9 || 10.0*d2 > d6 + 6.0 || d6 < 2.0 ||
            5.0*d2 + d6 < 4.75) &&
           (d2 > 10.0 || 3.0*d2 + 4.0 > 2.0*d4) &&
           (d2 < 2.5 || d4 > 3.0 || d2 < 12.0*d4 - 9.5) &&
           (d4 < 3.0 || d6 > 4.0 || 9.0*d6 > d4 + 6.0)) ||
          (d2 < 0.45 && (d4 < 0.7 || d4 + d6 < 2.7)));
}

double CraterCombinedObservable(const CraterDeposits& deposits) {
  if (!CraterCoincidence(deposits)) return -1.0;
  return deposits.kevPerUm[5] <= 20.0 ? deposits.kevPerUm[5]
                                      : deposits.kevPerUm[4];
}

bool CraterCombinedPublishedCut(const CraterDeposits& deposits) {
  if (!CraterCoincidence(deposits)) return false;
  return deposits.kevPerUm[5] <= 20.0 ? CraterD6PublishedCut(deposits)
                                      : CraterD5PublishedCut(deposits);
}

}  // namespace shield_sim
