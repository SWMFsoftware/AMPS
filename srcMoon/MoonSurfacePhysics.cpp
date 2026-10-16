#include "MoonSurfacePhysics.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Moon {
namespace SurfacePhysics {
namespace {

constexpr double kBoltzmannJk = 1.380649e-23;
constexpr double kPlanckJs = 6.62607015e-34;
constexpr double kGasConstantJmolK = 8.31446261815324;

bool PositiveFinite(double value) {
  return std::isfinite(value) && value > 0.0;
}

bool ValidateTable(const SpeciesParameters& parameters) {
  if (parameters.stickingTemperatureK.empty() ||
      parameters.stickingTemperatureK.size() !=
          parameters.stickingProbability.size()) return false;
  for (std::size_t index = 0;
       index < parameters.stickingTemperatureK.size(); ++index) {
    if (!PositiveFinite(parameters.stickingTemperatureK[index]) ||
        !std::isfinite(parameters.stickingProbability[index]) ||
        parameters.stickingProbability[index] < 0.0 ||
        parameters.stickingProbability[index] > 1.0 ||
        (index != 0 && parameters.stickingTemperatureK[index] <=
                           parameters.stickingTemperatureK[index - 1]))
      return false;
  }
  return true;
}

}  // namespace

double StickingProbability(const SpeciesParameters& parameters,
                           double temperatureK) {
  if (!PositiveFinite(temperatureK) || !ValidateTable(parameters))
    return std::numeric_limits<double>::quiet_NaN();
  if (temperatureK <= parameters.stickingTemperatureK.front())
    return parameters.stickingProbability.front();
  if (temperatureK >= parameters.stickingTemperatureK.back())
    return parameters.stickingProbability.back();
  const auto upper = std::upper_bound(parameters.stickingTemperatureK.begin(),
      parameters.stickingTemperatureK.end(), temperatureK);
  const std::size_t high =
      static_cast<std::size_t>(upper - parameters.stickingTemperatureK.begin());
  const std::size_t low = high - 1;
  const double fraction = (temperatureK-
      parameters.stickingTemperatureK[low]) /
      (parameters.stickingTemperatureK[high]-
       parameters.stickingTemperatureK[low]);
  return (1.0-fraction)*parameters.stickingProbability[low] +
      fraction*parameters.stickingProbability[high];
}

double ResidenceTimeS(const SpeciesParameters& parameters,
                      double temperatureK) {
  if (!PositiveFinite(temperatureK))
    return std::numeric_limits<double>::quiet_NaN();
  double logarithm = 0.0;
  switch (parameters.residenceLaw) {
    case ResidenceLaw::FirstOrderArrhenius:
      if (!PositiveFinite(parameters.attemptFrequencyHz) ||
          !PositiveFinite(parameters.activationEnergyJ))
        return std::numeric_limits<double>::quiet_NaN();
      logarithm = parameters.activationEnergyJ /
          (kBoltzmannJk*temperatureK) -
          std::log(parameters.attemptFrequencyHz);
      break;
    case ResidenceLaw::ModifiedArrhenius:
      if (!PositiveFinite(parameters.activationEnergyJmol))
        return std::numeric_limits<double>::quiet_NaN();
      // Kegerreis et al. (2017), Eq. 4: tau=h/(kT)*exp(Q/RT).
      logarithm = std::log(kPlanckJs/(kBoltzmannJk*temperatureK)) +
          parameters.activationEnergyJmol /
              (kGasConstantJmolK*temperatureK);
      break;
  }
  if (logarithm > std::log(std::numeric_limits<double>::max()))
    return std::numeric_limits<double>::infinity();
  if (logarithm < std::log(std::numeric_limits<double>::min())) return 0.0;
  return std::exp(logarithm);
}

double DesorptionProbability(const SpeciesParameters& parameters,
                             double temperatureK, double timeStepS) {
  if (!PositiveFinite(timeStepS))
    return std::numeric_limits<double>::quiet_NaN();
  const double residenceTime = ResidenceTimeS(parameters, temperatureK);
  if (std::isinf(residenceTime)) return 0.0;
  if (!PositiveFinite(residenceTime)) return
      residenceTime == 0.0 ? 1.0 : std::numeric_limits<double>::quiet_NaN();
  // -expm1(-x) is accurate for dt << tau and is algebraically identical to
  // the integrated first-order probability 1-exp(-dt/tau).
  return -std::expm1(-timeStepS/residenceTime);
}

double Budget::ClosureResidual() const {
  return initialExosphere + initialInventory + injected + convertedIn -
      convertedOut - escaped - numericalLoss - finalExosphere -
      finalInventory;
}

bool Inventory::Initialize(std::size_t cellCount, double initialAtomsPerCell,
                           std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  if (cellCount == 0 || !std::isfinite(initialAtomsPerCell) ||
      initialAtomsPerCell < 0.0) {
    *error = "surface inventory requires a positive cell count and a "
             "non-negative finite initial population";
    return false;
  }
  cells_.assign(cellCount, InventoryCell{});
  for (InventoryCell& cell : cells_) cell.adsorbedAtoms = initialAtomsPerCell;
  return true;
}

bool Inventory::Deposit(std::size_t cell, double atoms, Budget* budget,
                        std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  if (cell >= cells_.size() || !std::isfinite(atoms) || atoms < 0.0) {
    *error = "invalid surface-inventory deposition";
    return false;
  }
  cells_[cell].adsorbedAtoms += atoms;
  if (budget != nullptr) budget->adsorbed += atoms;
  return true;
}

bool Inventory::AdvanceCell(std::size_t cell,
                            const SpeciesParameters& parameters,
                            double temperatureK, double timeStepS,
                            double retentionTimescaleS,
                            double* releasedAtoms, Budget* budget,
                            std::string* error) {
  if (error == nullptr || releasedAtoms == nullptr) return false;
  error->clear();
  if (cell >= cells_.size() || !PositiveFinite(retentionTimescaleS)) {
    *error = "invalid surface-inventory cell or retention timescale";
    return false;
  }
  const double probability =
      DesorptionProbability(parameters, temperatureK, timeStepS);
  const double residence = ResidenceTimeS(parameters, temperatureK);
  if (!std::isfinite(probability) || probability < 0.0 || probability > 1.0 ||
      std::isnan(residence)) {
    *error = "invalid species residence/desorption law";
    return false;
  }
  InventoryCell& state = cells_[cell];
  *releasedAtoms = state.adsorbedAtoms * probability;
  state.adsorbedAtoms -= *releasedAtoms;
  // Thermal retention is a reported comparison of timescales.  It does not
  // override the continuous kinetic release above or create a permanent sink.
  state.thermalRetention = residence >= retentionTimescaleS;
  if (budget != nullptr) budget->desorbed += *releasedAtoms;
  return true;
}

bool Inventory::SetCellState(std::size_t cell, double adsorbedAtoms,
                             bool geometricPermanentShadow,
                             bool thermalRetention, std::string* error) {
  if (error == nullptr) return false;
  error->clear();
  if (cell >= cells_.size() || !std::isfinite(adsorbedAtoms) ||
      adsorbedAtoms < 0.0) {
    *error = "invalid production surface-inventory state";
    return false;
  }
  cells_[cell].adsorbedAtoms = adsorbedAtoms;
  cells_[cell].geometricPermanentShadow = geometricPermanentShadow;
  cells_[cell].thermalRetention = thermalRetention;
  return true;
}

bool Inventory::WriteCheckpoint(const std::string& path,
                                const std::string& speciesName,
                                double epochEt, std::string* error) const {
  if (error == nullptr) return false;
  error->clear();
  if (cells_.empty() || speciesName.empty() || !std::isfinite(epochEt)) {
    *error = "cannot checkpoint an uninitialized surface inventory";
    return false;
  }
  namespace fs = std::filesystem;
  std::error_code filesystemError;
  const fs::path target(path);
  if (!target.parent_path().empty())
    fs::create_directories(target.parent_path(), filesystemError);
  if (filesystemError) {
    *error = "cannot create checkpoint directory: " +
        filesystemError.message();
    return false;
  }
  const fs::path temporary = target.string()+".tmp";
  std::ofstream output(temporary, std::ios::trunc);
  if (!output) {
    *error = "cannot open surface checkpoint temporary file";
    return false;
  }
  output << "MOON_SURFACE_INVENTORY_V1\n";
  output << "species " << speciesName << "\n";
  output << std::setprecision(17) << "epoch_et_s " << epochEt << "\n";
  output << "cell_count " << cells_.size() << "\n";
  output << "cell,adsorbed_atoms,geometric_psr,thermal_retention\n";
  for (std::size_t index=0; index<cells_.size(); ++index)
    output << index << ',' << cells_[index].adsorbedAtoms << ','
           << (cells_[index].geometricPermanentShadow ? 1 : 0) << ','
           << (cells_[index].thermalRetention ? 1 : 0) << '\n';
  output.close();
  if (!output) {
    *error = "I/O failure while writing surface checkpoint";
    return false;
  }
  fs::rename(temporary, target, filesystemError);
  if (filesystemError) {
    *error = "cannot atomically install surface checkpoint: " +
        filesystemError.message();
    return false;
  }
  return true;
}

bool Inventory::ReadCheckpoint(const std::string& path,
                               const std::string& expectedSpeciesName,
                               double* epochEt, std::string* error) {
  if (error == nullptr || epochEt == nullptr) return false;
  error->clear();
  std::ifstream input(path);
  std::string line, species;
  std::size_t cellCount=0;
  if (!input || !std::getline(input,line) ||
      line!="MOON_SURFACE_INVENTORY_V1" ||
      !(input >> line >> species) || line!="species" ||
      species!=expectedSpeciesName ||
      !(input >> line >> *epochEt) || line!="epoch_et_s" ||
      !(input >> line >> cellCount) || line!="cell_count" ||
      !std::isfinite(*epochEt) || cellCount==0) {
    *error = "invalid or species-mismatched surface checkpoint header";
    return false;
  }
  std::getline(input,line);
  if (!std::getline(input,line) ||
      line!="cell,adsorbed_atoms,geometric_psr,thermal_retention") {
    *error = "invalid surface checkpoint column header";
    return false;
  }
  std::vector<InventoryCell> loaded(cellCount);
  for (std::size_t expected=0; expected<cellCount; ++expected) {
    if (!std::getline(input,line)) {
      *error = "truncated surface checkpoint";
      return false;
    }
    std::replace(line.begin(),line.end(),',',' ');
    std::istringstream record(line);
    std::size_t index=0;
    int geometric=0, thermal=0;
    if (!(record>>index>>loaded[expected].adsorbedAtoms>>geometric>>thermal) ||
        index!=expected || !std::isfinite(loaded[expected].adsorbedAtoms) ||
        loaded[expected].adsorbedAtoms<0.0 ||
        (geometric!=0 && geometric!=1) || (thermal!=0 && thermal!=1)) {
      *error = "invalid surface checkpoint record";
      return false;
    }
    loaded[expected].geometricPermanentShadow=geometric!=0;
    loaded[expected].thermalRetention=thermal!=0;
  }
  while (std::getline(input,line)) if (!line.empty()) {
    *error = "surface checkpoint contains trailing records";
    return false;
  }
  cells_=std::move(loaded);
  return true;
}

const char* SpeciesName(Species species) {
  return species==Species::Sodium ? "Na" : "Ar";
}

const char* ResidenceLawName(ResidenceLaw law) {
  return law==ResidenceLaw::FirstOrderArrhenius ?
      "first-order-arrhenius" : "modified-arrhenius";
}

}  // namespace SurfacePhysics
}  // namespace Moon
