#ifndef SRCMOON_MOON_SURFACE_PHYSICS_H
#define SRCMOON_MOON_SURFACE_PHYSICS_H

#include <cstddef>
#include <string>
#include <vector>

namespace Moon {
namespace SurfacePhysics {

enum class Species { Sodium, Argon };
enum class ResidenceLaw { FirstOrderArrhenius, ModifiedArrhenius };

struct SpeciesParameters {
  Species species = Species::Sodium;
  ResidenceLaw residenceLaw = ResidenceLaw::FirstOrderArrhenius;
  // FirstOrderArrhenius uses nu*exp(-E/kT).  ModifiedArrhenius uses the
  // transition-state prefactor kT/h and a molar activation energy E/(RT).
  double attemptFrequencyHz = 0.0;
  double activationEnergyJ = 0.0;
  double activationEnergyJmol = 0.0;
  std::vector<double> stickingTemperatureK;
  std::vector<double> stickingProbability;
  std::string citation;
};

double StickingProbability(const SpeciesParameters& parameters,
                           double temperatureK);
double ResidenceTimeS(const SpeciesParameters& parameters,
                      double temperatureK);
double DesorptionProbability(const SpeciesParameters& parameters,
                             double temperatureK, double timeStepS);

// Exact event accounting for one logical surface element.  Values are real
// particles/atoms represented by statistical weights, never model-particle
// counts.  Numerical loss is deliberately separate from physical escape.
struct Budget {
  double injected = 0.0;
  double adsorbed = 0.0;
  double desorbed = 0.0;
  double convertedIn = 0.0;
  double convertedOut = 0.0;
  double escaped = 0.0;
  double numericalLoss = 0.0;
  double initialExosphere = 0.0;
  double initialInventory = 0.0;
  double finalExosphere = 0.0;
  double finalInventory = 0.0;

  // Positive residual means more atoms entered the modeled system than are
  // present or left it. Adsorption/desorption are internal transfers and are
  // retained as diagnostics but cancel from the global closure equation.
  double ClosureResidual() const;
};

// Continuous adsorbed state used by I18--I20.  The cold-trap flag is a
// diagnostic classification; it never changes the kinetic equation into an
// irreversible sink.  Warming therefore releases inventory through the same
// temperature-dependent desorption law.
struct InventoryCell {
  double adsorbedAtoms = 0.0;
  bool geometricPermanentShadow = false;
  bool thermalRetention = false;
};

class Inventory {
 public:
  bool Initialize(std::size_t cellCount, double initialAtomsPerCell,
                  std::string* error);
  bool Deposit(std::size_t cell, double atoms, Budget* budget,
               std::string* error);
  bool AdvanceCell(std::size_t cell, const SpeciesParameters& parameters,
                   double temperatureK, double timeStepS,
                   double retentionTimescaleS, double* releasedAtoms,
                   Budget* budget, std::string* error);
  // Synchronize a production-owned logical surface array into the checkpoint
  // representation.  This is deliberately an explicit validated operation;
  // exposing mutable cells would allow callers to bypass non-negative and
  // finite inventory invariants.
  bool SetCellState(std::size_t cell, double adsorbedAtoms,
                    bool geometricPermanentShadow, bool thermalRetention,
                    std::string* error);
  const std::vector<InventoryCell>& cells() const { return cells_; }

  bool WriteCheckpoint(const std::string& path, const std::string& speciesName,
                       double epochEt, std::string* error) const;
  bool ReadCheckpoint(const std::string& path,
                      const std::string& expectedSpeciesName,
                      double* epochEt, std::string* error);

 private:
  std::vector<InventoryCell> cells_;
};

const char* SpeciesName(Species species);
const char* ResidenceLawName(ResidenceLaw law);

}  // namespace SurfacePhysics
}  // namespace Moon

#endif  // SRCMOON_MOON_SURFACE_PHYSICS_H
