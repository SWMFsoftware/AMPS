#ifndef AMPS_CHEMICAL_DATABASE_ELECTRON_IMPACT_H
#define AMPS_CHEMICAL_DATABASE_ELECTRON_IMPACT_H

#include <array>
#include <cstddef>
#include <string>

namespace AMPS { namespace ChemicalDatabase { namespace ElectronImpact {

// These are neutral ATOMS, not AMPS particle-buffer indices. The application
// adapter owns explicit neutral-to-ion allocation; this coefficient service
// never invents molecular channels, product kinematics or electron moments.
enum class NeutralSpecies { He, Na, Ne, Ar };

struct Coefficient {
  int atomicNumber=0;
  int electronCount=0;
  double ionizationEnergyEv=0.0;
  int p=0;
  double aCm3S=0.0;
  double x=0.0;
  double k=0.0;
  double temperatureMinEv=0.0;
  double temperatureMaxEv=0.0; // native Tmax [keV] multiplied by exactly1000
  std::size_t sourceLine=0;
};

// Exact case-sensitive species identifiers prevent a typo or an ion label
// from silently selecting a neutral channel. All failing pure-layer calls
// leave their outputs unchanged and provide an error instead of aborting.
bool NeutralSpeciesFromName(const std::string& name,NeutralSpecies* species,
                            std::string* error);
const char* NeutralSpeciesName(NeutralSpecies species);

// Unit conversion only: Te[eV]=k_B[exact SI J/K]*T[K]/e[exact SI J/eV].
// The Maxwellian closure is selected explicitly by the application, not by
// this converter. Finite positive, representable temperatures are required.
bool ElectronTemperatureEvFromK(double temperatureK,double* temperatureEv,
                                std::string* error);

class MaxwellianVoronovTable {
 public:
  // Parse the COMPLETE native406-row cfit.dat with fixed positions from
  // cfit.txt. Validate every atomic stage before committing the four neutral
  // rows. The caller separately verifies immutable-source/provenance hashes;
  // a structurally valid altered table is not automatically authoritative.
  bool LoadNativeTable(const std::string& path,std::string* error);
  bool GetCoefficient(NeutralSpecies species,Coefficient* coefficient,
                       std::string* error) const;
  bool RateCoefficientM3S(NeutralSpecies species,double electronTemperatureEv,
                          double* coefficientM3S,std::string* error) const;
  bool IonizationRateS(NeutralSpecies species,double electronTemperatureEv,
                        double electronDensityM3,double* rateS,
                        std::string* error) const;
  bool loaded() const { return loaded_; }
  std::size_t recordCount() const { return loaded_ ? 406 : 0; }

 private:
  std::array<Coefficient,4> neutral_{};
  bool loaded_=false;
};

} } } // namespace AMPS::ChemicalDatabase::ElectronImpact
#endif
