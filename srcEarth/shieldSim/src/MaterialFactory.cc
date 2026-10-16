#include "MaterialFactory.hh"

#include <G4Element.hh>
#include <G4Exception.hh>
#include <G4Isotope.hh>
#include <G4Material.hh>
#include <G4IonisParamMat.hh>
#include <G4NistManager.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>

#include <algorithm>
#include <cctype>
#include <cstdlib>

namespace shield_sim {

G4Material* MaterialFactory::Resolve(const std::string& idOrAlias, const std::string& role) {
  return Build(library_.Resolve(idOrAlias, role));
}

G4Material* MaterialFactory::Build(const MaterialRecord& record) {
  const auto cached = built_.find(record.contentHash);
  if (cached != built_.end()) return cached->second;
  if (record.status == MaterialStatus::Template)
    G4Exception("MaterialFactory::Build", "TemplateMaterial", FatalException,
                ("Cannot instantiate template material '" + record.id + "'").c_str());

  const auto fractions = library_.ElementalMassFractions(record);
  G4State state = kStateUndefined;
  if (record.phase == "solid") state = kStateSolid;
  else if (record.phase == "liquid") state = kStateLiquid;
  else if (record.phase == "gas") state = kStateGas;

  // Hash suffixes allow a reference record and one or more case-local
  // overrides of that record to coexist in Geant4's global material table.
  std::string safeId = record.id;
  std::replace_if(safeId.begin(), safeId.end(),
                  [](char c){ return !(std::isalnum(static_cast<unsigned char>(c)) || c == '_'); }, '_');
  const std::string hashSuffix = record.contentHash.substr(record.contentHash.size() - 12);
  const std::string name = "shieldSim_" + safeId + '_' + hashSuffix;
  auto* material = new G4Material(name, record.densityGPerCm3 * g/cm3,
                                  static_cast<G4int>(fractions.size()), state,
                                  record.temperatureK * kelvin,
                                  record.hasPressure ? record.pressurePa*pascal
                                                     : STP_Pressure);
  auto* nist = G4NistManager::Instance();
  for (const auto& element : fractions) {
    G4Element* g4Element=nullptr;
    if(element.first.rfind("isotope:",0)==0){
      const std::string isotopeId=element.first.substr(8);
      const std::size_t dash=isotopeId.rfind('-');
      const std::string symbol=isotopeId.substr(0,dash);
      const int massNumber=std::stoi(isotopeId.substr(dash+1));
      auto* natural=nist->FindOrBuildElement(symbol,true);
      if(!natural)
        G4Exception("MaterialFactory::Build","IsotopeElement",FatalException,
                    ("Geant4 cannot build isotope element '"+symbol+"'").c_str());
      const G4int atomicNumber=static_cast<G4int>(natural->GetZ());
      const std::string isotopeName=name+"_"+isotopeId;
      auto* isotope=new G4Isotope(isotopeName,atomicNumber,massNumber,massNumber*g/mole);
      g4Element=new G4Element(isotopeName,symbol,1);
      g4Element->AddIsotope(isotope,100.0*perCent);
    } else {
      g4Element = nist->FindOrBuildElement(element.first, true);
    }
    if (!g4Element)
      G4Exception("MaterialFactory::Build", "Element", FatalException,
                  ("Geant4 cannot build element '" + element.first + "'").c_str());
    material->AddElement(g4Element, element.second);
  }
  if(record.hasMeanExcitationEnergy)
    material->GetIonisation()->SetMeanExcitationEnergy(record.meanExcitationEnergyEV*eV);
  built_[record.contentHash] = material;
  return material;
}

} // namespace shield_sim
