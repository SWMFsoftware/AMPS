#ifndef SHIELDSIM_MATERIAL_FACTORY_HH
#define SHIELDSIM_MATERIAL_FACTORY_HH

#include "MaterialLibrary.hh"

#include <map>
#include <string>

class G4Material;

namespace shield_sim {

// Converts fully resolved material records to Geant4 materials.  The factory
// always uses the resolved density and elemental mass fractions; it does not
// silently fall back to a similarly named NIST material after an override.
class MaterialFactory {
public:
  explicit MaterialFactory(const MaterialLibrary& library) : library_(library) {}
  G4Material* Build(const MaterialRecord& record);
  G4Material* Resolve(const std::string& idOrAlias, const std::string& role);

private:
  const MaterialLibrary& library_;
  std::map<std::string, G4Material*> built_;
};

} // namespace shield_sim

#endif
