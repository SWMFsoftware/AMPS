#ifndef SHIELDSIM_MATERIAL_LIBRARY_HH
#define SHIELDSIM_MATERIAL_LIBRARY_HH

/*
 * Versioned material records used by Stage-1 case files.
 *
 * A MaterialLibrary is immutable after loading.  Case-local overrides create
 * new records with new IDs; they never modify the shared library entry.  Each
 * resolved record has a SHA-256 hash of its canonical, fully explicit JSON.
 */

#include "StrictJson.hh"

#include <map>
#include <string>
#include <vector>

namespace shield_sim {

enum class MaterialStatus { Reference, TradeStudy, Template, ProjectSpecific };

struct MaterialComponent {
  enum class Kind { Element, Isotope, Material };
  Kind kind = Kind::Element;
  std::string id;
  double value = 0.0;
};

struct MaterialRecord {
  std::string id;
  std::string revision;
  std::vector<std::string> aliases;
  std::string label;
  MaterialStatus status = MaterialStatus::Template;
  std::vector<std::string> roles;
  double densityGPerCm3 = 0.0;
  std::string phase;
  double temperatureK = 0.0;
  bool hasPressure = false;
  double pressurePa = 0.0;
  std::string compositionBasis;
  std::vector<MaterialComponent> components;
  std::string residualPolicy;
  std::string geant4Name;
  bool hasMeanExcitationEnergy = false;
  double meanExcitationEnergyEV = 0.0;
  std::string applicability;
  std::string parent;
  Json provenance;
  Json canonicalRecord;
  std::string contentHash;
};

class MaterialLibrary {
public:
  static MaterialLibrary Load(const std::string& path,
                              const Json::Array& caseOverrides = Json::Array{});

  const MaterialRecord& Resolve(const std::string& idOrAlias,
                                const std::string& requiredRole) const;
  const MaterialRecord& ById(const std::string& id) const;
  std::map<std::string, double> ElementalMassFractions(const MaterialRecord& record) const;

  const std::string& Id() const { return id_; }
  const std::string& Version() const { return version_; }
  const std::string& SourcePath() const { return sourcePath_; }
  const std::string& ContentHash() const { return contentHash_; }
  const Json& CanonicalDocument() const { return canonicalDocument_; }
  std::vector<std::string> MaterialIds() const;

private:
  static MaterialRecord ParseRecord(const Json& value, const std::string& context);
  static Json MergeOverride(const MaterialRecord& parent, const Json& overrideValue,
                            const std::string& context);
  void Insert(MaterialRecord record, const std::string& context);
  std::map<std::string, double> Flatten(const MaterialRecord& record,
                                       std::vector<std::string>& stack) const;

  std::string id_;
  std::string version_;
  std::string sourcePath_;
  std::string contentHash_;
  Json canonicalDocument_;
  std::map<std::string, MaterialRecord> records_;
  std::map<std::string, std::string> aliases_;
};

std::string MaterialStatusName(MaterialStatus status);

} // namespace shield_sim

#endif
