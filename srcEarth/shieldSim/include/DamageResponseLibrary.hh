#ifndef SHIELDSIM_DAMAGE_RESPONSE_LIBRARY_HH
#define SHIELDSIM_DAMAGE_RESPONSE_LIBRARY_HH

/*
 * Stage-5 displacement-damage response library.
 *
 * The library deliberately separates two response representations:
 *   - positive charged-particle NIEL nodes, interpolated log-log; and
 *   - ASTM neutron groups, evaluated as group constants on explicit edges.
 *
 * Neither representation extrapolates. Callers receive an explicit status so
 * an out-of-domain fluence contribution cannot be mistaken for zero damage.
 */

#include "StrictJson.hh"

#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace shield_sim {

enum class DamageResponseStatus { Available, UnsupportedSpecies, OutOfDomain };

std::string DamageResponseStatusName(DamageResponseStatus status);

struct DamageResponseValue {
  DamageResponseStatus status = DamageResponseStatus::UnsupportedSpecies;
  double nielMeVCm2PerG = 0.0;
  // Defined only for ASTM silicon neutron groups. It is the group damage
  // function divided by the prescribed 95 MeV-mb 1-MeV reference.
  double neutronHardnessFactor = 0.0;
  std::size_t intervalIndex = 0;
  std::string reason;
};

struct DifferentialFluenceBin {
  double lowerMeV = 0.0;
  double upperMeV = 0.0;
  double differentialFluencePerCm2PerMeV = 0.0;
};

struct DamageFoldResult {
  DamageResponseStatus status = DamageResponseStatus::Available;
  double fluencePerCm2 = 0.0;
  double dddMeVPerG = 0.0;
  // Nonzero only for a neutron fold. Proton and He-4 DDD are intentionally
  // never converted into an E722 neutron-equivalent fluence.
  double neutronEquivalentFluencePerCm2 = 0.0;
  std::string reason;
};

class DamageResponseLibrary {
public:
  static DamageResponseLibrary Load(const std::string& manifestPath);

  DamageResponseValue Evaluate(const std::string& species,double energyMeV) const;

  // Integrate a piecewise-constant differential fluence over explicit source
  // bin edges. Charged-particle log-log segments are integrated analytically;
  // neutron bins are intersected exactly with the E722 group edges.
  DamageFoldResult FoldPiecewiseConstantDifferential(
      const std::string& species,
      const std::vector<DifferentialFluenceBin>& bins) const;

  bool SupportsSpecies(const std::string& species) const;
  const std::string& LibraryId() const { return libraryId_; }
  const std::string& LibraryVersion() const { return libraryVersion_; }
  const std::string& TargetMaterialId() const { return targetMaterialId_; }
  const std::string& ManifestPath() const { return manifestPath_; }
  const std::string& ManifestHash() const { return manifestHash_; }
  double ReferenceNeutronNielMeVCm2PerG() const { return referenceNeutronNiel_; }
  const Json& CanonicalManifest() const { return canonicalManifest_; }

private:
  struct Table {
    std::string species;
    std::string representation;
    std::string path;
    std::string contentHash;
    double domainMinMeV = 0.0;
    double domainMaxMeV = 0.0;
    std::vector<double> energyMeV;
    std::vector<double> nielMeVCm2PerG;
    std::vector<double> lowerMeV;
    std::vector<double> upperMeV;
    std::vector<double> neutronHardness;
  };

  const Table* FindTable(const std::string& species) const;

  std::string libraryId_;
  std::string libraryVersion_;
  std::string targetMaterialId_;
  std::string manifestPath_;
  std::string manifestHash_;
  double referenceNeutronNiel_ = 0.0;
  Json canonicalManifest_;
  std::map<std::string,Table> tables_;
};

} // namespace shield_sim

#endif
