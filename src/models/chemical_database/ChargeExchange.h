#ifndef AMPS_CHEMICAL_DATABASE_CHARGE_EXCHANGE_H
#define AMPS_CHEMICAL_DATABASE_CHARGE_EXCHANGE_H

#include <cstddef>
#include <string>
#include <vector>

namespace AMPS { namespace ChemicalDatabase { namespace ChargeExchange {

// The ID explicitly selects its source/law. Na log-log interpolation is an
// approved named MODEL between the populated source knots, not extra measured
// data. Neither channel supplies a daughter state/angular distribution.
enum class ChannelId { ProtonHeOrnlA32, ProtonNaMorgan1985LogLog };

struct ReferenceFiles {
  std::string schema="srcMoon-total-capture-v1";
  std::string heliumPublicationPath,sodiumPublicationPath;
  std::string heliumPublicationSha256,sodiumPublicationSha256;
};
ReferenceFiles FrozenV1References(const std::string& researchDirectory);

struct ChannelMetadata {
  std::string id,reactants,products,source,sourceSha256,law;
  std::string energyConvention="target-rest kinetic energy per projectile u";
  double minimumEnergyEvPerU=0.0,maximumEnergyEvPerU=0.0;
  bool totalCapture=true;
  bool daughterKineticsAvailable=false;
  bool productionParticleExecutionEnabled=false;
};

class Registry {
 public:
  // Bind the original publication identities before publishing source-derived
  // constants. A file with a PDF-like name or a caller-updated hash is not a
  // newly reviewed v1 source. Failed load/query calls leave outputs unchanged.
  bool LoadPublishedReferences(const ReferenceFiles& references,std::string* error);
  bool GetMetadata(ChannelId id,ChannelMetadata* metadata,std::string* error) const;
  bool CrossSectionM2(ChannelId id,double energyEvPerU,
                      double* crossSectionM2,std::string* error) const;
  // This is ONLY a declared monoenergetic-beam rate n*sigma(E)*g. Energy and
  // relative speed are explicit inputs: no projectile mass, classical/high-
  // energy mapping, thermal bath or Maxwellian tail policy is silently chosen.
  bool ColdBeamRateS(ChannelId id,double energyEvPerU,double relativeSpeedMps,
      double projectileDensityM3,double* rateS,std::string* error) const;
  bool SodiumSourceKnots(std::vector<double>* energyEvPerU,
                          std::vector<double>* crossSectionM2,std::string* error) const;
  bool loaded() const {return loaded_;}
 private:
  bool loaded_=false;
  ReferenceFiles references_;
};

} } } // namespace AMPS::ChemicalDatabase::ChargeExchange
#endif
