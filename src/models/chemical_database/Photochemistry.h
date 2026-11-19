#ifndef AMPS_CHEMICAL_DATABASE_PHOTOCHEMISTRY_H
#define AMPS_CHEMICAL_DATABASE_PHOTOCHEMISTRY_H

#include "DataIO.h"

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace AMPS { namespace ChemicalDatabase {

namespace Photochemistry {

// PHIDRATES supplies TWO distinct solar spectra, not a continuous activity
// coordinate. Reject other enum values; never invent an interpolation law.
enum class SolarCase { Quiet=0, Active=1 };

struct SourceFiles {
  std::string schema="srcMoon-D04-phidrates-v1";
  std::string channelsPath,ratesPath,thresholdsPath,crossSectionsPath;
  std::string channelsSha256,ratesSha256,thresholdsSha256,crossSectionsSha256;
};

// Return the reviewed v1 identities; callers still supply actual file paths.
// A future dataset requires a separately reviewed schema/version, not changing
// an input hash to silently authorize different numbers under the same name.
SourceFiles FrozenV1Sources(const std::string& processedPhidratesDirectory);

struct Product {
  std::string species;       // E.g. "H", "H2", "OH+", "O", "e".
  std::string electronicState="not_supplied_by_source";
  unsigned multiplicity=1;   // The O+H+H branch has H multiplicity TWO.
  int charge=0;             // In units of the positive elementary charge.
  std::map<std::string,unsigned> nuclei; // Per one product molecule/atom.
};

struct TabulatedDomain {
  // Inventory of ARCHIVED samples, including duplicate wavelengths and zero
  // cross sections. These extrema do not license interpolation across gaps,
  // extrapolation, threshold repair, or use as a cross-section evaluator.
  double minimumWavelengthNm=0.0,maximumWavelengthNm=0.0;
  double minimumPhotonEnergyEv=0.0,maximumPhotonEnergyEv=0.0;
  std::size_t sampleCount=0;
};

struct Channel {
  std::string id,speciesKey,physicalParent,reactionType;
  std::string sourceProducts,sourceChannelLabel,sourceReference,sourceFile;
  std::string parentState="not_supplied_by_source";
  bool reactive=false; // Molecular total_photoabsorption rows are NOT events.
  bool experimentalOh=false;
  std::vector<Product> products;
  std::array<double,2> rate1AuS{{0.0,0.0}};
  std::array<double,2> meanExcessEnergyEv{{0.0,0.0}};
  std::array<std::string,2> rateSourceFile,rateSourceRecord;
  bool hasThreshold=false; // False is missing metadata, never threshold=0.
  double thresholdEnergyEv=0.0,thresholdWavelengthNm=0.0;
  TabulatedDomain tabulated;
  double rateOpticalDepth=0.0,rateReferenceDistanceAu=1.0;
  bool daughterKineticsAvailable=false; // Mean energy is not a joint PDF.
};

class Registry {
 public:
  // Native CSV parsing and source validation are transactional. Every failing
  // call preserves the previous registry/output. The complete reviewed v1
  // catalogue has 17 rows: four atomic + seven H2O + four OH_exp reactions,
  // and two non-reactive absorption-total rows. No missing family is accepted.
  bool LoadFrozenD04(const SourceFiles& sources,std::string* error);
  bool GetChannel(const std::string& id,Channel* channel,std::string* error) const;
  bool ChannelsForParent(const std::string& speciesKey,
                         std::vector<Channel>* channels,std::string* error) const;
  bool BranchFractions(const std::string& speciesKey,SolarCase solar,
      std::vector<std::pair<std::string,double>>* fractions,std::string* error) const;
  bool SelectChannel(const std::string& speciesKey,SolarCase solar,double unitUniform,
                     std::string* channelId,std::string* error) const;
  bool ChannelRate1AuS(const std::string& id,SolarCase solar,
                       double* rateS,std::string* error) const;
  // visiblePhotonFraction is a normalized GAS photon-visibility operator at
  // the actual particle position/epoch, supplied by the production geometry
  // adapter. It is not a facet-incidence cosine, surface albedo or temperature.
  // Only documented r^-2 scaling and this explicit [0,1] factor are applied.
  bool TotalDestructionRateS(const std::string& speciesKey,SolarCase solar,
      double distanceAu,double visiblePhotonFraction,double* rateS,std::string* error) const;
  bool loaded() const {return loaded_;}
  std::size_t catalogueSize() const {return loaded_ ? channels_.size() : 0;}
  std::size_t reactiveSize() const {return loaded_ ? 15 : 0;}
  const SourceFiles& sources() const {return sources_;}

 private:
  bool loaded_=false;
  SourceFiles sources_;
  std::vector<Channel> channels_;
  std::map<std::string,std::size_t> byId_;
};

} } } // namespace AMPS::ChemicalDatabase::Photochemistry
#endif
