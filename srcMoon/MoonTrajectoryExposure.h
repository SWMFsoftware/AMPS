#ifndef SRCMOON_MOON_TRAJECTORY_EXPOSURE_H
#define SRCMOON_MOON_TRAJECTORY_EXPOSURE_H
#include "MoonReactionController.h"
#include <memory>

namespace Moon { namespace Chemistry { namespace Exposure {
// Immutable frozen-environment provider, constructed from the actual normal
// Moon orbit matrices, sphere/LOLA boundary and BVH. The caller binds source
// rates to catalog channel IDs; no profile, species state or plasma value is
// inferred from a compiled species name. Its lifetime must cover all queries.
class FrozenProductionField {
 public:
  static bool Create(const Species::TrajectoryExposure::Controls& controls,
      bool earthOccultation,std::shared_ptr<const FrozenProductionField>* output,
      std::string* error);
  // Normal pre-buffer source installation uses explicit ET queries because
  // step-owned global Sun/Earth/epoch caches are not initialized yet. No
  // normal time/cache is mutated, and all transforms share this exact ET.
  static bool CreateAtEpoch(const Species::TrajectoryExposure::Controls& controls,
      bool earthOccultation,double ephemerisTimeS,
      std::shared_ptr<const FrozenProductionField>* output,std::string* error);
  bool Build(const Controller::ExposureRequest& request,
      const std::vector<Species::TrajectoryExposure::PhotonChannel>& photons,
      const std::vector<Species::TrajectoryExposure::ConstantChannel>& constants,
      Species::TrajectoryExposure::History* output,std::string* error) const;
  const Species::TrajectoryExposure::Field& geometry() const{return field_;}
 private:
  Species::TrajectoryExposure::Field field_;
  Species::TrajectoryExposure::Controls controls_;
};
}}}
#endif
