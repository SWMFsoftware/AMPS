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
  // VERIFICATION ONLY (I21 decision 1, approved 2026-11-26): the same frozen
  // field, but with a caller-prescribed static point-Sun position [m, frozen
  // LSO] instead of the SPICE Sun at this ET. Body transforms still come from
  // SPICE at the same ET. Earth occultation is not available, because a
  // prescribed Sun has no consistent Earth geometry. The Sun must lie outside
  // one solar radius of the Moon's centre, as for the normal field.
  static bool CreatePrescribedSunAtEpoch(const Species::TrajectoryExposure::Controls& controls,
      const Species::TrajectoryExposure::Vector& sunPositionM,double ephemerisTimeS,
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
