#ifndef SEP_CORONA_SWCME_SURFACE_SHOCK_H
#define SEP_CORONA_SWCME_SURFACE_SHOCK_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_coronal_cme/shock_provider.h"

#include <memory>
#include <vector>

namespace SEP { namespace CoronaSwcme {

struct GeometryPatchState {
  CoronalCME::SurfacePatch geometry;
  double normalSpeedMPerS = 0.0;
  AmbientPrimitive upstream;
};

// BG3D-3 independently qualifies the front and shock state, but its nested
// fixed-fraction surface is not the material contact used by BG3D-4.  Carrying
// this status in the epoch prevents a geometrically complete legacy reference
// from being mistaken for a cross-stage interface authority.
enum class ContactAuthorityQualification {
  FixedFractionReferenceUnqualified,
  SharedMaterialAuthority
};

struct SurfaceShockEpoch {
  EventKinematics event;
  CoronalCME::FixedOrientationEllipsoid front;
  CoronalCME::FixedOrientationEllipsoid contact;
  std::vector<GeometryPatchState> frontPatches;
  std::vector<CoronalCME::SurfacePatch> contactPatches;
  std::shared_ptr<const CoronalCME::ShockSurfaceSnapshot> shocks;
  ContactAuthorityQualification contactAuthority =
      ContactAuthorityQualification::FixedFractionReferenceUnqualified;
  std::uint64_t backgroundGeneration = 0;
  std::string eventIdentity;
};

// Background-only geometry/shock assembly. The provider owns no source rates:
// all particle-related fields are zero and the local shock exists solely as
// both-sided plasma/interface state for later regional material maps.
class SurfaceShockModel final {
 public:
  static Core::Result<std::shared_ptr<SurfaceShockModel>> Create(
      std::shared_ptr<const EventConfiguration> event,
      std::shared_ptr<const AmbientModel> ambient);

  Core::Result<std::shared_ptr<const SurfaceShockEpoch>> Prepare(
      double epochS,std::uint64_t backgroundGeneration,
      int polarCells,int azimuthCells);

  std::shared_ptr<const SurfaceShockEpoch> Current() const { return current_; }

 private:
  std::shared_ptr<const EventConfiguration> event_;
  std::shared_ptr<const AmbientModel> ambient_;
  CoronalCME::TransactionalShockProvider shockProvider_;
  std::shared_ptr<const SurfaceShockEpoch> current_;
};

// The contact is a nested, rear-aligned ellipsoid. The frozen configuration's
// contact fraction scales the rear-to-apex span and both lateral semiaxes; it
// is not an independent reset radius or a radial shell subtraction.
Core::Result<CoronalCME::EllipsoidKinematics> ContactKinematics(
    const EventConfiguration&,const CoronalCME::EllipsoidKinematics& front);

const char* Name(ContactAuthorityQualification) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
