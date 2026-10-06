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

// BG3D-3/Level A independently qualifies a prescribed front and local shock
// state, but its nested fixed-fraction surface is not the material contact used
// by the old BG3D-4 map and is not the selected Level B piston.  Carrying this
// status in the epoch prevents a geometrically complete legacy reference from
// being mistaken for a cross-stage interface authority.  Under the reviewed
// BG3D-4 replacement, a reconstructed ejecta body becomes the sole contact;
// its per-ray plasma evolution computes the Level B shock.  This prescribed
// front then remains a validation target, not a second production authority.
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

// Legacy/unqualified reference only.  The contact is a nested, rear-aligned
// ellipsoid whose fraction scales the rear-to-apex span and both lateral
// semiaxes; it is not an independent reset radius or a radial shell
// subtraction.  Production piston geometry must instead come from the one
// checksummed, differentiable ejecta/contact reconstruction and must pass
// star-shaped, non-grazing and ambient-coverage gates ray by ray.
Core::Result<CoronalCME::EllipsoidKinematics> ContactKinematics(
    const EventConfiguration&,const CoronalCME::EllipsoidKinematics& front);

const char* Name(ContactAuthorityQualification) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
