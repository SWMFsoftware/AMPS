#ifndef SEP_CORONA_SWCME_PISTON_SHEATH_H
#define SEP_CORONA_SWCME_PISTON_SHEATH_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_corona_swcme/piston_ambient.h"
#include "sep_corona_swcme/piston_contact.h"
#include "sep_corona_swcme/piston_solver.h"

#include <cstdint>
#include <memory>
#include <vector>

namespace SEP { namespace CoronaSwcme {

enum class PistonBackgroundRegion {
  // No Level-B tube exists on a ray excluded by contact/background support.
  UnsupportedFlank,
  // r<R_c belongs to the future BG3D-5 ejecta closure; never ambient-filled.
  EjectaNotOwned,
  // Resolved wave/shock transition ahead of the Q-centroid, or a sub-fast
  // disturbed column when no numerical shock exists.
  Compression,
  // Downstream material between the one contact authority and shock centroid.
  Sheath,
  // Either an evolved buffer cell below its disturbance threshold or a direct
  // ambient-authority sample beyond the finite material domain.
  Ambient
};

struct PistonRayReceipt {
  // All radii are heliocentric HCI metres at the committed epoch.  Receipts
  // describe owned material and diagnostics only; they do not transfer
  // ownership or permit mutation of the solver state.
  std::uint64_t rayId = 0;
  PistonRayDisposition disposition = PistonRayDisposition::NoIntersection;
  double solidAngleSr = 0.0;
  double contactRadiusM = 0.0;
  double outerRadiusM = 0.0;
  double disturbanceRadiusM = 0.0;
  std::size_t cells = 0;
  std::size_t ambientBufferCells = 0;
  PlanarShockState shock;
};

struct PistonBackgroundSample {
  // plasma carries the evolved material-cell state when radius lies inside
  // the tube; ambient carries the direct authority state outside it.  The
  // plasmaAvailable flag is false for unsupported and ejecta-not-owned
  // requests, preventing quiet ambient from masquerading as missing CME state.
  PistonBackgroundRegion region = PistonBackgroundRegion::UnsupportedFlank;
  bool plasmaAvailable = false;
  PlanarCellState plasma;
  AmbientPrimitive ambient;
  PistonRayState contact;
  PlanarShockState shock;
};

// One transactional collection of independent radial Level-B tubes.  The
// analytical ejecta body is the sole piston/contact authority.  The Level-A
// history is deliberately absent from this class; shock location and jump
// states come only from each evolved tube.
class PistonSheathModel final {
 public:
  ~PistonSheathModel();
  static Core::Result<std::unique_ptr<PistonSheathModel>> Create(
      std::shared_ptr<const EventConfiguration> event);

  // Multi-ray epoch commit: every supported solver is cloned, advanced and
  // buffer-checked first.  One failed ray or append leaves epoch, generation
  // and every committed inventory unchanged.
  Core::Status AdvanceTo(double epochS);
  Core::Result<std::vector<PistonRayReceipt>> Receipts() const;
  Core::Result<PistonBackgroundSample> QueryRay(
      std::uint64_t rayId,double radiusM) const;

  double EpochS() const noexcept { return epochS_; }
  std::uint64_t Generation() const noexcept { return generation_; }

 private:
  struct Ray;
  std::shared_ptr<const EventConfiguration> event_;
  std::shared_ptr<const AmbientModel> ambient_;
  std::shared_ptr<const PistonContactModel> contact_;
  PistonNumericsInput controls_;
  std::vector<std::unique_ptr<Ray>> rays_;
  double epochS_ = 0.0;
  std::uint64_t generation_ = 0;
};

const char* Name(PistonBackgroundRegion region) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
