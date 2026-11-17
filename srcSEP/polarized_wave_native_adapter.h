#ifndef SRCSEP_POLARIZED_WAVE_NATIVE_ADAPTER_H
#define SRCSEP_POLARIZED_WAVE_NATIVE_ADAPTER_H
#include "util/sep_polarized_wave_state.h"
#include "util/sep_turbulence_core.h"
#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_background_snapshot.h)

namespace SEP { namespace Turbulence { namespace Polarized { namespace PICAdapter {

// Current primitives consumed by the native wave adapter, independently read
// from native/decoded-checkpoint endpoints. This does not replace the enclosing
// checkpoint's validation of thermal and previous-epoch vertex slots.
struct NativeEndpoint {
  std::array<double,3> positionM{{0,0,0}},magneticFieldT{{0,0,0}},plasmaVelocityMPerS{{0,0,0}};
  double numberDensityM3=0;
};
struct NativePhysicalCell {
  std::uint64_t fieldLine=0,segment=0;
  double volumeM3=0,meanMagneticFieldT=0,massDensityKgM3=0;
  std::array<double,3> positionM{{0,0,0}},tangent{{0,0,0}},magneticFieldT{{0,0,0}},
      plasmaVelocityMPerS{{0,0,0}};
  NativeEndpoint begin,end;
};
Transport::Status ReadNativePhysicalCells(const Background::BackgroundSnapshot& background,
    std::vector<NativePhysicalCell>* cells);
Transport::Status ValidatePhysicalImage(const Snapshot& waves,
    const std::vector<NativePhysicalCell>& backgroundCells);
std::string NativeBackgroundIdentity(const Background::BackgroundSnapshot& background);
std::string NativePhysicalImageIdentity(const std::vector<NativePhysicalCell>& cells);

struct NativeDiagnostics {
  std::uint64_t cells=0,bins=0;
  double maximumBranchRelativeRoundoff=0;
  double maximumAbsoluteInitialRoundoffJ=0;
  bool mpiIdentityAgrees=false,readbackAgrees=false;
};
class NativeStateOwner;
class NativeRestorePlan {
 public:
  NativeRestorePlan();~NativeRestorePlan();
  NativeRestorePlan(NativeRestorePlan&&)noexcept;
  NativeRestorePlan& operator=(NativeRestorePlan&&)noexcept;
  NativeRestorePlan(const NativeRestorePlan&)=delete;
  NativeRestorePlan& operator=(const NativeRestorePlan&)=delete;
 private:
  friend class NativeStateOwner;
  struct Implementation;
  std::unique_ptr<Implementation> implementation_;
};

// Every native initialization is collective and quiescent. All four arrays
// remain in one immutable owner. Native integrated/two-branch arrays are only
// derived diagnostics after successful publication. Neither this owner nor a
// restored image is physically qualified for production evolution yet.
class NativeStateOwner {
 public:
  NativeStateOwner();
  Transport::Status InitializeFromIntegrated(const Configuration& configuration);
  std::shared_ptr<const Snapshot> Acquire()const;
  NativeDiagnostics Diagnostics()const;
  // Stage owns all allocations and validates independently decoded native
  // primitives, exact physical grid, provenance, geometry and four-array hash.
  // It never writes host storage. The outer checkpoint transaction must call
  // ConfirmCollectiveReady after every rank stages successfully and must keep
  // the model quiescent through CommitRestore. Commit allocates nothing on its
  // valid path. This component is not a complete checkpoint codec/transaction.
  Transport::Status StageRestore(const Snapshot& decoded,
      const Background::BackgroundSnapshot& restoredBackground,
      const std::vector<NativePhysicalCell>& restoredBackgroundCells,
      NativeRestorePlan* plan)const;
  Transport::Status ValidateReadyForCommit(const NativeRestorePlan& plan)const;
  // Collective, invoked on every rank after all ranks staged successfully.
  // Rechecks local readiness and actual restored identity agreement. The outer
  // transaction must remain quiescent from this gate through all host commits.
  Transport::Status ConfirmCollectiveReady(const NativeRestorePlan& plan)const;
  // An unready, owner/configuration/buffer-stale or wrong-owner plan returns
  // false before any write. A valid collectively qualified plan commits
  // without allocation and consumes itself. Endpoint data must remain frozen
  // by the enclosing quiescent transaction after the collective readiness gate.
  [[nodiscard]] bool CommitRestore(NativeRestorePlan* plan)noexcept;
 private:
  mutable std::mutex mutex_;
  std::shared_ptr<const Snapshot> current_;
  NativeDiagnostics diagnostics_;
  const std::uint64_t lifetimeId_;
};

// Call before PIC::FieldLine::Init activates segment offsets. The one resolved
// configured N determines derived 2*N diagnostic storage; no old live buffer
// is resized/remapped. The initial L/R policy is explicitly 0.5/0.5 in each
// branch as selected by the model owner, not a subsequent equality constraint.
Transport::Status ConfigureNativeDiagnosticsBeforeFieldLineInit();
Transport::Status InitializeNativePolarizedState();
bool HasNativePolarizedState();
std::shared_ptr<const Snapshot> AcquireNativePolarizedState();
NativeStateOwner& NativePolarizedStateOwner();
Transport::Status LegacyMutationPreflight(const char* operation);
// The explicit-owner overload permits a provider/fixture to query exactly the
// same immutable authority capability without publishing process-wide state.
Transport::Status LegacyMutationPreflight(const char* operation,const NativeStateOwner& authority);
// Legacy void APIs cannot return a typed capability. They call this once at
// entry, before cache changes, particle/wave writes or MPI operations.
void RequireLegacyMutationAllowed(const char* operation);
Transport::Status LegacyFixedSpectralPreflight(const Configuration& configuration,
    const char* operation);
Transport::Status CopyLegacyNativeSpectrum(const Configuration& configuration,
    const double* source,std::size_t sourceLength,std::vector<double>* destination);

} } } }
#endif
