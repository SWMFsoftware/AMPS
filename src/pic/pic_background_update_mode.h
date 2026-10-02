// Background update ownership is independent of the native field buffer ABI.
// Include after global.h so AMPS' host/device and managed-storage macros exist.
#ifndef AMPS_PIC_BACKGROUND_UPDATE_MODE_H
#define AMPS_PIC_BACKGROUND_UPDATE_MODE_H

namespace PIC {
namespace CPLR {
namespace DATAFILE {

// FileSchedule preserves the historical DATAFILE behavior: MULTIFILE controls
// loading, interpolation between file epochs, and termination at the last file.
// RuntimeProvider means the application publishes complete field snapshots at
// joined step boundaries. The DATAFILE offsets still describe the allocated
// center-node bytes, but there is no file schedule, next-file index or file EOF.
// SWCME, Parker and future runtime models share this policy; model selection
// belongs to the application's provider factory rather than this storage API.
enum class BackgroundUpdateMode {
  FileSchedule,
  RuntimeProvider
};

// The definition in pic_datafile.cpp defaults to FileSchedule, preserving other
// AMPS applications. Set the policy on every rank after storage initialization
// and before the first publication/read/time step. Do not change it while any
// host thread or device kernel is consuming a snapshot. Managed storage matches
// CurrDataFileOffset: the host/device field getter must see the same ownership.
extern BackgroundUpdateMode _TARGET_DEVICE_ _CUDA_MANAGED_
    BackgroundUpdatePolicy;

_TARGET_HOST_ _TARGET_DEVICE_
inline bool UsesFileSchedule() {
  return BackgroundUpdatePolicy == BackgroundUpdateMode::FileSchedule;
}

} // namespace DATAFILE
} // namespace CPLR
} // namespace PIC

#endif // AMPS_PIC_BACKGROUND_UPDATE_MODE_H
