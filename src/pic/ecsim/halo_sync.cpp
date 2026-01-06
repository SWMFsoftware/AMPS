// =====================================================================================
// pic_field_solver_ecsim_halo_sync.cpp
//
// Implementation of ECSIM halo sync wrappers.
// Requires:
//   - PIC::Parallel::cNodeHaloSyncManager, eNodeType, eHaloOp
//   - PIC::Parallel::SyncNodeHalo_DomainBoundaryLayer()
//   - PIC::CPLR::DATAFILE::Offset::<...>.RelativeOffset constants for the stored vectors
//
// =====================================================================================

#include "../pic.h"
#include "halo_sync.h"

namespace PIC {
namespace FieldSolver {
namespace Electromagnetic {
namespace ECSIM {

// These are declared extern in the header; define them in your ECSIM module.
// int CurrentEOffset = ...; etc.

void SyncE() {
  using PIC::Parallel::cNodeHaloSyncManager;
  using PIC::Parallel::eNodeType;
  using PIC::Parallel::eHaloOp;

  // ---- Sync Current E (Ex,Ey,Ez) ----
  {
    cNodeHaloSyncManager m;
    m.NodeType = eNodeType::Corner;
    m.RelativeOffsetBytesFromAssociated =
        PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset + CurrentEOffset;
    m.nDoubles = 3;
    m.Op = eHaloOp::Replace;

    // Unique tagBase for this exchange
    PIC::Parallel::SyncNodeHalo_DomainBoundaryLayer(m, /*tagBase=*/42000);
  }

  // ---- Sync Half-step E (Ex,Ey,Ez) ----
  {
    cNodeHaloSyncManager m;
    m.NodeType = eNodeType::Corner;
    m.RelativeOffsetBytesFromAssociated =
        PIC::CPLR::DATAFILE::Offset::ElectricField.RelativeOffset + OffsetE_HalfTimeStep;
    m.nDoubles = 3;
    m.Op = eHaloOp::Replace;

    PIC::Parallel::SyncNodeHalo_DomainBoundaryLayer(m, /*tagBase=*/42010);
  }
}

void SyncB() {
  using PIC::Parallel::cNodeHaloSyncManager;
  using PIC::Parallel::eNodeType;
  using PIC::Parallel::eHaloOp;

  // ---- Sync Current B (Bx,By,Bz) ----
  {
    cNodeHaloSyncManager m;
    m.NodeType = eNodeType::Center;
    m.RelativeOffsetBytesFromAssociated =
        PIC::CPLR::DATAFILE::Offset::MagneticField.RelativeOffset + CurrentBOffset; // adjust if needed
    m.nDoubles = 3;
    m.Op = eHaloOp::Replace;

    PIC::Parallel::SyncNodeHalo_DomainBoundaryLayer(m, /*tagBase=*/42100);
  }
}



} // namespace ECSIM
} // namespace Electromagnetic
} // namespace FieldSolver
} // namespace PIC

