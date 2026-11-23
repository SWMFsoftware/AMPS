#ifndef SRCSEP_UTIL_SEP_BACKEND_ADMISSION_H
#define SRCSEP_UTIL_SEP_BACKEND_ADMISSION_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_status.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_material_background.h)

namespace SEP { namespace Background { namespace Admission {

// This inventory names existing authorities, not substitute background laws.
// It is intentionally independent of mover choice: selecting another solver
// cannot manufacture material identities, temporal coverage or a source law.
enum class Backend { ReducedSseParker, Swcme1D, CoronalBundle, FullCoronaSwcme };

// New production callers bind the exact canonical profile and physical region.
// Keep the original backend overload below as its source-compatible inventory
// API; neither overload admits an unqualified family by capability assertion.
inline Core::Status RequireMovingTransport(Material::Profile profile,
    Material::Region region=Material::Region::ParkerExterior) {
  return Material::Require(profile,Material::Use::MaterialTransport,region);
}

inline Core::Status RequireMovingTransport(Backend backend) {
  // Success here is only the backend-level prerequisite. The reduced adapter
  // still checks its actual profile, immutable interval, frame/composition,
  // supported roots, explicit source and absorbing-boundary contract. It must
  // never interpret this inventory as permission to use PFSS/downstream state.
  switch(backend) {
    case Backend::ReducedSseParker:
      return Core::Status::Success();
    case Backend::Swcme1D:
      // The maintained SHOCK_ONLY profile uses constant radial U and the
      // Leblanc r^-2+r^-4+r^-6 density. Consequently r^2 rho U is not constant:
      // simply moving vertices with U does not produce a source-free material
      // tube. FULL_ICME also lacks a complete exported regional thermal state.
      // Neither reconstructing density from another provider nor adding U to
      // relative particle displacement supplies the missing core contract.
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
          "SWCME1D moving SEP transport is unsupported: the maintained Leblanc-density/constant-wind profile has no conservative material-state evolution contract; FULL_ICME additionally lacks exported regional thermodynamics. Diagnostic core queries remain available.");
    case Backend::CoronalBundle:
      // Bundle history describes connection/observer coverage; it is NOT a
      // sequence of plasma states. The current serialized neutral ABI has one primitive/geometry epoch,
      // no material vertex lineage or time-dependent grad(U), and no explicit
      // species-density basis. Extending historyEndS cannot cure these gaps.
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
          "Coronal-CME moving SEP transport is unsupported: the neutral single-epoch bundle has no material trajectory/vertex identity, time-dependent primitive and velocity-gradient coverage, or explicit species-density contract. Snapshot exchange remains available.");
    case Backend::FullCoronaSwcme:
      // This is a core qualification prerequisite, not an adapter TODO which
      // may be bypassed by painting local RH limits throughout a CME volume.
      return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
          "Full corona-SWCME moving SEP transport is unsupported: BG3D-4 regional closure and BG3D-5/6 material state/handoff plus native lifecycle qualification are incomplete. Local ambient/shock queries do not qualify a downstream volume.");
  }
  return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
      "Unknown moving SEP background backend; no provider fallback is permitted");
}

} } }
#endif
