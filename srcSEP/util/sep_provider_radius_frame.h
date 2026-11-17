#ifndef SRCSEP_UTIL_SEP_PROVIDER_RADIUS_FRAME_H
#define SRCSEP_UTIL_SEP_PROVIDER_RADIUS_FRAME_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_background_snapshot.h)

namespace SEP { namespace Background { namespace RadiusFrame {
// One immutable semantic identity for the provider-defined coefficient radius.
// Analytic origin VALUES remain owned by Initialization/DeclaredNativeWind;
// SWCME and reduced models already declare a global heliocentric frame. No
// origin follows from a raw SWMF coordinate array or a local-evolution copy.
// Include this contract in background fingerprints and coefficient provenance
// so same-wire-layout checkpoints cannot resume the old global-zero analytic
// coefficient interpretation under an indistinguishable physical identity.
inline const char* ContractIdentity(Provider provider) {
  switch(provider) {
    case Provider::Analytic:return "provider-radius-v1:analytic-frozen-declared-origin";
    case Provider::Swcme:return "provider-radius-v1:swcme-global-heliocentric-origin-zero";
    case Provider::ReducedShock:return "provider-radius-v1:reduced-global-heliocentric-origin-zero";
    case Provider::Swmf:return "provider-radius-v1:swmf-origin-unresolved";
    case Provider::LocalEvolution:return "provider-radius-v1:local-handoff-origin-unresolved";
  }
  return "provider-radius-v1:unknown-provider-origin-unresolved";
}
} } }
#endif
