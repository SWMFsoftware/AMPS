#ifndef SRCSEP_UTIL_SEP_ANALYTIC_RADIAL_BACKGROUND_H
#define SRCSEP_UTIL_SEP_ANALYTIC_RADIAL_BACKGROUND_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)

#include <array>

namespace SEP { namespace Background { namespace AnalyticRadial {

// An explicit stationary, constant-speed radial velocity declaration. Position
// and origin use the same Cartesian SI frame. This value is neither a new
// background authority nor a density/thermal closure: the frozen snapshot and
// its existing initialization configuration retain that ownership.
struct Wind {
  std::array<double,3> originM{{0,0,0}};
  double speedMPerS=0;
};

struct Sample {
  double radiusM=0;
  std::array<double,3> velocityMPerS{{0,0,0}};
  // G_ij = partial U_i / partial x_j [s^-1], row-major Cartesian components.
  std::array<double,9> gradientPerS{{0,0,0,0,0,0,0,0,0}};
};

// U=v*rhat and G=(v/r)(I-rhat*rhat), with r=x-origin. Every output is
// staged and checked; an undefined origin or unrepresentable derivative is an
// explicit error, leaving the caller's sample unchanged. No field or speed
// is inferred from storage values, magnetic direction or a provider fallback.
Transport::Status Evaluate(const Wind& wind,const double positionM[3],Sample* output);

// Compare native endpoint storage with the declared field, using a numerical
// relative tolerance scaled only by its positive speed. A disagreement never
// causes a replacement of imported/native velocity or a new gradient model.
Transport::Status ValidateVelocity(const Wind& wind,const double positionM[3],
    const double velocityMPerS[3],double relativeTolerance);

// Locale-independent, 17-digit binary64 declaration text for the existing
// background fingerprint owner. This formats explicit values only; it creates
// neither another hash/configuration owner nor a physical validation bypass.
// Mover selection is deliberately absent from this physical declaration.
std::string CanonicalDeclaration(const Wind& wind);

// Native-only adapter, defined beside the startup owner in main_lib.cpp. It
// reads the active initialization's origin/speed once per request; absence of
// runtime input selects the already documented legacy zero-origin/400 km/s
// declaration. It owns no mutable cache or second configuration fingerprint.
Wind DeclaredNativeWind();

// Native-only collective startup transaction: every MPI rank must call it,
// including nonanalytic providers (which collectively take the no-op path).
// Applies only to analytic generated domains;
// it validates all existing B/U first, then installs the existing n~r^-2 and
// identical previous B/U/n at the initial fixed-background epoch. It does not
// supply thermal values, overwrite imported primitives or publish a snapshot.
// All ranks agree on successful preparation before any datum is written; a
// local defect returns the same failing-rank/code/reason on every rank.
Transport::Status InitializeFixedNativeHistory();

} } }
#endif
