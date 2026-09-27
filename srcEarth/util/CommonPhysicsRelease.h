#ifndef _SRC_EARTH_UTIL_COMMON_PHYSICS_RELEASE_H_
#define _SRC_EARTH_UTIL_COMMON_PHYSICS_RELEASE_H_

//======================================================================================
// CommonPhysicsRelease.h -- Roadmap Step 12 release identity
//======================================================================================
//
// Standalone and SWMF-coupled executables are allowed to differ in *where* the frozen
// magnetic field comes from.  They are not allowed to silently select different
// particle equations, units, characteristic integration, or product semantics.  This
// dependency-free header gives both builds one compile-time identity for that shared
// physics contract.  The tag is written to manifests and Tecplot AUXDATA, then checked
// by release_validation/run_release.py before cross-path results can be accepted.
//
// Changing any common particle-physics convention requires changing the tag and
// re-running the complete U/I/C/F/O validation matrix.  A build-system or source-tree
// SHA-256 is recorded separately by the Step-12 release manifest; the human-readable
// tag must never be used as a substitute for that immutable content digest.
//
// Phase 1 is intentionally limited to an instantaneous or quasi-static magnetic
// snapshot.  A dynamic E/B characteristic, long-duration trapping, local acceleration,
// and physical loss modelling are not capabilities implied by this release identity.
//======================================================================================

namespace Earth {
namespace CommonPhysicsRelease {

// Internal linkage keeps this C++11-compatible for the dependency-free unit targets;
// every translation unit still receives the identical literal selected at build time.
static constexpr const char* kTag =
    "sep-in-geospace-phase1-static-characteristics-v1";
static constexpr const char* kScope =
    "INSTANTANEOUS_QUASI_STATIC_MAGNETIC";

} // namespace CommonPhysicsRelease
} // namespace Earth

#endif // _SRC_EARTH_UTIL_COMMON_PHYSICS_RELEASE_H_
