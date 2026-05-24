#ifndef _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_
#define _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_

//======================================================================================
// Mode3DForwardSWMF.h
//======================================================================================
//
// PURPOSE
// -------
// Thin SWMF-coupled initialization bridge for the 3d_forward energetic-particle
// branch.
//
// This directory intentionally contains only the AMPS-facing initialization hook
// below.  The actual forward-particle implementation remains in srcEarth/3d_forward:
// particle injection, time-step/weight calculation, trajectory tracking, density
// sampling, sphere-flux sampling, and mover selection are all reused from there.
//
// The only intended physical difference from standalone 3d_forward is the source of
// the background fields.  In an SWMF-coupled build, the 3d_forward movers obtain B and
// E from the SWMF coupler through PIC::CPLR rather than from standalone analytic or
// DATAFILE-populated fields.  Because the normal standalone argument line is not
// available in a coupled SWMF run, amps_init() reads the conventional local input
// file AMPS_PARAM.in.
//
// LIFECYCLE
// ---------
// main_lib.cpp::amps_init() calls Earth::Mode3DForwardSWMF::amps_init() when
// _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_.  The bridge performs the coupled
// initialization directly by calling existing srcEarth/3d_forward helpers and
// callbacks; it does not depend on a separate coupled-runtime wrapper.
//
//======================================================================================

namespace Earth {
namespace Mode3DForwardSWMF {

// Initialize the 3d_forward runtime for SWMF-coupled execution.
void amps_init();

} // namespace Mode3DForwardSWMF
} // namespace Earth

#endif // _SRC_EARTH_3D_FORWARD_SWMF_MODE3DFORWARDSWMF_H_
