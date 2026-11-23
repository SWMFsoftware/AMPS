#ifndef SRCSEP_REDUCED_SHOCK_BOUNDARY_H
#define SRCSEP_REDUCED_SHOCK_BOUNDARY_H
#include "reduced_shock_particles.h"
#include "util/sep_absorbing_bridge.h"
namespace SEP { namespace ReducedParticles { namespace Native {
struct BoundaryInterval {
  bool active=false;
  double durationS=0,particleArcM=0,frontBeginArcM=0,frontEndArcM=0;
  int upstreamSign=1;
};
Transport::Status AbsorbInitialReturn(const Transport::PICAdapter::ParticleContext&,
    double momentumKgMPerS,bool* removed);
// Predictable refinement happens BEFORE the particle Gaussian/event draws.
// Only an internal interval is shortened: the frozen physical host endpoint,
// its sampling exposure and each birth's remaining duration stay unchanged.
Transport::Status PrepareBoundary(const Transport::PICAdapter::ParticleContext&,
    double elapsedS,double* numericalIntervalS,BoundaryInterval*,
    Transport::BoundaryIntervalControl* history=nullptr);
// The Euler Brownian interval uses its exact conditional bridge law. A focused
// deterministic interval uses its accepted linear streaming interpolant. MFP
// invokes this for EVERY between-event interval, before any endpoint scatter.
Transport::AbsorbingBridgeResult CheckBoundary(const BoundaryInterval&,
    double beginDisplacementM,double endDisplacementM,double elapsedBeginS,
    double elapsedEndS,double kappaTimesDurationM2,Transport::RandomStream&);
} } }
#endif
