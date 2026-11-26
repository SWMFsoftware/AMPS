#ifndef SRCSEP_REDUCED_PARTICLE_OPTIONS_H
#define SRCSEP_REDUCED_PARTICLE_OPTIONS_H

#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>

namespace SEP { namespace ReducedParticles {

// This application-owned contract selects a consumer of the existing reduced
// background. It does not change the core event's equations or manufacture a
// downstream volume. Keep return and source-population semantics independent.
enum class ReturnPolicy { Delete, Pass };
enum class Normalization { NetUpstreamEscape, GrossEmission };
struct Configuration {
  bool enabled=false;
  ReturnPolicy returnPolicy=ReturnPolicy::Delete;
  Normalization normalization=Normalization::NetUpstreamEscape;
  // L_ref is a physical signed NORMAL distance [m], not a segment fraction,
  // diffusion length or numerical placement epsilon. Zero means unspecified;
  // enabled runs must supply a strictly positive value without a default.
  double releaseDistanceM=0;
  // Alternative L_ref law: L_ref(t) = f * R_curv(t), with R_curv the radius of
  // the front's generating sphere at the birth epoch (the heliocentric shock
  // radius for the v1.3 sphere). Zero means unspecified. Exactly one of the
  // two release laws is given for an enabled source.
  double releaseDistanceShockRadiusFraction=0;
  // The existing source:shock block owns total kinetic-energy bounds [J] and
  // net release fraction. These fields never borrow a SWCME plasma/compression.
  double minimumEnergyJ=0,maximumEnergyJ=0,releaseFraction=0;
  std::uint64_t samplesPerStep=0,campaignSeed=0;
  // Numerical control on curved moving-front linearization. Source times use
  // stratified Monte Carlo integration over the full physical window; neither
  // operation changes AMPS dt or adds a sampling opportunity.
  double boundaryToleranceM=1.0;
  // Zero selects the production core DSA law. A positive SI energy is an
  // explicitly labelled manufactured monoenergetic verification source; it
  // cannot be mistaken for a physical production default or calibration.
  double verificationMonoenergeticEnergyJ=0;
  // Opt-in importance sampling of the production DSA momentum law: proposals
  // are log-uniform in p on [pmin,pmax] and each carries the exact weight
  // f(p)/g(p)=f(p) p ln(pmax/pmin), whose proposal mean is one. The physical
  // number, energy and spectrum are unbiased; only Monte Carlo variance moves
  // from the flux-dominant low energies toward the sparse high-energy tail.
  bool logUniformMomentumSampling=false;
  // Opt-in numerical representation of the absorbing shock's material
  // coordinate during one AMPS step: a piecewise-linear table in time whose
  // nodes are exact core roots and whose every interval is verified, at its
  // quarter, middle and three-quarter epochs, against exact roots to within
  // boundaryToleranceM (the same tolerance that already accepts the linear
  // front of each particle interval). Off: every query solves an exact root.
  bool boundaryRootTable=false;
};
inline const char* Name(ReturnPolicy p){return p==ReturnPolicy::Delete?"delete":p==ReturnPolicy::Pass?"pass":"invalid";}
inline const char* Name(Normalization n){return n==Normalization::NetUpstreamEscape?"net-upstream-escape":n==Normalization::GrossEmission?"gross-emission":"invalid";}
inline Transport::Status Validate(const Configuration& c) {
  using namespace Transport;
  if((c.returnPolicy!=ReturnPolicy::Delete&&c.returnPolicy!=ReturnPolicy::Pass)||
      (c.normalization!=Normalization::NetUpstreamEscape&&c.normalization!=Normalization::GrossEmission))
    return Status::Error(StatusCode::InvalidArgument,"unknown reduced particle source/return policy");
  if(c.returnPolicy==ReturnPolicy::Pass)return Status::Error(StatusCode::UnsupportedConfiguration,
      "Shock return policy 'pass' is not implemented yet for reduced SSE/Parker transport: no downstream transport volume is available. Select 'delete'.");
  if(c.normalization==Normalization::GrossEmission)return Status::Error(StatusCode::UnsupportedConfiguration,
      "gross-emission normalization is not implemented for the finite-reference source: an accelerator-to-reference escape law is required; select net-upstream-escape");
  if(!std::isfinite(c.releaseDistanceShockRadiusFraction)||c.releaseDistanceShockRadiusFraction<0||
      c.releaseDistanceShockRadiusFraction>10)
    return Status::Error(StatusCode::InvalidArgument,"release distance shock-radius fraction must lie in (0,10]");
  if(c.releaseDistanceM>0&&c.releaseDistanceShockRadiusFraction>0)
    return Status::Error(StatusCode::InvalidArgument,"give either a fixed release distance or a shock-radius fraction, not both");
  if(!std::isfinite(c.releaseDistanceM)||c.releaseDistanceM<0||!std::isfinite(c.boundaryToleranceM)||
      c.boundaryToleranceM<=0||!std::isfinite(c.verificationMonoenergeticEnergyJ)||c.verificationMonoenergeticEnergyJ<0)
    return Status::Error(StatusCode::InvalidArgument,"invalid reduced source release distance or numerical controls");
  if(!c.enabled)return Status::Ok();
  if(c.logUniformMomentumSampling&&c.verificationMonoenergeticEnergyJ>0)
    return Status::Error(StatusCode::InvalidArgument,
        "log-uniform momentum sampling applies to the production DSA law, not to a monoenergetic verification source");
  if(!(c.releaseDistanceM>0||c.releaseDistanceShockRadiusFraction>0)||!(c.minimumEnergyJ>0)||c.maximumEnergyJ<c.minimumEnergyJ||
      !std::isfinite(c.minimumEnergyJ)||!std::isfinite(c.maximumEnergyJ)||
      !std::isfinite(c.releaseFraction)||c.releaseFraction<0||c.releaseFraction>1||c.samplesPerStep==0)
    return Status::Error(StatusCode::InvalidArgument,
        "particle-enabled reduced transport requires an explicit positive L_ref [m] or L_ref/R_sh fraction, valid source energy bounds [J], release fraction in [0,1], and a positive sample count");
  return Status::Ok();
}
inline std::string Canonical(const Configuration& c) {
  std::ostringstream out;out<<std::setprecision(17)<<c.enabled<<' '<<static_cast<int>(c.returnPolicy)<<' '
      <<static_cast<int>(c.normalization)<<' '<<c.releaseDistanceM<<' '<<c.minimumEnergyJ<<' '
      <<c.maximumEnergyJ<<' '<<c.releaseFraction<<' '<<c.samplesPerStep<<' '<<c.campaignSeed<<' '
      <<c.boundaryToleranceM<<' '<<c.verificationMonoenergeticEnergyJ;
  if(c.logUniformMomentumSampling)out<<" log-uniform-momentum";
  if(c.releaseDistanceShockRadiusFraction>0)out<<" release-shock-radius-fraction "<<c.releaseDistanceShockRadiusFraction;
  if(c.boundaryRootTable)out<<" boundary-root-table";
  return out.str();
}
// L_ref [m] at a birth epoch whose front generating sphere has radius R [m].
inline double ReleaseDistanceM(const Configuration& c,double generatingRadiusM) {
  return c.releaseDistanceShockRadiusFraction>0?c.releaseDistanceShockRadiusFraction*generatingRadiusM:
      c.releaseDistanceM;
}
} }
#endif
