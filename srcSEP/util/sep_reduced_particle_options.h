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
  if(!std::isfinite(c.releaseDistanceM)||c.releaseDistanceM<0||!std::isfinite(c.boundaryToleranceM)||
      c.boundaryToleranceM<=0||!std::isfinite(c.verificationMonoenergeticEnergyJ)||c.verificationMonoenergeticEnergyJ<0)
    return Status::Error(StatusCode::InvalidArgument,"invalid reduced source release distance or numerical controls");
  if(!c.enabled)return Status::Ok();
  if(!(c.releaseDistanceM>0)||!(c.minimumEnergyJ>0)||c.maximumEnergyJ<c.minimumEnergyJ||
      !std::isfinite(c.minimumEnergyJ)||!std::isfinite(c.maximumEnergyJ)||
      !std::isfinite(c.releaseFraction)||c.releaseFraction<0||c.releaseFraction>1||c.samplesPerStep==0)
    return Status::Error(StatusCode::InvalidArgument,
        "particle-enabled reduced transport requires explicit positive L_ref [m], valid source energy bounds [J], release fraction in [0,1], and a positive sample count");
  return Status::Ok();
}
inline std::string Canonical(const Configuration& c) {
  std::ostringstream out;out<<std::setprecision(17)<<c.enabled<<' '<<static_cast<int>(c.returnPolicy)<<' '
      <<static_cast<int>(c.normalization)<<' '<<c.releaseDistanceM<<' '<<c.minimumEnergyJ<<' '
      <<c.maximumEnergyJ<<' '<<c.releaseFraction<<' '<<c.samplesPerStep<<' '<<c.campaignSeed<<' '
      <<c.boundaryToleranceM<<' '<<c.verificationMonoenergeticEnergyJ;return out.str();
}
} }
#endif
