#ifndef _HARTLEY2_DAUGHTER_PACKET_PLANNER_H_
#define _HARTLEY2_DAUGHTER_PACKET_PLANNER_H_

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

// Pure Stage-6 transaction planner.  It deliberately knows nothing about
// AMPS particles, MPI, meshes, reaction kinematics, or RNG state.  Production
// chemistry can therefore prove the complete K-products-by-N-packets
// allocation and every per-product correction sum before it changes a parent
// particle or diagnostics.
namespace Comet {
namespace DaughterPacketPlanner {

struct cPlan {
  std::size_t ConcreteProductCount=0;
  std::size_t PacketCount=1;
  std::size_t TotalDaughterCount=0;
  std::size_t RequiredAdditionalParticles=0;

  // Packet-major storage: Correction(packet,product).  Packet-major ordering
  // matches the production commit, where one independent complete recoil
  // solution is emitted before the next packet is sampled.
  std::vector<double> Corrections;

  double Correction(std::size_t packet,std::size_t product) const {
    if (packet>=PacketCount || product>=ConcreteProductCount)
      throw std::out_of_range("daughter-packet correction index is invalid");
    return Corrections[packet*ConcreteProductCount+product];
  }
};

inline cPlan BuildPlan(const std::vector<double> &unsplitCorrections,
    std::size_t packetCount,std::size_t maximumPacketCount,
    std::size_t maximumDaughterCount,
    std::size_t availableAdditionalParticles=
        std::numeric_limits<std::size_t>::max()) {
  cPlan result;
  result.ConcreteProductCount=unsplitCorrections.size();
  result.PacketCount=packetCount;

  if (packetCount<1)
    throw std::invalid_argument("daughter packet count must be positive");
  if (maximumPacketCount<1 || packetCount>maximumPacketCount)
    throw std::length_error("daughter packet count exceeds its hard maximum");
  if (result.ConcreteProductCount==0 && packetCount!=1)
    throw std::invalid_argument(
        "a product-free sink cannot be split into daughter packets");
  if (result.ConcreteProductCount!=0 &&
      packetCount>std::numeric_limits<std::size_t>::max()/
          result.ConcreteProductCount)
    throw std::length_error("daughter packet transaction count overflows");

  result.TotalDaughterCount=result.ConcreteProductCount*packetCount;
  if (result.TotalDaughterCount>maximumDaughterCount)
    throw std::length_error(
        "daughter packet transaction exceeds its hard product limit");
  result.RequiredAdditionalParticles=result.TotalDaughterCount==0 ? 0 :
      result.TotalDaughterCount-1;
  if (result.RequiredAdditionalParticles>availableAdditionalParticles)
    throw std::length_error(
        "insufficient particle capacity for daughter packet transaction");

  result.Corrections.assign(result.TotalDaughterCount,0.0);
  for (std::size_t product=0;product<result.ConcreteProductCount;product++) {
    const double correction=unsplitCorrections[product];
    if (!(correction>0.0) || !std::isfinite(correction))
      throw std::invalid_argument(
          "unsplit daughter correction must be positive and finite");

    const double equal=correction/static_cast<double>(packetCount);
    if (!(equal>0.0) || !std::isfinite(equal))
      throw std::runtime_error("daughter packet correction underflowed");
    double assigned=0.0;
    for (std::size_t packet=0;packet+1<packetCount;packet++) {
      result.Corrections[packet*result.ConcreteProductCount+product]=equal;
      assigned+=equal;
    }
    const double finalCorrection=correction-assigned;
    if (!(finalCorrection>0.0) || !std::isfinite(finalCorrection) ||
        assigned+finalCorrection!=correction)
      throw std::runtime_error(
          "daughter packet corrections do not exactly reconstruct the input");
    result.Corrections[(packetCount-1)*result.ConcreteProductCount+product]=
        finalCorrection;
  }
  return result;
}

} // namespace DaughterPacketPlanner
} // namespace Comet

#endif
