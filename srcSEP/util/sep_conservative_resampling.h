#ifndef SRCSEP_UTIL_SEP_CONSERVATIVE_RESAMPLING_H
#define SRCSEP_UTIL_SEP_CONSERVATIVE_RESAMPLING_H
#include "sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include "pic_conservative_resampling.h"
#include <cstdint>
#include <vector>

namespace SEP { namespace PopulationControl { namespace Conservative {
struct Particle {
  std::uint64_t identity=0;
  int species=-1;
  double weight=0.0,localFraction=0.0,parallelMPerS=0.0,normalMPerS=0.0;
  double massKg=0.0,chargeC=0.0;
  // A reduced shock source retains its per-proton kinetic energy at birth
  // [J] even after cooling/heating changes the current momentum. Weight-only
  // resampling must preserve the represented birth-energy measure as well
  // as current energy; otherwise deleting computational representatives would
  // appear as an unrecorded physical source/loss in the return ledger. All
  // candidates must either carry this history or explicitly lack it. Zero is
  // not a substitute for missing source history in an accelerated population.
  bool hasBirthEnergy=false;
  double birthKineticEnergyJ=0.0;
};
struct UnpreservedMoments {
  double parallelCurrentCMPerS=0.0;
  double parallelStressJ=0.0;
  double perpendicularStressJ=0.0;
};
struct Endpoints {
  PIC::ConservativeResampling::Endpoints weights;
  UnpreservedMoments before,positive,negative;
};
Transport::Status EvaluateConstraints(const std::vector<Particle>& particles,
    double speedOfLightMPerS,std::vector<std::vector<double>>* rows,
    std::vector<double>* weights);
// SI-scaled final transaction check, including every deleted input at weight0.
double MaximumRelativeResidual(const std::vector<std::vector<double>>& rows,
    const std::vector<double>& before,const std::vector<double>& after);
UnpreservedMoments EvaluateUnpreserved(const std::vector<Particle>& particles,
    const std::vector<double>& weights,double speedOfLightMPerS=299792458.0);
// Same species/charge/mass in one segment; retain all phase points/identities.
// Exact declared rows: number, relativistic kinetic energy [J], parallel
// momentum [kg m/s], mu, mu^2, local segment fraction. Charge follows number.
// Reduced-source populations additionally conserve birth kinetic energy [J].
// Current and stress moments are diagnostics, not secretly conserved by mu^2.
Endpoints Prepare(const std::vector<Particle>& particles,double tolerance,
                  double speedOfLightMPerS=299792458.0);
Transport::Status Select(const Endpoints& endpoints,double openUniform,
                         std::vector<double>* weights,
                         UnpreservedMoments* unpreserved=nullptr);
}}} // namespace SEP::PopulationControl::Conservative
#endif
