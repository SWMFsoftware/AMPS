#ifndef SEP_CORONA_SWCME_AMBIENT_STATE_H
#define SEP_CORONA_SWCME_AMBIENT_STATE_H

#include "sep_corona_swcme/cme_event.h"
#include "sep_coronal_cme/pfss_harmonics.h"
#include "sep_coronal_cme/plasma_eos.h"

#include <array>
#include <memory>
#include <vector>

namespace SEP { namespace CoronaSwcme {

enum class AmbientRegion { PfssOpen, PfssClosed, ParkerExterior };
enum class DifferenceStencil { CentralSecondOrder, ForwardSecondOrder,
                               BackwardSecondOrder };

struct AmbientPrimitive {
  // SI state in the event's inertial HCI Cartesian frame. Pressure is total
  // thermal species pressure; temperatures retain the configured proton/
  // electron partition instead of being reconstructed from an implicit mu.
  CoronalCME::Vec3 magneticFieldT;
  CoronalCME::Vec3 velocityMPerS;
  CoronalCME::PlasmaState plasma;
  double protonTemperatureK = 0.0;
  double electronTemperatureK = 0.0;
  double electronPressurePa = 0.0;
  AmbientRegion region = AmbientRegion::PfssClosed;
  int magneticSector = 0; // zero only for closed PFSS loops
};

struct AmbientState {
  AmbientPrimitive primitive;
  // Row-major component/coordinate tensors in T/m and 1/s. Each derivative
  // stencil stays on one topology/sector branch; interface-crossing finite
  // differences fail instead of smoothing a physical discontinuity.
  std::array<double,9> gradientB{};
  std::array<double,9> gradientU{};
  std::array<DifferenceStencil,3> stencils{};
  CoronalCME::Vec3 gradientLogBPerM;
  double epochS = 0.0;
  std::uint64_t generation = 0;
  std::string eventIdentity;
  bool derivativesValid = false;
};

// Minimal immutable construction record for the canonical ambient authority.
// The full-CME EventConfiguration remains a supported caller below, but a
// reduced front must not manufacture inactive contact/sheath/ejecta inputs
// merely to obtain rho, p, U and B.  All quantities are SI and the identity
// binds the composition, PFSS/Parker inputs and radial/time support supplied by
// the reduced event resolver.
struct AmbientDefinition {
  Composition composition;
  AmbientInput ambient;
  EventSupport support;
  std::string coordinateFrame;
  std::string physicsFingerprint;
};

// One immutable ambient authority from the physical Sun through declared
// >1-AU coverage. The low corona uses maintained PFSS/topology/closed-plasma
// kernels; open tubes and the exterior share one isothermal mass-flux/Parker
// continuation. There are no field/density floors or fallback providers.
class AmbientModel final {
 public:
  static Core::Result<std::shared_ptr<const AmbientModel>> Create(
      std::shared_ptr<const EventConfiguration> event);
  static Core::Result<std::shared_ptr<const AmbientModel>> Create(
      const AmbientDefinition& definition);

  Core::Result<AmbientPrimitive> Evaluate(
      CoronalCME::Vec3 positionM,double epochS) const;
  Core::Result<AmbientState> EvaluateWithDerivatives(
      CoronalCME::Vec3 positionM,double epochS,std::uint64_t generation) const;

  const EventConfiguration& Event() const { return *event_; }
  const AmbientDefinition& Definition() const { return definition_; }
  const std::string& Identity() const { return definition_.physicsFingerprint; }
  double MassFluxPerSteradianKgPerS() const { return massFluxPerSr_; }

 private:
  static Core::Result<std::shared_ptr<const AmbientModel>> CreateResolved(
      const AmbientDefinition& definition,
      std::shared_ptr<const EventConfiguration> event);
  std::shared_ptr<const EventConfiguration> event_;
  AmbientDefinition definition_;
  CoronalCME::PfssHarmonics pfss_;
  std::vector<CoronalCME::IonSpecies> ions_;
  std::vector<double> logRadius_,logSpeed_,slope_,winding_;
  double massFluxPerSr_=0.0,soundSquared_=0.0,criticalRadiusM_=0.0;

  double WindSpeed(double radiusM) const;
  double IntegrateWinding(double beginM,double endM) const;
  double Winding(double radiusM) const;
  Core::Result<std::pair<CoronalCME::FieldLineTopology,CoronalCME::Vec3>> Trace(
      CoronalCME::Vec3 positionM) const;
};

const char* Name(AmbientRegion region) noexcept;
bool ValidAmbientRegionSector(AmbientRegion region,int sector) noexcept;

} } // namespace SEP::CoronaSwcme

#endif
