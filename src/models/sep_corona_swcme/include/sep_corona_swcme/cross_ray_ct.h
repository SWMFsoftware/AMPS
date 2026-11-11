#ifndef SEP_CORONA_SWCME_CROSS_RAY_CT_H
#define SEP_CORONA_SWCME_CROSS_RAY_CT_H

#include "sep_corona_swcme/cross_ray_mesh.h"

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

struct EdgeEmfInterval {
  // Values are time-averaged oriented line integrals integral E.dl, volts,
  // on the mesh's canonical curved edges over [startS,endS]. A production
  // integrator must supply its actual qualified temporal quadrature; this
  // kernel does not invent a midpoint/hold or reconstruct E from particle data.
  std::string meshFingerprint,eventFingerprint;
  std::uint64_t fieldGeneration=0;
  double startS=0,endS=0;
  std::vector<double> valuesVolt;
};

// Immutable STATIC-mesh constrained-transport state. One face flux is shared
// by both cells and one edge circulation by all incident faces. The exact
// integer boundary-of-boundary identity makes the discrete divergence update
// zero in real arithmetic. Raw floating-point divergence and explicit
// arithmetic uncertainty are retained separately; neither is cleaned/floored.
// This does not yet couple magnetic energy to the fluid, reconstruct cell B,
// supply ideal-MHD edge E, or handle moving/remapped geometry. Those remain
// required before regional or native qualification.
class StaticMagneticFlux final {
 public:
  static Core::Result<std::shared_ptr<const StaticMagneticFlux>> FromFaceFlux(
      std::shared_ptr<const StaticSphericalMesh> mesh,
      std::vector<double> faceFluxTm2);
  // The edge data are integral A.dl in T m^2, with B=curl A. This is an
  // initial-condition interface, not permission to silently project an
  // inconsistent measured/background B onto a different field.
  static Core::Result<std::shared_ptr<const StaticMagneticFlux>> FromEdgePotential(
      std::shared_ptr<const StaticSphericalMesh> mesh,
      const std::vector<double>& edgePotentialTm2);
  Core::Result<std::shared_ptr<const StaticMagneticFlux>> Advance(
      const EdgeEmfInterval& interval) const;

  const std::shared_ptr<const StaticSphericalMesh>& Mesh() const noexcept {return mesh_;}
  const std::vector<double>& FaceFluxTm2() const noexcept {return flux_;}
  const std::vector<double>& NumericalFluxBoundsTm2() const noexcept {return uncertainty_;}
  const std::vector<double>& RawDivergenceTPerM() const noexcept {return divergence_;}
  double EpochS() const noexcept {return epochS_;}
  std::uint64_t Generation() const noexcept {return generation_;}
  std::string Fingerprint() const;

 private:
  // Only a validated factory may create an authority; an empty public default
  // object would bypass mesh/epoch/positivity checks before a later query.
  StaticMagneticFlux()=default;
  static Core::Result<std::shared_ptr<const StaticMagneticFlux>> Publish(
      std::shared_ptr<const StaticSphericalMesh> mesh,
      std::vector<double> flux,std::vector<double> uncertainty,
      double epoch,std::uint64_t generation);
  std::shared_ptr<const StaticSphericalMesh> mesh_;
  std::vector<double> flux_,uncertainty_,divergence_;
  double epochS_=0;
  std::uint64_t generation_=0;
};

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif
