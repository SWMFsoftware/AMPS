#ifndef SEP_CORONA_SWCME_CROSS_RAY_FLUID_H
#define SEP_CORONA_SWCME_CROSS_RAY_FLUID_H

#include "sep_corona_swcme/cross_ray_reconstruction.h"
#include "sep_corona_swcme/cross_ray_transport.h"

namespace SEP { namespace CoronaSwcme { namespace CrossRay {

struct FluidDensity {
  double massKgM3=0;
  CoronalCME::Vec3 momentumKgM2S;
  // Total HCI energy includes internal, kinetic and magnetic energy. B has
  // exactly one authority (the shared CT face flux), not a second cell array.
  double totalEnergyJPerM3=0;
};

struct IntegratedFluidFlux {
  // Time-averaged integral of normal flux over the ACTUAL curved face, in
  // the mesh's canonical left-to-right orientation: kg/s, newtons and watts.
  // A point flux times |integral n dA| is not generally this integral.
  double massKgPerS=0;
  CoronalCME::Vec3 momentumN;
  double totalEnergyW=0;
};

struct IntegratedVolumeSource {
  // Explicit volume integrals of prescribed HCI force and total power. No
  // mass source is supported. Force work, heating and electromagnetic source
  // power must be included by the caller's independently qualified closure;
  // this kernel neither invents missing work nor resets thermal energy.
  CoronalCME::Vec3 forceN;
  double totalPowerW=0;
};

struct FluidFluxInterval {
  EdgeEmfInterval magnetic;
  std::vector<IntegratedFluidFlux> faces;
  std::vector<IntegratedVolumeSource> sources;
};

struct ConservedInventory {
  double massKg=0;
  CoronalCME::Vec3 momentumKgMPerS;
  double totalEnergyJ=0;
};

struct FluidStepLedger {
  ConservedInventory initial,final,outwardBoundaryTransfer,volumeAddition;
  // final-initial+outwardBoundaryTransfer-volumeAddition, with the raw
  // floating-point residual retained in SI, never used to correct a cell.
  ConservedInventory residual;
};

// Atomic finite-volume/CT epoch coupling on a STATIC mesh. This is the
// conservative reference-update kernel, not a production MHD time integrator.
// The user selected analytical material transport without an MHD runtime; these
// operators independently check conservation/energy, not select that runtime.
// Qualified
// face/time quadrature, CFL selection, ideal upwind edge EMF, physical domain
// boundaries, regional source laws and moving piston geometry remain caller
// responsibilities and are not inferred by this API.
//
// Every shared face is evaluated once by the caller and applied with opposite
// signs here. The new CT field is reconstructed only after every fluid flux
// and source has been accumulated. Pressure is recovered from the conserved
// total energy using that new field; nonpositive density/thermal energy
// rejects the WHOLE epoch. No face-by-face partial state can be published.
class StaticFluidState final {
 public:
  static Core::Result<std::shared_ptr<const StaticFluidState>> Create(
      std::shared_ptr<const StaticMagneticFlux> magnetic,
      std::vector<FluidDensity> cells,Thermodynamics gas);
  Core::Result<std::shared_ptr<const StaticFluidState>> Advance(
      const FluidFluxInterval& interval) const;

  const std::shared_ptr<const StaticMagneticFlux>& Magnetic() const noexcept {return magnetic_;}
  const std::vector<FluidDensity>& Cells() const noexcept {return cells_;}
  const std::vector<Primitive>& Primitives() const noexcept {return primitive_;}
  const Thermodynamics& Gas() const noexcept {return gas_;}
  const FluidStepLedger& StepLedger() const noexcept {return ledger_;}
  std::string Fingerprint() const;

 private:
  // Only a validated factory may create an authority; an empty public default
  // object would bypass mesh/epoch/positivity checks before a later query.
  StaticFluidState()=default;
  std::shared_ptr<const StaticMagneticFlux> magnetic_;
  std::vector<FluidDensity> cells_;
  std::vector<Primitive> primitive_;
  Thermodynamics gas_;
  FluidStepLedger ledger_;
};

} } } // namespace SEP::CoronaSwcme::CrossRay

#endif
