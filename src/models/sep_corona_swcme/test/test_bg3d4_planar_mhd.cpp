#include "sep_corona_swcme/piston_solver.h"
#include "sep_coronal_cme/constants.h"
#include "sep_coronal_cme/mhd_jump_solver.h"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using SEP::CoronaSwcme::PlanarPistonInput;
using SEP::CoronaSwcme::PlanarPistonSolver;
using SEP::CoronalCME::KinematicValue;
using SEP::Core::Result;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

double Relative(double value,double reference) {
  return std::abs(value-reference)/std::max(std::abs(reference),1e-300);
}

struct Row {
  double position=0,speed=0,compression=0,flux=0,energy=0;
};

Row Run(int cells) {
  constexpr double pistonSpeed=0.5,end=0.9;
  PlanarPistonInput in;
  in.startS=0;in.endS=end;in.leftPositionM=0;in.columnLengthM=3;
  in.areaM2=1;in.initialDensityKgM3=1;in.initialPressurePa=1e-8;
  in.initialVelocityMPerS=0;in.initialTransverseMagneticFieldT=1;
  // gamma=2 makes deposited shock heat and transverse magnetic pressure part
  // of one barotropic gamma=2 analogue.  A gamma=5/3 gas heated by Q would be
  // a different finite-beta MHD closure and must use the canonical RH oracle.
  in.magneticPermeabilityNPerA2=1;in.gammaAdiabatic=2;
  in.cells=cells;in.cfl=0.25;in.quadraticViscosity=2;
  in.linearViscosity=0;in.linearViscosityActivation=0;
  in.shockThreshold=0.02;
  auto solver=PlanarPistonSolver::Create(in,[](double time) {
    return Result<KinematicValue>::Success(
        {pistonSpeed*time,pistonSpeed,0});
  });
  Require(solver.ok(),solver.status.message);
  const auto advanced=solver.value->AdvanceTo(end);
  Require(advanced.ok(),advanced.message);
  const auto shock=solver.value->DetectShock();
  Require(shock.ok()&&shock.value.present&&shock.value.statesAvailable,
      "perpendicular piston lacks a detected shock/two-sided state");
  const auto energy=solver.value->EnergyLedger();
  Require(energy.ok(),energy.status.message);

  // In the p->0 perpendicular limit, transverse magnetic energy has an
  // effective gamma=2 because B_t/rho is frozen.  This independent reduction
  // gives the exact piston-shock speed without calling the production solver.
  constexpr double alfvenSpeed=1;
  const double exactSpeed=0.75*pistonSpeed+
      std::sqrt(std::pow(0.75*pistonSpeed,2)+alfvenSpeed*alfvenSpeed);
  const double exactCompression=exactSpeed/(exactSpeed-pistonSpeed);
  const double energyScale=std::abs(energy.value.currentEnergyJ)+
      std::abs(energy.value.initialEnergyJ)+std::abs(energy.value.pistonWorkJ)+
      std::abs(energy.value.outerBoundaryWorkJ);
  return {Relative(shock.value.radiusM,exactSpeed*end),
      Relative(shock.value.speedMPerS,exactSpeed),
      Relative(shock.value.compressionRatio,exactCompression),
      Relative(shock.value.downstreamTransverseMagneticFieldT/
          shock.value.upstreamTransverseMagneticFieldT,
          shock.value.compressionRatio),
      std::abs(energy.value.residualJ)/energyScale};
}

void TestColdPerpendicularPiston() {
  const std::array<int,3> cells={400,800,1600};
  std::array<Row,3> row;
  for(std::size_t i=0;i<row.size();++i)row[i]=Run(cells[i]);
  std::cout<<"[CSWC0623-EVIDENCE] position_errors="<<row[0].position<<','
           <<row[1].position<<','<<row[2].position<<" speed_errors="
           <<row[0].speed<<','<<row[1].speed<<','<<row[2].speed
           <<" compression_errors="<<row[0].compression<<','
           <<row[1].compression<<','<<row[2].compression
           <<" frozen_flux_errors="<<row[0].flux<<','<<row[1].flux<<','
           <<row[2].flux<<" energy_errors="<<row[0].energy<<','
           <<row[1].energy<<','<<row[2].energy<<'\n';
  Require(row[2].position<0.02&&row[2].speed<0.02&&
      row[2].compression<0.02&&row[2].flux<2e-12&&row[2].energy<2e-4,
      "cold perpendicular MHD piston exceeds preregistered accuracy");
  Require(row[2].position<row[1].position&&row[1].position<row[0].position,
      "perpendicular MHD shock position does not converge");
}

struct FiniteBetaReference {
  double shockSpeed=0,compression=0,downstreamPressure=0;
};

FiniteBetaReference CanonicalFiniteBeta(double pistonSpeed) {
  const double mu0=SEP::CoronalCME::Constants::kVacuumPermeabilityHPerM;
  SEP::CoronalCME::MhdPrimitiveState upstream;
  upstream.massDensityKgM3=1;upstream.pressurePa=0.1;
  upstream.velocityMPerS={0,0,0};
  upstream.magneticFieldT={0,std::sqrt(mu0),0};
  const SEP::CoronalCME::Vec3 normal={1,0,0};
  const auto characteristics=SEP::CoronalCME::EvaluateMhdCharacteristics(
      upstream,normal,5.0/3.0);
  Require(characteristics.ok(),characteristics.status.message);
  double low=1.001*characteristics.value.fastSpeedMPerS,high=4;
  auto velocity=[&](double speed) {
    const auto jump=SEP::CoronalCME::SolveObliqueFastShock(
        upstream,normal,speed,5.0/3.0,1e-10);
    Require(jump.ok(),jump.status.message);
    return jump;
  };
  Require(velocity(low).value.downstream.velocityMPerS.x<pistonSpeed&&
      velocity(high).value.downstream.velocityMPerS.x>pistonSpeed,
      "canonical finite-beta piston root is not bracketed");
  for(int i=0;i<90;++i) {
    const double middle=0.5*(low+high);
    if(velocity(middle).value.downstream.velocityMPerS.x<pistonSpeed)
      low=middle;else high=middle;
  }
  const auto exact=velocity(0.5*(low+high));
  return {0.5*(low+high),exact.value.compressionRatio,
      exact.value.downstream.pressurePa};
}

Row RunFiniteBeta(int cells,const FiniteBetaReference& exact) {
  constexpr double pistonSpeed=0.5,end=0.9;
  const double mu0=SEP::CoronalCME::Constants::kVacuumPermeabilityHPerM;
  PlanarPistonInput in;
  in.startS=0;in.endS=end;in.leftPositionM=0;in.columnLengthM=3;
  in.areaM2=1;in.initialDensityKgM3=1;in.initialPressurePa=0.1;
  in.initialVelocityMPerS=0;in.initialTransverseMagneticFieldT=std::sqrt(mu0);
  in.magneticPermeabilityNPerA2=mu0;in.gammaAdiabatic=5.0/3.0;
  in.cells=cells;in.cfl=0.25;in.quadraticViscosity=2;
  in.linearViscosity=0;in.shockThreshold=0.02;
  auto solver=PlanarPistonSolver::Create(in,[](double time) {
    return Result<KinematicValue>::Success(
        {pistonSpeed*time,pistonSpeed,0});
  });
  Require(solver.ok(),solver.status.message);
  Require(solver.value->AdvanceTo(end).ok(),
      "finite-beta numerical piston advance failed");
  const auto shock=solver.value->DetectShock();
  Require(shock.ok()&&shock.value.present&&shock.value.statesAvailable,
      "finite-beta piston lacks two-sided shock states");
  const auto energy=solver.value->EnergyLedger();
  Require(energy.ok(),energy.status.message);
  const double scale=std::abs(energy.value.currentEnergyJ)+
      std::abs(energy.value.initialEnergyJ)+std::abs(energy.value.pistonWorkJ)+
      std::abs(energy.value.outerBoundaryWorkJ);
  return {Relative(shock.value.radiusM,exact.shockSpeed*end),
      Relative(shock.value.speedMPerS,exact.shockSpeed),
      Relative(shock.value.compressionRatio,exact.compression),
      Relative(shock.value.downstreamPressurePa,exact.downstreamPressure),
      std::abs(energy.value.residualJ)/scale};
}

void TestFiniteBetaCanonicalRh() {
  constexpr double pistonSpeed=0.5;
  const auto exact=CanonicalFiniteBeta(pistonSpeed);
  const std::array<int,3> cells={400,800,1600};
  std::array<Row,3> row;
  for(std::size_t i=0;i<row.size();++i)
    row[i]=RunFiniteBeta(cells[i],exact);
  std::cout<<"[CSWC0623-FINITE-BETA] exact_speed="<<exact.shockSpeed
           <<" exact_compression="<<exact.compression
           <<" exact_pressure="<<exact.downstreamPressure
           <<" position_errors="<<row[0].position<<','<<row[1].position<<','
           <<row[2].position<<" speed_errors="<<row[0].speed<<','
           <<row[1].speed<<','<<row[2].speed<<" compression_errors="
           <<row[0].compression<<','<<row[1].compression<<','
           <<row[2].compression<<" pressure_errors="<<row[0].flux<<','
           <<row[1].flux<<','<<row[2].flux<<" energy_errors="
           <<row[0].energy<<','<<row[1].energy<<','<<row[2].energy<<'\n';
  Require(row[2].position<0.02&&row[2].speed<0.02&&
      row[2].compression<0.02&&row[2].flux<0.03&&row[2].energy<2e-4,
      "finite-beta piston disagrees with the canonical RH oracle");
  Require(row[2].position<row[1].position&&row[1].position<row[0].position,
      "finite-beta shock position does not converge");
}

} // namespace

int main() {
  try {
    TestColdPerpendicularPiston();
    TestFiniteBetaCanonicalRh();
    std::cout<<"[CSWC0623] PASS cold perpendicular-MHD piston and flux\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-4 planar MHD FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
