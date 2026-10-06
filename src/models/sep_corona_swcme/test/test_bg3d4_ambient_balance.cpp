#include "sep_corona_swcme/piston_ambient.h"
#include "sep_corona_swcme/piston_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronaSwcme::AmbientRegion;
using SEP::CoronaSwcme::PistonAmbientProjection;
using SEP::CoronaSwcme::PistonAmbientSourceTable;
using SEP::CoronaSwcme::PistonAmbientTrajectory;
using SEP::CoronaSwcme::PistonInitialState;
using SEP::CoronaSwcme::PistonTubeGeometry;
using SEP::CoronaSwcme::PlanarPistonInput;
using SEP::CoronaSwcme::PlanarPistonSolver;
using SEP::CoronalCME::KinematicValue;
using SEP::CoronalCME::Vec3;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read fixture: "+path);
  std::ostringstream bytes;bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"cannot finish reading fixture: "+path);
  return bytes.str();
}

SEP::Core::Result<std::string> ReadFile(const std::string& path) {
  try {
    return SEP::Core::Result<std::string>::Success(FileBytes(path));
  } catch(const std::exception& error) {
    return SEP::Core::Result<std::string>::Failure(
        SEP::Core::StatusCode::DataIntegrityFailure,error.what());
  }
}

std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration> Event() {
  const auto event=SEP::CoronaSwcme::ResolveEventConfiguration(
      FileBytes("examples/bg3d4_piston/event.conf"),ReadFile);
  Require(event.ok(),"piston event rejected: "+event.status.message);
  return event.value;
}

double ShellCentroid(double left,double right) {
  return 0.75*(std::pow(right,4)-std::pow(left,4))/
      (std::pow(right,3)-std::pow(left,3));
}

struct Errors {
  double density=0,pressure=0,velocity=0,magnetic=0,energy=0;
  double maxForceRatio=0,maxHeating=0;
};

struct ParkerErrors {
  double density=0,pressure=0,velocity=0,magnetic=0,energy=0;
  double radialFlux=0,maximumHeating=0,trajectoryVelocity=0;
};

double TestSourceLedger() {
  PlanarPistonInput in;
  in.geometry=PistonTubeGeometry::Planar;in.startS=0;in.endS=2;
  in.leftPositionM=0;in.columnLengthM=1;in.areaM2=1;
  in.initialDensityKgM3=1;in.initialPressurePa=1;
  in.initialVelocityMPerS=0;in.initialTransverseMagneticFieldT=0.2;
  in.initialTransverseMagneticField2T=-0.1;
  in.magneticPermeabilityNPerA2=1;in.gammaAdiabatic=5.0/3.0;
  in.cells=64;in.cfl=0.2;in.quadraticViscosity=2;
  in.shockThreshold=0.01;
  const auto fixed=[](double) {
    return SEP::Core::Result<KinematicValue>::Success({0,0,0});
  };
  const auto source=[](double,double) {
    SEP::CoronaSwcme::PistonVolumeSource value;
    value.heatingWPerM3=0.01;
    value.transverseInvariantRate={0.02,-0.01};
    return SEP::Core::Result<SEP::CoronaSwcme::PistonVolumeSource>::Success(value);
  };
  auto solver=PlanarPistonSolver::Create(in,fixed,source);
  Require(solver.ok(),"source-ledger solver construction failed");
  const auto advanced=solver.value->AdvanceTo(2);
  Require(advanced.ok(),"source-ledger advance failed: "+advanced.message);
  const auto cells=solver.value->Cells();
  const auto ledger=solver.value->EnergyLedger();
  Require(cells.ok()&&ledger.ok(),"source-ledger state is unavailable");
  const double expectedE=1/(in.gammaAdiabatic-1)+0.02;
  for(const auto& cell:cells.value)Require(
      std::abs(cell.specificInternalEnergyJPerKg-expectedE)<2e-13&&
      std::abs(cell.transverseMagneticFieldT-0.24)<2e-13&&
      std::abs(cell.transverseMagneticField2T+0.12)<2e-13&&
      std::abs(cell.velocityMPerS)<2e-13,
      "uniform heating/induction source does not recover its exact solution");
  Require(std::abs(ledger.value.volumeHeatingJ-0.02)<2e-13&&
      ledger.value.magneticSourceWorkJ>0,
      "heating or magnetic-source work is absent from the ledger");
  return std::abs(ledger.value.residualJ)/ledger.value.initialEnergyJ;
}

Errors Run(int cells) {
  const auto event=Event();
  const auto ambient=AmbientModel::Create(event);
  Require(ambient.ok(),"ambient construction failed: "+ambient.status.message);
  // The equatorial dipole ray is a closed PFSS branch over this fixture.  Its
  // radial flow is exactly zero while the transverse field, gas pressure and
  // density are strongly nonuniform, so it tests magnetic hoop and thermal/
  // magnetic pressure balance without hiding errors behind uniform data.
  const auto projection=PistonAmbientProjection::Create(
      ambient.value,Vec3{1,0,0},event->support.startS,2e-5);
  Require(projection.ok(),"ambient projection failed: "+projection.status.message);
  const double solar=event->support.solarRadiusM;
  const double inner=1.15*solar,length=0.2*solar,end=200;
  const auto innerState=projection.value->Evaluate(inner);
  Require(innerState.ok()&&innerState.value.region==AmbientRegion::PfssClosed,
      "well-balance ray is not on the declared closed branch");
  const auto sourceTable=PistonAmbientSourceTable::Create(projection.value,
      inner,inner+length,2*cells+1);
  Require(sourceTable.ok(),"ambient source table failed: "+
      sourceTable.status.message);

  PlanarPistonInput in;
  in.geometry=PistonTubeGeometry::RadialSpherical;
  in.startS=0;in.endS=end;in.leftPositionM=inner;in.columnLengthM=length;
  in.solidAngleSr=0.1;in.initialDensityKgM3=innerState.value.densityKgM3;
  in.initialPressurePa=innerState.value.pressurePa;
  in.initialVelocityMPerS=0;
  in.initialTransverseMagneticFieldT=
      innerState.value.transverseMagneticFieldT[0];
  in.initialTransverseMagneticField2T=
      innerState.value.transverseMagneticFieldT[1];
  in.gammaAdiabatic=event->composition.gammaAdiabatic;in.cells=cells;
  in.cfl=0.2;in.quadraticViscosity=2;in.linearViscosity=0;
  in.shockThreshold=0.01;

  PistonInitialState state;
  state.nodeVelocityMPerS.assign(cells+1,0);
  state.cellDensityKgM3.resize(cells);state.cellPressurePa.resize(cells);
  state.cellTransverseMagneticFieldT.resize(cells);
  state.cellTransverseMagneticField2T.resize(cells);
  for(int i=0;i<cells;++i) {
    const double left=inner+length*i/cells,right=inner+length*(i+1)/cells;
    const auto sample=projection.value->Evaluate(ShellCentroid(left,right));
    Require(sample.ok()&&sample.value.region==AmbientRegion::PfssClosed,
        "initial ambient cells cross a topology branch");
    state.cellDensityKgM3[i]=sample.value.densityKgM3;
    state.cellPressurePa[i]=sample.value.pressurePa;
    state.cellTransverseMagneticFieldT[i]=
        sample.value.transverseMagneticFieldT[0];
    state.cellTransverseMagneticField2T[i]=
        sample.value.transverseMagneticFieldT[1];
  }
  const auto piston=[inner](double) {
    return SEP::Core::Result<KinematicValue>::Success({inner,0,0});
  };
  const auto source=[table=sourceTable.value](double radius,double) {
    return table->Evaluate(radius);
  };
  auto solver=PlanarPistonSolver::CreateInitialized(in,piston,std::move(state),
      source);
  Require(solver.ok(),"well-balance solver creation failed: "+solver.status.message);
  const auto advanced=solver.value->AdvanceTo(end);
  Require(advanced.ok(),"well-balance advance failed: "+advanced.message);
  const auto output=solver.value->Cells();
  const auto ledger=solver.value->EnergyLedger();
  Require(output.ok()&&ledger.ok(),"well-balance output unavailable");
  Errors error;
  for(int cellIndex=0;cellIndex<cells;++cellIndex) {
    const auto& cell=output.value[cellIndex];
    const auto exact=projection.value->Evaluate(cell.centerM);
    Require(exact.ok(),"final ambient reference left support");
    error.density=std::max(error.density,std::abs(cell.densityKgM3-
        exact.value.densityKgM3)/exact.value.densityKgM3);
    error.pressure=std::max(error.pressure,std::abs(cell.pressurePa-
        exact.value.pressurePa)/exact.value.pressurePa);
    error.velocity=std::max(error.velocity,std::abs(cell.velocityMPerS)/
        exact.value.canonicalFastSpeedMPerS);
    const double bExact=std::hypot(exact.value.transverseMagneticFieldT[0],
        exact.value.transverseMagneticFieldT[1]);
    const double bNumerical=std::hypot(cell.transverseMagneticFieldT,
        cell.transverseMagneticField2T);
    error.magnetic=std::max(error.magnetic,
        std::abs(bNumerical-bExact)/bExact);
    const double gravity=std::abs(exact.value.gravityAccelerationMPerS2);
    error.maxForceRatio=std::max(error.maxForceRatio,
        std::abs(exact.value.ambientMaintainingAccelerationMPerS2)/gravity);
    error.maxHeating=std::max(error.maxHeating,
        std::abs(exact.value.source.heatingWPerM3));
  }
  error.energy=std::abs(ledger.value.residualJ)/
      std::max(std::abs(ledger.value.initialEnergyJ),1e-300);
  error.maxForceRatio=sourceTable.value->MaximumAmbientForceOverGravity();
  error.maxHeating=sourceTable.value->MaximumAbsoluteHeatingWPerM3();
  return error;
}

ParkerErrors RunParker(int cells) {
  constexpr double au=149597870700.0,duration=10000;
  const auto event=Event();
  const auto ambient=AmbientModel::Create(event);
  Require(ambient.ok(),"Parker ambient construction failed");
  const Vec3 ray={std::sqrt(0.5),0,std::sqrt(0.5)};
  const auto projection=PistonAmbientProjection::Create(
      ambient.value,ray,event->support.startS,2e-5);
  Require(projection.ok(),"Parker projection failed");
  const double inner=0.9*au,length=0.05*au,outer=inner+length;
  const auto innerPath=PistonAmbientTrajectory::Create(projection.value,
      inner,0,duration,500);
  const auto outerPath=PistonAmbientTrajectory::Create(projection.value,
      outer,0,duration,500);
  Require(innerPath.ok()&&outerPath.ok(),"Parker material paths failed");
  const auto finalInner=innerPath.value->Evaluate(duration),
      finalOuter=outerPath.value->Evaluate(duration);
  Require(finalInner.ok()&&finalOuter.ok(),"Parker paths lack final support");
  const auto sourceTable=PistonAmbientSourceTable::Create(projection.value,
      inner,finalOuter.value.value,2*cells+1);
  Require(sourceTable.ok(),"Parker source table failed");
  const auto first=projection.value->Evaluate(inner),last=projection.value->Evaluate(outer);
  Require(first.ok()&&last.ok()&&first.value.region==
      SEP::CoronaSwcme::AmbientRegion::ParkerExterior,
      "Parker well-balance interval is not on the exterior branch");

  PlanarPistonInput in;
  in.geometry=PistonTubeGeometry::RadialSpherical;
  in.startS=0;in.endS=duration;in.leftPositionM=inner;
  in.columnLengthM=length;in.solidAngleSr=0.1;
  in.initialDensityKgM3=first.value.densityKgM3;
  in.initialPressurePa=first.value.pressurePa;
  in.initialVelocityMPerS=last.value.radialVelocityMPerS;
  in.initialTransverseMagneticFieldT=first.value.transverseMagneticFieldT[0];
  in.initialTransverseMagneticField2T=first.value.transverseMagneticFieldT[1];
  in.initialRadialMagneticFieldT=first.value.radialMagneticFieldT;
  in.gammaAdiabatic=event->composition.gammaAdiabatic;in.cells=cells;
  in.cfl=0.2;in.quadraticViscosity=2;in.shockThreshold=0.01;
  PistonInitialState state;
  state.nodeVelocityMPerS.resize(cells+1);
  state.cellDensityKgM3.resize(cells);state.cellPressurePa.resize(cells);
  state.cellTransverseMagneticFieldT.resize(cells);
  state.cellTransverseMagneticField2T.resize(cells);
  state.cellRadialMagneticFieldT.resize(cells);
  std::vector<double> initialRadialFlux(cells);
  for(int node=0;node<=cells;++node) {
    const auto sample=projection.value->Evaluate(inner+length*node/cells);
    Require(sample.ok(),"Parker node initialization failed");
    state.nodeVelocityMPerS[node]=sample.value.radialVelocityMPerS;
  }
  for(int cell=0;cell<cells;++cell) {
    const double leftRadius=inner+length*cell/cells;
    const double rightRadius=inner+length*(cell+1)/cells;
    const auto sample=projection.value->Evaluate(
        ShellCentroid(leftRadius,rightRadius));
    Require(sample.ok(),"Parker cell initialization failed");
    state.cellDensityKgM3[cell]=sample.value.densityKgM3;
    state.cellPressurePa[cell]=sample.value.pressurePa;
    state.cellTransverseMagneticFieldT[cell]=
        sample.value.transverseMagneticFieldT[0];
    state.cellTransverseMagneticField2T[cell]=
        sample.value.transverseMagneticFieldT[1];
    state.cellRadialMagneticFieldT[cell]=sample.value.radialMagneticFieldT;
    initialRadialFlux[cell]=sample.value.radialMagneticFieldT*
        sample.value.radiusM*sample.value.radiusM;
  }
  const auto piston=[path=innerPath.value](double time) {
    return path->Evaluate(time);
  };
  const auto outerBoundary=[path=outerPath.value](double time) {
    return path->Evaluate(time);
  };
  const auto source=[table=sourceTable.value](double radius,double) {
    return table->Evaluate(radius);
  };
  auto solver=PlanarPistonSolver::CreateInitialized(in,piston,std::move(state),
      source,outerBoundary);
  Require(solver.ok(),"Parker solver creation failed: "+solver.status.message);
  const auto advanced=solver.value->AdvanceTo(duration);
  Require(advanced.ok(),"Parker advance failed: "+advanced.message);
  const auto output=solver.value->Cells();
  const auto ledger=solver.value->EnergyLedger();
  Require(output.ok()&&ledger.ok(),"Parker output is unavailable");
  ParkerErrors error;
  for(int cellIndex=0;cellIndex<cells;++cellIndex) {
    const auto& cell=output.value[cellIndex];
    const auto exact=projection.value->Evaluate(cell.centerM);
    Require(exact.ok(),"Parker final reference failed");
    error.density=std::max(error.density,std::abs(cell.densityKgM3-
        exact.value.densityKgM3)/exact.value.densityKgM3);
    error.pressure=std::max(error.pressure,std::abs(cell.pressurePa-
        exact.value.pressurePa)/exact.value.pressurePa);
    error.velocity=std::max(error.velocity,std::abs(cell.velocityMPerS-
        exact.value.radialVelocityMPerS)/exact.value.canonicalFastSpeedMPerS);
    const double exactB=std::sqrt(
        exact.value.transverseMagneticFieldT[0]*
            exact.value.transverseMagneticFieldT[0]+
        exact.value.transverseMagneticFieldT[1]*
            exact.value.transverseMagneticFieldT[1]+
        exact.value.radialMagneticFieldT*exact.value.radialMagneticFieldT);
    const double numericalB=std::sqrt(
        cell.transverseMagneticFieldT*cell.transverseMagneticFieldT+
        cell.transverseMagneticField2T*cell.transverseMagneticField2T+
        cell.radialMagneticFieldT*cell.radialMagneticFieldT);
    error.magnetic=std::max(error.magnetic,
        std::abs(numericalB-exactB)/exactB);
    error.radialFlux=std::max(error.radialFlux,std::abs(
        cell.radialMagneticFieldT*cell.centerM*cell.centerM-
        initialRadialFlux[cellIndex])/
        std::max(std::abs(initialRadialFlux[cellIndex]),1e-300));
  }
  error.energy=std::abs(ledger.value.residualJ)/ledger.value.initialEnergyJ;
  error.maximumHeating=sourceTable.value->MaximumAbsoluteHeatingWPerM3();
  const auto middle=innerPath.value->Evaluate(0.5*duration);
  Require(middle.ok(),"Parker trajectory midpoint left support");
  const auto middleState=projection.value->Evaluate(middle.value.value);
  Require(middleState.ok(),"Parker trajectory midpoint failed");
  error.trajectoryVelocity=std::abs(middle.value.firstDerivative-
      middleState.value.radialVelocityMPerS)/middleState.value.radialVelocityMPerS;
  return error;
}

} // namespace

int main() {
  try {
    const double sourceLedger=TestSourceLedger();
    const std::array<int,3> cells={64,128,256};
    std::array<Errors,3> e;
    std::array<ParkerErrors,3> parker;
    for(int i=0;i<3;++i)e[i]=Run(cells[i]);
    for(int i=0;i<3;++i)parker[i]=RunParker(cells[i]);
    std::cout<<"[CSWC0625-EVIDENCE] density="<<e[0].density<<','<<e[1].density
             <<','<<e[2].density<<" pressure="<<e[0].pressure<<','
             <<e[1].pressure<<','<<e[2].pressure<<" velocity="
             <<e[0].velocity<<','<<e[1].velocity<<','<<e[2].velocity
             <<" magnetic="<<e[0].magnetic<<','<<e[1].magnetic<<','
             <<e[2].magnetic<<" energy="<<e[0].energy<<','<<e[1].energy
             <<','<<e[2].energy<<" max_famb_over_g="<<e[2].maxForceRatio
             <<" max_heating_w_m3="<<e[2].maxHeating
             <<" source_ledger_residual="<<sourceLedger<<'\n';
    std::cout<<"[CSWC0625-PARKER] density="<<parker[0].density<<','
             <<parker[1].density<<','<<parker[2].density<<" pressure="
             <<parker[0].pressure<<','<<parker[1].pressure<<','
             <<parker[2].pressure<<" velocity="<<parker[0].velocity<<','
             <<parker[1].velocity<<','<<parker[2].velocity<<" magnetic="
             <<parker[0].magnetic<<','<<parker[1].magnetic<<','
             <<parker[2].magnetic<<" energy="<<parker[0].energy<<','
             <<parker[1].energy<<','<<parker[2].energy
             <<" max_heating_w_m3="<<parker[2].maximumHeating
             <<" radial_flux_error="<<parker[2].radialFlux
             <<" trajectory_velocity_error="<<parker[2].trajectoryVelocity<<'\n';
    // These gates were selected before seeing the numerical results.  The
    // finite-volume equilibrium is not made exact with a mesh-dependent
    // cancelling source; it must converge under physical-profile sampling.
    Require(e[2].density<2e-4&&e[2].pressure<2e-4&&
        e[2].velocity<2e-4&&e[2].magnetic<2e-4&&e[2].energy<2e-6,
        "projected ambient exceeds frozen well-balance bounds");
    Require(sourceLedger<2e-12,
        "independent heating/induction source ledger does not close");
    Require(e[2].density<e[1].density&&e[1].density<e[0].density&&
        e[2].pressure<e[1].pressure&&e[1].pressure<e[0].pressure&&
        e[2].velocity<e[1].velocity&&e[1].velocity<e[0].velocity&&
        e[2].magnetic<e[1].magnetic&&e[1].magnetic<e[0].magnetic,
        "projected ambient errors do not decrease at three resolutions");
    Require(e[2].maxForceRatio>0&&std::isfinite(e[2].maxForceRatio)&&
        e[2].maxHeating<1e-30,
        "ambient source magnitude ledger is incomplete or closed-ray heating is nonzero");
    Require(parker[2].density<5e-4&&parker[2].pressure<5e-4&&
        parker[2].velocity<5e-4&&parker[2].magnetic<5e-4&&
        parker[2].energy<5e-6&&parker[2].maximumHeating>0&&
        parker[2].radialFlux<2e-13&&
        parker[2].trajectoryVelocity<2e-8,
        "advecting Parker ambient exceeds preregistered well-balance bounds");
    Require(parker[2].density<parker[1].density&&
        parker[1].density<parker[0].density&&
        parker[2].pressure<parker[1].pressure&&
        parker[1].pressure<parker[0].pressure&&
        parker[2].velocity<parker[1].velocity&&
        parker[1].velocity<parker[0].velocity&&
        parker[2].magnetic<parker[1].magnetic&&
        parker[1].magnetic<parker[0].magnetic&&
        parker[2].energy<parker[1].energy&&
        parker[1].energy<parker[0].energy,
        "advecting Parker errors do not decrease at three resolutions");
    std::cout<<"[CSWC0625] PASS projected coronal ambient well balance\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"CSWC0625 FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
