#include "sep_corona_swcme/piston_solver.h"

#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::CoronaSwcme::PistonTubeGeometry;
using SEP::CoronaSwcme::PistonInitialState;
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

struct TaylorReference {
  double pistonToShockRadius = 0.0;
  double pistonPressureOverRhoVs2 = 0.0;
  std::vector<double> xi,u,density,pressure;
};

struct TaylorState { double u=0,density=0,pressure=0; };

TaylorState EvaluateTaylor(const TaylorReference& reference,double xi) {
  if(xi>=1)return {reference.u.front(),reference.density.front(),
      reference.pressure.front()};
  if(xi<=reference.pistonToShockRadius)return {reference.u.back(),
      reference.density.back(),reference.pressure.back()};
  for(std::size_t i=1;i<reference.xi.size();++i)if(reference.xi[i]<=xi) {
    const double w=(reference.xi[i-1]-xi)/
        (reference.xi[i-1]-reference.xi[i]);
    return {(1-w)*reference.u[i-1]+w*reference.u[i],
        (1-w)*reference.density[i-1]+w*reference.density[i],
        (1-w)*reference.pressure[i-1]+w*reference.pressure[i]};
  }
  throw std::runtime_error("Taylor interpolation left similarity support");
}

// Independent integration of the constant-speed strong spherical-piston
// similarity equations (Taylor 1946).  With xi=r/(V_s t), u=V_s U,
// rho=rho_0 G and p=rho_0 V_s^2 P, the smooth post-shock ODE is
//
//   G'/G = -2 U (U-xi) / [xi ((U-xi)^2-c^2)],
//   U'   =  2 U c^2       / [xi ((U-xi)^2-c^2)],
//   P'/P = gamma G'/G,        c^2=gamma P/G.
//
// Strong-shock data are imposed at xi=1.  The material piston is the first
// inward point where U=xi; thus xi_p=U_p/V_s and P(xi_p) is the face pressure.
TaylorReference IntegrateTaylor(double gamma) {
  struct State { double u,g,p; };
  auto derivative=[&](double xi,const State& s) {
    const double w=s.u-xi,c2=gamma*s.p/s.g;
    const double logG=-2*s.u*w/(xi*(w*w-c2));
    return State{2*s.u*c2/(xi*(w*w-c2)),s.g*logG,
        gamma*s.p*logG};
  };
  auto add=[](State a,double factor,State b) {
    return State{a.u+factor*b.u,a.g+factor*b.g,a.p+factor*b.p};
  };
  double xi=1;
  State state={2/(gamma+1),(gamma+1)/(gamma-1),2/(gamma+1)};
  TaylorReference result;
  result.xi.push_back(xi);result.u.push_back(state.u);
  result.density.push_back(state.g);result.pressure.push_back(state.p);
  constexpr double step=-2e-5;
  double previousXi=xi;State previous=state;
  for(int iteration=0;iteration<50000;++iteration) {
    previousXi=xi;previous=state;
    const State k1=derivative(xi,state);
    const State k2=derivative(xi+0.5*step,add(state,0.5*step,k1));
    const State k3=derivative(xi+0.5*step,add(state,0.5*step,k2));
    const State k4=derivative(xi+step,add(state,step,k3));
    state={state.u+step*(k1.u+2*k2.u+2*k3.u+k4.u)/6,
        state.g+step*(k1.g+2*k2.g+2*k3.g+k4.g)/6,
        state.p+step*(k1.p+2*k2.p+2*k3.p+k4.p)/6};
    xi+=step;
    if(state.u-xi>=0) {
      const double before=previous.u-previousXi,after=state.u-xi;
      const double weight=-before/(after-before);
      result.pistonToShockRadius=previousXi+weight*(xi-previousXi);
      result.pistonPressureOverRhoVs2=
          previous.p+weight*(state.p-previous.p);
      result.xi.push_back(result.pistonToShockRadius);
      result.u.push_back(previous.u+weight*(state.u-previous.u));
      result.density.push_back(previous.g+weight*(state.g-previous.g));
      result.pressure.push_back(result.pistonPressureOverRhoVs2);
      return result;
    }
    result.xi.push_back(xi);result.u.push_back(state.u);
    result.density.push_back(state.g);result.pressure.push_back(state.p);
  }
  throw std::runtime_error("independent Taylor ODE did not reach the piston");
}

struct Row { double radius=0,speed=0,pressure=0,energy=0; };

Row Run(int cells,const TaylorReference& reference) {
  constexpr double gamma=5.0/3.0,pistonSpeed=1,start=0.2,end=0.8;
  const double shockSpeed=pistonSpeed/reference.pistonToShockRadius;
  PlanarPistonInput in;
  in.geometry=PistonTubeGeometry::RadialSpherical;
  in.startS=start;in.endS=end;in.leftPositionM=pistonSpeed*start;
  in.columnLengthM=3-in.leftPositionM;in.solidAngleSr=1;
  in.initialDensityKgM3=1;in.initialPressurePa=1e-3;
  in.initialVelocityMPerS=0;in.gammaAdiabatic=gamma;
  in.cells=cells;in.cfl=0.2;in.quadraticViscosity=2;
  in.linearViscosity=0;in.shockThreshold=0.02;
  const auto history=[](double time) {
    return Result<KinematicValue>::Success(
        {pistonSpeed*time,pistonSpeed,0});
  };
  PistonInitialState initial;
  initial.nodeVelocityMPerS.resize(cells+1);
  initial.cellDensityKgM3.resize(cells);
  initial.cellPressurePa.resize(cells);
  initial.cellTransverseMagneticFieldT.assign(cells,0);
  const double shockAtStart=shockSpeed*start;
  for(int node=0;node<=cells;++node) {
    const double radius=in.leftPositionM+in.columnLengthM*node/cells;
    initial.nodeVelocityMPerS[node]=radius<=shockAtStart?
        shockSpeed*EvaluateTaylor(reference,radius/shockAtStart).u:0;
  }
  initial.nodeVelocityMPerS.front()=pistonSpeed;
  initial.nodeVelocityMPerS.back()=0;
  for(int cell=0;cell<cells;++cell) {
    const double left=in.leftPositionM+in.columnLengthM*cell/cells;
    const double right=in.leftPositionM+in.columnLengthM*(cell+1)/cells;
    const double radius=0.75*(std::pow(right,4)-std::pow(left,4))/
        (std::pow(right,3)-std::pow(left,3));
    if(radius<=shockAtStart) {
      const auto exact=EvaluateTaylor(reference,radius/shockAtStart);
      initial.cellDensityKgM3[cell]=exact.density;
      initial.cellPressurePa[cell]=shockSpeed*shockSpeed*exact.pressure;
    } else {
      initial.cellDensityKgM3[cell]=1;
      initial.cellPressurePa[cell]=in.initialPressurePa;
    }
  }
  auto solver=PlanarPistonSolver::CreateInitialized(in,history,
      std::move(initial));
  Require(solver.ok(),solver.status.message);
  const auto advanced=solver.value->AdvanceTo(end);
  Require(advanced.ok(),advanced.message);
  const auto shock=solver.value->DetectShock();
  Require(shock.ok()&&shock.value.present&&shock.value.statesAvailable,
      "spherical piston lacks detected shock/two-sided states");
  const auto state=solver.value->Cells();
  Require(state.ok()&&state.value.size()>=3,
      "spherical piston cells unavailable");
  double facePressure=0;
  for(int i=0;i<3;++i)facePressure+=state.value[i].pressurePa/3;
  const double exactPressure=shockSpeed*shockSpeed*
      reference.pistonPressureOverRhoVs2;
  const auto energy=solver.value->EnergyLedger();
  Require(energy.ok(),energy.status.message);
  const double scale=std::abs(energy.value.currentEnergyJ)+
      std::abs(energy.value.initialEnergyJ)+std::abs(energy.value.pistonWorkJ)+
      std::abs(energy.value.outerBoundaryWorkJ);
  return {Relative(shock.value.radiusM,shockSpeed*end),
      Relative(shock.value.speedMPerS,shockSpeed),
      Relative(facePressure,exactPressure),
      std::abs(energy.value.residualJ)/scale};
}

void TestTaylorSphere() {
  const auto reference=IntegrateTaylor(5.0/3.0);
  const std::array<int,3> cells={150,300,600};
  std::array<Row,3> row;
  for(std::size_t i=0;i<row.size();++i)row[i]=Run(cells[i],reference);
  std::cout<<"[CSWC0624-EVIDENCE] piston_shock_ratio="
           <<reference.pistonToShockRadius<<" piston_pressure="
           <<reference.pistonPressureOverRhoVs2<<" radius_errors="
           <<row[0].radius<<','<<row[1].radius<<','<<row[2].radius
           <<" speed_errors="<<row[0].speed<<','<<row[1].speed<<','
           <<row[2].speed<<" pressure_errors="<<row[0].pressure<<','
           <<row[1].pressure<<','<<row[2].pressure<<" energy_errors="
           <<row[0].energy<<','<<row[1].energy<<','<<row[2].energy<<'\n';
  Require(row[2].radius<0.02&&row[2].speed<0.02&&row[2].pressure<0.05&&
      row[2].energy<2e-4,"Taylor spherical piston exceeds frozen accuracy");
  Require(row[2].radius<row[1].radius&&row[1].radius<row[0].radius&&
      row[2].pressure<row[1].pressure&&row[1].pressure<row[0].pressure,
      "Taylor spherical piston does not converge under mass refinement");
}

struct MagneticRow { double energy=0,frozen=0; };

MagneticRow RunMagneticEnergy(int cells) {
  PlanarPistonInput in;
  in.geometry=PistonTubeGeometry::RadialSpherical;
  in.startS=0;in.endS=2;in.leftPositionM=10;in.columnLengthM=10;
  in.solidAngleSr=0.1;in.initialDensityKgM3=1;in.initialPressurePa=1;
  in.initialVelocityMPerS=0;in.initialTransverseMagneticFieldT=0.3;
  in.initialTransverseMagneticField2T=0.4;
  in.magneticPermeabilityNPerA2=1;in.gammaAdiabatic=5.0/3.0;
  in.cells=cells;in.cfl=0.2;in.quadraticViscosity=2;
  in.shockThreshold=0.01;
  const auto piston=[](double time) {
    // A modest constant-speed boundary launches a curved MHD compression.
    // The test does not claim an exact state solution; its independent
    // identities are frozen flux and total energy against boundary work.
    return Result<KinematicValue>::Success({10+0.2*time,0.2,0});
  };
  auto solver=PlanarPistonSolver::Create(in,piston);
  Require(solver.ok(),"spherical magnetic solver creation failed");
  const auto advanced=solver.value->AdvanceTo(2);
  Require(advanced.ok(),"spherical magnetic advance failed: "+advanced.message);
  const auto state=solver.value->Cells();
  const auto ledger=solver.value->EnergyLedger();
  Require(state.ok()&&ledger.ok(),"spherical magnetic output unavailable");
  MagneticRow out;
  for(int i=0;i<cells;++i) {
    const double left=10.0+10.0*i/cells,right=10.0+10.0*(i+1)/cells;
    const double initialRadius=0.75*(std::pow(right,4)-std::pow(left,4))/
        (std::pow(right,3)-std::pow(left,3));
    const double expected1=0.3/initialRadius,expected2=0.4/initialRadius;
    const double actual1=state.value[i].transverseMagneticFieldT/
        (state.value[i].densityKgM3*state.value[i].centerM);
    const double actual2=state.value[i].transverseMagneticField2T/
        (state.value[i].densityKgM3*state.value[i].centerM);
    out.frozen=std::max(out.frozen,std::max(Relative(actual1,expected1),
        Relative(actual2,expected2)));
  }
  const double scale=std::abs(ledger.value.currentEnergyJ)+
      std::abs(ledger.value.initialEnergyJ)+std::abs(ledger.value.pistonWorkJ)+
      std::abs(ledger.value.outerBoundaryWorkJ);
  out.energy=std::abs(ledger.value.residualJ)/scale;
  return out;
}

void TestSphericalMagneticEnergy() {
  const std::array<int,3> cells={100,200,400};
  std::array<MagneticRow,3> row;
  for(int i=0;i<3;++i)row[i]=RunMagneticEnergy(cells[i]);
  std::cout<<"[CSWC0624-MAGNETIC] energy_errors="<<row[0].energy<<','
           <<row[1].energy<<','<<row[2].energy<<" frozen_errors="
           <<row[0].frozen<<','<<row[1].frozen<<','<<row[2].frozen<<'\n';
  Require(row[2].energy<2e-4&&row[2].energy<row[1].energy&&
      row[1].energy<row[0].energy&&
      std::max({row[0].frozen,row[1].frozen,row[2].frozen})<2e-13,
      "spherical magnetic energy or frozen-flux identity is unqualified");
}

} // namespace

int main() {
  try {
    TestTaylorSphere();
    TestSphericalMagneticEnergy();
    std::cout<<"[CSWC0624] PASS Taylor sphere, curved MHD energy and volumes\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-4 spherical piston FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
