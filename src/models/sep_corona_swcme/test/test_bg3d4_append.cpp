#include "sep_corona_swcme/piston_ambient.h"
#include "sep_corona_swcme/piston_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronaSwcme::PistonAmbientProjection;
using SEP::CoronaSwcme::PistonAmbientSourceTable;
using SEP::CoronaSwcme::PistonAmbientTrajectory;
using SEP::CoronaSwcme::PistonAppendState;
using SEP::CoronaSwcme::PistonInitialState;
using SEP::CoronaSwcme::PistonTubeGeometry;
using SEP::CoronaSwcme::PlanarPistonInput;
using SEP::CoronaSwcme::PlanarPistonSolver;
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
  Require(event.ok(),event.status.message);return event.value;
}

double ShellCentroid(double left,double right) {
  return 0.75*(std::pow(right,4)-std::pow(left,4))/
      (std::pow(right,3)-std::pow(left,3));
}

PistonInitialState AmbientState(
    const std::shared_ptr<const PistonAmbientProjection>& projection,
    const std::vector<double>& nodes) {
  PistonInitialState state;
  const std::size_t cells=nodes.size()-1;
  state.nodeVelocityMPerS.resize(nodes.size());
  state.cellDensityKgM3.resize(cells);state.cellPressurePa.resize(cells);
  state.cellTransverseMagneticFieldT.resize(cells);
  state.cellTransverseMagneticField2T.resize(cells);
  state.cellRadialMagneticFieldT.resize(cells);
  for(std::size_t i=0;i<nodes.size();++i) {
    const auto sample=projection->Evaluate(nodes[i]);
    Require(sample.ok(),"ambient append node left projection support");
    state.nodeVelocityMPerS[i]=sample.value.radialVelocityMPerS;
  }
  for(std::size_t i=0;i<cells;++i) {
    const auto sample=projection->Evaluate(ShellCentroid(nodes[i],nodes[i+1]));
    Require(sample.ok(),"ambient append cell left projection support");
    state.cellDensityKgM3[i]=sample.value.densityKgM3;
    state.cellPressurePa[i]=sample.value.pressurePa;
    state.cellTransverseMagneticFieldT[i]=
        sample.value.transverseMagneticFieldT[0];
    state.cellTransverseMagneticField2T[i]=
        sample.value.transverseMagneticFieldT[1];
    state.cellRadialMagneticFieldT[i]=sample.value.radialMagneticFieldT;
  }
  return state;
}

struct Row {
  double density=0,pressure=0,velocity=0,magnetic=0,energy=0,mass=0;
  double appendedEnergy=0;
};

Row Run(int initialCells) {
  constexpr double au=149597870700.0,appendTime=5000,end=10000;
  const double inner0=0.9*au,outer0=0.95*au,newOuter0=1.0*au;
  const auto event=Event();
  const auto ambient=AmbientModel::Create(event);
  Require(ambient.ok(),"append ambient creation failed");
  const auto projection=PistonAmbientProjection::Create(ambient.value,
      Vec3{std::sqrt(0.5),0,std::sqrt(0.5)},0,2e-5);
  Require(projection.ok(),"append ambient projection failed");
  const auto innerPath=PistonAmbientTrajectory::Create(projection.value,
      inner0,0,end,500);
  const auto outerPath=PistonAmbientTrajectory::Create(projection.value,
      outer0,0,end,500);
  const auto newOuterPath=PistonAmbientTrajectory::Create(projection.value,
      newOuter0,0,end,500);
  Require(innerPath.ok()&&outerPath.ok()&&newOuterPath.ok(),
      "append boundary trajectory creation failed");
  const auto finalOuter=newOuterPath.value->Evaluate(end);
  Require(finalOuter.ok(),"new outer trajectory lacks final state");
  const auto sourceTable=PistonAmbientSourceTable::Create(projection.value,
      inner0,finalOuter.value.value,4*initialCells+1);
  Require(sourceTable.ok(),"append source table failed");

  std::vector<double> initialNodes(initialCells+1);
  for(int i=0;i<=initialCells;++i)
    initialNodes[i]=inner0+(outer0-inner0)*i/initialCells;
  auto initial=AmbientState(projection.value,initialNodes);
  const auto first=projection.value->Evaluate(inner0),last=projection.value->Evaluate(outer0);
  Require(first.ok()&&last.ok(),"append endpoints lack ambient state");
  PlanarPistonInput in;
  in.geometry=PistonTubeGeometry::RadialSpherical;in.startS=0;in.endS=end;
  in.leftPositionM=inner0;in.columnLengthM=outer0-inner0;
  in.solidAngleSr=0.1;in.initialDensityKgM3=first.value.densityKgM3;
  in.initialPressurePa=first.value.pressurePa;
  in.initialVelocityMPerS=last.value.radialVelocityMPerS;
  in.initialTransverseMagneticFieldT=first.value.transverseMagneticFieldT[0];
  in.initialTransverseMagneticField2T=first.value.transverseMagneticFieldT[1];
  in.initialRadialMagneticFieldT=first.value.radialMagneticFieldT;
  in.gammaAdiabatic=event->composition.gammaAdiabatic;in.cells=initialCells;
  in.cfl=0.2;in.quadraticViscosity=2;in.shockThreshold=0.01;
  const auto piston=[path=innerPath.value](double time) {return path->Evaluate(time);};
  const auto oldOuter=[path=outerPath.value](double time) {return path->Evaluate(time);};
  const auto source=[table=sourceTable.value](double radius,double) {
    return table->Evaluate(radius);
  };
  auto solver=PlanarPistonSolver::CreateInitialized(in,piston,std::move(initial),
      source,oldOuter);
  Require(solver.ok(),"append solver creation failed: "+solver.status.message);
  Require(solver.value->AdvanceTo(appendTime).ok(),"pre-append advance failed");

  const auto oldBoundary=outerPath.value->Evaluate(appendTime);
  const auto newBoundary=newOuterPath.value->Evaluate(appendTime);
  Require(oldBoundary.ok()&&newBoundary.ok(),"append epoch boundary unavailable");
  std::vector<double> appendNodes(initialCells+1);
  for(int i=0;i<=initialCells;++i)appendNodes[i]=oldBoundary.value.value+
      (newBoundary.value.value-oldBoundary.value.value)*i/initialCells;
  PistonAppendState valid;
  valid.nodePositionM=appendNodes;
  valid.state=AmbientState(projection.value,appendNodes);
  valid.outerBoundary=[path=newOuterPath.value](double time) {
    return path->Evaluate(time);
  };

  // Deliberately corrupt only the shared node.  The rejection must preserve
  // cell count, positions, energy ledger and the old outer-boundary authority.
  PistonAppendState invalid=valid;
  invalid.nodePositionM.front()+=1000;
  const std::size_t beforeCount=solver.value->CellCount();
  const auto beforeNodes=solver.value->NodePositionsM();
  const auto beforeLedger=solver.value->EnergyLedger();
  const auto rejected=solver.value->AppendAmbient(std::move(invalid));
  Require(!rejected.ok()&&solver.value->CellCount()==beforeCount&&
      solver.value->NodePositionsM()==beforeNodes,
      "failed append changed committed material inventory");
  const auto afterRejectedLedger=solver.value->EnergyLedger();
  Require(beforeLedger.ok()&&afterRejectedLedger.ok()&&
      beforeLedger.value.currentEnergyJ==afterRejectedLedger.value.currentEnergyJ&&
      afterRejectedLedger.value.appendedEnergyJ==0,
      "failed append changed committed energy identity");

  const auto accepted=solver.value->AppendAmbient(std::move(valid));
  Require(accepted.ok()&&solver.value->CellCount()==2*beforeCount,
      "valid ambient append failed: "+accepted.message);
  const auto appendedCells=solver.value->Cells();
  const auto appendedLedger=solver.value->EnergyLedger();
  Require(appendedCells.ok()&&appendedLedger.ok()&&
      appendedLedger.value.appendedEnergyJ>0,
      "accepted append lacks material/energy receipt");
  const double massAfterAppend=std::accumulate(appendedCells.value.begin(),
      appendedCells.value.end(),0.0,[](double sum,const auto& cell) {
        return sum+cell.massKg;
      });
  Require(solver.value->AdvanceTo(end).ok(),"post-append advance failed");
  const auto output=solver.value->Cells();
  const auto ledger=solver.value->EnergyLedger();
  Require(output.ok()&&ledger.ok(),"post-append output unavailable");
  Row row;double finalMass=0;
  for(const auto& cell:output.value) {
    finalMass+=cell.massKg;
    const auto exact=projection.value->Evaluate(cell.centerM);
    Require(exact.ok(),"post-append reference left support");
    row.density=std::max(row.density,std::abs(cell.densityKgM3-
        exact.value.densityKgM3)/exact.value.densityKgM3);
    row.pressure=std::max(row.pressure,std::abs(cell.pressurePa-
        exact.value.pressurePa)/exact.value.pressurePa);
    row.velocity=std::max(row.velocity,std::abs(cell.velocityMPerS-
        exact.value.radialVelocityMPerS)/exact.value.canonicalFastSpeedMPerS);
    const double exactB=std::sqrt(
        exact.value.radialMagneticFieldT*exact.value.radialMagneticFieldT+
        exact.value.transverseMagneticFieldT[0]*
            exact.value.transverseMagneticFieldT[0]+
        exact.value.transverseMagneticFieldT[1]*
            exact.value.transverseMagneticFieldT[1]);
    const double actualB=std::sqrt(
        cell.radialMagneticFieldT*cell.radialMagneticFieldT+
        cell.transverseMagneticFieldT*cell.transverseMagneticFieldT+
        cell.transverseMagneticField2T*cell.transverseMagneticField2T);
    row.magnetic=std::max(row.magnetic,std::abs(actualB-exactB)/exactB);
  }
  row.mass=std::abs(finalMass-massAfterAppend)/massAfterAppend;
  row.energy=std::abs(ledger.value.residualJ)/ledger.value.initialEnergyJ;
  row.appendedEnergy=ledger.value.appendedEnergyJ;
  return row;
}

} // namespace

int main() {
  try {
    const std::array<int,3> cells={32,64,128};
    std::array<Row,3> row;
    for(int i=0;i<3;++i)row[i]=Run(cells[i]);
    std::cout<<"[CSWC0626-EVIDENCE] density="<<row[0].density<<','
             <<row[1].density<<','<<row[2].density<<" pressure="
             <<row[0].pressure<<','<<row[1].pressure<<','<<row[2].pressure
             <<" velocity="<<row[0].velocity<<','<<row[1].velocity<<','
             <<row[2].velocity<<" magnetic="<<row[0].magnetic<<','
             <<row[1].magnetic<<','<<row[2].magnetic<<" energy="
             <<row[0].energy<<','<<row[1].energy<<','<<row[2].energy
             <<" mass="<<row[0].mass<<','<<row[1].mass<<','<<row[2].mass
             <<" appended_energy_j="<<row[2].appendedEnergy<<'\n';
    Require(row[2].density<5e-4&&row[2].pressure<5e-4&&
        row[2].velocity<5e-4&&row[2].magnetic<5e-4&&
        row[2].energy<5e-6&&row[2].mass<2e-15&&row[2].appendedEnergy>0,
        "ambient append exceeds preregistered conservation bounds");
    Require(row[2].density<row[1].density&&row[1].density<row[0].density&&
        row[2].pressure<row[1].pressure&&row[1].pressure<row[0].pressure&&
        row[2].velocity<row[1].velocity&&row[1].velocity<row[0].velocity&&
        row[2].magnetic<row[1].magnetic&&row[1].magnetic<row[0].magnetic&&
        row[2].energy<row[1].energy&&row[1].energy<row[0].energy,
        "ambient append does not converge under material refinement");
    std::cout<<"[CSWC0626] PASS transactional ambient append and ledgers\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"CSWC0626 FAIL: "<<error.what()<<'\n';return 1;
  }
}
