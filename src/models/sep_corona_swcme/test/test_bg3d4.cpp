#include "sep_corona_swcme/sheath_model.h"
#include "sep_coronal_cme/configuration_parser.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace {

using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronaSwcme::EventConfiguration;
using SEP::CoronaSwcme::SheathCellDisposition;
using SEP::CoronaSwcme::SheathMaterialLabel;
using SEP::CoronaSwcme::ShockFedSheathModel;
using SEP::CoronaSwcme::SurfaceShockModel;
using SEP::CoronalCME::Cross;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

bool Close(double a,double b,double tolerance) {
  return std::abs(a-b)<=tolerance*std::max({std::abs(a),std::abs(b),1.0});
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read sheath fixture: "+path);
  std::ostringstream bytes;bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"failed reading sheath fixture: "+path);
  return bytes.str();
}

std::shared_ptr<const EventConfiguration> FrozenEvent() {
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(
      FileBytes("examples/bg3d1/event.conf"),[](const std::string& path) {
        try {
          return SEP::Core::Result<std::string>::Success(FileBytes(path));
        } catch(const std::exception& error) {
          return SEP::Core::Result<std::string>::Failure(
              SEP::Core::StatusCode::DataIntegrityFailure,error.what());
        }
      });
  Require(result.ok(),result.status.message);return result.value;
}

std::shared_ptr<const EventConfiguration> EventWithContactFraction(
    const std::string& fraction) {
  std::string asset=FileBytes("examples/bg3d1/ejecta.asset");
  const std::string assignment="contact_apex_fraction=0.99";
  const std::size_t position=asset.find(assignment);
  Require(position!=std::string::npos,"contact fraction fixture assignment is absent");
  const std::string oldHash=SEP::CoronalCME::ComputeContentChecksum(asset);
  asset.replace(position,assignment.size(),"contact_apex_fraction="+fraction);
  const std::string newHash=SEP::CoronalCME::ComputeContentChecksum(asset);
  std::string configuration=FileBytes("examples/bg3d1/event.conf");
  const std::size_t digest=configuration.find(oldHash);
  Require(digest!=std::string::npos,"regional asset digest is absent from fixture");
  configuration.replace(digest,oldHash.size(),newHash);
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(configuration,
      [&](const std::string& path) {
        try {
          return SEP::Core::Result<std::string>::Success(
              path=="examples/bg3d1/ejecta.asset"?asset:FileBytes(path));
        } catch(const std::exception& error) {
          return SEP::Core::Result<std::string>::Failure(
              SEP::Core::StatusCode::DataIntegrityFailure,error.what());
        }
      });
  Require(result.ok(),result.status.message);return result.value;
}

struct Models {
  std::shared_ptr<const AmbientModel> ambient;
  std::shared_ptr<SurfaceShockModel> shock;
  std::shared_ptr<ShockFedSheathModel> sheath;
};

Models CreateModels(const std::shared_ptr<const EventConfiguration>& event) {
  const auto ambient=AmbientModel::Create(event);Require(ambient.ok(),ambient.status.message);
  const auto shock=SurfaceShockModel::Create(event,ambient.value);
  Require(shock.ok(),shock.status.message);
  const auto sheath=ShockFedSheathModel::Create(event,ambient.value);
  Require(sheath.ok(),sheath.status.message);
  return {ambient.value,shock.value,sheath.value};
}

double Determinant(const std::array<Vec3,3>& columns) {
  return Dot(columns[0],Cross(columns[1],columns[2]));
}

Vec3 Apply(const std::array<Vec3,3>& columns,Vec3 value) {
  return value.x*columns[0]+value.y*columns[1]+value.z*columns[2];
}

Vec3 Fourth(Vec3 minus2,Vec3 minus,Vec3 plus,Vec3 plus2,double step) {
  return (minus2-8*minus+8*plus-plus2)/(12*step);
}

double Fourth(double minus2,double minus,double plus,double plus2,double step) {
  return (minus2-8*minus+8*plus-plus2)/(12*step);
}

std::array<Vec3,3> SpatialGradient(const std::array<Vec3,3>& positions,
    const std::array<Vec3,3>& values) {
  const double determinant=Determinant(positions);
  Require(std::abs(determinant)>1e-30,"material-coordinate gradient is singular");
  const auto derivative=[&](Vec3 direction) {
    const double c0=Dot(direction,Cross(positions[1],positions[2]))/determinant;
    const double c1=Dot(direction,Cross(positions[2],positions[0]))/determinant;
    const double c2=Dot(direction,Cross(positions[0],positions[1]))/determinant;
    return c0*values[0]+c1*values[1]+c2*values[2];
  };
  return {derivative({1,0,0}),derivative({0,1,0}),derivative({0,0,1})};
}

void TestShockLimitAndCauchyMap() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  const double crossing=10000;
  const auto surface=models.shock->Prepare(crossing,31,12,24);
  Require(surface.ok(),surface.status.message);
  std::size_t index=0;
  while(index<surface.value->shocks->patches.size()&&
      !surface.value->shocks->patches[index].fast)++index;
  Require(index<surface.value->shocks->patches.size(),"fixture has no fast birth patch");
  const auto& geometry=surface.value->frontPatches[index].geometry;
  const auto& downstream=surface.value->shocks->patches[index].jump.downstream;
  const SheathMaterialLabel label{geometry.polarParameterRad,
      geometry.azimuthParameterRad,crossing,geometry.physicalId};
  const auto boundary=models.sheath->Evaluate(label,crossing);
  Require(boundary.ok(),boundary.status.message);
  Require(boundary.value.jacobian==1&&
      Norm(boundary.value.positionM-geometry.centerM)<=1e-12*Norm(geometry.centerM)&&
      Close(boundary.value.primitive.massDensityKgM3,downstream.massDensityKgM3,1e-14)&&
      Close(boundary.value.primitive.pressurePa,downstream.pressurePa,1e-14)&&
      Norm(boundary.value.primitive.velocityMPerS-downstream.velocityMPerS)<=1e-10&&
      Norm(boundary.value.primitive.magneticFieldT-downstream.magneticFieldT)<=1e-18,
      "shock-fed map does not recover the canonical downstream boundary state");

  const auto evolved=models.sheath->Evaluate(label,crossing+10);
  Require(evolved.ok(),evolved.status.message);
  const double determinant=Determinant(evolved.value.deformationColumns);
  Require(Close(determinant,evolved.value.jacobian,2e-7)&&
      evolved.value.jacobian>=event->regional.minimumJacobian,
      "deformation determinant and reported positive Jacobian disagree");
  const Vec3 expectedB=Apply(evolved.value.deformationColumns,
      downstream.magneticFieldT)/evolved.value.jacobian;
  Require(Norm(expectedB-evolved.value.primitive.magneticFieldT)<=
      2e-12*std::max(Norm(expectedB),1e-30)&&
      Close(evolved.value.primitive.massDensityKgM3*evolved.value.jacobian,
          downstream.massDensityKgM3,2e-12)&&
      Close(evolved.value.primitive.pressurePa*
          std::pow(evolved.value.jacobian,event->composition.gammaAdiabatic),
          downstream.pressurePa,2e-12),
      "Cauchy mass/flux or adiabatic material identity failed");
}

void TestInventoryAndRefinement() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  // The early, low-coronal interval exercises solar clipping, fast/sub-fast
  // admission, contact exits and the full closed angular parameterization.
  const double epoch=40;
  const auto coarse=models.sheath->PrepareInventory(epoch,41,4,8,2);
  Require(coarse.ok(),coarse.status.message);
  const auto medium=models.sheath->PrepareInventory(epoch,42,4,8,4);
  Require(medium.ok(),medium.status.message);
  const auto fine=models.sheath->PrepareInventory(epoch,43,4,8,8);
  Require(fine.ok(),fine.status.message);
  const auto reference=models.sheath->PrepareInventory(epoch,44,4,8,16);
  Require(reference.ok(),reference.status.message);
  const double e0=std::abs(coarse.value->admittedMassKg-reference.value->admittedMassKg);
  const double e1=std::abs(medium.value->admittedMassKg-reference.value->admittedMassKg);
  const double e2=std::abs(fine.value->admittedMassKg-reference.value->admittedMassKg);
  Require(e2<e1&&e1<e0,"admission-time mass quadrature does not converge");
  for(const auto* inventory:{coarse.value.get(),medium.value.get(),fine.value.get(),
      reference.value.get()}) {
    Require(inventory->maximumAdmissionMassResidual<=1e-9&&
        inventory->maximumInventoryMassResidual<=1e-12&&
        inventory->minimumJacobian>=event->regional.minimumJacobian&&
        inventory->fastAreaTimeM2S>0&&inventory->subfastAreaTimeM2S>0,
        "sheath admission/inventory/positive-J gate failed");
    std::set<std::tuple<std::uint64_t,double>> identities;
    for(const auto& cell:inventory->cells)Require(identities.insert(
        {cell.label.shockPatchLineage,cell.label.crossingTimeS}).second,
        "duplicate admitted material cell");
    const double closed=inventory->retainedMassKg+inventory->contactExitMassKg+
        inventory->solarExitMassKg+inventory->outerExitMassKg;
    Require(Close(closed,inventory->admittedMassKg,1e-13),
        "contact/solar/outer disposition ledger does not close");
  }
  const auto committed=models.sheath->Current();
  const auto failed=models.sheath->PrepareInventory(epoch,0,4,8,8);
  Require(!failed.ok()&&models.sheath->Current()==committed,
      "invalid inventory changed the committed sheath epoch");
  const auto exits=models.sheath->PrepareInventory(200,45,4,8,16);
  Require(exits.ok(),"first-exit inventory failed: "+exits.status.message);
  Require(exits.value->contactExitMassKg+exits.value->solarExitMassKg+
      exits.value->outerExitMassKg>0&&
      exits.value->maximumInventoryMassResidual<=1e-12,
      "longer inventory did not exercise and close explicit boundary exits");
  std::cout<<"[CMBGU05-LIMIT] permeable_contact_mass_kg="
           <<exits.value->contactExitMassKg
           <<" material_contact_qualified=false\n";
  const auto incompatibleEvent=EventWithContactFraction("0.82");
  auto incompatible=CreateModels(incompatibleEvent);
  const auto rejected=incompatible.sheath->PrepareInventory(200,46,4,8,16);
  Require(!rejected.ok()&&!incompatible.sheath->Current(),
      "incompatible thick contact hid a ballistic-map fold or published partial state");
  std::cout<<"[CMBGU05-EVIDENCE] mass_kg="<<reference.value->admittedMassKg
           <<" time_errors="<<e0<<','<<e1<<','<<e2
           <<" min_j="<<reference.value->minimumJacobian
           <<" admission_res="<<reference.value->maximumAdmissionMassResidual
           <<" inventory_res="<<reference.value->maximumInventoryMassResidual
           <<" retained="<<reference.value->retainedMassKg
           <<" contact_exit="<<reference.value->contactExitMassKg
           <<" solar_exit="<<reference.value->solarExitMassKg
           <<" outer_exit="<<reference.value->outerExitMassKg
           <<" long_contact_exit="<<exits.value->contactExitMassKg
           <<" long_solar_exit="<<exits.value->solarExitMassKg
           <<" long_outer_exit="<<exits.value->outerExitMassKg<<'\n';
}

void TestMaterialDifferentialIdentities() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  const double crossing=10000,current=10001;
  const auto surface=models.shock->Prepare(crossing,51,12,24);
  Require(surface.ok(),surface.status.message);
  const double angleStep=5e-4,timeLabelStep=0.5;
  SheathMaterialLabel center;
  bool found=false;
  for(std::size_t i=0;i<surface.value->frontPatches.size()&&!found;++i) {
    if(!surface.value->shocks->patches[i].fast)continue;
    const auto& p=surface.value->frontPatches[i].geometry;
    SheathMaterialLabel trial{p.polarParameterRad,p.azimuthParameterRad,
        crossing,p.physicalId};
    const SheathMaterialLabel neighbors[]={
      {trial.polarRad-2*angleStep,trial.azimuthRad,crossing,p.physicalId},
      {trial.polarRad-angleStep,trial.azimuthRad,crossing,p.physicalId},
      {trial.polarRad+angleStep,trial.azimuthRad,crossing,p.physicalId},
      {trial.polarRad+2*angleStep,trial.azimuthRad,crossing,p.physicalId},
      {trial.polarRad,trial.azimuthRad-2*angleStep,crossing,p.physicalId},
      {trial.polarRad,trial.azimuthRad-angleStep,crossing,p.physicalId},
      {trial.polarRad,trial.azimuthRad+angleStep,crossing,p.physicalId},
      {trial.polarRad,trial.azimuthRad+2*angleStep,crossing,p.physicalId},
      {trial.polarRad,trial.azimuthRad,crossing-2*timeLabelStep,p.physicalId},
      {trial.polarRad,trial.azimuthRad,crossing-timeLabelStep,p.physicalId},
      {trial.polarRad,trial.azimuthRad,crossing+timeLabelStep,p.physicalId},
      {trial.polarRad,trial.azimuthRad,crossing+2*timeLabelStep,p.physicalId}};
    bool valid=true;for(const auto& label:neighbors)
      valid=valid&&models.sheath->Evaluate(label,current).ok();
    if(valid){center=trial;found=true;}
  }
  Require(found,"no smooth interior fast patch is available for material derivatives");
  const auto sample=[&](SheathMaterialLabel label,double time) {
    const auto state=models.sheath->Evaluate(label,time);
    Require(state.ok(),state.status.message);return state.value;
  };
  const auto thetaMinus2=sample({center.polarRad-2*angleStep,center.azimuthRad,
      crossing,center.shockPatchLineage},current);
  const auto thetaMinus=sample({center.polarRad-angleStep,center.azimuthRad,
      crossing,center.shockPatchLineage},current);
  const auto thetaPlus=sample({center.polarRad+angleStep,center.azimuthRad,
      crossing,center.shockPatchLineage},current);
  const auto thetaPlus2=sample({center.polarRad+2*angleStep,center.azimuthRad,
      crossing,center.shockPatchLineage},current);
  const auto phiMinus2=sample({center.polarRad,center.azimuthRad-2*angleStep,
      crossing,center.shockPatchLineage},current);
  const auto phiMinus=sample({center.polarRad,center.azimuthRad-angleStep,
      crossing,center.shockPatchLineage},current);
  const auto phiPlus=sample({center.polarRad,center.azimuthRad+angleStep,
      crossing,center.shockPatchLineage},current);
  const auto phiPlus2=sample({center.polarRad,center.azimuthRad+2*angleStep,
      crossing,center.shockPatchLineage},current);
  const auto tauMinus2=sample({center.polarRad,center.azimuthRad,
      crossing-2*timeLabelStep,center.shockPatchLineage},current);
  const auto tauMinus=sample({center.polarRad,center.azimuthRad,
      crossing-timeLabelStep,center.shockPatchLineage},current);
  const auto tauPlus=sample({center.polarRad,center.azimuthRad,
      crossing+timeLabelStep,center.shockPatchLineage},current);
  const auto tauPlus2=sample({center.polarRad,center.azimuthRad,
      crossing+2*timeLabelStep,center.shockPatchLineage},current);
  const auto central=sample(center,current);
  const std::array<Vec3,3> dx={
      Fourth(thetaMinus2.positionM,thetaMinus.positionM,thetaPlus.positionM,
          thetaPlus2.positionM,angleStep),
      Fourth(phiMinus2.positionM,phiMinus.positionM,phiPlus.positionM,
          phiPlus2.positionM,angleStep),
      Fourth(tauMinus2.positionM,tauMinus.positionM,tauPlus.positionM,
          tauPlus2.positionM,timeLabelStep)};
  const std::array<Vec3,3> du={
      Fourth(thetaMinus2.primitive.velocityMPerS,thetaMinus.primitive.velocityMPerS,
          thetaPlus.primitive.velocityMPerS,thetaPlus2.primitive.velocityMPerS,angleStep),
      Fourth(phiMinus2.primitive.velocityMPerS,phiMinus.primitive.velocityMPerS,
          phiPlus.primitive.velocityMPerS,phiPlus2.primitive.velocityMPerS,angleStep),
      Fourth(tauMinus2.primitive.velocityMPerS,tauMinus.primitive.velocityMPerS,
          tauPlus.primitive.velocityMPerS,tauPlus2.primitive.velocityMPerS,timeLabelStep)};
  const auto gradientU=SpatialGradient(dx,du);
  const double divergenceU=gradientU[0].x+gradientU[1].y+gradientU[2].z;
  const double evolutionStep=0.1;
  const auto before2=sample(center,current-2*evolutionStep);
  const auto before=sample(center,current-evolutionStep);
  const auto after=sample(center,current+evolutionStep);
  const auto after2=sample(center,current+2*evolutionStep);
  const double densityRate=Fourth(before2.primitive.massDensityKgM3,
      before.primitive.massDensityKgM3,after.primitive.massDensityKgM3,
      after2.primitive.massDensityKgM3,evolutionStep);
  const double pressureRate=Fourth(before2.primitive.pressurePa,
      before.primitive.pressurePa,after.primitive.pressurePa,
      after2.primitive.pressurePa,evolutionStep);
  const Vec3 specificBefore2=before2.primitive.magneticFieldT/
      before2.primitive.massDensityKgM3;
  const Vec3 specificBefore=before.primitive.magneticFieldT/
      before.primitive.massDensityKgM3;
  const Vec3 specificAfter=after.primitive.magneticFieldT/
      after.primitive.massDensityKgM3;
  const Vec3 specificAfter2=after2.primitive.magneticFieldT/
      after2.primitive.massDensityKgM3;
  const Vec3 specific=central.primitive.magneticFieldT/
      central.primitive.massDensityKgM3;
  const Vec3 specificRate=Fourth(specificBefore2,specificBefore,specificAfter,
      specificAfter2,evolutionStep);
  const Vec3 stretching=Apply(gradientU,specific);
  const double massResidual=std::abs(densityRate+
      central.primitive.massDensityKgM3*divergenceU)/
      std::max(std::abs(densityRate)+
          std::abs(central.primitive.massDensityKgM3*divergenceU),1e-300);
  const double thermalResidual=std::abs(pressureRate+
      event->composition.gammaAdiabatic*central.primitive.pressurePa*divergenceU)/
      std::max(std::abs(pressureRate)+std::abs(event->composition.gammaAdiabatic*
          central.primitive.pressurePa*divergenceU),1e-300);
  const double inductionResidual=Norm(specificRate-stretching)/
      std::max(Norm(specificRate)+Norm(stretching),1e-300);
  std::cout<<"[CMBGU05-DIFFERENTIAL] mass="<<massResidual
           <<" thermal="<<thermalResidual<<" induction="<<inductionResidual<<'\n';
  Require(massResidual<=1e-7&&thermalResidual<=1e-7&&inductionResidual<=1e-6,
      "independent Lagrangian mass/adiabatic/induction residual exceeds budget");

  // A second oracle stays entirely on one material characteristic and
  // refines only its time stencil.  Fdot*F^-1 is the velocity gradient of a
  // regular map, so it checks the continuity, adiabatic and frozen-flux laws
  // without reusing the spatial-gradient construction above.
  std::array<double,3> refined{};
  int level=0;
  for(const double step:{0.4,0.2,0.1}) {
    const auto m2=sample(center,current-2*step),m1=sample(center,current-step);
    const auto p1=sample(center,current+step),p2=sample(center,current+2*step);
    std::array<Vec3,3> fDot;
    for(int column=0;column<3;++column)fDot[column]=Fourth(
        m2.deformationColumns[column],m1.deformationColumns[column],
        p1.deformationColumns[column],p2.deformationColumns[column],step);
    const auto velocityGradient=SpatialGradient(central.deformationColumns,fDot);
    const double jRate=Fourth(m2.jacobian,m1.jacobian,p1.jacobian,p2.jacobian,step);
    const double div=jRate/central.jacobian;
    const double rhoRate=Fourth(m2.primitive.massDensityKgM3,
        m1.primitive.massDensityKgM3,p1.primitive.massDensityKgM3,
        p2.primitive.massDensityKgM3,step);
    const double pRate=Fourth(m2.primitive.pressurePa,m1.primitive.pressurePa,
        p1.primitive.pressurePa,p2.primitive.pressurePa,step);
    const Vec3 qM2=m2.primitive.magneticFieldT/m2.primitive.massDensityKgM3;
    const Vec3 qM1=m1.primitive.magneticFieldT/m1.primitive.massDensityKgM3;
    const Vec3 qP1=p1.primitive.magneticFieldT/p1.primitive.massDensityKgM3;
    const Vec3 qP2=p2.primitive.magneticFieldT/p2.primitive.massDensityKgM3;
    const Vec3 qRate=Fourth(qM2,qM1,qP1,qP2,step);
    const Vec3 stretch=Apply(velocityGradient,specific);
    const double rMass=std::abs(rhoRate+central.primitive.massDensityKgM3*div)/
        std::max(std::abs(rhoRate)+
            std::abs(central.primitive.massDensityKgM3*div),1e-300);
    const double rThermal=std::abs(pRate+event->composition.gammaAdiabatic*
        central.primitive.pressurePa*div)/std::max(std::abs(pRate)+
        std::abs(event->composition.gammaAdiabatic*
            central.primitive.pressurePa*div),1e-300);
    const double rInduction=Norm(qRate-stretch)/
        std::max(Norm(qRate)+Norm(stretch),1e-300);
    refined[level++]=std::max({rMass,rThermal,rInduction});
  }
  std::cout<<"[CMBGU05-REFINEMENT] material_residuals="<<refined[0]<<','
           <<refined[1]<<','<<refined[2]<<'\n';
  Require(refined[2]<=1e-7&&refined[1]<=1e-7&&refined[0]<=1e-7,
      "three-level material time refinement exceeds the frozen smooth budget");
}

void TestSpatialInverse() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  const auto inventory=models.sheath->PrepareInventory(40,61,6,12,16);
  Require(inventory.ok(),inventory.status.message);
  const auto surface=models.shock->Prepare(39,62,12,24);
  Require(surface.ok(),surface.status.message);
  SheathMaterialLabel label;bool found=false;
  for(std::size_t i=0;i<surface.value->frontPatches.size();++i)
    if(surface.value->shocks->patches[i].fast) {
      const auto& patch=surface.value->frontPatches[i].geometry;
      label={patch.polarParameterRad+2e-4,patch.azimuthParameterRad+3e-4,
          39,patch.physicalId};
      if(models.sheath->Evaluate(label,40).ok()){found=true;break;}
    }
  Require(found,"no retained off-grid sheath point is available for inverse test");
  const auto direct=models.sheath->Evaluate(label,40);
  const auto inverse=models.sheath->EvaluateAtPosition(direct.value.positionM,40,0.1,24);
  Require(inverse.ok(),inverse.status.message);
  Require(Norm(inverse.value.positionM-direct.value.positionM)<=0.1&&
      Close(inverse.value.primitive.massDensityKgM3,
          direct.value.primitive.massDensityKgM3,2e-8)&&
      Close(inverse.value.primitive.pressurePa,direct.value.primitive.pressurePa,2e-8)&&
      Norm(inverse.value.primitive.magneticFieldT-direct.value.primitive.magneticFieldT)<=
          2e-8*std::max(Norm(direct.value.primitive.magneticFieldT),1e-30),
      "off-grid spatial inverse does not recover the material sheath state");
  const Vec3 upstreamPoint=surface.value->frontPatches.front().geometry.centerM+
      1e8*surface.value->frontPatches.front().geometry.outwardNormal;
  Require(!models.sheath->EvaluateAtPosition(upstreamPoint,40,0.1,12).ok(),
      "spatial sheath inverse silently filled an unsupported ambient point");
}

} // namespace

int main() {
  try {
    TestShockLimitAndCauchyMap();
    TestInventoryAndRefinement();
    TestMaterialDifferentialIdentities();
    TestSpatialInverse();
    std::cout<<"[CMBGU05] PASS shock-fed spatial sheath map and inventory\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-4 FAIL: "<<error.what()<<'\n';return 1;
  }
}
