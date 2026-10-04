#include "sep_corona_swcme/sheath_model.h"
#include "sep_corona_swcme/sheath_diagnostics.h"
#include "sep_coronal_cme/configuration_parser.h"
#include "bg3d4_reference_fixtures.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

namespace {

constexpr double kPi=3.141592653589793238462643383279502884;

using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronaSwcme::EventConfiguration;
using SEP::CoronaSwcme::ContactAuthorityQualification;
using SEP::CoronaSwcme::IntegrateSheathResiduals;
using SEP::CoronaSwcme::SheathCellDisposition;
using SEP::CoronaSwcme::SheathContactState;
using SEP::CoronaSwcme::SheathDiagnosticOptions;
using SEP::CoronaSwcme::SheathMaterialLabel;
using SEP::CoronaSwcme::ShockFedSheathModel;
using SEP::CoronaSwcme::SurfaceShockModel;
using SEP::CoronaSwcme::WeightedSheathDiagnosticPoint;
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

std::shared_ptr<const EventConfiguration> EventWithHistory(
    const std::string& history) {
  std::string configuration=FileBytes("examples/bg3d1/event.conf");
  const std::string original=FileBytes("examples/bg3d1/history.csv");
  const std::string oldHash=SEP::CoronalCME::ComputeContentChecksum(original);
  const std::string newHash=SEP::CoronalCME::ComputeContentChecksum(history);
  const std::size_t digest=configuration.find(oldHash);
  Require(digest!=std::string::npos,"history asset digest is absent from fixture");
  configuration.replace(digest,oldHash.size(),newHash);
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(configuration,
      [&](const std::string& path) {
        try {
          return SEP::Core::Result<std::string>::Success(
              path=="examples/bg3d1/history.csv"?history:FileBytes(path));
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

Vec3 SpatialScalarGradient(const std::array<Vec3,3>& positions,
    const std::array<double,3>& values) {
  const double determinant=Determinant(positions);
  Require(std::abs(determinant)>1e-30,"material-coordinate scalar gradient is singular");
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
      Close(boundary.value.referenceVolumeDensityM3PerRad2S,
          boundary.value.currentVolumeDensityM3PerRad2S,1e-14)&&
      Close(Determinant(boundary.value.referenceDerivativeColumns),
          boundary.value.referenceVolumeDensityM3PerRad2S,1e-14)&&
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

void TestIndependentPlanarAndCurvedReferences() {
  const BG3D4Reference::FinitePlanarShockPiston planar;
  const double w1=planar.UpstreamInflow(),w2=planar.DownstreamInflow();
  const double mass1=planar.upstream.massDensityKgM3*w1;
  const double mass2=planar.downstream.massDensityKgM3*w2;
  const double momentum1=planar.upstream.massDensityKgM3*w1*w1+
      planar.upstream.pressurePa;
  const double momentum2=planar.downstream.massDensityKgM3*w2*w2+
      planar.downstream.pressurePa;
  const auto shockEnergyFlux=[&](const SEP::CoronalCME::MhdPrimitiveState& state,
      double inflow) {
    const double magnetic=Dot(state.magneticFieldT,state.magneticFieldT)/
        (2*BG3D4Reference::kMu0);
    const double energy=0.5*state.massDensityKgM3*inflow*inflow+
        state.pressurePa/(planar.gamma-1)+magnetic;
    return (energy+state.pressurePa+magnetic)*inflow-
        inflow*state.magneticFieldT.x*state.magneticFieldT.x/
            BG3D4Reference::kMu0;
  };
  const double energy1=shockEnergyFlux(planar.upstream,w1);
  const double energy2=shockEnergyFlux(planar.downstream,w2);
  const double time=7.0;
  Require(Close(mass1,mass2,1e-14)&&Close(momentum1,momentum2,1e-14)&&
      Close(energy1,energy2,1e-14),
      "finite planar reference does not satisfy its exact MHD shock jump");
  Require(Close(planar.TotalMassKg(time),
      planar.InitialMassKg()+planar.AdmittedMassKg(time),1e-14)&&
      Close(planar.InitialMap(-planar.initialThicknessM,time),
          planar.ContactPosition(time),1e-14)&&
      Close(planar.AdmittedMap(time,time),planar.ShockPosition(time),1e-14)&&
      Close(planar.EnergyRateW(),planar.MovingBoundaryEnergyPowerW(),1e-14),
      "finite planar inventory/map or independently evaluated energy budget failed");
  Require(Dot(planar.downstream.velocityMPerS-planar.pistonSide.velocityMPerS,
      Vec3{-1,0,0})==0&&planar.downstream.pressurePa==planar.pistonSide.pressurePa&&
      Norm(planar.downstream.magneticFieldT-planar.pistonSide.magneticFieldT)==0,
      "planar transverse-field contact boundary states are incompatible");

  const BG3D4Reference::HomologousCurvedShell curved;
  const double curvedTime=5.0,lambda=curved.Lambda(curvedTime);
  const double exact=curved.VolumeM3(curvedTime);
  const double unweighted=curved.UnweightedInnerAreaColumnM3(curvedTime);
  Require(std::abs(unweighted-exact)>0.1*exact,
      "curved-volume negative control cannot distinguish an unweighted column");
  std::array<double,3> volumeErrors{};
  int level=0;
  for(const int radialCells:{8,16,32}) {
    const int polarCells=2*radialCells,azimuthCells=4*radialCells;
    const double ds=(curved.outerRadiusM-curved.innerRadiusM)/radialCells;
    const double dtheta=kPi/polarCells,dphi=2*kPi/azimuthCells;
    double volume=0;
    for(int ir=0;ir<radialCells;++ir)for(int it=0;it<polarCells;++it)
      for(int ip=0;ip<azimuthCells;++ip) {
        const double s=curved.innerRadiusM+(ir+0.5)*ds;
        const double theta=(it+0.5)*dtheta;
        volume+=lambda*lambda*lambda*s*s*std::sin(theta)*ds*dtheta*dphi;
      }
    volumeErrors[level++]=std::abs(volume-exact)/exact;
  }
  const double order01=std::log(volumeErrors[0]/volumeErrors[1])/std::log(2.0);
  const double order12=std::log(volumeErrors[1]/volumeErrors[2])/std::log(2.0);
  const double jacobian=curved.Jacobian(curvedTime);
  const Vec3 field=curved.MagneticField(curvedTime);
  Require(Close(curved.Density(curvedTime)*jacobian,curved.rho0KgM3,2e-14)&&
      Close(curved.Pressure(curvedTime)*std::pow(jacobian,curved.gamma),
          curved.pressure0Pa,2e-14)&&
      Norm(field-curved.magnetic0T/(lambda*lambda))<=1e-14&&
      volumeErrors[2]<=1e-3&&order01>=1.8&&order12>=1.8,
      "curved reference violates its map/state identity or frozen quadrature criteria");
  std::array<double,3> energyErrors{};level=0;
  for(const double h:{0.4,0.2,0.1}) {
    const double derivative=Fourth(curved.TotalEnergyJ(curvedTime-2*h),
        curved.TotalEnergyJ(curvedTime-h),curved.TotalEnergyJ(curvedTime+h),
        curved.TotalEnergyJ(curvedTime+2*h),h);
    energyErrors[level++]=std::abs(derivative-curved.BoundaryWorkW(curvedTime))/
        std::max(std::abs(curved.BoundaryWorkW(curvedTime)),1e-30);
  }
  Require(energyErrors[2]<=1e-7&&energyErrors[1]<energyErrors[0]&&
      energyErrors[2]<energyErrors[1],
      "curved reference energy/boundary-work identity does not converge");
  std::cout<<"[CMBGU05-PLANAR-REFERENCE] initial_mass_kg="
           <<planar.InitialMassKg()<<" admitted_mass_kg="<<planar.AdmittedMassKg(time)
           <<" energy_rate_w="<<planar.EnergyRateW()<<'\n';
  std::cout<<"[CMBGU05-CURVED-REFERENCE] exact_volume_m3="<<exact
           <<" unweighted_volume_m3="<<unweighted<<" volume_errors="
           <<volumeErrors[0]<<','<<volumeErrors[1]<<','<<volumeErrors[2]
           <<" orders="<<order01<<','<<order12<<" energy_errors="
           <<energyErrors[0]<<','<<energyErrors[1]<<','<<energyErrors[2]<<'\n';
}

void TestCrossStageContactAuthorityRegression() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  const double time=20;
  const auto surface=models.shock->Prepare(time,34,12,24);
  Require(surface.ok(),surface.status.message);
  Require(surface.value->contactAuthority==
      ContactAuthorityQualification::FixedFractionReferenceUnqualified,
      "fixed-fraction BG3D-3 surface was incorrectly advertised as a shared contact");
  double maximumRelativeMismatch=0;std::size_t compared=0;
  for(const auto& patch:surface.value->contactPatches) {
    const auto material=models.sheath->EvaluateContact(
        patch.polarParameterRad,patch.azimuthParameterRad,time);
    if(!material.ok())continue;
    const Vec3 fixed=surface.value->contact.Point(patch.polarParameterRad,
        patch.azimuthParameterRad);
    maximumRelativeMismatch=std::max(maximumRelativeMismatch,
        Norm(fixed-material.value.positionM)/
            std::max({Norm(fixed),Norm(material.value.positionM),1.0}));
    Require(material.value.eventIdentity==surface.value->eventIdentity,
        "contact candidates do not share the event authority");
    ++compared;
  }
  Require(compared>0&&maximumRelativeMismatch>1e-10,
      "cross-stage regression no longer detects the unqualified contact mismatch");
  std::cout<<"[CMBGU05-CONTACT-AUTHORITY] status="
           <<SEP::CoronaSwcme::Name(surface.value->contactAuthority)
           <<" compared="<<compared
           <<" max_relative_position_mismatch="<<maximumRelativeMismatch<<'\n';
}

void TestMaterialContactAndStartup() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  const auto startup=models.sheath->PrepareInventory(0,35,4,8,4);
  Require(startup.ok()&&startup.value->zeroVolumeStartup&&startup.value->cells.empty()&&
      startup.value->admittedMassKg==0&&startup.value->retainedMassKg==0&&
      startup.value->minimumJacobian==1,
      "zero-volume startup attempted to invert or floor a degenerate sheath");
  Require(!models.sheath->EvaluateAtPosition({1,2,3},0,1).ok(),
      "zero-volume startup exposed a spatial sheath state");

  // Locate a globally supported contact chart, then differentiate its physical
  // position independently.  The finite-difference boundary velocity is not
  // taken from the returned contact velocity, so the leakage check is not a
  // self-comparison.
  const auto surface=models.shock->Prepare(20,36,12,24);
  Require(surface.ok(),surface.status.message);
  SheathContactState contact;bool found=false;
  for(const auto& patch:surface.value->frontPatches) {
    const auto candidate=models.sheath->EvaluateContact(
        patch.geometry.polarParameterRad,patch.geometry.azimuthParameterRad,20);
    if(candidate.ok()){contact=candidate.value;found=true;break;}
  }
  Require(found,"no globally initialized material-contact chart is supported");
  std::array<double,3> leakage{};int level=0;
  for(const double h:{0.1,0.05,0.025}) {
    const auto before2=models.sheath->EvaluateContact(contact.label.polarRad,
        contact.label.azimuthRad,20-2*h);
    const auto before=models.sheath->EvaluateContact(contact.label.polarRad,
        contact.label.azimuthRad,20-h);
    const auto after=models.sheath->EvaluateContact(contact.label.polarRad,
        contact.label.azimuthRad,20+h);
    const auto after2=models.sheath->EvaluateContact(contact.label.polarRad,
        contact.label.azimuthRad,20+2*h);
    Require(before2.ok()&&before.ok()&&after.ok()&&after2.ok(),
        "contact finite-difference stencil is unsupported");
    const Vec3 boundaryVelocity=Fourth(before2.value.positionM,before.value.positionM,
        after.value.positionM,after2.value.positionM,h);
    leakage[level++]=std::abs(contact.primitive.massDensityKgM3*
        Dot(contact.primitive.velocityMPerS-boundaryVelocity,contact.outwardNormal));
  }
  const double bound=event->regional.contactFluxAbsoluteToleranceKgM2S+
      event->regional.contactFluxRelativeTolerance*
          event->regional.contactFluxReferenceKgM2S;
  std::cout<<"[CMBGU05-CONTACT] flux_bound="<<bound
           <<" fd_fluxes="<<leakage[0]<<','<<leakage[1]<<','<<leakage[2]<<'\n';
  Require(contact.relativeMassFluxKgM2S==0&&
      leakage[1]<=bound&&leakage[2]<=bound&&leakage[1]<leakage[0],
      "material-contact pointwise dimensional leakage exceeds its frozen bound");
  const double magneticPressure=Dot(contact.primitive.magneticFieldT,
      contact.primitive.magneticFieldT)/(2*1.25663706212e-6);
  const Vec3 expectedTraction=(contact.primitive.pressurePa+magneticPressure)*
      contact.outwardNormal-(contact.normalMagneticFieldT/1.25663706212e-6)*
          contact.primitive.magneticFieldT;
  Require(Norm(expectedTraction-contact.tractionPa)<=
      1e-12*std::max(Norm(expectedTraction),1e-30),
      "contact output omitted the tensor magnetic traction");

  // Absolute leakage cannot cancel between opposite signs.  Integrate an
  // independently differentiated boundary velocity over the supported
  // contact charts and time; compare against the same area-time measure times
  // the frozen dimensional-plus-relative flux bound.  Unsupported charts are
  // counted rather than silently assigned zero state.
  std::array<double,3> absoluteLeakageKg{},reportedLeakageKg{},allowedLeakageKg{};
  std::array<double,3> derivativeUncertaintyKg{},maximumVelocityError{};
  std::array<int,3> supported{},unsupported{};
  bool anyIndependentLevelUnqualified=false;
  for(int levelIndex=0;levelIndex<3;++levelIndex) {
    const int polarCells=2+levelIndex,azimuthCells=2*polarCells;
    const int timeCells=2+levelIndex;
    const double dTheta=kPi/polarCells,dPhi=2*kPi/azimuthCells;
    const double dt=40.0/timeCells;
    for(int it=0;it<timeCells;++it) {
      const double time=(it+0.5)*dt;
      for(int i=0;i<polarCells;++i)for(int j=0;j<azimuthCells;++j) {
        const double theta=(i+0.5)*dTheta,phi=(j+0.5)*dPhi;
        const auto center=models.sheath->EvaluateContact(theta,phi,time);
        if(!center.ok()){++unsupported[levelIndex];continue;}
        const auto derivative=[&](double h) {
          const auto m2=models.sheath->EvaluateContact(theta,phi,time-2*h);
          const auto m1=models.sheath->EvaluateContact(theta,phi,time-h);
          const auto p1=models.sheath->EvaluateContact(theta,phi,time+h);
          const auto p2=models.sheath->EvaluateContact(theta,phi,time+2*h);
          Require(m2.ok()&&m1.ok()&&p1.ok()&&p2.ok(),
            "absolute contact-leakage stencil left supported contact history");
          return Fourth(m2.value.positionM,m1.value.positionM,
              p1.value.positionM,p2.value.positionM,h);
        };
        // The position is O(1e9 m), so extrapolating already small fourth-order
        // differences eventually amplifies roundoff.  The independently
        // established h=0.05 s stencil is used directly; its difference from
        // h=0.1 s is retained as numerical uncertainty rather than folded into
        // the physical leakage.
        const Vec3 coarseVelocity=derivative(0.1);
        const Vec3 fineVelocity=derivative(0.05);
        const Vec3 vc=fineVelocity;
        const double derivativeFlux=std::abs(center.value.primitive.massDensityKgM3*
            Dot(center.value.primitive.velocityMPerS-vc,
                center.value.outwardNormal));
        const double reportedFlux=std::abs(center.value.relativeMassFluxKgM2S);
        const double areaTime=center.value.areaDensityM2PerRad2*dTheta*dPhi*dt;
        // Qualification uses the independently differentiated boundary speed.
        // The API-reported value is retained separately because it is zero by
        // construction when contact velocity is copied from material velocity.
        absoluteLeakageKg[levelIndex]+=derivativeFlux*areaTime;
        reportedLeakageKg[levelIndex]+=reportedFlux*areaTime;
        derivativeUncertaintyKg[levelIndex]+=
            center.value.primitive.massDensityKgM3*
            std::abs(Dot(coarseVelocity-fineVelocity,
                center.value.outwardNormal))*areaTime;
        maximumVelocityError[levelIndex]=std::max(maximumVelocityError[levelIndex],
            std::abs(Dot(center.value.velocityMPerS-vc,
                center.value.outwardNormal))/
                std::max(std::abs(center.value.normalSpeedMPerS),1.0));
        allowedLeakageKg[levelIndex]+=bound*areaTime;
        ++supported[levelIndex];
      }
    }
    const bool independentlyQualified=
        absoluteLeakageKg[levelIndex]<=allowedLeakageKg[levelIndex]&&
        maximumVelocityError[levelIndex]<=
            event->regional.contactNormalVelocityNumericalTolerance;
    anyIndependentLevelUnqualified=anyIndependentLevelUnqualified||
        !independentlyQualified;
    std::cout<<"[CMBGU05-CONTACT-INTEGRAL-L"<<levelIndex<<"] absolute_kg="
             <<absoluteLeakageKg[levelIndex]<<" allowed_kg="
             <<allowedLeakageKg[levelIndex]<<" supported="<<supported[levelIndex]
             <<" unsupported="<<unsupported[levelIndex]
             <<" reported_kg="<<reportedLeakageKg[levelIndex]
             <<" derivative_uncertainty_kg="<<derivativeUncertaintyKg[levelIndex]
             <<" max_velocity_error="<<maximumVelocityError[levelIndex]
             <<" independently_qualified="
             <<(independentlyQualified?"true":"false")<<'\n';
    Require(supported[levelIndex]>0&&unsupported[levelIndex]>0&&
        std::isfinite(absoluteLeakageKg[levelIndex])&&
        std::isfinite(derivativeUncertaintyKg[levelIndex]),
        "absolute material-contact integration did not cover its typed domain");
  }
  std::cout<<"[CMBGU05-CONTACT-INTEGRAL] absolute_kg="
           <<absoluteLeakageKg[0]<<','<<absoluteLeakageKg[1]<<','
           <<absoluteLeakageKg[2]<<" allowed_kg="<<allowedLeakageKg[0]<<','
           <<allowedLeakageKg[1]<<','<<allowedLeakageKg[2]
           <<" reported_kg="<<reportedLeakageKg[0]<<','<<reportedLeakageKg[1]
           <<','<<reportedLeakageKg[2]
           <<" derivative_uncertainty_kg="<<derivativeUncertaintyKg[0]<<','
           <<derivativeUncertaintyKg[1]<<','<<derivativeUncertaintyKg[2]
           <<" supported="<<supported[0]<<','<<supported[1]<<','<<supported[2]
           <<" unsupported="<<unsupported[0]<<','<<unsupported[1]<<','
           <<unsupported[2]<<'\n';
  Require(anyIndependentLevelUnqualified,
      "experimental contact was unexpectedly qualified by its distributed derivative");
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
    const double closed=inventory->initialMassKg+inventory->retainedMassKg+
        inventory->solarExitMassKg+inventory->outerExitMassKg;
    Require(std::abs(closed-inventory->admittedMassKg)<=
        event->regional.inventoryMassAbsoluteToleranceKg+
        1e-12*std::max(closed,inventory->admittedMassKg),
        "initial/admission/solar/outer inventory ledger does not close");
  }
  const auto committed=models.sheath->Current();
  const auto failed=models.sheath->PrepareInventory(epoch,0,4,8,8);
  Require(!failed.ok()&&models.sheath->Current()==committed,
      "invalid inventory changed the committed sheath epoch");
  // The relaxing map must remain regular through the transition without the
  // old contact-exit deletion.  Every admitted cell is retained unless it
  // reaches a genuine solar or outer-domain boundary.
  const auto longCandidate=models.sheath->PrepareInventory(200,45,4,8,16);
  Require(longCandidate.ok(),"relaxing map failed through handoff: "+
      longCandidate.status.message);
  Require(longCandidate.value->minimumJacobian>=event->regional.minimumJacobian&&
      longCandidate.value->retainedMassKg>0&&
      longCandidate.value->maximumInventoryMassResidual<=1e-12,
      "handoff-epoch relaxing map failed its material inventory gates");
  const auto outerCandidate=models.sheath->PrepareInventory(200000,47,2,4,8);
  Require(outerCandidate.ok(),"relaxing map failed at the configured >1-AU endpoint: "+
      outerCandidate.status.message);
  Require(outerCandidate.value->minimumJacobian>=event->regional.minimumJacobian&&
      outerCandidate.value->retainedMassKg>0,
      "outer relaxing map lacks positive-J retained material");
  const auto incompatibleEvent=EventWithContactFraction("0.82");
  auto incompatible=CreateModels(incompatibleEvent);
  const auto fractionIndependent=incompatible.sheath->PrepareInventory(epoch,46,4,8,16);
  Require(fractionIndependent.ok()&&Close(fractionIndependent.value->admittedMassKg,
      reference.value->admittedMassKg,1e-13),
      "BG3D-4 still uses the BG3D-5 ejecta fraction as a sheath contact");
  std::cout<<"[CMBGU05-EVIDENCE] mass_kg="<<reference.value->admittedMassKg
           <<" time_errors="<<e0<<','<<e1<<','<<e2
           <<" min_j="<<reference.value->minimumJacobian
           <<" admission_res="<<reference.value->maximumAdmissionMassResidual
           <<" inventory_res="<<reference.value->maximumInventoryMassResidual
           <<" retained="<<reference.value->retainedMassKg
           <<" solar_exit="<<reference.value->solarExitMassKg
           <<" outer_exit="<<reference.value->outerExitMassKg
           <<" handoff_min_j="<<longCandidate.value->minimumJacobian
           <<" outer_min_j="<<outerCandidate.value->minimumJacobian
           <<" outer_mass_kg="<<outerCandidate.value->admittedMassKg<<'\n';
}

void TestRejectedSubfastEpochRetainsCommittedMaterial() {
  // This analytical history starts with an admissible fast launch and slows
  // continuously to a sub-fast front before handoff.  It is a failure-policy
  // fixture, not a compression-wave solution: the unsupported candidate must
  // be rejected without changing the last complete material epoch.
  const std::string history=
      "time_s,center_m,center_rate_m_s,radial_m,radial_rate_m_s,"
      "lateral1_m,lateral1_rate_m_s,lateral2_m,lateral2_rate_m_s\n"
      "0,800000000,500000,200000000,300000,300000000,200000,250000000,180000\n"
      "100,825500000,10000,215500000,10000,310250000,5000,259250000,5000\n"
      "200,826500000,10000,216500000,10000,310750000,5000,259750000,5000\n";
  const auto event=EventWithHistory(history);auto models=CreateModels(event);
  const auto committedResult=models.sheath->PrepareInventory(20,71,4,8,8);
  Require(committedResult.ok(),"fast-launch inventory failed: "+
      committedResult.status.message);
  const auto committed=models.sheath->Current();
  const auto retained=std::find_if(committed->cells.begin(),committed->cells.end(),
      [](const auto& cell) {return cell.disposition==SheathCellDisposition::Retained;});
  Require(retained!=committed->cells.end(),
      "fast-launch fixture has no committed retained material");
  const auto before=models.sheath->QueryCommitted(retained->label);
  Require(before.ok(),before.status.message);

  double rejectedEpoch=0;SEP::Core::Status rejection;
  for(const double candidateEpoch:{60.0,80.0,95.0}) {
    const auto candidate=models.sheath->PrepareInventory(candidateEpoch,72,4,8,8);
    if(!candidate.ok()&&candidate.status.message.find("sub-fast")!=std::string::npos) {
      rejectedEpoch=candidateEpoch;rejection=candidate.status;break;
    }
  }
  Require(rejectedEpoch>0&&models.sheath->Current()==committed,
      "sub-fast candidate was not rejected transactionally");
  const auto after=models.sheath->QueryCommitted(retained->label);
  Require(after.ok()&&Norm(after.value.positionM-before.value.positionM)==0&&
      after.value.jacobian==before.value.jacobian&&
      after.value.eventIdentity==before.value.eventIdentity,
      "rejected current shock admission erased or recomputed committed material");
  std::cout<<"[CMBGU05-SUBFAST-TRANSACTION] rejected_epoch_s="<<rejectedEpoch
           <<" retained_epoch_s="<<committed->epochS
           <<" retained_cells="<<committed->cells.size()
           <<" status="<<rejection.message<<'\n';
}

void TestDistributedProductionResiduals() {
  const auto event=FrozenEvent();auto models=CreateModels(event);
  const double crossing=10000,epoch=10001;
  const auto surface=models.shock->Prepare(crossing,81,12,24);
  Require(surface.ok(),surface.status.message);
  std::vector<WeightedSheathDiagnosticPoint> points;
  const double angleStep=2e-4,labelStep=0.25;
  const double dtheta=kPi/12,dphi=2*kPi/24;
  for(const auto& patch:surface.value->frontPatches) {
    if(points.size()>=8)break;
    const auto shock=std::find_if(surface.value->shocks->patches.begin(),
        surface.value->shocks->patches.end(),[&](const auto& value) {
          return value.stableId==patch.geometry.physicalId;
        });
    if(shock==surface.value->shocks->patches.end()||!shock->fast||
        patch.geometry.polarParameterRad<=3*angleStep||
        patch.geometry.polarParameterRad>=kPi-3*angleStep||
        patch.geometry.azimuthParameterRad<=3*angleStep||
        patch.geometry.azimuthParameterRad>=2*kPi-3*angleStep)continue;
    SheathMaterialLabel label{patch.geometry.polarParameterRad,
        patch.geometry.azimuthParameterRad,crossing,patch.geometry.physicalId};
    const auto state=models.sheath->Evaluate(label,epoch);
    if(!state.ok())continue;
    points.push_back({label,state.value.currentVolumeDensityM3PerRad2S*
        dtheta*dphi*(2*labelStep)});
  }
  Require(points.size()>=4,
      "production map has too few distributed smooth residual points");
  std::array<double,3> integratedForce{},p99{},energy{},work{};
  std::array<std::size_t,3> supported{},unsupported{};
  int level=0;
  for(const double evolutionStep:{0.2,0.1,0.05}) {
    SheathDiagnosticOptions options;
    options.angularStepRad=angleStep;
    options.admissionTimeStepS=labelStep;
    options.evolutionTimeStepS=evolutionStep;
    const auto report=IntegrateSheathResiduals(*models.sheath,points,epoch,
        event->composition.gammaAdiabatic,options);
    Require(report.ok(),report.status.message);
    Require(report.value.samples.size()>=4&&report.value.sampledVolumeM3>0&&
        std::isfinite(report.value.integratedMomentumRatio)&&
        std::isfinite(report.value.volumeWeightedLocalMomentumP99)&&
        std::isfinite(report.value.absoluteEnergyResidualW)&&
        std::isfinite(report.value.forceWorkRatio),
        "distributed production force/energy report is incomplete or nonfinite");
    integratedForce[level]=report.value.integratedMomentumRatio;
    p99[level]=report.value.volumeWeightedLocalMomentumP99;
    energy[level]=report.value.absoluteEnergyResidualW;
    work[level]=report.value.forceWorkRatio;
    supported[level]=report.value.samples.size();
    unsupported[level]=report.value.unsupportedSamples;
    std::cout<<"[CMBGU05-DISTRIBUTED-L"<<level<<"] samples="
             <<report.value.samples.size()<<" unsupported="
             <<report.value.unsupportedSamples<<" volume_m3="
             <<report.value.sampledVolumeM3<<" integrated_force_ratio="
             <<report.value.integratedMomentumRatio<<" local_p99="
             <<report.value.volumeWeightedLocalMomentumP99
             <<" signed_energy_residual_w="<<report.value.signedEnergyResidualW
             <<" absolute_energy_residual_w="<<report.value.absoluteEnergyResidualW
             <<" signed_force_work_w="<<report.value.signedResidualForceWorkW
             <<" absolute_force_work_w="<<report.value.absoluteResidualForceWorkW
             <<" inertia_work_w="<<report.value.signedInertiaWorkW
             <<" pressure_gradient_work_w="
             <<report.value.signedPressureGradientWorkW
             <<" lorentz_work_w="<<report.value.signedLorentzWorkW
             <<" gravity_work_w="<<report.value.signedGravityWorkW
             <<" force_work_ratio="<<report.value.forceWorkRatio<<'\n';
    ++level;
  }
  // These are convergence observations, not qualification thresholds.  The
  // empirical map remains unqualified even if a broad configured force cap is
  // met.  Requiring only finite, nonzero distributed measures prevents the old
  // single-point diagnostic from being silently substituted for this budget.
  Require(integratedForce[2]>0&&p99[2]>0&&energy[2]>0&&work[2]>0&&
      supported[0]==supported[1]&&supported[1]==supported[2],
      "production residual distribution did not expose a physical model residual");
  std::cout<<"[CMBGU05-DISTRIBUTED-REFINEMENT] integrated_force_ratio="
           <<integratedForce[0]<<','<<integratedForce[1]<<','<<integratedForce[2]
           <<" local_p99="<<p99[0]<<','<<p99[1]<<','<<p99[2]
           <<" absolute_energy_w="<<energy[0]<<','<<energy[1]<<','<<energy[2]
           <<" work_ratio="<<work[0]<<','<<work[1]<<','<<work[2]<<'\n';
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
  const std::array<Vec3,3> db={
      Fourth(thetaMinus2.primitive.magneticFieldT,thetaMinus.primitive.magneticFieldT,
          thetaPlus.primitive.magneticFieldT,thetaPlus2.primitive.magneticFieldT,angleStep),
      Fourth(phiMinus2.primitive.magneticFieldT,phiMinus.primitive.magneticFieldT,
          phiPlus.primitive.magneticFieldT,phiPlus2.primitive.magneticFieldT,angleStep),
      Fourth(tauMinus2.primitive.magneticFieldT,tauMinus.primitive.magneticFieldT,
          tauPlus.primitive.magneticFieldT,tauPlus2.primitive.magneticFieldT,timeLabelStep)};
  const std::array<double,3> dp={
      Fourth(thetaMinus2.primitive.pressurePa,thetaMinus.primitive.pressurePa,
          thetaPlus.primitive.pressurePa,thetaPlus2.primitive.pressurePa,angleStep),
      Fourth(phiMinus2.primitive.pressurePa,phiMinus.primitive.pressurePa,
          phiPlus.primitive.pressurePa,phiPlus2.primitive.pressurePa,angleStep),
      Fourth(tauMinus2.primitive.pressurePa,tauMinus.primitive.pressurePa,
          tauPlus.primitive.pressurePa,tauPlus2.primitive.pressurePa,timeLabelStep)};
  const auto gradientU=SpatialGradient(dx,du);
  const auto gradientB=SpatialGradient(dx,db);
  const Vec3 gradientP=SpatialScalarGradient(dx,dp);
  const double divergenceU=gradientU[0].x+gradientU[1].y+gradientU[2].z;
  const double divergenceB=gradientB[0].x+gradientB[1].y+gradientB[2].z;
  const double evolutionStep=0.05;
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
  const Vec3 acceleration=Fourth(before2.primitive.velocityMPerS,
      before.primitive.velocityMPerS,after.primitive.velocityMPerS,
      after2.primitive.velocityMPerS,evolutionStep);
  const Vec3 curlB{gradientB[1].z-gradientB[2].y,
      gradientB[2].x-gradientB[0].z,gradientB[0].y-gradientB[1].x};
  const Vec3 gravity=(-1.32712440018e20/std::pow(Norm(central.positionM),3))*
      central.positionM;
  const Vec3 inertia=central.primitive.massDensityKgM3*acceleration;
  const Vec3 lorentz=(1/1.25663706212e-6)*Cross(curlB,
      central.primitive.magneticFieldT);
  const Vec3 gravityForce=central.primitive.massDensityKgM3*gravity;
  const Vec3 momentumResidual=inertia+gradientP-lorentz-gravityForce;
  const double momentumScale=Norm(inertia)+Norm(gradientP)+Norm(lorentz)+
      Norm(gravityForce);
  const double momentumRatio=Norm(momentumResidual)/std::max(momentumScale,1e-300);
  const double divergenceScale=Norm(central.primitive.magneticFieldT)/
      std::max(Norm(central.positionM),1.0);
  const double divergenceResidual=std::abs(divergenceB)/
      std::max(divergenceScale,1e-300);
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
           <<" thermal="<<thermalResidual<<" induction="<<inductionResidual
           <<" div_b="<<divergenceResidual
           <<" momentum_pa_m="<<Norm(momentumResidual)
           <<" momentum_ratio="<<momentumRatio<<'\n';
  Require(massResidual<=1e-7&&thermalResidual<=1e-7&&inductionResidual<=1e-6,
      "independent Lagrangian mass/adiabatic/induction residual exceeds budget");
  Require(divergenceResidual<=1e-6,
      "shock-born magnetic reference is not spatially flux compatible");
  Require(momentumRatio<=event->regional.maximumLocalForceRatioP99,
      "prescribed sheath sustaining-force ratio exceeds the frozen local budget");

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
    TestIndependentPlanarAndCurvedReferences();
    TestCrossStageContactAuthorityRegression();
    TestMaterialContactAndStartup();
    TestInventoryAndRefinement();
    TestRejectedSubfastEpochRetainsCommittedMaterial();
    TestMaterialDifferentialIdentities();
    TestDistributedProductionResiduals();
    TestSpatialInverse();
    std::cout<<"[CMBGU05] PASS shock-fed spatial sheath map and inventory\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-4 FAIL: "<<error.what()<<'\n';return 1;
  }
}
