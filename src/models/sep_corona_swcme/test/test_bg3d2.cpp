#include "sep_corona_swcme/ambient_state.h"
#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronalCME::Dot;
using SEP::CoronalCME::Norm;
using SEP::CoronalCME::Vec3;

constexpr double kPi=3.14159265358979323846;
constexpr double kAuM=149597870700.0;

void Require(bool condition,const std::string& message) {
  if(!condition)throw std::runtime_error(message);
}

std::string FileBytes(const std::string& path) {
  std::ifstream input(path,std::ios::binary);
  Require(static_cast<bool>(input),"cannot read ambient fixture: "+path);
  std::ostringstream bytes;
  bytes<<input.rdbuf();
  Require(input.good()||input.eof(),"failed reading ambient fixture: "+path);
  return bytes.str();
}

std::shared_ptr<const SEP::CoronaSwcme::EventConfiguration> Event() {
  const auto result=SEP::CoronaSwcme::ResolveEventConfiguration(
      FileBytes("examples/bg3d1/event.conf"),[](const std::string& path) {
        try {
          return SEP::Core::Result<std::string>::Success(FileBytes(path));
        } catch(const std::exception& error) {
          return SEP::Core::Result<std::string>::Failure(
              SEP::Core::StatusCode::DataIntegrityFailure,error.what());
        }
      });
  Require(result.ok(),result.status.message);
  return result.value;
}

Vec3 Direction(double theta,double phi) {
  return {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)};
}

void TestRegionsAndInvariants(const AmbientModel& model) {
  const double solar=model.Event().support.solarRadiusM;
  const double epoch=1000;
  const Vec3 north=Direction(kPi/4,0);
  const auto outer=model.Evaluate(kAuM*north,epoch);
  Require(outer.ok(),"1-AU ambient sample failed: "+outer.status.message);
  Require(outer.value.region==SEP::CoronaSwcme::AmbientRegion::ParkerExterior&&
      outer.value.magneticSector==1,"northern exterior region/sector is wrong");
  const double radialSpeed=Dot(outer.value.velocityMPerS,north);
  const double flux=outer.value.plasma.massDensityKgM3*radialSpeed*kAuM*kAuM;
  Require(std::abs(flux-model.MassFluxPerSteradianKgPerS())<=
      2e-12*model.MassFluxPerSteradianKgPerS(),"exterior mass flux is not invariant");
  const Vec3 azimuthScaled={-north.y,north.x,0};
  const double sineTheta=Norm(azimuthScaled);
  const Vec3 azimuth=azimuthScaled/sineTheta;
  const double radialField=Dot(outer.value.magneticFieldT,north);
  const double azimuthalField=Dot(outer.value.magneticFieldT,azimuth);
  const auto& a=model.Event().ambient;
  const double expectedRatio=-a.rotationRateRadPerS*
      (kAuM-a.sourceSurfaceRadiusM*a.sourceSurfaceRadiusM/kAuM)*sineTheta/radialSpeed;
  Require(std::abs(azimuthalField/radialField-expectedRatio)<=1e-12*
      std::max(1.0,std::abs(expectedRatio)),"Parker winding identity failed");
  Require(std::abs(Dot(outer.value.velocityMPerS,azimuth)-
      a.rotationRateRadPerS*a.sourceSurfaceRadiusM*a.sourceSurfaceRadiusM*
      sineTheta/kAuM)<1e-9,
      "angular-momentum velocity identity failed");
  Require(outer.value.plasma.pressurePa>outer.value.electronPressurePa&&
      outer.value.electronPressurePa>0&&outer.value.protonTemperatureK>0,
      "ambient EOS/temperature state is incomplete");
  const double expectedElectronPressure=outer.value.plasma.electronNumberDensityM3*
      SEP::CoronalCME::Constants::kBoltzmannJPerK*outer.value.electronTemperatureK;
  Require(std::abs(outer.value.electronPressurePa-expectedElectronPressure)<=
      2e-15*expectedElectronPressure,"independent electron EOS identity failed");

  const auto southern=model.Evaluate(kAuM*Direction(3*kPi/4,0),epoch);
  Require(southern.ok()&&southern.value.magneticSector==-1,
      "southern signed magnetic sector is wrong");
  const auto open=model.Evaluate(1.2*solar*Vec3{0,0,1},epoch);
  Require(open.ok()&&open.value.region==SEP::CoronaSwcme::AmbientRegion::PfssOpen,
      "polar low-coronal open branch missing");
  const auto closed=model.Evaluate(1.2*solar*Vec3{1,0,0},epoch);
  Require(closed.ok()&&closed.value.region==SEP::CoronaSwcme::AmbientRegion::PfssClosed&&
      closed.value.magneticSector==0,"equatorial closed-corona branch missing");
  const Vec3 expectedCorotation={-a.rotationRateRadPerS*0,
      a.rotationRateRadPerS*1.2*solar,0};
  Require(Norm(closed.value.velocityMPerS-expectedCorotation)<1e-8*
      std::max(1.0,Norm(expectedCorotation)),"closed corona is not rigidly corotating");

  Require(!model.Evaluate((solar-1)*north,epoch).ok(),
      "sample below first plasma radius accepted");
  Require(!model.Evaluate((model.Event().support.coverageRadiusM+1)*north,epoch).ok(),
      "sample beyond outer coverage accepted");
  Require(!model.Evaluate(kAuM*north,-1).ok(),"epoch extrapolation accepted");
  const auto nullSample=model.Evaluate(kAuM*Vec3{1,0,0},epoch);
  Require(!nullSample.ok()&&nullSample.status.code==SEP::Core::StatusCode::NumericalFailure,
      "dipole current-sheet null was not a typed numerical disposition");
}

void TestDerivatives(const AmbientModel& model) {
  const Vec3 point=kAuM*Direction(kPi/3,0.4);
  const auto state=model.EvaluateWithDerivatives(point,1200,7);
  Require(state.ok()&&state.value.derivativesValid&&state.value.generation==7&&
      state.value.eventIdentity==model.Event().physicsFingerprint,
      "ambient derivative identity/readiness failed");
  Require(!model.EvaluateWithDerivatives(point,1200,0).ok(),
      "reserved background generation accepted");
  const double divergence=state.value.gradientB[0]+state.value.gradientB[4]+
      state.value.gradientB[8];
  const double normalized=kAuM*std::abs(divergence)/Norm(state.value.primitive.magneticFieldT);
  Require(normalized<1e-7,"Parker exterior divergence exceeds BG3D-0 target");
  const Vec3 radial=point/Norm(point);
  Vec3 radialDerivative;
  const double er[]={radial.x,radial.y,radial.z};
  for(int component=0;component<3;++component) {
    double value=0;
    for(int coordinate=0;coordinate<3;++coordinate)
      value+=state.value.gradientU[3*component+coordinate]*er[coordinate];
    if(component==0)radialDerivative.x=value;
    else if(component==1)radialDerivative.y=value;
    else radialDerivative.z=value;
  }
  const double speed=Dot(state.value.primitive.velocityMPerS,radial);
  const double derivative=Dot(radialDerivative,radial);
  const double soundSquared=state.value.primitive.plasma.pressurePa/
      state.value.primitive.plasma.massDensityKgM3;
  const double lhs=(speed-soundSquared/speed)*derivative;
  const double rhs=2*soundSquared/Norm(point)-
      SEP::CoronalCME::Constants::kSolarGravitationalParameterM3PerS2/
      (Norm(point)*Norm(point));
  Require(std::abs(lhs-rhs)<=2e-6*std::max(std::abs(lhs),std::abs(rhs)),
      "isothermal Parker momentum identity failed");

  // A separate five-point derivative oracle does not reuse the provider's
  // three-point stencil. Compare all tensor components away from interfaces.
  const double step=2e-4*kAuM;
  for(int axis=0;axis<3;++axis) {
    Vec3 delta;
    if(axis==0)delta.x=step;else if(axis==1)delta.y=step;else delta.z=step;
    const auto m2=model.Evaluate(point-2*delta,1200),m1=model.Evaluate(point-delta,1200);
    const auto p1=model.Evaluate(point+delta,1200),p2=model.Evaluate(point+2*delta,1200);
    Require(m2.ok()&&m1.ok()&&p1.ok()&&p2.ok(),
        "five-point derivative oracle left branch: "+m2.status.message+" | "+
        m1.status.message+" | "+p1.status.message+" | "+p2.status.message);
    const Vec3 dB=(m2.value.magneticFieldT-8*m1.value.magneticFieldT+
        8*p1.value.magneticFieldT-p2.value.magneticFieldT)/(12*step);
    const Vec3 dU=(m2.value.velocityMPerS-8*m1.value.velocityMPerS+
        8*p1.value.velocityMPerS-p2.value.velocityMPerS)/(12*step);
    const double expectedB[]={dB.x,dB.y,dB.z},expectedU[]={dU.x,dU.y,dU.z};
    for(int component=0;component<3;++component) {
      const int index=3*component+axis;
      Require(std::abs(state.value.gradientB[index]-expectedB[component])<=
          2e-6*std::max(1e-30,std::abs(expectedB[component])),
          "magnetic derivative differs from five-point oracle");
      Require(std::abs(state.value.gradientU[index]-expectedU[component])<=
          2e-6*std::max(1e-20,std::abs(expectedU[component])),
          "velocity derivative differs from five-point oracle");
    }
  }
}

double FluxResidual(const AmbientModel& model,int polar,int azimuth) {
  double signedFlux=0,unsignedFlux=0;
  const double radius=1.1*kAuM;
  for(int i=0;i<polar;++i) {
    const double theta=kPi*(i+0.5)/polar;
    const double area=radius*radius*std::sin(theta)*(kPi/polar)*(2*kPi/azimuth);
    for(int j=0;j<azimuth;++j) {
      const Vec3 radial=Direction(theta,2*kPi*(j+0.5)/azimuth);
      const auto state=model.Evaluate(radius*radial,5000);
      Require(state.ok(),"sphere flux quadrature encountered invalid ambient sample");
      const double contribution=Dot(state.value.magneticFieldT,radial)*area;
      signedFlux+=contribution;
      unsignedFlux+=std::abs(contribution);
    }
  }
  return std::abs(signedFlux)/unsignedFlux;
}

void TestMagneticFluxConvergence(const AmbientModel& model) {
  const double coarse=FluxResidual(model,12,24);
  const double medium=FluxResidual(model,24,48);
  const double fine=FluxResidual(model,48,96);
  std::cout<<"[CMBGU03-EVIDENCE] flux_residuals="<<coarse<<','<<medium<<','<<fine
           <<" mass_flux_kg_s_sr="<<model.MassFluxPerSteradianKgPerS()<<'\n';
  // The axisymmetric dipole quadrature is already at roundoff on the coarsest
  // grid, so monotone convergence is neither expected nor meaningful here.
  Require(std::max({coarse,medium,fine})<1e-12,
      "signed exterior magnetic flux does not close at three resolutions");
}

} // namespace

int main() {
  try {
    const auto created=AmbientModel::Create(Event());
    Require(created.ok(),created.status.message);
    TestRegionsAndInvariants(*created.value);
    TestDerivatives(*created.value);
    TestMagneticFluxConvergence(*created.value);
    std::cout<<"[CMBGU03] PASS ambient plasma/IMF/EOS/coverage/derivatives\n";
    return 0;
  } catch(const std::exception& error) {
    std::cerr<<"BG3D-2 FAIL: "<<error.what()<<'\n';
    return 1;
  }
}
