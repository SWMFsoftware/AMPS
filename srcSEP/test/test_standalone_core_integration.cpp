#include "reduced_shock_background_adapter.h"
#include "provider.h"
#include "sep_reduced_geometry.h"
#include "sep_background_snapshot.h"
#include "sep_parker_core.h"
#include "sep_focused_transport_core.h"
#include "sep_focused_transport_mfp_core.h"
#include "sep_reduced_source.h"
#include "sep_absorbing_bridge.h"
#include <sep_coronal_cme/particle_source.h>
#include <sep_coronal_cme/mhd_jump_solver.h>

#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>

namespace RS=SEP::ReducedShock;
namespace SF=SEP::CoronaSwcme::ShockFront;
namespace T=SEP::Transport;
namespace {
struct ZeroRandom final:T::RandomStream {
  unsigned draws=0;
  double UniformOpen01() override {++draws;return 0.5;}
  double Normal01() override {++draws;return 0;}
};
struct Spatial final:T::SpatialDiffusionProvider {
  double kappa=4,derivative=0;
  T::SpatialDiffusionSample Evaluate(double,double) const override {
    T::SpatialDiffusionSample out;out.status=T::Status::Ok();out.kappaParallelM2PerS=kappa;
    out.dKappaParallelDsMPerS=derivative;out.provenance="manufactured-SI-diffusion";return out;
  }
};
struct Pitch final:T::PitchAngleDiffusionProvider {
  double rate=0;
  T::PitchAngleDiffusionSample Evaluate(double,double,double mu) const override {
    T::PitchAngleDiffusionSample out;out.status=T::Status::Ok();
    out.dMuMuPerS=rate*(1-mu*mu)/2;out.dDmuMuDmuPerS=-rate*mu;
    out.valueState=T::CoefficientPhysics::ValueState::Finite;
    out.provenance="manufactured-isotropic-pitch-diffusion";out.turbulenceStateIdentity="manufactured-fixed-frame";return out;
  }
};
struct Mfp final:T::MeanFreePathProvider {
  double rate=4,speed=10;
  T::MeanFreePathSample Evaluate(double,double,double) const override {
    T::MeanFreePathSample out;out.status=T::Status::Ok();out.lambdaParallelM=speed/rate;
    out.valueState=T::CoefficientPhysics::ValueState::Finite;
    out.nuPlusPerS=out.nuMinusPerS=rate/2;out.hasBranchResolvedRates=true;
    out.provenance="manufactured-isotropic-reset";out.turbulenceStateIdentity="manufactured-fixed-frame";return out;
  }
};
void Require(bool condition,const std::string& why) {
  if(!condition)throw std::runtime_error(why);
}
unsigned referenceChecks=0;
double maximumReferenceBoundFraction=0;
void Near(double actual,double expected,double relative=2e-12) {
  ++referenceChecks;
  if(std::isfinite(actual)&&std::isfinite(expected))maximumReferenceBoundFraction=std::max(
      maximumReferenceBoundFraction,std::abs(actual-expected)/(relative*std::max(1e-30,std::abs(expected))));
  Require(std::isfinite(actual)&&std::isfinite(expected)&&
      std::abs(actual-expected)<=relative*std::max(1e-30,std::abs(expected)),
      "independent/reference value mismatch: "+std::to_string(actual)+" versus "+std::to_string(expected));
}
std::shared_ptr<SF::Provider> Direct(const std::string& path) {
  const auto read=[](const std::string& file) {
    std::ifstream input(file);return std::string(std::istreambuf_iterator<char>(input),{});
  };
  const auto slash=path.find_last_of('/');
  const std::string directory=slash==std::string::npos?".":path.substr(0,slash);
  const auto configuration=SF::ResolveConfiguration(read(path),[&](const std::string& name) {
    return SEP::Core::Result<std::string>::Success(read(directory+"/"+name));
  });
  Require(configuration.ok(),configuration.status.message);
  const auto provider=SF::Provider::Create(configuration.value);
  Require(provider.ok(),provider.status.message);return provider.value;
}
}
int main(int argc,char** argv) {
  if(argc!=3){std::cerr<<"usage: test_standalone_core_integration EVENT ALPHA_EVENT\n";return 2;}
  int failures=0;
  const auto run=[&](const char* id,const std::function<void()>& test) {
    const auto start=std::chrono::steady_clock::now();referenceChecks=0;maximumReferenceBoundFraction=0;
    try {test();std::cout<<'['<<id<<"] PASS SI/HCI core integration contract\n";}
    catch(const std::exception& error){++failures;std::cout<<'['<<id<<"] FAIL "<<error.what()<<'\n';}
    std::cout<<std::setprecision(17)<<"METRIC "<<id<<" duration_s="
        <<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
        <<" reference_checks="<<referenceChecks<<" maximum_reference_bound_fraction="<<maximumReferenceBoundFraction<<'\n';
  };
  std::string error;
  run("SINT01",[&] {
    Require(RS::Configure(argv[1],&error),error);Require(RS::Prepare(0,&error),error);
    const auto view=RS::CaptureInterval(0,600,RS::EventIdentity());
    Require(view.ok(),view.status.message);
    const RS::Position point{{4e10,1e10,3e10}};
    const auto actual=view.value.AmbientAt(point,30);
    const auto expected=Direct(argv[1])->QueryAmbient({point[0],point[1],point[2]},30,1);
    Require(actual.ok(),actual.status.message);Require(expected.ok(),expected.status.message);
    Near(actual.value.magneticFieldT[0],expected.value.primitive.magneticFieldT.x);
    Near(actual.value.velocityMPerS[1],expected.value.primitive.velocityMPerS.y);
    Near(actual.value.massDensityKgM3,expected.value.primitive.plasma.massDensityKgM3);
    Require(actual.value.epochS==30&&actual.value.generation==1&&
        actual.value.eventIdentity==RS::EventIdentity(),"query provenance lost");
  });
  run("SINT02",[&] {
    Require(RS::Configure(argv[2],&error),error);Require(RS::Prepare(0,&error),error);
    const auto view=RS::CaptureInterval(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    const auto sample=view.value.AmbientAt({{0,0,4e10}},0);Require(sample.ok(),sample.status.message);
    // Charge neutrality of explicitly declared H+ and He++: ne=np+2*nalpha.
    // The fixture changes composition before evaluation; no observed residual
    // or production output is used to choose this independent expectation.
    Near(sample.value.electronNumberDensityM3,1.2*sample.value.protonNumberDensityM3);
    Near(sample.value.alphaNumberDensityM3,0.1*sample.value.protonNumberDensityM3);
    Near(sample.value.numberDensityM3,sample.value.electronNumberDensityM3);
    Require(sample.value.pressurePa>0&&sample.value.protonTemperatureK>0,"thermal state absent");
  });
  run("SINT03",[&] {
    const auto view=RS::CaptureInterval(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    const RS::Position point{{4e10,1e10,3e10}};
    const auto actual=view.value.AmbientAt(point,100);
    const auto expected=Direct(argv[2])->QueryAmbient({point[0],point[1],point[2]},100,1);
    Require(actual.ok()&&expected.ok(),"derivative query failed");
    for(int i=0;i<9;++i) {
      Near(actual.value.velocityGradientPerS[i],expected.value.gradientU[i]);
      Near(actual.value.magneticGradientTPerM[i],expected.value.gradientB[i]);
    }
    Require(actual.value.derivativesValid,"derivative validity lost");
  });
  run("SINT04",[&] {
    const auto view=RS::CaptureInterval(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    const auto before=RS::FrontEpoch();
    Require(!RS::CaptureInterval(0,601,RS::EventIdentity()).ok(),"uncovered interval accepted");
    Require(!RS::CaptureInterval(0,1,"wrong-event").ok(),"foreign identity accepted");
    Require(!view.value.AmbientAt({{0,0,4e10}},601).ok(),"future sample accepted");
    Require(!view.value.AmbientAt({{0,0,0}},30).ok(),"invalid domain accepted");
    Require(!view.value.RequirePhysicalDownstream().ok(),"downstream volume invented");
    Require(RS::FrontEpoch()==before&&RS::Metadata().epochS==0,"read query mutated committed epoch");
    const auto old=view.value.AmbientAt({{0,0,4e10}},30);Require(old.ok(),old.status.message);
    Require(RS::Configure(argv[1],&error),error);
    const auto retained=view.value.AmbientAt({{0,0,4e10}},30);Require(retained.ok(),retained.status.message);
    Near(retained.value.protonNumberDensityM3,old.value.protonNumberDensityM3);
    Require(retained.value.eventIdentity==old.value.eventIdentity,"retained view rebound to new provider");
    auto& store=SEP::Background::SnapshotStore::Instance();store.ResetForTests();
    const auto snapshot=[](double t,std::uint64_t generation) {
      return SEP::Background::BackgroundSnapshot(SEP::Background::Provider::ReducedShock,
          SEP::Background::Ownership::ModelOwned,t,t,t+1,generation,"integration-contract","immutable core",t,t);
    };
    store.Publish(snapshot(0,1));const auto first=store.Current();
    auto prepared=store.PreparePublication(snapshot(1,2));
    Require(store.Current()==first,"metadata preparation published early");
    {
      auto guard=store.BeginParticleRead(0);bool rejected=false;
      try {store.CommitPublication(&prepared);}catch(const std::logic_error&){rejected=true;}
      Require(rejected&&store.Current()==first,"metadata commit overlapped active reads");
    }
    store.CommitPublication(&prepared);Require(store.Current()->epoch_seconds()==1,"prepared metadata did not commit");
    bool rejected=false;try {store.CommitPublication(&prepared);}catch(const std::logic_error&){rejected=true;}
    Require(rejected,"single-use metadata publication was reused");
    auto stale=store.PreparePublication(snapshot(2,3));store.Publish(snapshot(3,4));
    rejected=false;try {store.CommitPublication(&stale);}catch(const std::logic_error&){rejected=true;}
    Require(rejected&&store.Current()->epoch_seconds()==3,"stale metadata replaced another writer");
    store.ResetForTests();
    Require(RS::Prepare(0,&error),error);const auto originalEpoch=RS::FrontEpoch();
    auto cadence=RS::PrepareCadencePublication(600.25);Require(cadence.ok(),cadence.status.message);
    Require(RS::Metadata().epochS==0&&RS::FrontEpoch()==originalEpoch&&cadence.value.Generation()==2,
        "cadence staging changed the committed provider");
    Require(RS::CommitCadencePublication(&cadence.value).ok(),"cadence candidate did not commit");
    Require(RS::Metadata().epochS==600&&!RS::CommitCadencePublication(&cadence.value).ok(),"cadence token was reused or incorrectly timed");
    auto foreign=RS::PrepareCadencePublication(1200);Require(foreign.ok(),foreign.status.message);
    Require(RS::Configure(argv[1],&error)&&RS::Prepare(0,&error),error);
    Require(!RS::CommitCadencePublication(&foreign.value).ok()&&RS::Metadata().epochS==0,"foreign cadence owner committed after reconfiguration");
  });
  run("SINT05",[&] {
    // Odd endpoint mantissa and two half-steps expose a one-ULP RK stage
    // overshoot. A strict closed provider must never be queried beyond its
    // declared endpoint, even though the constant-velocity oracle is simple.
    const double begin=std::nextafter(std::nextafter(12000.125,INFINITY),INFINITY);
    const double end=std::nextafter(12300.0,INFINITY);
    const auto closed=RS::IntegrateTrajectory([&](const RS::Position&,double t) {
      return t>=begin&&t<=end?SEP::Core::Result<RS::Position>::Success({{2,0,0}}):
          SEP::Core::Result<RS::Position>::Failure(SEP::Core::StatusCode::OutOfDomain,"RK stage escaped a closed interval");
    },{{1,0,0}},begin,end,{1e-9,400,1,1000});
    Require(closed.ok(),closed.status.message);Near(closed.value[0],1+2*(end-begin));
    RS::GeometryControls controls{1e-9,0.1,1,10000};
    const auto translated=RS::IntegrateTrajectory([](const RS::Position&,double) {
      return SEP::Core::Result<RS::Position>::Success({{2,3,-4}});
    },{{1,2,3}},0,0.75,controls);
    Require(translated.ok(),translated.status.message);
    Near(translated.value[0],2.5);Near(translated.value[1],4.25);
    // The exact z endpoint is zero: use an absolute SI roundoff bound there,
    // since relative error against zero has no mathematical meaning.
    Require(std::abs(translated.value[2])<1e-12,"zero translation endpoint exceeds roundoff bound [m]");
    const auto rotated=RS::IntegrateTrajectory([](const RS::Position& p,double) {
      return SEP::Core::Result<RS::Position>::Success({{-p[1],p[0],0}});
    },{{1,0,0}},0,0.75,controls);
    Require(rotated.ok(),rotated.status.message);
    Require(std::hypot(rotated.value[0]-std::cos(0.75),rotated.value[1]-std::sin(0.75))<1e-9,
        "manufactured rotating trajectory fails its exact reference");
    const auto reversed=RS::IntegrateTrajectory([](const RS::Position& p,double) {
      return SEP::Core::Result<RS::Position>::Success({{-p[1],p[0],0}});
    },rotated.value,0.75,0,controls);
    Require(reversed.ok()&&std::hypot(reversed.value[0]-1,reversed.value[1])<2e-9,
        "inverse material map fails its independent round-trip bound");
  });
  run("SINT06",[&] {
    Require(RS::Configure(argv[1],&error),error);Require(RS::Prepare(0,&error),error);
    const auto view=RS::CaptureInterval(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    const RS::Position initial{{4e10,1e10,3e10}};
    RS::GeometryControls controls{0.01,60,1e8,10000};
    const auto moved=RS::MoveVertex(view.value,initial,0,600,controls);Require(moved.ok(),moved.status.message);
    const auto direct=Direct(argv[1]);
    const double r0=std::hypot(std::hypot(initial[0],initial[1]),initial[2]);
    // Independent scalar quadrature/inversion of dr/dt=Ur(r), followed by
    // angular quadrature of dphi/dt=Omega*Rs^2/r^2. No RK trajectory routine
    // or production displacement is used to construct the expected point.
    const auto integrate=[&](double radius,bool angle) {
      constexpr int n=512;double sum=0;
      for(int i=0;i<=n;++i) {
        const double r=r0+(radius-r0)*i/n;
        const auto state=direct->Ambient().Evaluate({0,0,r},0);Require(state.ok(),state.status.message);
        const double u=state.value.velocityMPerS.z;
        const double value=angle?1/(r*r*u):1/u;
        sum+=(i==0||i==n?1:(i%2?4:2))*value;
      }
      return sum*(radius-r0)/(3*n);
    };
    double lo=r0,hi=r0+600*2e6;
    for(int i=0;i<52;++i) {const double middle=(lo+hi)/2;if(integrate(middle,false)<600)lo=middle;else hi=middle;}
    const double r=(lo+hi)/2;
    const auto& ambient=direct->Event().ambient.ambient;
    const double phi=ambient.rotationRateRadPerS*ambient.sourceSurfaceRadiusM*
        ambient.sourceSurfaceRadiusM*integrate(r,true);
    RS::Position reference{{r/r0*(std::cos(phi)*initial[0]-std::sin(phi)*initial[1]),
        r/r0*(std::sin(phi)*initial[0]+std::cos(phi)*initial[1]),r/r0*initial[2]}};
    double norm=0;for(int d=0;d<3;++d)norm=std::hypot(norm,moved.value[d]-reference[d]);
    Require(norm<0.2,"core trajectory differs from independent radial/azimuth quadrature by "+std::to_string(norm)+" m");
    const auto back=RS::MoveVertex(view.value,moved.value,600,0,controls);Require(back.ok(),back.status.message);
    norm=0;for(int d=0;d<3;++d)norm=std::hypot(norm,back.value[d]-initial[d]);
    Require(norm<0.04,"core inverse material map exceeds its declared accuracy");
    // Independent ideal-MHD material identities grade compatibility of the
    // provider's time evolution and its actual U, rather than just comparing
    // two adapters to the same point sample. Central differencing follows the
    // physical worldline at t=300 s. The normalized residual allowance 1e-5
    // covers the declared core Cartesian derivative stencil and RK error.
    const auto middle=RS::MoveVertex(view.value,initial,0,300,controls);
    Require(middle.ok(),middle.status.message);
    const auto minus=RS::MoveVertex(view.value,middle.value,300,299.5,controls),
        plus=RS::MoveVertex(view.value,middle.value,300,300.5,controls);
    Require(minus.ok()&&plus.ok(),"material derivative stencil trajectory failed");
    const auto c=view.value.AmbientAt(middle.value,300),a=view.value.PrimitiveAt(minus.value,299.5),
        b=view.value.PrimitiveAt(plus.value,300.5);
    Require(c.ok()&&a.ok()&&b.ok(),"material MHD identity state unavailable");
    const auto& g=c.value.velocityGradientPerS;const double theta=g[0]+g[4]+g[8];
    const double continuity=(b.value.massDensityKgM3-a.value.massDensityKgM3)+c.value.massDensityKgM3*theta;
    Require(std::abs(continuity)<1e-5*c.value.massDensityKgM3*std::abs(theta),"Parker material flow violates mass continuity");
    double induction=0,scale=0;
    for(int i=0;i<3;++i) {
      double expected=-c.value.magneticFieldT[i]*theta;
      for(int j=0;j<3;++j)expected+=g[3*i+j]*c.value.magneticFieldT[j];
      induction=std::hypot(induction,b.value.magneticFieldT[i]-a.value.magneticFieldT[i]-expected);
      scale=std::hypot(scale,c.value.magneticFieldT[i]*theta);
    }
    Require(induction<1e-5*scale,"Parker material flow violates ideal magnetic induction");
    const auto crossing=RS::CapturePhysicalStep(599.5,600.25,RS::EventIdentity());
    Require(crossing.ok()&&crossing.value.PrimitiveAt(initial,600.25).ok(),
        "continuous prescribed core cannot service the unchanged full step across a cache cadence");
    Require(RS::Metadata().epochS==0,"full-step query advanced the committed front");
  });
  run("SINT07",[&] {
    const auto view=RS::CaptureInterval(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    RS::GeometryControls controls{0.01,60,1e7,10000};
    for(double polarity:{1.,-1.}) {
      const RS::Position initial{{4e10,1e10,polarity*3e10}};
      const auto line=RS::BuildInitialLine(view.value,initial,1e9,65,0,controls);
      Require(line.ok(),line.status.message);
      double worst=0;
      for(std::size_t i=1;i<line.value.size();++i) {
        RS::Position midpoint,delta;
        double length=0;
        for(int d=0;d<3;++d){midpoint[d]=(line.value[i][d]+line.value[i-1][d])/2;delta[d]=line.value[i][d]-line.value[i-1][d];length=std::hypot(length,delta[d]);}
        const auto state=view.value.PrimitiveAt(midpoint,0);Require(state.ok(),state.status.message);
        double b=0;for(double value:state.value.magneticFieldT)b=std::hypot(b,value);
        double mismatch=0;for(int d=0;d<3;++d)mismatch=std::hypot(mismatch,delta[d]/length-state.value.magneticSector*state.value.magneticFieldT[d]/b);
        worst=std::max(worst,mismatch);
      }
      Require(worst<2e-6,"core magnetic trace does not resolve the local field tangent");
      Require(std::hypot(std::hypot(line.value.back()[0],line.value.back()[1]),line.value.back()[2])>
          std::hypot(std::hypot(initial[0],initial[1]),initial[2]),"negative magnetic sector reversed material ordering");
    }
  });
  run("SINT07R",[&] {
    const auto view=RS::CapturePhysicalStep(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    RS::GeometryControls controls{0.001,30,1e7,100000};const RS::Position initial{{4e10,1e10,3e10}};
    const auto measure=[&](int count) {
      const auto traced=RS::BuildInitialLine(view.value,initial,1e9,count,0,controls);Require(traced.ok(),traced.status.message);
      std::vector<RS::Position> moved;
      for(const auto& x:traced.value){const auto y=RS::MoveVertex(view.value,x,0,600,controls);Require(y.ok(),y.status.message);moved.push_back(y.value);}
      double maximum=0;
      for(std::size_t i=1;i<moved.size();++i) {
        RS::Position midpoint{},chord{};double length=0;
        for(int d=0;d<3;++d){midpoint[d]=(moved[i-1][d]+moved[i][d])/2;chord[d]=moved[i][d]-moved[i-1][d];length=std::hypot(length,chord[d]);}
        const auto field=view.value.PrimitiveAt(midpoint,600);Require(field.ok(),field.status.message);
        double magnitude=0;for(double value:field.value.magneticFieldT)magnitude=std::hypot(magnitude,value);
        double error=0;for(int d=0;d<3;++d)error=std::hypot(error,chord[d]/length-field.value.magneticSector*field.value.magneticFieldT[d]/magnitude);
        maximum=std::max(maximum,error);
      }
      return maximum;
    };
    const double coarse=measure(9),fine=measure(17);
    std::cout<<"METRIC SINT07R coarse_tangent="<<coarse<<" fine_tangent="<<fine<<" fixed_dt_s=600\n";
    Require(fine<2e-6&&coarse>3*fine,"advected Parker line lost B tangency or second-order chord refinement");
  });
  run("SINT10",[&] {
    const double epoch=12000;
    const auto view=RS::CapturePhysicalStep(epoch,epoch,RS::EventIdentity());Require(view.ok(),view.status.message);
    const auto direct=Direct(argv[1]);const auto trajectory=direct->Trajectory(epoch);Require(trajectory.ok(),trajectory.status.message);
    const double r=trajectory.value.apexRadiusM;
    const auto connected=RS::ConnectLine(view.value,{{{0,0,0.5*r}},{{0,0,1.5*r}}},epoch,0.001);
    Require(connected.ok(),connected.status.message);Require(connected.value.size()==1,"axis front connection missing or duplicated");
    const auto& hit=connected.value[0];Near(hit.fraction,0.5);Near(hit.front.positionM[2],r);
    Near(hit.front.normalSpeedMPerS,trajectory.value.apexSpeedMPerS);
    Require(hit.front.upstreamValid&&hit.front.acceptedShock&&hit.front.downstreamValid,"accepted core jump was lost");
    const auto ambient=view.value.PrimitiveAt(hit.front.positionM,epoch);Require(ambient.ok(),ambient.status.message);
    Near(hit.front.upstream.protonNumberDensityM3,ambient.value.protonNumberDensityM3);
    const auto exact=direct->EvaluateFrontPoint({0,0,r},epoch);Require(exact.ok(),exact.status.message);
    Near(hit.front.compressionRatio,exact.value.jump.compressionRatio);
  });
  run("SINT11",[&] {
    const double epoch=12000;const auto view=RS::CapturePhysicalStep(epoch,epoch,RS::EventIdentity());Require(view.ok(),view.status.message);
    const auto sphere=view.value.SphereAt(epoch);Require(sphere.ok(),sphere.status.message);
    const double apex=sphere.value.centerM[2]+sphere.value.radiusM;
    const auto connected=RS::ConnectLine(view.value,{{{0,0,0}},{{0,0,2*apex}}},epoch,0.001);
    Require(connected.ok(),connected.status.message);Require(connected.value.size()==1,"generating sphere rear root became a shock");
    Near(connected.value[0].front.positionM[2],apex);
    const auto absent=RS::ConnectLine(view.value,{{{0,0,-2*apex}},{{0,0,0}}},epoch,0.001);
    Require(absent.ok()&&absent.value.empty(),"unsupported backward ray became a shock");
  });
  run("SINT12",[&] {
    const double epoch=12000;const auto view=RS::CapturePhysicalStep(epoch,epoch,RS::EventIdentity());Require(view.ok(),view.status.message);
    const auto sphere=view.value.SphereAt(epoch);Require(sphere.ok(),sphere.status.message);
    const double a=sphere.value.radiusM,c=sphere.value.centerM[2],r=a+c;
    const auto endpoint=RS::ConnectLine(view.value,{{{0,0,0.5*r}},{{0,0,r}},{{0,0,1.5*r}}},epoch,0.001);
    Require(endpoint.ok(),endpoint.status.message);
    Require(endpoint.value.size()==1&&endpoint.value[0].segment==1&&endpoint.value[0].fraction==0,"endpoint ownership differs from floor(S)");
    const auto both=RS::ConnectLine(view.value,{{{-2*a,0,c}},{{2*a,0,c}}},epoch,0.001);
    Require(both.ok(),both.status.message);Require(both.value.size()==2,"genuine distinct cap intersections were coalesced");
    Near(both.value[0].fraction,0.25);Near(both.value[1].fraction,0.75);
    const auto tangent=RS::ConnectLine(view.value,{{{-a,0,r}},{{a,0,r}}},epoch,0.001);
    Require(tangent.ok(),tangent.status.message);
    Require(tangent.value.size()==1&&tangent.value[0].tangent,"apex tangency is missing or double-counted");
  });
  run("SINT13",[&] {
    const auto view=RS::CapturePhysicalStep(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    const auto before=RS::FrontEpoch();const auto r0=Direct(argv[1])->Trajectory(0).value.apexRadiusM;
    const std::vector<RS::Position> departure{{{0,0,0.95*r0}},{{0,0,1.5*r0}}};
    RS::GeometryControls controls{0.001,60,1e7,100000};
    const auto hit=RS::ConnectMaterialLine(view.value,departure,600,0.001,controls);
    Require(hit.ok(),hit.status.message);Require(hit.value.size()==1,"moving material line lost its front intersection");
    const auto a=RS::MoveVertex(view.value,departure[0],0,600,controls),b=RS::MoveVertex(view.value,departure[1],0,600,controls);
    Require(a.ok()&&b.ok(),"material endpoints unavailable");
    const auto& connection=hit.value[0];const double position=a.value[2]+connection.fraction*(b.value[2]-a.value[2]);
    Near(position,Direct(argv[1])->Trajectory(600).value.apexRadiusM);
    Require(connection.front.epochS==600&&RS::FrontEpoch()==before&&RS::Metadata().epochS==0,"time query mutated or mislabeled committed provider state");
    const double radius600=Direct(argv[1])->Trajectory(600).value.apexRadiusM;
    const std::vector<RS::Position> exitLine{{{0,0,r0-1e7}},{{0,0,r0+1e7}}};
    const auto exitBefore=RS::ConnectMaterialLine(view.value,exitLine,0,0.001,controls),
        exitAfter=RS::ConnectMaterialLine(view.value,exitLine,600,0.001,controls);
    Require(exitBefore.ok()&&exitAfter.ok()&&exitBefore.value.size()==1&&exitAfter.value.empty(),
        "material connection loss during the physical interval was not resolved");
    std::vector<RS::Position> entryLine;
    for(double offset:{-1e7,1e7}){const auto inverse=RS::MoveVertex(view.value,{{0,0,radius600+offset}},600,0,controls);
      Require(inverse.ok(),inverse.status.message);entryLine.push_back(inverse.value);}
    const auto entryBefore=RS::ConnectMaterialLine(view.value,entryLine,0,0.001,controls),
        entryAfter=RS::ConnectMaterialLine(view.value,entryLine,600,0.001,controls);
    Require(entryBefore.ok()&&entryAfter.ok()&&entryBefore.value.empty()&&entryAfter.value.size()==1,
        "material connection entry during the physical interval was not resolved");
  });
  run("SINT14",[&] {
    Spatial spatial;spatial.derivative=1;ZeroRandom random;
    T::ParkerBackground background(0,0.3);background.logAreaGradientPerM=0.25;
    const auto increment=T::AdvanceParker({0,2},background,10,0.75,spatial,random);
    Require(increment.status.ok(),increment.status.message);
    // Forward physical-number particles represent N=4*pi*p^2*A*f per ds dp.
    // Rewriting (1/A) d_s(A*kappa*d_s f) in forward conservative form gives
    // drift kappa'+kappa*dln(A)/ds, independently of the S-to-arc conversion.
    Near(increment.deterministicDisplacementM,(1+4*0.25)*0.75);
    Near(increment.state.momentumKgMPerS,2*std::exp(-0.3*0.75/3));
    spatial.derivative=0;background.velocityDivergencePerS=0;
    constexpr std::size_t n=32768;double sum=0,squares=0;
    for(std::size_t i=0;i<n;++i) {
      T::KeyedRandomStream stream(4289,i+1,101,0);
      const auto value=T::AdvanceParker({0,2},background,10,0.75,spatial,stream);
      Require(value.status.ok(),value.status.message);sum+=value.displacementM;squares+=value.displacementM*value.displacementM;
    }
    const double mean=sum/n,variance=squares/n-mean*mean,expectedMean=0.75,expectedVariance=6;
    // Fixed predeclared 6-sigma normal/chi-square moment bounds. This exact
    // Gaussian reference has no fitted acceptance parameter or observed data.
    Require(std::abs(mean-expectedMean)<=6*std::sqrt(expectedVariance/n),"Parker metric Gaussian mean failed");
    Require(std::abs(variance-expectedVariance)<=6*expectedVariance*std::sqrt(2./(n-1)),"Parker metric Gaussian variance failed");
    const auto draws=random.draws;background.logAreaGradientPerM=std::numeric_limits<double>::quiet_NaN();
    Require(!T::AdvanceParker({0,2},background,10,0.75,spatial,random).status.ok()&&random.draws==draws,"invalid metric consumed randomness or was accepted");
    std::cout<<"METRIC SINT14 N="<<n<<" mean_m="<<mean<<" variance_m2="<<variance<<" expected_mean_m=0.75 expected_variance_m2=6\n";
  });
  run("SINT15",[&] {
    T::FocusedTransportBackground background(-0.1,0,0,0.3);background.fieldAlignedStrainPerS=0.07;
    background.equationMode=T::FocusedEquationMode::FullGyrotropic;
    Near(T::FocusedPitchDriftPerS(0,10,background),0.5);
    Near(T::FocusedPitchDriftPerS(0.5,10,background),0.391875);
    Near(T::FocusedLogMomentumRatePerS(0.5,background),-0.10375);
    // Three-point Gauss integration is exact for the quadratic cooling rate.
    const double x=std::sqrt(3./5);
    const double average=(5*T::FocusedLogMomentumRatePerS(-x,background)+
        8*T::FocusedLogMomentumRatePerS(0,background)+5*T::FocusedLogMomentumRatePerS(x,background))/18;
    Near(average,-0.3/3);
    background.velocityDivergencePerS=background.fieldAlignedStrainPerS=0;
    Pitch pitch;ZeroRandom random;const auto p=T::MomentumFromSpeed(10,1,1000);Require(p.status.ok(),p.status.message);
    const auto evolve=[&](int n) {
      T::FocusedTransportState state(0,p.value,0.2);
      for(int i=0;i<n;++i) {
        const auto increment=T::AdvanceFocusedTransportDmumu(state,background,1,1000,0.75/n,pitch,random,nullptr);
        Require(increment.status.ok(),increment.status.message);state=increment.state;
      }
      return state;
    };
    const auto coarse=evolve(75),fine=evolve(150);
    const double angle=std::atanh(0.2),expected=std::tanh(angle+0.375);
    const double path=20*std::log(std::cosh(angle+0.375)/std::cosh(angle));
    Require(std::abs(fine.mu-expected)<3e-7&&std::abs(fine.arcLengthM-path)<1e-5,"focused analytical focusing/path failed");
    Require(std::abs(coarse.mu-expected)>3*std::abs(fine.mu-expected),"focused deterministic refinement lost second-order convergence");
    // Exact simultaneous compression and streaming: theta=3 sigma makes
    // gyrotropic pitch drift vanish, p=p0 exp(-sigma t), while relativistic
    // v(p) changes. This detects final-speed-for-the-whole-step integration.
    background.dLnAbsBdsPerM=0;background.velocityDivergencePerS=0.6;
    background.fieldAlignedStrainPerS=0.2;
    const double duration=0.75,mu=0.6,finishMomentum=p.value*std::exp(-0.2*duration);
    const double exactPath=mu*1000/0.2*(std::asinh(p.value/1000)-std::asinh(finishMomentum/1000));
    const auto coolingPath=[&](int intervals) {
      T::FocusedTransportState state(0,p.value,mu);
      for(int i=0;i<intervals;++i){const auto step=T::AdvanceFocusedTransportDmumu(state,background,1,1000,duration/intervals,pitch,random,nullptr);
        Require(step.status.ok(),step.status.message);state=step.state;}
      Near(state.momentumKgMPerS,finishMomentum);Near(state.mu,mu);return std::abs(state.arcLengthM-exactPath);
    };
    const double coolingCoarse=coolingPath(8),coolingFine=coolingPath(16);
    Require(coolingFine<4e-5&&coolingCoarse>3*coolingFine,"focused cooling/streaming lost second-order time centring");
  });
  run("SINT16",[&] {
    constexpr double v=10,nu=4,total=12.5,c=1000;
    constexpr std::size_t n=2048;
    const double kappa=v*v/(3*nu),exactVariance=2*kappa*(total-(1-std::exp(-nu*total))/nu);
    const auto momentum=T::MomentumFromSpeed(v,1,c);Require(momentum.status.ok(),momentum.status.message);
    Pitch pitch;pitch.rate=nu;Mfp mfp;Spatial spatial;spatial.kappa=kappa;
    T::FocusedTransportBackground diffusion;diffusion.equationMode=T::FocusedEquationMode::FullGyrotropic;
    T::FocusedTransportMfpBackground events;events.equationMode=T::FocusedEquationMode::FullGyrotropic;
    double sums[3]={},squares[3]={},fourth[3]={};
    for(std::size_t i=0;i<n;++i) {
      const double mu=-1+2*(i+0.5)/n; // exact stratified isotropic initial directions
      T::KeyedRandomStream rd(5297,i+1,102,0),re(5297,i+1,103,0),rp(5297,i+1,104,0);
      T::FocusedTransportState state(0,momentum.value,mu);
      constexpr int steps=5000;
      for(int j=0;j<steps;++j) {
        const auto next=T::AdvanceFocusedTransportDmumu(state,diffusion,1,c,total/steps,pitch,rd,nullptr);
        Require(next.status.ok(),next.status.message);state=next.state;
      }
      const auto event=T::AdvanceFocusedTransportMfp({0,momentum.value,mu},events,1,c,total,0.25,mfp,re,nullptr);
      Require(event.status.ok(),event.status.message);
      const auto parker=T::AdvanceParker({0,momentum.value},{},v,total,spatial,rp);Require(parker.status.ok(),parker.status.message);
      const double positions[3]={state.arcLengthM,event.state.arcLengthM,parker.state.arcLengthM};
      for(int model=0;model<3;++model){const double x=positions[model];sums[model]+=x;squares[model]+=x*x;fourth[model]+=x*x*x*x;}
    }
    for(int model=0;model<3;++model) {
      const double mean=sums[model]/n,second=squares[model]/n,expected=model==2?2*kappa*total:exactVariance;
      // Six-standard-error ensemble criteria are fixed before evaluation.
      // Both focused processes have <mu(t)mu(0)>=(1/3)exp(-nu*t), giving
      // the stated finite-time variance independently of their different
      // higher angular harmonics. Parker is its long-time diffusion limit.
      const double error=6*std::sqrt(std::max(0.,fourth[model]/n-second*second)/n);
      Require(std::abs(mean)<6*std::sqrt(expected/n),"common-limit mean failed");
      Require(std::abs(second-expected)<=error+0.01*expected,"common-limit second moment failed");
      std::cout<<"METRIC SINT16 model="<<model<<" N="<<n<<" second_m2="<<second<<" expected_m2="<<expected<<" sampling_6se_m2="<<error<<" declared_discretization_fraction=0.01\n";
    }
  });
  run("SINT17",[&] {
    namespace RP=SEP::ReducedParticles;
    RP::Configuration c;Require(RP::Validate(c).ok(),"background-only default changed");
    c.enabled=true;Require(!RP::Validate(c).ok(),"missing L_ref was invented");
    c.releaseDistanceM=1e7;c.minimumEnergyJ=1e-14;c.maximumEnergyJ=2e-14;
    c.releaseFraction=0.001;c.samplesPerStep=8;
    Require(RP::Validate(c).ok(),"explicit manufactured source input rejected");
    c.returnPolicy=RP::ReturnPolicy::Pass;Require(!RP::Validate(c).ok(),"pass policy silently accepted");
    c.returnPolicy=RP::ReturnPolicy::Delete;c.normalization=RP::Normalization::GrossEmission;
    Require(!RP::Validate(c).ok(),"gross source silently reinterpreted as net");
  });
  run("SINT18",[&] {
    for(double compression:{2.,3.,4.}) {
      const auto spectrum=SEP::CoronalCME::MomentumSpectrum::LocalCompressionDsa(2,10,compression);
      Require(spectrum.ok(),spectrum.status.message);
      const double exponent=2-3*compression/(compression-1);
      Near(spectrum.value.NumberExponent(),exponent);
      for(double u:{0.1,0.5,0.9}) {
        const double p=spectrum.value.InverseCdf(u),s=exponent+1;
        const double independent=(std::pow(p,s)-std::pow(2.,s))/(std::pow(10.,s)-std::pow(2.,s));
        Near(independent,u);
      }
    }
  });
  run("SINT19",[&] {
    const auto view=RS::CapturePhysicalStep(7200,7201,RS::EventIdentity());Require(view.ok(),view.status.message);
    constexpr double distance=1e7;
    for(double epoch:{7200.,7200.37,7201.}) {
      const auto sphere=view.value.SphereAt(epoch);Require(sphere.ok(),sphere.status.message);
      // An oblique chord is independent of the core's tracing algorithm. Its
      // analytic offset-sphere root verifies that normal displacement is not
      // mistaken for an intersection on the same material field line.
      const double z=sphere.value.centerM[2]+sphere.value.radiusM;
      std::vector<RS::Position> vertices{{{0,0,z-1e8}},{{1e7,0,z+1e8}}};
      const auto points=SEP::ReducedParticles::ReferenceIntersections(view.value,vertices,epoch,distance,1e-3);
      Require(points.ok(),points.status.message);Require(points.value.size()==1,"reference root missing or rear retained");
      const auto& p=points.value[0];
      for(int d=0;d<3;++d)Near(p.xM[d],vertices[0][d]+p.fraction*(vertices[1][d]-vertices[0][d]),1e-12);
      double r=0;for(int d=0;d<3;++d)r=std::hypot(r,p.xM[d]-sphere.value.centerM[d]);
      Require(std::abs(r-sphere.value.radiusM-distance)<1e-3,"normal distance incorrect");
      Require(p.normalFoot.epochS==epoch,"release used a different birth epoch");
    }
  });
  run("SINT20",[&] {
    RS::FrontSample front;front.acceptedShock=front.upstreamValid=true;
    front.inflowMPerS=5;front.outwardNormal={{0,0,1}};
    front.upstream.magneticFieldT={{3,0,4}};
    front.upstream.protonNumberDensityM3=7;front.upstream.electronNumberDensityM3=70;
    const auto rate=SEP::ReducedParticles::NetProtonRate(front,8,0.2);
    Require(rate.ok(),rate.status.message);Near(rate.value,0.2*7*5*8/4);
    front.upstream.electronNumberDensityM3*=10;
    Near(SEP::ReducedParticles::NetProtonRate(front,8,0.2).value,rate.value);
    Near(SEP::ReducedParticles::NetProtonRate(front,8,0).value,0);
    front.upstream.magneticFieldT={{1,0,0}};
    Require(!SEP::ReducedParticles::NetProtonRate(front,8,0.2).ok(),"tangency silently floored");
  });
  run("SINT21",[&] {
    const double speed=10,critical=0.12345;
    const auto distribution=SEP::CoronalCME::BuildFocusedEscapeDistribution(speed,1,1,-speed*critical,1024,0);
    Require(distribution.ok(),distribution.status.message);
    Near(distribution.value.normalizationMPerS,speed*(1-critical)*(1-critical)/4,1e-12);
    for(double u:{1e-12,0.1,0.5,0.9}) {
      const auto mu=SEP::CoronalCME::SampleFocusedMu(distribution.value,u);Require(mu.ok(),mu.status.message);
      Require(mu.value>critical,"sample has negative/non-outward normal flux");
      Near(mu.value,critical+(1-critical)*std::sqrt(u),1e-10);
    }
  });
  run("SINT22",[&] {
    constexpr unsigned n=65536;unsigned hits=0;
    T::KeyedRandomStream random(773,1,105,0);
    for(unsigned i=0;i<n;++i){const auto hit=T::AbsorbingLinearBridge(1,1,1,random);Require(hit.status.ok(),hit.status.message);hits+=hit.hit;
      if(hit.hit)Require(hit.fraction>0&&hit.fraction<=1,"invalid first-hit fraction");}
    const double p=std::exp(-1.),bound=6*std::sqrt(n*p*(1-p));
    Require(std::abs(hits-n*p)<=bound,"hidden Brownian hits disagree with independent bridge probability");
    std::cout<<"METRIC SINT22 N="<<n<<" hits="<<hits<<" expected="<<n*p<<" six_sigma="<<bound<<'\n';
  });
  run("SINT23",[&] {
    constexpr unsigned n=65536;unsigned hits=0,early=0;
    T::KeyedRandomStream endpoint(997,1,106,0),boundary(997,1,107,0);
    for(unsigned i=0;i<n;++i){const auto hit=T::AbsorbingLinearBridge(2,2+std::sqrt(2.)*endpoint.Normal01(),1,boundary);
      Require(hit.status.ok(),hit.status.message);hits+=hit.hit;early+=hit.hit&&hit.fraction<=0.25;}
    const double p=std::erfc(1.),q=std::erfc(2.);
    Require(std::abs(hits-n*p)<=6*std::sqrt(n*p*(1-p)),"unconditional first passage is wrong");
    Require(std::abs(early-n*q)<=6*std::sqrt(n*q*(1-q)),"first-hit time law is wrong");
    ZeroRandom none;const auto deterministic=T::AbsorbingLinearBridge(2,-2,0,none);
    Require(deterministic.status.ok()&&deterministic.hit&&deterministic.fraction==0.5&&none.draws==0,"deterministic crossing drew noise");
    Require(!T::AbsorbingLinearBridge(1,1,-1,none).status.ok()&&none.draws==0,"invalid bridge consumed randomness");
    std::cout<<"METRIC SINT23 N="<<n<<" hits="<<hits<<" expected="<<n*p<<" early_hits="<<early<<" expected_early="<<n*q<<'\n';
    unsigned movingHits=0,movingEarly=0;
    T::KeyedRandomStream movingEndpoint(20261122,2,108,0),movingBoundary(20261122,2,109,0);
    for(unsigned i=0;i<n;++i) {
      const auto hit=T::AbsorbingLinearBridge(2,1+std::sqrt(2.)*movingEndpoint.Normal01(),1,movingBoundary);
      Require(hit.status.ok(),hit.status.message);movingHits+=hit.hit;movingEarly+=hit.hit&&hit.fraction<=0.25;
    }
    // A boundary moving toward an initially stationary Brownian particle at
    // 1 m/s is drift a=-1 in signed distance. The independent inverse-Gaussian
    // first-passage CDF is Phi(-(d+a*t)/sqrt(2*k*t)) plus its image term.
    const auto movingCdf=[](double t) {
      return 0.5*std::erfc((2-t)/std::sqrt(4*t))+
          0.5*std::exp(2.)*std::erfc((2+t)/std::sqrt(4*t));
    };
    for(const auto pair:{std::make_pair(movingHits,movingCdf(1)),std::make_pair(movingEarly,movingCdf(0.25))})
      Require(std::abs(pair.first-n*pair.second)<=6*std::sqrt(n*pair.second*(1-pair.second)),
          "moving-boundary first-passage distribution is wrong");
    std::cout<<"METRIC SINT23 moving_front_m_per_s=1 moving_hits="<<movingHits<<" expected_moving="<<n*movingCdf(1)
        <<" moving_early_hits="<<movingEarly<<" expected_moving_early="<<n*movingCdf(0.25)<<'\n';
  });
  run("SINT24",[&] {
    Mfp mfp;ZeroRandom random;
    const auto p0=T::MomentumFromSpeed(10,1,10000);Require(p0.status.ok(),p0.status.message);
    T::FocusedTransportMfpState initial(1,p0.value,-1);initial.remainingOpticalDepth=2;
    unsigned checks=0;
    const auto boundary=[&](const T::FocusedTransportMfpState& a,
        const T::FocusedTransportMfpState& b,double elapsed,double elapsedEnd) {
      ++checks;T::MfpBoundaryDecision d;d.status=T::Status::Ok();
      Require(elapsed==0&&elapsedEnd>elapsed,"boundary was queried after the first event");
      if(a.arcLengthM>0&&b.arcLengthM<=0){d.absorb=true;d.fraction=a.arcLengthM/(a.arcLengthM-b.arcLengthM);}
      return d;
    };
    const auto killed=T::AdvanceFocusedTransportMfp(initial,{},1,10000,1,1,mfp,random,nullptr,boundary);
    Require(killed.status.ok(),killed.status.message);
    Require(killed.absorbed&&checks==1&&killed.diagnostics.scatteringEvents==0&&random.draws==0,
        "MFP scattered/drew a new clock after the absorbing hit");
    Near(killed.elapsedS,0.1);Require(std::abs(killed.state.arcLengthM)<1e-12,"MFP first-hit coordinate is wrong");
    Near(killed.state.momentumKgMPerS,p0.value);
    Near(killed.state.remainingOpticalDepth,1.6);
  });
  run("SINT26",[&] {
    const auto view=RS::CapturePhysicalStep(0,600,RS::EventIdentity());Require(view.ok(),view.status.message);
    for(double radius:{1e10,4e10,149597870700.})for(double epoch:{0.,30.,599.25}) {
      const auto pole=view.value.PrimitiveAt({{0,0,radius}},epoch);
      const auto tilted=view.value.PrimitiveAt({{0.6*radius,0,0.8*radius}},epoch);
      Require(pole.ok(),pole.status.message);Require(tilted.ok(),tilted.status.message);
      // Independent axial-dipole and radial-wind identities at one radius:
      // B_r proportional cos(theta), B_z=B_r*cos(theta), U_z=U_r*cos(theta).
      // Toroidal components have zero z projection. The off-axis query retains
      // the full winding calculation and is not an axis-specialized oracle.
      Near(pole.value.magneticFieldT[2],tilted.value.magneticFieldT[2]/0.64);
      Near(pole.value.velocityMPerS[2],tilted.value.velocityMPerS[2]/0.8);
      Near(pole.value.protonNumberDensityM3,tilted.value.protonNumberDensityM3);
    }
  });
  run("SINT27",[&] {
    namespace C=SEP::CoronalCME;
    // Independent gas-dynamic RH identities for exactly parallel magnetic
    // fields. B_n is continuous and cancels from the normal momentum/energy
    // jump. Vary gamma, magnetic sign and strength, Mach number and a common
    // tangential flow; the ordinary production residual/entropy/fast-branch
    // checks must still accept the physical compressive state. A nonzero B_t
    // is deliberately not approximated as parallel by the implementation.
    for(double gamma:{4./3.,5./3.})for(double mach:{1.01,2.,5.})
      for(double field:{-1e-9,0.,1e-9})for(double tangent:{0.,2e4}) {
        C::MhdPrimitiveState upstream{1e-20,2e-10,{4e5,tangent,-1.5*tangent},{field,0,0}};
        const double sound=std::sqrt(gamma*upstream.pressurePa/upstream.massDensityKgM3);
        const double speed=upstream.velocityMPerS.x+mach*sound;
        const auto solved=C::SolveObliqueFastShock(upstream,{1,0,0},speed,gamma,1e-9);
        Require(solved.ok(),solved.status.message);
        const double compression=(gamma+1)*mach*mach/((gamma-1)*mach*mach+2);
        Near(solved.value.compressionRatio,compression);
        Near(solved.value.downstream.massDensityKgM3/upstream.massDensityKgM3,compression);
        Near(solved.value.downstream.pressurePa/upstream.pressurePa,
            (2*gamma*mach*mach-(gamma-1))/(gamma+1));
        Near(solved.value.downstream.velocityMPerS.x,speed-mach*sound/compression);
        Near(solved.value.downstream.velocityMPerS.y,upstream.velocityMPerS.y);
        Near(solved.value.downstream.velocityMPerS.z,upstream.velocityMPerS.z);
        Require(solved.value.downstream.magneticFieldT.x==field&&
                solved.value.downstream.magneticFieldT.y==0&&solved.value.downstream.magneticFieldT.z==0,
                "parallel RH optimization changed magnetic flux");
        Require(solved.value.residuals.maximum<=1e-9,"parallel RH residual gate was bypassed");
      }
  });
  run("SINT28",[&] {
    Mfp mfp;const auto momentum=T::MomentumFromSpeed(10,1,10000);
    Require(momentum.status.ok(),momentum.status.message);
    // These are internal numerical durations, not changes to AMPS physical
    // dt. A scheduled scatter leaves a non-dyadic remainder. Completing that
    // remainder must report the caller's EXACT closed endpoint; recomputing
    // elapsed+(total-elapsed) can exceed it by an ULP and caused the native
    // moving-boundary guard to reject the longer application campaign.
    unsigned wrong=0;mfp.rate=0.01;
    T::KeyedRandomStream clockCases(20261122,1,112,0);
    for(unsigned i=1;i<=4096;++i) {
      const double duration=1+39*clockCases.UniformOpen01();
      const double eventTime=duration*(0.05+0.9*clockCases.UniformOpen01());
      T::FocusedTransportMfpState initial(1,momentum.value,1);
      initial.remainingOpticalDepth=mfp.rate*eventTime;ZeroRandom random;
      const auto advanced=T::AdvanceFocusedTransportMfp(initial,{},1,10000,
          duration,duration,mfp,random,nullptr);
      Require(advanced.status.ok(),advanced.status.message);
      Require(!advanced.absorbed&&advanced.diagnostics.scatteringEvents==1,
          "closed-interval fixture did not exercise a scatter and a remainder");
      wrong+=advanced.elapsedS!=duration;
    }
    std::cout<<"METRIC SINT28 intervals=4096 inexact_closed_endpoints="<<wrong<<'\n';
    Require(wrong==0,"MFP final numerical time differs from its exact closed endpoint");
    mfp.rate=4;
    for(double fraction:{0.,0.25,0.75,std::nextafter(1.,0.),1.}) {
      T::FocusedTransportMfpState initial(1,momentum.value,1);initial.remainingOpticalDepth=0.1;
      ZeroRandom random;unsigned calls=0;double hitBegin=0,hitEnd=0;
      const double duration=0.133;
      const auto boundary=[&](const T::FocusedTransportMfpState&,
          const T::FocusedTransportMfpState&,double begin,double end) {
        ++calls;T::MfpBoundaryDecision decision;decision.status=T::Status::Ok();
        if(calls==2){decision.absorb=true;decision.fraction=fraction;hitBegin=begin;hitEnd=end;}
        return decision;
      };
      const auto killed=T::AdvanceFocusedTransportMfp(initial,{},1,10000,
          duration,duration,mfp,random,nullptr,boundary);
      Require(killed.status.ok(),killed.status.message);
      Require(killed.absorbed&&calls==2&&random.draws==3&&
              killed.elapsedS>=hitBegin&&killed.elapsedS<=hitEnd,
              "absorbed MFP clock escaped its explicit interval or consumed later noise");
      if(fraction==0)Require(killed.elapsedS==hitBegin,"zero-fraction hit changed its beginning");
      if(fraction==1)Require(killed.elapsedS==duration,"unit-fraction hit lost its exact endpoint");
    }
  });
  run("SINT29",[&] {
    unsigned rejected=0;double maximumError=0;
    // The final microsecond-scale remainder is much shorter than the other
    // RK intervals. Constant flow has an independent exact trajectory even at
    // a heliospheric coordinate where position arithmetic resolves ~1e-5 m.
    // The absolute endpoint budget stays 1 mm; distributing an infinitesimal
    // fresh budget to that tail must not create thousands of useless retries.
    for(int sign:{-1,1})for(int i=1;i<=32;++i) {
      const double begin=25200,end=begin+sign*(300+i*1e-6);
      RS::GeometryControls controls;controls.positionToleranceM=0.001;
      controls.maximumTimeSubstepS=37.5;controls.maximumAttempts=128;
      const auto moved=RS::IntegrateTrajectory([&](const RS::Position&,double t) {
        Require(t>=std::min(begin,end)&&t<=std::max(begin,end),"trajectory queried beyond closed time support");
        return SEP::Core::Result<RS::Position>::Success({{0,0,6e5}});
      },{{0,0,1e11}},begin,end,controls);
      if(!moved.ok()){++rejected;continue;}
      const double error=std::abs(moved.value[2]-(1e11+6e5*(end-begin)));
      maximumError=std::max(maximumError,error);
      Require(error<=0.001,"trajectory exceeded its unchanged global endpoint tolerance");
    }
    std::cout<<"METRIC SINT29 trajectories=64 rejected="<<rejected<<" maximum_endpoint_error_m="<<maximumError<<'\n';
    Require(rejected==0,"constant-flow trajectory exhausted work/progress at a short final remainder");
  });
  run("SINT30",[&] {
    // Exact failing inputs captured from the complete FTE-MFP application at
    // iteration 84, identically on one and four MPI ranks. This is a material
    // vertex query inside a 300-s physical step, not a shorter particle/PIC
    // step. The old controller accepted 49989 almost-zero intervals and
    // rejected 50011 more, exhausting its unchanged 100000-attempt bound.
    const double begin=25200,end=25237.503973302857;
    const RS::Position initial{{0,0,149804795566.1752}};
    const RS::GeometryControls controls{0.0085036544248850987,37.5,1,100000};
    const auto view=RS::CapturePhysicalStep(begin,25500,RS::EventIdentity());
    Require(view.ok(),view.status.message);
    const auto moved=RS::MoveVertex(view.value,initial,begin,end,controls);
    Require(moved.ok(),moved.status.message);

    // Independent trajectory oracle: fixed-step explicit midpoint, rather
    // than production RK4/step doubling. Keep the accumulator and time grid
    // in long double so repeated additions at 1 AU do not dominate a mm-scale
    // endpoint comparison. Only the common physical core velocity is shared.
    // Three declared grids must converge before grading the production map;
    // no tolerance is inferred from the observed production residual.
    const auto reference=[&](unsigned count) {
      long double z=initial[2];
      const long double h=(static_cast<long double>(end)-begin)/count;
      for(unsigned i=0;i<count;++i) {
        const long double t=begin+i*h;
        const auto a=view.value.PrimitiveAt({{0,0,static_cast<double>(z)}},static_cast<double>(t));
        Require(a.ok(),a.status.message);
        const auto b=view.value.PrimitiveAt({{0,0,static_cast<double>(z+h*a.value.velocityMPerS[2]/2)}},
            static_cast<double>(t+h/2));
        Require(b.ok(),b.status.message);z+=h*b.value.velocityMPerS[2];
      }
      return z;
    };
    const auto coarse=reference(8192),middle=reference(16384),fine=reference(32768);
    Require(std::abs(coarse-middle)<=controls.positionToleranceM/64&&
        std::abs(middle-fine)<=controls.positionToleranceM/64,
        "independent midpoint trajectory oracle did not converge at its declared grids");
    const double residual=static_cast<double>(std::abs(static_cast<long double>(moved.value[2])-fine));
    Require(residual<=controls.positionToleranceM,"captured material trajectory exceeds its unchanged SI endpoint tolerance");
    Require(moved.value[0]==0&&moved.value[1]==0,"polar material trajectory left its invariant axis");
    std::cout<<"METRIC SINT30 endpoint_error_m="<<residual<<" position_tolerance_m="<<controls.positionToleranceM
        <<" midpoint_refinement_difference_m="<<static_cast<double>(std::abs(middle-fine))<<'\n';
  });
  run("SINT31",[&] {
    namespace C=SEP::CoronalCME;namespace A=SEP::CoronaSwcme;
    auto definition=Direct(argv[1])->Ambient().Definition();
    // Non-axisymmetric harmonics deliberately distinguish the PFSS pole
    // evaluator's azimuthal limits. Signed Cartesian zeros affect atan2 even
    // though they describe the same mathematical axis. An exact cache must
    // preserve those supplied inputs and the existing pole convention.
    definition.ambient.harmonics="1,0,0.0001,0;1,1,0.00002,-0.00003;2,1,-0.00001,0.000015";
    definition.physicsFingerprint="manufactured-axial-source-cache-identity";
    const std::vector<C::HarmonicCoefficient> harmonics{{1,0,0.0001,0},
        {1,1,0.00002,-0.00003},{2,1,-0.00001,0.000015}};
    const auto uncached=C::PfssHarmonics::Create(definition.support.solarRadiusM,
        definition.ambient.sourceSurfaceRadiusM,harmonics);
    Require(uncached.ok(),uncached.status.message);
    const auto model=A::AmbientModel::Create(definition);Require(model.ok(),model.status.message);
    const double inside=definition.ambient.sourceSurfaceRadiusM*(1-64*std::numeric_limits<double>::epsilon());
    unsigned checks=0;
    for(int sign:{-1,1})for(double xzero:{0.,-0.})for(double yzero:{0.,-0.})
      for(double radius:{1e10,4e10,149597870700.})for(double epoch:{0.,30.,599.25}) {
        const C::Vec3 p{xzero,yzero,sign*radius};const double r=C::Norm(p);
        const double angle=-definition.ambient.rotationRateRadPerS*epoch;
        const C::Vec3 body{std::cos(angle)*p.x-std::sin(angle)*p.y,
            std::sin(angle)*p.x+std::cos(angle)*p.y,p.z};
        const auto direction=C::Unit(body),foot=(inside/C::Norm(direction))*direction;
        const auto source=uncached.value.EvaluateCartesian(foot);Require(source.ok(),source.status.message);
        const double br=C::Dot(source.value,C::Unit(foot))*
            std::pow(definition.ambient.sourceSurfaceRadiusM/r,2);
        const auto expected=br*(p/r);
        const auto actual=model.value->Evaluate(p,epoch);
        Require(actual.ok(),actual.status.message);
        Require(actual.value.magneticFieldT.x==expected.x&&actual.value.magneticFieldT.y==expected.y&&
            actual.value.magneticFieldT.z==expected.z,"axial ambient differs from the uncached canonical PFSS source query");
        ++checks;
      }
    // Cache presence must not bypass the ordinary radial magnetic-null gate.
    // This manufactured threshold rejects the outer point; it does not floor
    // its field or make this synthetic model a production configuration.
    definition.ambient.minimumMagneticFieldT=1e-3;
    const auto rejecting=A::AmbientModel::Create(definition);Require(rejecting.ok(),rejecting.status.message);
    Require(!rejecting.value->Evaluate({0,0,149597870700.},30).ok(),"axial source cache bypassed the physical magnetic-null rejection");
    std::cout<<"METRIC SINT31 exact_uncached_pfss_checks="<<checks<<'\n';
  });
  run("SINT32",[&] {
    // A manufactured material-coordinate front s(t)=t^2 has an independent
    // maximum linear-chord error h^2/4. Keep the spatial bound 1 mm, the
    // quarter/mid checks and every accepted h identical in cold/warm controls.
    // Repeatedly starting at the rejected 8-s proposal costs eight checks;
    // remembering the last CHECKED h must avoid repeating that work.
    const auto exercise=[&](bool warm) {
      T::BoundaryIntervalControl history;unsigned calls=0;double begin=0;
      for(unsigned step=0;step<128;++step) {
        double accepted=-1;
        const auto status=T::RefineBoundaryInterval(begin,8,1e-12,[&](double h,bool* resolved) {
          ++calls;*resolved=true;
          for(double fraction:{0.25,0.5,0.75}) {
            const double t=begin+fraction*h;
            const double chord=begin*begin+fraction*((begin+h)*(begin+h)-begin*begin);
            if(std::abs(t*t-chord)>1e-3)*resolved=false;
          }
          return T::Status::Ok();
        },&accepted,warm?&history:nullptr);
        Require(status.ok(),status.message);
        Require(accepted==0.0625&&accepted*accepted/4<=1e-3,
            "boundary hint changed the independent accepted curvature interval");
        begin+=accepted;
      }
      Require(begin==8,"boundary history changed the complete manufactured physical endpoint");
      return calls;
    };
    const unsigned cold=exercise(false),warm=exercise(true);
    std::cout<<"METRIC SINT32 cold_checks="<<cold<<" warm_checks="<<warm<<" accepted_intervals=128 spatial_bound_m=0.001\n";
    Require(cold==1024&&warm<=200,"predictable boundary refinement repeats cold coarse queries");
    T::BoundaryIntervalControl failed;failed.previousAcceptedS=1;failed.successfulSinceGrowth=2;
    double accepted=-1;
    const auto failure=T::RefineBoundaryInterval(0,1,1e-30,[](double,bool* resolved) {
      *resolved=false;return T::Status::Ok();
    },&accepted,&failed);
    Require(!failure.ok()&&accepted==-1&&failed.previousAcceptedS==1&&failed.successfulSinceGrowth==2,
        "failed boundary refinement published interval/history state");
  });
  return failures?1:0;
}
