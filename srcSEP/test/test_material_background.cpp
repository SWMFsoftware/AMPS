#include "sep_material_background.h"
#include "reduced_material_background_adapter.h"
#include "swcme1d_adapter.h"
#include "swcme1d_material_adapter.h"
#include "swcme_material_upstream.hpp"
#include "provider.h"
#include <sep_coronal_cme/constants.h>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <set>
#include <stdexcept>

namespace M=SEP::Background::Material;
namespace RS=SEP::ReducedShock;
namespace SF=SEP::CoronaSwcme::ShockFront;
namespace C=SEP::CoronalCME::Constants;
using SEP::Core::StatusCode;

// Test-only measurements report actual assertions and reference residuals.
// They never select physics, change the existing comparison bound, or replace
// physical event time with wall time. Reset them for each isolated case.
std::size_t assertionCount=0,nearCount=0;
double maximumRelativeError=0;
constexpr double referenceRelativeTolerance=2e-12;
void Require(bool condition,const std::string& why) {
  ++assertionCount;if(!condition)throw std::runtime_error(why);
}
void Near(double a,double b,double relative=referenceRelativeTolerance) {
  ++nearCount;
  const double residual=std::isfinite(a)&&std::isfinite(b)?
      std::abs(a-b)/std::max(1e-30,std::abs(b)):std::numeric_limits<double>::infinity();
  maximumRelativeError=std::max(maximumRelativeError,residual);
  Require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=relative*std::max(1e-30,std::abs(b)),
          "independent SI reference differs");
}
M::Descriptor SoftwareDescriptor() {
  // This deliberately manufactured SI record exercises software contracts,
  // not a physical CME, source calibration or observational validation.
  M::Descriptor d;d.profile=M::Profile::ReducedParkerV1;d.revision=1;
  d.capabilities=M::MaterialRequirements;d.backgroundGeneration=3;
  d.geometryGeneration=8;d.geometryGenerationPresent=true;d.snapshotGeneration=9;
  d.eventIdentity="software-event";d.physicalConfigurationIdentity="software-physics";
  d.nativeConfigurationIdentity="software-native";d.materialIdentity="software-line";
  d.compositionIdentity="software-composition";d.eosIdentity="software-EOS";
  d.frame="HCI";d.units="SI";d.densityBasis=M::DensityBasis::CoreNonrelativisticFluid;
  d.beginS=0;d.endS=600;return d;
}
struct SoftwareOwner final:M::Owner {
  M::Descriptor descriptor=SoftwareDescriptor();
  mutable unsigned samples=0,moves=0;
  int fault=0;
  const M::Descriptor& Identity() const override {return descriptor;}
  SEP::Core::Result<M::State> Sample(const M::Query& q) const override {
    ++samples;M::State s;s.eventIdentity=descriptor.eventIdentity;
    s.physicalConfigurationIdentity=descriptor.physicalConfigurationIdentity;
    s.backgroundGeneration=descriptor.backgroundGeneration;s.geometryGeneration=descriptor.geometryGeneration;
    s.positionM=q.positionM;s.epochS=q.epochS;s.region=q.region;s.side=q.side;
    s.interfaceIdentity=q.interfaceIdentity;s.velocityMPerS={{2,3,4}};s.magneticFieldT={{0,0,2}};
    s.massDensityKgM3=4;s.pressurePa=6;s.gammaAdiabatic=5.0/3;
    s.species=std::vector<M::Species>{{"electron",4,10},{"proton",2,20},{"alpha",1,30}};
    // All nine entries are distinct, so scalar divergence or a transposed
    // tensor cannot accidentally satisfy the component-by-component oracle.
    s.velocityGradientPerS=M::Gradient{{{1,2,3,4,5,6,7,8,9}},"software-branch",true};
    s.magneticGradientTPerM=M::Gradient{{{9,8,7,6,5,4,3,2,1}},"software-branch",true};
    if(fault==1)s.eventIdentity="foreign";
    if(fault==2)s.epochS+=1;
    if(fault==3)++s.backgroundGeneration;
    if(fault==4)s.positionM[0]+=1;
    if(fault==5)s.region=M::Region::Sheath;
    if(fault==6)s.side=M::Side::DownstreamLimit;
    if(fault==7)s.physicalConfigurationIdentity="foreign-physics";
    if(fault==8)++s.geometryGeneration;
    if(fault==9)s.velocityGradientPerS.reset();
    if(fault==10)s.pressurePa=std::numeric_limits<double>::quiet_NaN();
    if(fault==11)s.species->push_back(s.species->front());
    if(fault==12)s.velocityGradientPerS->values[4]=std::numeric_limits<double>::quiet_NaN();
    if(fault==13)s.velocityGradientPerS->branchPreserved=false;
    if(fault==14)s.species->at(0).temperatureK.reset();
    if(fault==15)s.massDensityKgM3.reset();
    if(fault==16)s.species.reset();
    // Exercise an owner that forgot to populate mandatory numerical fields,
    // independently of the optional-data tests. Missing values must never
    // become a physical stationary flow, magnetic null or zero tensor.
    if(fault==17)s.velocityMPerS=M::State{}.velocityMPerS;
    if(fault==18)s.magneticFieldT=M::State{}.magneticFieldT;
    if(fault==19)s.velocityGradientPerS->values=M::Gradient{}.values;
    return SEP::Core::Result<M::State>::Success(std::move(s));
  }
  SEP::Core::Result<M::Motion> Move(const M::Vector& p,double begin,double end) const override {
    ++moves;M::Motion m;m.positionM=p;m.positionM[0]+=2*(end-begin);m.velocityMPerS={{2,0,0}};
    return SEP::Core::Result<M::Motion>::Success(m);
  }
};
M::Query Query(double time=300) {
  M::Query q;q.positionM={{1,2,3}};q.epochS=time;q.region=M::Region::ParkerExterior;return q;
}
std::shared_ptr<const M::Interval> Capture(const std::shared_ptr<SoftwareOwner>& owner) {
  const auto value=M::Interval::Create(owner);Require(value.ok(),value.status.message);return value.value;
}
std::shared_ptr<SF::Provider> Direct(const std::string& path) {
  const auto read=[](const std::string& p) {
    std::ifstream in(p);return std::string(std::istreambuf_iterator<char>(in),{});
  };
  const auto configuration=SF::ResolveConfiguration(read(path),[&](const std::string& relative) {
    const auto position=path.find_last_of('/');
    return SEP::Core::Result<std::string>::Success(read(path.substr(0,position+1)+relative));
  });
  Require(configuration.ok(),configuration.status.message);
  const auto provider=SF::Provider::Create(configuration.value);Require(provider.ok(),provider.status.message);
  return provider.value;
}
RS::IntervalView CoreView(const std::string& event,double end=1800) {
  std::string error;Require(RS::Configure(event,&error),error);Require(RS::Prepare(0,&error),error);
  const auto interval=RS::CapturePhysicalStep(0,end,RS::EventIdentity());
  Require(interval.ok(),interval.status.message);return interval.value;
}
const RS::GeometryControls controls{0.001,30,1e8,100000};
std::shared_ptr<const M::Interval> CoreInterval(const RS::IntervalView& view) {
  const auto interval=RS::MaterialAdapter::Capture(view,8,9,"native-core-test",controls);
  Require(interval.ok(),interval.status.message);return interval.value;
}

int main(int argc,char** argv) {
  if(argc!=3)return 2;
  unsigned failures=0;
  const auto run=[&](const char* id,const std::function<void()>& body) {
    assertionCount=nearCount=0;maximumRelativeError=0;
    const auto started=std::chrono::steady_clock::now();
    try {body();std::cout<<'['<<id<<"] PASS\n";}
    catch(const std::exception& e){++failures;std::cout<<'['<<id<<"] FAIL "<<e.what()<<'\n';}
    const double elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    std::cout<<std::setprecision(17)<<'['<<id<<"] METRICS assertions="<<assertionCount
        <<" near_calls="<<nearCount<<" maximum_relative_error="<<maximumRelativeError
        <<" relative_tolerance="<<referenceRelativeTolerance<<" seconds="<<elapsed<<'\n';
  };
  run("MBC01",[] {
    const char* expected[]={"reduced-sse-parker-exterior-v1","swcme1d-shock-only-legacy-v1",
      "swcme1d-full-icme-diagnostic-v1","swcme3d-sphere-shock-only-legacy-v1",
      "swcme3d-ellipsoid-shock-only-legacy-v1","swcme3d-sse-shock-only-legacy-v1",
      "swcme3d-sphere-full-icme-diagnostic-v1","swcme3d-ellipsoid-full-icme-diagnostic-v1",
      "swcme3d-sse-full-icme-diagnostic-v1","coronal-bundle-snapshot-v1",
      "coronal-ellipsoid-live-planned-v1","full-corona-cross-ray-planned-v1",
      "swcme1d-adiabatic-outer-wind-v1"};
    constexpr unsigned count=sizeof(expected)/sizeof(expected[0]);
    std::set<std::string> keys;unsigned admitted=0;
    for(unsigned i=1;i<=count;++i) {
      const auto profile=static_cast<M::Profile>(i);const auto info=M::Lookup(profile);
      Require(info.ok()&&keys.insert(info.value.key).second,"missing/duplicate profile identity");
      Require(std::string(info.value.key)==expected[i-1],"frozen profile name/serialized ID changed");
      if(M::Require(profile,M::Use::MaterialTransport,M::Region::ParkerExterior).ok())++admitted;
    }
    Require(keys.size()==count&&admitted==1,"profile family admission broadened");
    Require(!M::Lookup(static_cast<M::Profile>(999)).ok(),"unknown profile defaulted");
  });
  run("MBC02",[] {
    // Production selection comes from the actual frozen canonical mode, not
    // an application family label or a manually asserted capability flag.
    SEP::SW1DAdapter::ConfigurationRequest request;
    Require(SEP::SW1DAdapter::Configure(request).ok(),"canonical default configure failed");
    Require(SEP::SW1DAdapter::MaterialProfile()==M::Profile::Swcme1dShockOnlyV1,"SHOCK_ONLY profile mapping changed");
    request.command_line_assignments={{"shock.region_mode","full_icme","MBC02",0},
        {"shock.acceleration_mode","resolved_compression","MBC02",0}};
    Require(SEP::SW1DAdapter::Configure(request).ok(),"canonical diagnostic configure failed");
    Require(SEP::SW1DAdapter::MaterialProfile()==M::Profile::Swcme1dFullIcmeDiagnosticV1,"FULL_ICME profile mapping changed");
    for(unsigned i=2;i<=13;++i) {
      Require(!M::Require(static_cast<M::Profile>(i),M::Use::MaterialTransport,M::Region::ParkerExterior).ok(),"legacy/planned profile admitted");
      Require(!M::Require(static_cast<M::Profile>(i),M::Use::PhysicalDownstream).ok(),"local diagnostic became volume support");
    }
    Require(!M::Require(M::Profile::ReducedParkerV1,M::Use::MaterialTransport,M::Region::CoronalOpen).ok(),"PFSS inherited Parker admission");
    Require(!M::Require(M::Profile::ReducedParkerV1,M::Use::SnapshotExchange,M::Region::Unspecified,2).ok(),"unknown revision admitted");
    Require(!M::Require(M::Profile::ReducedParkerV1,static_cast<M::Use>(99)).ok(),"unknown operation admitted");

    // Exercise the real typed Stage 12 owner separately from both legacy
    // modes. Numerical physics has an independent oracle in the SWCME core
    // suite. Here exact comparisons detect changed units, swapped Cartesian
    // derivatives, foreign lineage, and invented thermodynamic export.
    namespace W=swcme::material1d;
    const double au=swcme::constants::AU_M;
    W::Config wind;
    wind.inner_radius_m=.30*au;wind.outer_radius_m=1.20*au;wind.reference_radius_m=au;
    wind.wind_m_s=400000;wind.reference_electron_density_m3=6e6;
    wind.reference_electron_temperature_k=120000;wind.reference_proton_temperature_k=120000;
    wind.reference_field_t=5e-9;wind.reference_sin_colatitude=1;wind.polarity=-1;
    wind.rotation_axis={{0,0,1}};wind.rotation_rad_s=2.86533e-6;
    wind.solar_gravity_m3_s2=C::kSolarGravitationalParameterM3PerS2;wind.gamma=5.0/3.0;
    W::EventConfig config;
    config.identity="MBC02 controlled SWCME1D material event";config.begin_s=0;config.end_s=600;
    config.front.mode=swcme::kinematics::Mode::Ballistic;
    config.front.r0_m=.9*au;config.front.V0_m_s=1100000;config.front.Vsw_m_s=wind.wind_m_s;
    const auto event=std::make_shared<const W::Event>(wind,config);
    auto foreign=config;foreign.identity.append("\0foreign",8);
    Require(W::Event(wind,foreign).fingerprint()!=event->fingerprint(),"typed event fingerprint truncated its identity");
    Require(!SEP::SW1DAdapter::CaptureUpstreamMaterial(nullptr,4,0).ok(),"missing canonical event captured");
    Require(!SEP::SW1DAdapter::CaptureUpstreamMaterial(event,0,0).ok(),"unpublished background generation captured");
    const auto captured=SEP::SW1DAdapter::CaptureUpstreamMaterial(event,4,0);
    Require(captured.ok(),captured.status.message);
    const auto& interval=*captured.value;
    const auto& descriptor=interval.Identity();
    Require(descriptor.profile==M::Profile::Swcme1dAdiabaticOuterWindV1&&
        descriptor.eventIdentity==config.identity&&descriptor.backgroundGeneration==4&&
        descriptor.geometryGenerationPresent&&descriptor.geometryGeneration==0&&
        descriptor.beginS==config.begin_s&&descriptor.endS==config.end_s,
        "controlled upstream owner lost its actual event or geometry lineage");
    auto q=Query();q.positionM={{.6*au,.5*au,.4*au}};q.required=M::MaterialRequirements;
    W::Sample direct;
    Require(event->sample(q.positionM,q.epochS,direct).ok(),"control point left canonical coverage");
    const auto sample=interval.Sample(q);Require(sample.ok(),sample.status.message);
    const auto& s=sample.value;
    Require(s.velocityMPerS==direct.velocity_m_s&&s.magneticFieldT==direct.magnetic_field_t&&
        s.massDensityKgM3==direct.mass_density_kg_m3&&s.pressurePa==direct.pressure_pa&&
        s.gammaAdiabatic==wind.gamma&&s.velocityGradientPerS->values==direct.velocity_gradient_s_inv&&
        s.magneticGradientTPerM->values==direct.magnetic_gradient_t_m_inv,
        "neutral Cartesian state differs from the canonical upstream state");
    Require(s.species&&s.species->size()==3&&s.species->at(0).identity=="electron"&&
        s.species->at(0).numberDensityM3==direct.electron_density_m3&&
        s.species->at(0).temperatureK==direct.electron_temperature_k&&
        s.species->at(1).identity=="proton"&&s.species->at(1).numberDensityM3==direct.proton_density_m3&&
        s.species->at(1).temperatureK==direct.proton_temperature_k&&s.species->at(2).identity=="alpha"&&
        s.species->at(2).numberDensityM3==0&&!s.species->at(2).temperatureK,
        "upstream composition or absence of alpha temperature was lost");
    W::Motion moved;
    Require(event->advance_fixed(q.positionM,300,600,300,moved).ok(),"canonical fixed motion failed");
    const auto neutralMotion=interval.Move(q.positionM,300,600);
    Require(neutralMotion.ok()&&neutralMotion.value.positionM==moved.position_m&&
        neutralMotion.value.velocityMPerS==moved.velocity_m_s&&!neutralMotion.value.deformation,
        "neutral trajectory differs from canonical constant-step motion");
    W::Shock shock;
    Require(event->shock_at(q.positionM,q.epochS,shock).ok()&&shock.jump.has_shock,"fast control front failed");
    const auto surface=interval.Shock(shock.position_m,q.epochS);Require(surface.ok(),surface.status.message);
    Require(surface.value.upstream&&surface.value.downstream&&
        surface.value.upstream->pressurePa==shock.ambient.pressure_pa&&
        surface.value.downstream->pressurePa==shock.jump.downstream.pressure_Pa&&
        surface.value.downstream->gammaAdiabatic==wind.gamma&&
        !surface.value.downstream->species&&!surface.value.downstream->velocityGradientPerS,
        "canonical RH limit lost bulk state or invented downstream species/gradient");
    Require(!interval.Shock(q.positionM,q.epochS).ok(),"off-front ambient was labeled a shock");
    q.region=M::Region::Sheath;
    Require(!interval.Sample(q).ok(),"local RH state became sheath volume");
    q.region=M::Region::ParkerExterior;q.epochS=std::nextafter(600.,1000.);
    Require(!interval.Sample(q).ok(),"neutral export extrapolated the physical interval");
    q.epochS=300;q.positionM={{.25*au,0,0}};
    Require(!interval.Sample(q).ok()&&!interval.Move(q.positionM,300,600).ok(),"neutral export crossed radial coverage");
    const auto native=interval.RequireStep(300,600,300,M::Region::ParkerExterior);
    Require(!native.ok()&&native.code==StatusCode::UnsupportedCapability,"upstream-only build admitted native motion");
    namespace B=SEP::Background;
    const B::BackgroundSnapshot snapshot(B::Provider::ReducedShock,B::Ownership::ModelOwned,
        0,0,600,4,"unrelated-native-owner","MBC02 admission control");
    const auto binding=interval.BindSnapshot(snapshot);
    Require(!binding.ok()&&binding.code==StatusCode::UnsupportedCapability,"upstream owner bound a native snapshot");
  });
  run("MBC03",[] {
    const auto d=SoftwareDescriptor();const auto bytes=M::Encode(d);Require(bytes.ok(),bytes.status.message);
    const auto restored=M::Decode(bytes.value);
    Require(restored.ok()&&M::Payload(restored.value)==M::Payload(d),"descriptor lost exact identity/precision");
    Require(!M::Interval::Create(nullptr).ok(),"metadata reconstructed a missing physical owner");
    // Generation zero is PIC's real registered initial material state. The
    // explicit presence bit must survive serialization; absent lineage still
    // rejects, without renumbering the container or changing any timestep.
    auto initial=d;initial.geometryGeneration=0;
    Require(M::Validate(initial).ok(),"registered initial generation zero rejected");
    const auto initialBytes=M::Encode(initial);Require(initialBytes.ok(),initialBytes.status.message);
    const auto initialRestored=M::Decode(initialBytes.value);
    Require(initialRestored.ok()&&initialRestored.value.geometryGenerationPresent&&initialRestored.value.geometryGeneration==0,
        "initial generation presence lost in metadata codec");
  });
  run("MBC04",[] {
    const auto bytes=M::Encode(SoftwareDescriptor()).value;
    for(auto damaged:{bytes+"foreign",bytes.substr(0,bytes.size()/2),std::string("SEP_MATERIAL_DESCRIPTOR 2\n")})Require(!M::Decode(damaged).ok(),"corrupt/foreign descriptor decoded");
    auto changed=bytes;changed[changed.find("software-event")]='X';Require(!M::Decode(changed).ok(),"changed identity passed digest");
    auto d=SoftwareDescriptor();d.capabilities|=M::Bit(M::Capability::DownstreamVolume);
    Require(!M::Encode(d).ok(),"forged capability encoded");
  });
  run("MBC05",[] {
    for(unsigned fault=0;fault<9;++fault) {
      auto owner=std::make_shared<SoftwareOwner>();
      if(fault==0)owner->descriptor.frame="corotating";
      if(fault==1)owner->descriptor.units="km";
      if(fault==2)owner->descriptor.densityBasis=M::DensityBasis::Unspecified;
      if(fault==3)owner->descriptor.materialIdentity.clear();
      if(fault==4)owner->descriptor.eosIdentity.clear();
      if(fault==5)owner->descriptor.compositionIdentity.clear();
      if(fault==6)owner->descriptor.backgroundGeneration=0;
      if(fault==7)owner->descriptor.snapshotGeneration=0;
      if(fault==8)owner->descriptor.geometryGenerationPresent=false;
      Require(!M::Interval::Create(owner).ok(),"incomplete descriptor accepted");
    }
  });
  run("MBC06",[] {
    auto owner=std::make_shared<SoftwareOwner>();const auto interval=Capture(owner);
    Require(interval->Sample(Query(0)).ok()&&interval->Sample(Query(600)).ok(),"closed interval endpoint lost");
    const auto calls=owner->samples;
    Require(!interval->Sample(Query(std::nextafter(600.,1000.))).ok()&&!interval->Sample(Query(-1)).ok(),"query extrapolated");
    Require(owner->samples==calls,"uncovered query reached physical owner");
  });
  run("MBC07",[] {
    auto owner=std::make_shared<SoftwareOwner>();owner->fault=9;
    auto q=Query();q.required|=M::Bit(M::Capability::VelocityGradient);
    const auto result=Capture(owner)->Sample(q);
    Require(!result.ok()&&result.status.code==StatusCode::UnsupportedCapability,"missing full gradient became div(U)/zero");
    owner->fault=16;q.required=M::Bit(M::Capability::Species);
    Require(!Capture(owner)->Sample(q).ok(),"missing species became pure hydrogen");
    owner->fault=14;q.required=M::Bit(M::Capability::Species)|M::Bit(M::Capability::Thermodynamics);
    Require(!Capture(owner)->Sample(q).ok(),"missing populated temperature became a physical zero");
  });
  run("MBC08",[] {
    for(int fault:{10,11,12,13,15,17,18,19}) {
      auto owner=std::make_shared<SoftwareOwner>();owner->fault=fault;
      Require(!Capture(owner)->Sample(Query()).ok(),"invalid or unpopulated primitive/derivative survived");
    }
  });
  run("MBC09",[] {
    for(int fault=1;fault<=8;++fault) {
      auto owner=std::make_shared<SoftwareOwner>();owner->fault=fault;
      const auto result=Capture(owner)->Sample(Query());
      Require(!result.ok()&&result.status.code==StatusCode::DataIntegrityFailure,"foreign sample realization accepted");
    }
  });
  run("MBC10",[] {
    auto owner=std::make_shared<SoftwareOwner>();const auto interval=Capture(owner);
    auto q=Query();q.required=M::Bit(M::Capability::VelocityGradient)|M::Bit(M::Capability::MagneticGradient);
    const auto state=interval->Sample(q);Require(state.ok(),state.status.message);
    for(unsigned i=0;i<9;++i)Require(state.value.velocityGradientPerS->values[i]==i+1&&
        state.value.magneticGradientTPerM->values[i]==9-i,"tensor components were transposed/replaced by scalar divergence");
    q.side=M::Side::UpstreamLimit;
    Require(!interval->Sample(q).ok(),"one-sided query lost interface identity");
    q.interfaceIdentity="explicit-manufactured-interface";Require(interval->Sample(q).ok(),"explicit one-sided software fixture rejected");
    q.side=static_cast<M::Side>(99);Require(!interval->Sample(q).ok(),"unknown side accepted");
    q=Query();q.region=static_cast<M::Region>(99);Require(!interval->Sample(q).ok(),"unknown region accepted");
  });
  run("MBC11",[] {
    const auto interval=Capture(std::make_shared<SoftwareOwner>());
    Require(interval->RequireStep(0,300,300,M::Region::ParkerExterior).ok(),"full fixed step rejected");
    for(double end:{299.,301.,601.})Require(!interval->RequireStep(0,end,300,M::Region::ParkerExterior).ok(),"shortened/uncovered AMPS dt accepted");
    Require(!interval->RequireStep(0,300,0,M::Region::ParkerExterior).ok(),"missing fixed dt accepted");
  });
  run("MBC12",[] {
    namespace B=SEP::Background;auto& store=B::SnapshotStore::Instance();store.ResetForTests();
    const auto interval=Capture(std::make_shared<SoftwareOwner>());
    const B::BackgroundSnapshot good(B::Provider::ReducedShock,B::Ownership::ModelOwned,0,0,600,9,"software-native","software contract");
    Require(interval->BindSnapshot(good).ok(),"distinct geometry/cache/native generations conflated");
    auto sparse=std::make_shared<SoftwareOwner>();sparse->descriptor.capabilities=M::Bit(M::Capability::Primitives);
    Require(!Capture(sparse)->BindSnapshot(good).ok(),"profile label replaced required native material capabilities");
    store.Publish(good);const auto saved=store.Current();
    const B::BackgroundSnapshot foreign(B::Provider::Analytic,B::Ownership::ModelOwned,0,0,600,9,"software-native","foreign provider");
    Require(!interval->BindSnapshot(foreign).ok()&&store.Current()==saved,"foreign provider touched existing publication");
    {auto reader=store.BeginParticleRead(0);bool rejected=false;
      try{store.Publish(good);}catch(const std::exception&){rejected=true;}
      Require(rejected&&store.Current()==saved,"native particle read phase lost its publication guard");}
    store.ResetForTests();
  });
  run("MBC13",[] {
    auto owner=std::make_shared<SoftwareOwner>();const auto interval=Capture(owner);
    owner->descriptor.physicalConfigurationIdentity="changed-after-capture";
    const auto result=interval->Sample(Query());
    Require(!result.ok()&&result.status.code==StatusCode::DataIntegrityFailure&&owner->samples==0,"mutated owner rebound captured realization");
  });
  run("MBC14",[] {
    auto owner=std::make_shared<SoftwareOwner>();const auto interval=Capture(owner);
    Require(!interval->Inverse({{1,2,3}},0).ok()&&!interval->Contact({{1,2,3}},0).ok()&&!interval->Shock({{1,2,3}},0).ok(),
        "absent inverse/contact/front acquired a default");
    auto q=Query();q.required=M::Bit(M::Capability::DownstreamVolume);
    Require(!interval->Sample(q).ok()&&owner->samples==0,"downstream capability reached unsupported owner");
  });
  run("MBC15",[&] {
    const auto view=CoreView(argv[1]);const auto interval=CoreInterval(view);const auto direct=Direct(argv[1]);
    M::Query q;q.positionM={{0,0,20*C::kSolarRadiusM}};q.epochS=600;q.region=M::Region::ParkerExterior;q.required=M::MaterialRequirements&~M::Bit(M::Capability::Trajectory)&~M::Bit(M::Capability::DenseTime);
    const auto state=interval->Sample(q);Require(state.ok(),state.status.message);
    const auto expected=direct->QueryAmbient({q.positionM[0],q.positionM[1],q.positionM[2]},q.epochS,view.Generation());
    Require(expected.ok(),expected.status.message);
    Require(*state.value.massDensityKgM3==expected.value.primitive.plasma.massDensityKgM3&&
        *state.value.pressurePa==expected.value.primitive.plasma.pressurePa,"neutral export changed canonical SI primitives");
    Require(state.value.velocityGradientPerS->values==expected.value.gradientU&&
        state.value.magneticGradientTPerM->values==expected.value.gradientB,"neutral tensor export differs from actual core");
    // Compare every inertial Cartesian component with the independent core
    // owner at the identical SI point/time. A hidden rotation, plasma-frame
    // subtraction or unit scaling must not pass merely because rho and p do.
    const auto& primitive=expected.value.primitive;
    Require(state.value.velocityMPerS==M::Vector{{primitive.velocityMPerS.x,primitive.velocityMPerS.y,primitive.velocityMPerS.z}}&&
        state.value.magneticFieldT==M::Vector{{primitive.magneticFieldT.x,primitive.magneticFieldT.y,primitive.magneticFieldT.z}},
        "neutral Cartesian frame or SI U/B components differ from canonical owner");
    const auto& species=*state.value.species;
    Require(species[0].identity=="electron"&&species[1].identity=="proton"&&species[2].identity=="alpha",
        "explicit species identities differ from the declared composition");
    Near(*species[0].temperatureK,primitive.electronTemperatureK);
    Near(*species[1].temperatureK,primitive.protonTemperatureK);
    Near(*species[2].temperatureK,direct->Event().ambient.composition.alphaTemperatureK);
  });
  run("MBC16",[&] {
    const auto old=CoreInterval(CoreView(argv[1]));auto q=Query();q.positionM={{0,0,20*C::kSolarRadiusM}};
    const auto before=old->Sample(q);Require(before.ok(),before.status.message);
    const auto other=CoreInterval(CoreView(argv[2]));
    Require(other->Identity().eventIdentity!=old->Identity().eventIdentity,"different physical composition did not change identity");
    const auto after=old->Sample(q);Require(after.ok()&&M::Payload(old->Identity())!=M::Payload(other->Identity())&&
        after.value.species->at(2).numberDensityM3==before.value.species->at(2).numberDensityM3,"reconfiguration changed an in-flight physical owner");
  });
  run("MBC17",[&] {
    const auto interval=CoreInterval(CoreView(argv[2]));auto q=Query();q.positionM={{0,0,20*C::kSolarRadiusM}};
    q.required=M::Bit(M::Capability::Species)|M::Bit(M::Capability::Thermodynamics);
    const auto state=interval->Sample(q);Require(state.ok(),state.status.message);const auto& s=*state.value.species;
    const double ne=s[0].numberDensityM3,np=s[1].numberDensityM3,na=s[2].numberDensityM3;
    Near(ne,np+2*na);Near(na,0.1*np);Require(ne!=np,"multispecies fixture did not distinguish density basis");
    Near(*state.value.massDensityKgM3,np*C::kProtonMassKg+na*C::kAlphaMassKg);
    Near(*state.value.pressurePa,C::kBoltzmannJPerK*(ne**s[0].temperatureK+np**s[1].temperatureK+na**s[2].temperatureK));
    Require(interval->Identity().densityBasis==M::DensityBasis::CoreNonrelativisticFluid,"core density was relabeled as exact relativistic proper density");
  });
  run("MBC18",[&] {
    const auto view=CoreView(argv[1]);const auto interval=CoreInterval(view);
    const M::Vector x{{0,0,20*C::kSolarRadiusM}};const auto motion=interval->Move(x,0,300);
    const auto reference=RS::MoveVertex(view,x,0,300,controls);
    Require(motion.ok()&&reference.ok()&&motion.value.positionM==reference.value,"trajectory changed the existing numerical method");
    // The configured isothermal Parker wind accelerates with radius. An
    // r0+U(r0)*dt oracle would incorrectly freeze that real wind. Independently
    // solve integral(r0,r1) dr/U(r)=300 s by composite Simpson quadrature and
    // bisection, using a separate canonical primitive owner, not the production
    // RK trajectory. Two quadrature resolutions must agree under the same
    // pre-existing relative comparison bound before judging the actual move.
    const auto direct=Direct(argv[1]);
    const auto flightTime=[&](double radius,int panels) {
      const double h=(radius-x[2])/panels;double sum=0;
      for(int i=0;i<=panels;++i) {
        const auto state=direct->Ambient().Evaluate({0,0,x[2]+i*h},0);
        Require(state.ok(),state.status.message);
        const double speed=state.value.velocityMPerS.z;
        Require(std::isfinite(speed)&&speed>0,"independent radial oracle requires outward wind");
        sum+=(i==0||i==panels?1:i%2?4:2)/speed;
      }
      return h*sum/3;
    };
    const auto independentRadius=[&](int panels) {
      double left=x[2],right=x[2]+1e8;
      for(int i=0;flightTime(right,panels)<300;++i) {
        Require(i<16,"independent trajectory bracket failed");right=x[2]+2*(right-x[2]);
      }
      for(int i=0;i<60;++i) {
        const double middle=left+(right-left)/2;
        if(middle==left||middle==right)break;
        if(flightTime(middle,panels)<300)left=middle;else right=middle;
      }
      return left+(right-left)/2;
    };
    const double coarse=independentRadius(128),fine=independentRadius(256);
    Near(coarse,fine);Near(motion.value.positionM[2],fine);
    Require(motion.value.positionM[0]==0&&motion.value.positionM[1]==0,
        "polar radial material trajectory acquired transverse motion");
    Require(!interval->Move(x,0,1801).ok(),"trajectory extrapolated beyond captured interval");
    // The existing event includes a diagnostic PFSS corona inside its source
    // surface. Test two ACTUAL valid coronal primitives (open polar and closed
    // equatorial), not a manufactured invalid coordinate. The neutral Parker
    // state and trajectory APIs must agree that material evolution there is
    // unavailable, including a zero-duration query that executes no RK stages.
    for(const M::Vector point:{M::Vector{{0,0,2*C::kSolarRadiusM}},
                               M::Vector{{2*C::kSolarRadiusM,0,0}}}) {
      const auto primitive=view.PrimitiveAt(point,0);
      Require(primitive.ok()&&!primitive.value.parkerExterior,
          "region control must exercise an actual unqualified coronal primitive");
      M::Query query;query.positionM=point;query.epochS=0;query.region=M::Region::ParkerExterior;
      Require(!interval->Sample(query).ok(),"PFSS sample acquired Parker material authority");
      const auto coronalMove=interval->Move(point,0,1),zeroMove=interval->Move(point,0,0);
      Require(!coronalMove.ok()&&coronalMove.status.code==StatusCode::UnsupportedCapability,
          "PFSS trajectory acquired Parker material authority or rejected for an unrelated reason");
      Require(!zeroMove.ok()&&zeroMove.status.code==StatusCode::UnsupportedCapability,
          "zero-duration PFSS query acquired material authority or rejected for an unrelated reason");
    }
    // Reverse-time geometry queries are supported by the same RK controller.
    // Start just outside the DECLARED source surface, then independently
    // confirm that the legacy diagnostic trajectory enters the PFSS corona.
    // The material wrapper must reject an internal-stage region crossing;
    // it may not freeze the velocity, clip the step or accept a fallback.
    const M::Vector edge{{0,0,direct->Event().ambient.ambient.sourceSurfaceRadiusM+
        1e-4*C::kSolarRadiusM}};
    const auto edgeState=view.PrimitiveAt(edge,300);
    Require(edgeState.ok()&&edgeState.value.parkerExterior,"reverse control does not start in Parker exterior");
    const auto legacyReverse=RS::MoveVertex(view,edge,300,0,controls);
    Require(legacyReverse.ok(),legacyReverse.status.message);
    const auto coronalEnd=view.PrimitiveAt(legacyReverse.value,0);
    Require(coronalEnd.ok()&&!coronalEnd.value.parkerExterior,
        "reverse control did not actually enter an unqualified coronal region");
    const auto reverse=interval->Move(edge,300,0);
    Require(!reverse.ok()&&reverse.status.code==StatusCode::UnsupportedCapability,
        "reverse-time path crossed an unqualified material region or failed for an unrelated reason");
    // The adapter's public constructor must not allow a valid-looking foreign
    // descriptor to relabel this captured physical owner in the motion API.
    auto foreignDescriptor=interval->Identity();foreignDescriptor.physicalConfigurationIdentity="foreign-physics";
    const auto foreign=M::Interval::Create(std::make_shared<const RS::MaterialAdapter::Owner>(
        view,foreignDescriptor,controls));
    Require(foreign.ok(),foreign.status.message);
    const auto foreignMove=foreign.value->Move(x,0,300);
    Require(!foreignMove.ok()&&foreignMove.status.code==StatusCode::DataIntegrityFailure,
        "foreign material trajectory descriptor rebound a physical owner");
  });
  run("MBC19",[&] {
    const auto view=CoreView(argv[1]);const auto interval=CoreInterval(view);
    const auto point=view.LeadingPointOnRay({{0,0,1}},1800);Require(point.ok(),point.status.message);
    const auto surface=interval->Shock(point.value,1800);Require(surface.ok(),surface.status.message);
    const auto direct=view.FrontAt(point.value,1800);Require(direct.ok(),direct.status.message);
    Require(surface.value.normalSpeedMPerS==direct.value.normalSpeedMPerS&&surface.value.downstream.has_value()==direct.value.downstreamValid,"actual RH surface limits changed");
    Require(!M::Require(interval->Identity().profile,M::Use::PhysicalDownstream).ok(),"local RH primitive became a downstream volume");
    if(surface.value.downstream)Require(!surface.value.downstream->species&&!surface.value.downstream->velocityGradientPerS,"unavailable downstream heating/derivative was fabricated");
  });
  run("MBC20",[&] {
    namespace B=SEP::Background;auto& store=B::SnapshotStore::Instance();store.ResetForTests();
    const auto view=CoreView(argv[1]);const auto interval=CoreInterval(view);
    const B::BackgroundSnapshot good(B::Provider::ReducedShock,B::Ownership::ModelOwned,0,0,1800,9,"native-core-test","actual core contract");
    Require(interval->BindSnapshot(good).ok(),"actual neutral/native binding failed");store.Publish(good);const auto before=store.Current();
    const auto foreignRun=RS::MaterialAdapter::Capture(view,8,9,"foreign-native-configuration",controls);
    Require(foreignRun.ok()&&!foreignRun.value->BindSnapshot(*store.Current()).ok()&&store.Current()==before,
        "reused actual native snapshot accepted a foreign configuration");
    const auto bytes=M::Encode(interval->Identity());Require(bytes.ok(),bytes.status.message);
    const auto decoded=M::Decode(bytes.value);Require(decoded.ok(),"actual descriptor did not restore exactly");
    auto corrupted=decoded.value;corrupted.eventIdentity="foreign-event";
    const auto rejected=view.MaterialStateAt(corrupted,Query());
    Require(!rejected.ok()&&store.Current()==before,"late foreign realization mutated native publication");
    Require(interval->Sample(M::Query{{{0,0,20*C::kSolarRadiusM}},300,M::Region::ParkerExterior,M::Side::Smooth,"",M::Bit(M::Capability::Primitives)}).ok(),"rejection destroyed valid core owner");
    // Stage a newer cadence without publishing it. A low geometry generation
    // deliberately makes the real cache generation control the native join,
    // exposing accidental attachment to the old live cache. Physical queries
    // remain on the same immutable event and full covered time interval.
    const auto staged=RS::PrepareCadencePublication(1200);Require(staged.ok(),staged.status.message);
    const auto stagedView=staged.value.PhysicalInterval(1200,1500);Require(stagedView.ok(),stagedView.status.message);
    Require(stagedView.value.Generation()>view.Generation(),"staged cadence did not change cache generation");
    const auto next=RS::MaterialAdapter::Capture(stagedView.value,1,staged.value.Generation(),"native-core-test",controls);
    Require(next.ok(),next.status.message);
    const B::BackgroundSnapshot candidate(B::Provider::ReducedShock,B::Ownership::ModelOwned,
        1200,1200,1500,staged.value.Generation(),"native-core-test","unpublished cadence candidate");
    Require(next.value->BindSnapshot(candidate).ok(),"staged cadence owner/native generation join failed");
    Require(store.Current()==before,"staging a new physical view published it early");
    store.ResetForTests();
  });
  return failures?1:0;
}
