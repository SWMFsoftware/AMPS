// Independent field/jump/gradient oracles and the real provider publication
// contracts. These are portable prerequisites; they do not claim MPI evidence.
#include "sep3d_test_registry.h"
#include "configuration_io.h"
#include "background_factory.h"
#include "background_snapshot.h"
#include "bg_swcme.h"
#include "bg_parker.h"
#include "runtime_adapters.h"
#include "swcme3d_input.hpp"
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
namespace {
namespace R=SEP3D::RuntimeModel;namespace B=SEP3D::Background;namespace C=SEP3D::Core;
using Result=SEP3D::Testing::Result;
Result Finish(bool ok,const std::string& message) { Result r;
  r.status=ok?SEP3D::Testing::Status::Pass:SEP3D::Testing::Status::Fail;r.message=message;return r; }
bool Near(double a,double b,double tolerance=1e-9) {
  // Relative comparisons also work for SI magnetic derivatives far below one;
  // the tiny floor prevents a zero reference from becoming a zero tolerance.
  return std::fabs(a-b)<=tolerance*std::max(1e-30,std::max(std::fabs(a),std::fabs(b)));
}
struct Fixture {
  // Load the real source-free deck and resolve its canonical assignments
  // independently of factory creation. Each callback gets fresh model state;
  // no fixture installs an AMPS mesh or simulates an MPI halo exchange.
  R::RunConfiguration3DOptions options;
  std::shared_ptr<const R::RunConfiguration3D> configuration;
  std::shared_ptr<B::BackgroundProvider> provider;
  swcme::input3d::ResolvedConfiguration canonical;
  C::Status Setup(const char* deck="examples/sep3d_swcme_mesh_background_20rs_1au.in") {
    auto s=R::LoadConfigurationFile(deck,&options);if(!s.ok())return s;
    s=R::RunConfiguration3D::Create(options,&configuration);if(!s.ok())return s;
    std::vector<swcme::input3d::Assignment> assignments;
    for(const auto& raw:options.swcmeAssignments){swcme::input3d::Assignment a;a.key=raw.key;a.value=raw.value;assignments.push_back(a);}
    const auto resolved=swcme::input3d::Resolve(assignments);
    if(!resolved.ok())return C::Status::Error(resolved.status.message);
    canonical=resolved.configuration;
    s=R::CreateBackgroundProvider(*configuration,&provider);if(s.ok())s=provider->Prepare(0);
    return s;
  }
};
Result Configuration() {
  // Exercise the public parser/frozen factory, not just option construction.
  // Mutations probe double-acceleration/source conflicts and output preservation
  // when a canonical shock/background fingerprint no longer agrees.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  if(f.options.background!=R::BackgroundAuthority::Swcme||f.options.source.enabled||
     f.canonical.model.region_mode!=swcme::regions::Mode::FullICME||
     f.canonical.model.shock_acceleration_mode!=swcme::acceleration::Mode::ResolvedCompression)
    return Finish(false,"resolved mesh example selects the wrong physical mode");
  std::ifstream file("examples/sep3d_swcme_mesh_background_20rs_1au.in");
  std::string deck((std::istreambuf_iterator<char>(file)),{});
  auto reject=[&](const std::string& a,const std::string& b) {
    auto text=deck;const auto at=text.find(a);if(at==std::string::npos)return false;text.replace(at,a.size(),b);
    R::RunConfiguration3DOptions unused;return !R::ParseConfigurationText(text,&unused).ok();
  };
  const bool good=reject("shock.acceleration_mode = resolved_compression","shock.acceleration_mode = source") &&
     reject("provider = swcme","provider = analytic-parker") &&
     reject("enabled = false\n", "enabled = true\n");
  auto changed=f.options;changed.swcmeConfigurationFingerprint="foreign";
  std::shared_ptr<const R::RunConfiguration3D> frozen;
  s=R::RunConfiguration3D::Create(changed,&frozen);
  std::shared_ptr<B::BackgroundProvider> original=f.provider;
  if(!s.ok()||R::CreateBackgroundProvider(*frozen,&original).ok()||original!=f.provider)
    return Finish(false,"background factory accepted a different shock identity or modified its output on rejection");
  return Finish(good,"schema-4 FULL_ICME publishes resolved fields; mode/source conflicts and foreign canonical identity are rejected");
}
Result Fields() {
  // Query the canonical RH state as a coupling oracle: this proves the adapter
  // exports exact vector endpoints and total pressure without re-solving MHD
  // with its own formula. Independent shock-physics verification belongs to
  // SWCME's canonical tests, rather than being claimed by this endpoint check.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  swcme3d::Model oracle(f.canonical.model);const auto step=oracle.prepare_step(0);
  const double direction[3]={1,0,0};swcme3d::LocalShockState shock;
  if(!oracle.shock_state_direction(step,direction,shock)||!shock.has_shock)
    return Finish(false,"canonical launch has no MHD shock");
  const auto boundary=swcme::regions::make_boundaries(step.r_sh_m,step.region_config);
  const double inner=step.r_sh_m-0.5*boundary.smooth_shock_width_m;
  const auto downstream=f.provider->Evaluate({inner,0,0});
  const auto upstream=f.provider->Evaluate({step.r_sh_m+boundary.smooth_shock_width_m,0,0});
  const double b2[3]={downstream.B.x,downstream.B.y,downstream.B.z};
  const double u2[3]={downstream.U.x,downstream.U.y,downstream.U.z};
  bool good=downstream.valid&&upstream.valid&&Near(downstream.numberDensityM3,shock.downstream_n_m3)&&
      Near(downstream.pressurePa,shock.downstream.pressure_Pa);
  for(int k=0;k<3;++k)good=good&&Near(b2[k],shock.downstream.magnetic_T[k])&&Near(u2[k],shock.downstream.velocity_m_s[k]);
  const auto ambient=swcme::solarwind::thermodynamic_state(step.common.solar_wind,downstream.numberDensityM3);
  good=good&&Near(downstream.temperatureK,step.common.solar_wind.T_K*downstream.pressurePa/ambient.pressure_Pa)&&
      downstream.temperatureK>upstream.temperatureK && downstream.B.Norm()>upstream.B.Norm();
  // The 1.00..1.05 Rs continuation is explicit, not a swallowed model error.
  const auto innerSun=f.provider->Evaluate({1.02*C::Const::R_sun,0,0});
  const auto handoff1=f.provider->Evaluate({(1.05-1e-8)*C::Const::R_sun,0,0});
  const auto handoff2=f.provider->Evaluate({(1.05+1e-8)*C::Const::R_sun,0,0});
  good=good&&f.provider->Evaluate({C::Const::R_sun,0,0}).valid&&innerSun.valid&&handoff1.valid&&handoff2.valid&&Near(handoff1.B.x,handoff2.B.x,1e-6)&&
      Near(handoff1.numberDensityM3,handoff2.numberDensityM3,1e-6);
  if (!good) {
    std::ostringstream detail;detail.precision(17);
    detail << "RH/inner-shell mismatch: n=" << downstream.numberDensityM3 << '/' << shock.downstream_n_m3
           << " p=" << downstream.pressurePa << '/' << shock.downstream.pressure_Pa
           << " T=" << downstream.temperatureK << '/' << upstream.temperatureK
           << " B=" << downstream.B.Norm() << '/' << upstream.B.Norm()
           << " inner=" << innerSun.valid << ':' << innerSun.status.message
           << " handoff=" << handoff1.valid << ':' << handoff2.valid
           << " Br=" << handoff1.B.x << '/' << handoff2.B.x
           << " ne=" << handoff1.numberDensityM3 << '/' << handoff2.numberDensityM3;
    return Finish(false,detail.str());
  }
  return Finish(true,"mesh fields recover exact vector RH downstream and heated pressure; inner ambient handoff is continuous");
}
Result Evolution() {
  // The 24-Rs cell begins ahead of the launch front and is reached later.
  // Compare immutable snapshots at fixed coordinates so changing U/B proves
  // actual field evolution, not merely a moved front or incremented metadata.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  B::BackgroundSnapshotBuilder builder;
  const std::vector<C::Vec3> positions={{24*C::Const::R_sun,0,0},{C::Const::AU,0,0}};
  std::shared_ptr<const B::BackgroundSnapshot> first,second;
  s=builder.Build(*f.provider,positions,&first);if(!s.ok())return Finish(false,s.message);
  s=f.provider->Prepare(4000);if(s.ok())s=builder.Build(*f.provider,positions,&second);
  if(!s.ok())return Finish(false,s.message);
  const auto& a=first->samples()[0];const auto& b=second->samples()[0];
  const bool good=first->metadata().epochS==0 &&second->metadata().epochS==4000 &&
      second->metadata().generation>first->metadata().generation &&a.generation==1 &&
      b.generation==second->metadata().generation && b.U.x>a.U.x &&b.B.Norm()>a.B.Norm()&&
      first->metadata().configurationFingerprint==second->metadata().configurationFingerprint;
  return Finish(good,"one prepared runtime advances plasma and IMF across a fixed cell while the prior snapshot stays immutable");
}
Result Gradients() {
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  swcme3d::Model oracle(f.canonical.model);const auto step=oracle.prepare_step(0);
  const C::Vec3 at(step.r_sh_m,0,0);const auto sample=f.provider->Evaluate(at);
  if(!sample.valid)return Finish(false,sample.status.message);
  // Independent smaller stencil queries the canonical vectors directly, not
  // the application's derivative helper. Includes the finite shock layer.
  const double h=1e-5*at.Norm();C::Tensor3 gradU,gradB;
  for(int j=0;j<3;++j) {
    C::Vec3 d;j==0?d.x=h:(j==1?d.y=h:d.z=h);
    const auto a=at+d,b=at-d;double x[2]={a.x,b.x},y[2]={a.y,b.y},z[2]={a.z,b.z};
    double n[2],ux[2],uy[2],uz[2],bx[2],by[2],bz[2];
    if(!oracle.evaluate_cartesian_with_B_checked(step,x,y,z,n,ux,uy,uz,bx,by,bz,2).ok())return Finish(false,"gradient oracle failed");
    const double* us[3]={ux,uy,uz};const double* bs[3]={bx,by,bz};
    for(int i=0;i<3;++i){gradU(i,j)=(us[i][0]-us[i][1])/(2*h);gradB(i,j)=(bs[i][0]-bs[i][1])/(2*h);}
  }
  const auto c=f.provider->Capabilities();bool good=c.hasGradB&&c.hasGradU&&!c.hasAnalyticGradB&&!c.hasAnalyticDivU;
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    good=good&&Near(gradU(i,j),sample.gradU(i,j),5e-4)&&Near(gradB(i,j),sample.gradB(i,j),5e-4);
  good=good&&sample.divU<0&&Near(sample.divU,gradU.Trace(),5e-4)&&
      Near(sample.fieldAlignedStrain,gradU.DoubleContract(sample.bHat,sample.bHat),5e-4);
  // SHOCK_ONLY is a separate limiting case: the front moves but sampled U is
  // constant radial wind, whose exact Cartesian divergence is 2U/r.
  Fixture ambient;s=ambient.Setup("examples/sep3d_swcme_20rs_1au.in");
  if(!s.ok())return Finish(false,s.message);
  auto o=ambient.options;o.background=R::BackgroundAuthority::Swcme;
  std::shared_ptr<const R::RunConfiguration3D> config;
  if(!R::RunConfiguration3D::Create(o,&config).ok()||!R::CreateBackgroundProvider(*config,&ambient.provider).ok()||!ambient.provider->Prepare(0).ok())
    return Finish(false,"SHOCK_ONLY mesh provider rejected");
  const auto a=ambient.provider->Evaluate({C::Const::AU,0,0});
  good=good&&Near(a.divU,2*400000/C::Const::AU,1e-6);
  return Finish(good,"full vector gradients resolve compression and match independent Cartesian stencil; SHOCK_ONLY retains exact radial divergence");
}
Result Rejections() {
  // Sentinel outputs and pointer identity expose accidental writes on failure.
  // A mixed valid/inside-Sun batch must preserve the bad point while returning
  // valid neighbours; a mesh snapshot rejects that same candidate atomically.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  const auto metadata=*f.provider->PreparedMetadata();
  if(f.provider->Prepare(-1).ok()||f.provider->Prepare(std::numeric_limits<double>::quiet_NaN()).ok()||
      f.provider->Prepare(432001).ok()||f.provider->PreparedMetadata()->generation!=metadata.generation)
    return Finish(false,"failed preparation changed the published state");
  double x[3]={C::Const::AU,0,0.3*C::Const::AU},y[3]={0,0,0},z[3]={0,0,0};
  B::BackgroundSample outputs[3];outputs[1].pressurePa=777;C::Status statuses[3];
  if(f.provider->EvaluateBatchDetailed(x,y,z,3,outputs,statuses).ok()||
      !statuses[0].ok()||statuses[1].ok()||!statuses[2].ok()||outputs[1].pressurePa!=777)
    return Finish(false,"detailed batch lost valid points or changed a rejected output");
  B::BackgroundSnapshotBuilder builder;std::shared_ptr<const B::BackgroundSnapshot> snapshot;
  if(!builder.Build(*f.provider,{{C::Const::AU,0,0}},&snapshot).ok())return Finish(false,"baseline build failed");
  const auto prior=snapshot;
  if(builder.Build(*f.provider,{{0,0,0}},&snapshot).ok()||snapshot!=prior)return Finish(false,"bad candidate replaced active snapshot");
  // An empty local snapshot carries the prepared metadata and is necessary
  // for a rank with no owner cells. Global coverage remains a native MPI gate.
  if(!builder.Build(*f.provider,{},&snapshot).ok()||!snapshot->samples().empty()||snapshot->metadata().generation!=metadata.generation)
    return Finish(false,"zero-owned-cell MPI rank cannot join publication");
  B::ParkerConfiguration parker;B::AnalyticParkerProvider restored(parker);
  if(!restored.Prepare(0).ok()||!restored.RestorePreparedGeneration(12).ok()||
      !restored.Prepare(60).ok()||restored.PreparedMetadata()->generation!=13)
    return Finish(false,"restored analytic provider reused an older generation on update");
  return Finish(true,"invalid epochs/cells preserve old state; empty owner ranks join; restored Parker advances its generation");
}

Result PolarSnapshots() {
  // The native failure occurred during initial owner-cell sampling, before
  // stepping or halo exchange. Recreate that provider/builder route across the
  // full sphere instead of checking only the equatorial CME apex. Include the
  // two-step native epochs and a later state where the parallel gas-dynamic
  // branch is admissible. Sampling the finite layer/sheath/ejecta also executes
  // every provider's seven-point Cartesian derivative stencil at these angles.
  Fixture f;auto status=f.Setup();if(!status.ok())return Finish(false,status.message);
  swcme3d::Model oracle(f.canonical.model);
  std::shared_ptr<const B::BackgroundSnapshot> snapshot;
  B::BackgroundSnapshotBuilder builder;
  std::size_t checked=0;
  for (double epoch:{0.,60.,120.,4000.}) {
    status=f.provider->Prepare(epoch);if(!status.ok())return Finish(false,status.message);
    const auto step=oracle.prepare_step(epoch);
    const auto boundaries=swcme::regions::make_boundaries(step.r_sh_m,step.region_config);
    std::vector<C::Vec3> points;
    for (double thetaDegrees:{0.,0.01,0.25,0.5,1.,2.,5.,30.,90.,150.,175.,178.,179.,179.5,179.75,179.99,180.}) {
      const double theta=thetaDegrees*C::Const::kPi/180;
      for (int azimuth=0;azimuth<8;++azimuth) {
        const double phi=azimuth*C::Const::kPi/4;
        const C::Vec3 direction(std::sin(theta)*std::cos(phi),
            std::sin(theta)*std::sin(phi),std::cos(theta));
        for (double radius:{step.r_sh_m,(step.r_sh_m+boundaries.R_le_m)/2,
                             (boundaries.R_le_m+boundaries.R_te_m)/2})
          points.push_back(direction*radius);
      }
    }
    status=builder.Build(*f.provider,points,&snapshot);
    if(!status.ok())return Finish(false,"polar snapshot epoch="+std::to_string(epoch)+": "+status.message);
    // Successful validation must expose complete plasma/transport state at
    // the requested epoch, not silently substitute an untagged ambient field.
    for (const auto& sample:snapshot->samples()) {
      if (!sample.valid || snapshot->metadata().epochS!=epoch ||
          sample.generation!=snapshot->metadata().generation || !(sample.pressurePa>0) ||
          !(sample.numberDensityM3>0) || !std::isfinite(sample.divU))
        return Finish(false,"polar snapshot contains incomplete or wrong-epoch state");
      ++checked;
    }
  }
  return Finish(checked==1632,"1632 polar/oblique region samples and their derivative stencils publish at 0/60/120/4000 s");
}

Result ParallelShockLimit() {
  // Independent RH oracle for the switch-on limit, not another call to the
  // production scalar residual. Choose rho, Bn and p directly: M_An=1.7 and
  // low beta put the gas-dynamic root below the downstream Alfven speed, so
  // that root cannot be the fast branch. The parallel fast limit instead has
  // r=M_An^2 and a nonzero downstream transverse field. Exercise polarity and
  // a tangential Galilean boost; no particular transverse azimuth is imposed
  // by the physical RH equations in this degenerate limit.
  const double mu0=4e-7*C::Const::kPi, gamma=5./3., rho=1e-20, Bn=5e-8;
  const double U=1.7*Bn/std::sqrt(mu0*rho), r=1.7*1.7;
  const double p1=0.005*Bn*Bn/mu0;
  const double p2=(gamma*r-(gamma-1))*p1+
      (gamma-1)*rho*U*U*(r-1)*(r-1)/(2*r);
  const double Bt2=2*mu0*(rho*U*U*(1-1/r)+p1-p2);
  for (double polarity:{-1.,1.}) for(double boost:{0.,12345.}) {
    swcme::shock::PrimitiveState up;
    up.rho_kg_m3=rho;up.pressure_Pa=p1;
    up.velocity_m_s={{400000,boost,-2*boost}};up.magnetic_T={{polarity*Bn,0,0}};
    const auto jump=swcme::shock::solve_ideal_mhd_fast_shock(up,{{1,0,0}},400000+U,gamma);
    const auto& down=jump.downstream;
    if (!jump.solver_converged || !jump.has_shock || !jump.evolutionary_fast_branch ||
        !Near(jump.compression,r,1e-9) || !Near(down.pressure_Pa,p2,1e-9) ||
        !Near(down.magnetic_T[0],polarity*Bn,1e-9) ||
        !Near(down.magnetic_T[1]*down.magnetic_T[1]+down.magnetic_T[2]*down.magnetic_T[2],Bt2,1e-9))
      return Finish(false,"parallel RH limit mismatch: "+std::string(swcme::shock::solve_status_name(jump.status)));
    // Verify tangential electric/momentum jumps from supplied/final primitives
    // and independent component formulas, outside production diagnostics.
    const double un2=down.velocity_m_s[0]-(400000+U);
    for(int j=1;j<3;++j) {
      const double du=down.velocity_m_s[j]-up.velocity_m_s[j];
      const double induction=un2*down.magnetic_T[j]-polarity*Bn*du;
      const double momentum=-rho*U*du-polarity*Bn*down.magnetic_T[j]/mu0;
      if (std::fabs(induction)>1e-9*U*std::sqrt(Bn*Bn+Bt2) ||
          std::fabs(momentum)>1e-9*rho*U*U)
        return Finish(false,"parallel limit violates independent conserved tangential fluxes");
    }
    if (jump.mass_residual>1e-9 || jump.momentum_residual>1e-8 ||
        jump.energy_residual>1e-8 || jump.electric_residual>1e-8)
      return Finish(false,"parallel limit failed unchanged acceptance thresholds");
  }
  return Finish(true,"parallel switch-on r/p/B agree with analytic RH; polarity and tangential boosts conserve fluxes");
}
// Future-model proof: a differently named provider joins the same factory,
// builder and Runtime adapter without a new switch in AMPS cell publication.
class FutureProvider final:public B::BackgroundProvider {
  // Reuse Parker only as this deterministic extension's physics fixture.
  // Relabel metadata with RuntimeModel provenance so the factory/publication
  // route is exercised without pretending to implement another physical model.
  B::AnalyticParkerProvider inner_;B::SnapshotMetadata metadata_;bool ready_=false;
 public:
  explicit FutureProvider(const B::ParkerConfiguration& p):inner_(p){}
  const char* CanonicalName()const override{return "future-model-test";}
  C::Status Validate()const override{return inner_.Validate();}
  C::Status Prepare(double t)override{auto s=inner_.Prepare(t);if(s.ok()){metadata_=*inner_.PreparedMetadata();metadata_.provider=B::ProviderKind::RuntimeModel;metadata_.providerIdentity=CanonicalName();ready_=true;}return s;}
  const B::SnapshotMetadata* PreparedMetadata()const override{return ready_?&metadata_:nullptr;}
  B::BackgroundSample Evaluate(const C::Vec3& x)const override{return inner_.Evaluate(x);}
  std::string ResolvedManifest()const override{return inner_.ResolvedManifest();}
  B::ProviderCapabilities Capabilities()const override{return inner_.Capabilities();}
};
Result Extension() {
  // Registration is process-global: this deterministic callback runs once in
  // the test registry and explicitly verifies duplicate/built-in protection.
  // The actual Runtime/StandaloneAdapter consumes the constructed snapshot.
  Fixture f;auto s=f.Setup();if(!s.ok())return Finish(false,s.message);
  const auto factory=[](const R::RunConfiguration3D&,std::shared_ptr<B::BackgroundProvider>* out){
    B::ParkerConfiguration p;p.validityCadenceS=60;out->reset(new FutureProvider(p));return C::Status::OK();};
  s=R::RegisterBackgroundModel("test-future-background",factory);
  if(!s.ok()||R::RegisterBackgroundModel("test-future-background",factory).ok()||R::RegisterBackgroundModel("swcme",factory).ok())
    return Finish(false,"factory registration/duplicate protection failed");
  auto options=f.options;options.background=R::BackgroundAuthority::RuntimeModel;options.backgroundModelId="test-future-background";
  std::shared_ptr<const R::RunConfiguration3D> config;
  s=R::RunConfiguration3D::Create(options,&config);if(!s.ok())return Finish(false,s.message);
  std::shared_ptr<B::BackgroundProvider> provider;
  s=R::CreateBackgroundProvider(*config,&provider);if(s.ok())s=provider->Prepare(0);
  if(!s.ok())return Finish(false,s.message);
  std::shared_ptr<const B::BackgroundSnapshot> snapshot;B::BackgroundSnapshotBuilder builder;
  s=builder.Build(*provider,{{C::Const::AU,0,0}},&snapshot);if(!s.ok())return Finish(false,s.message);
  R::Runtime runtime;R::MeshBinding mesh;mesh.layout=config->storage_layout();
  R::StandaloneAdapter adapter;
  s=runtime.Configure(config);if(s.ok())s=runtime.BindMesh(mesh);
  if(s.ok())s=adapter.Initialize(&runtime);
  if(s.ok())s=adapter.PublishSnapshot(&runtime,*snapshot);
  if(!s.ok())return Finish(false,s.message);
  return Finish(runtime.active_snapshot()&&runtime.active_snapshot()->authority==R::BackgroundAuthority::RuntimeModel,
      "a registered future model publishes through the existing provider/snapshot/Runtime contract; built-ins cannot be replaced");
}
}
std::vector<SEP3D::Testing::Descriptor> RegisterSwcmeBackgroundTests() {
  // Stable IDs shared with the Python application catalog. No case requires
  // AMPS initialization; native readback/cadence/ghost evidence is SWBGAMPS.
  using D=SEP3D::Testing::Descriptor;
  auto make=[](const char* id,const char* name,SEP3D::Testing::TestCallback callback){D d;
    d.id=id;d.name=name;d.group="SWBG3D";d.description="SWCME runtime mesh background prerequisite";
    d.initialization=SEP3D::Testing::InitializationLevel::None;d.supportedBuildModes="standalone-no-AMPS";
    d.runtime=SEP3D::Testing::RuntimeClass::Routine;d.seedPolicy="deterministic";
    d.stateIsolation="fresh model/provider per callback";d.callback=std::move(callback);return d;};
  return {make("SWBG3D01","SWCME background input/identity",Configuration),make("SWBG3D02","MHD vector and heated plasma closure",Fields),
    make("SWBG3D03","Moving mesh background snapshots",Evolution),make("SWBG3D04","Full vector compression derivatives",Gradients),
    make("SWBG3D05","Rejected and empty-rank snapshots",Rejections),make("SWBG3D06","Future model registration/publication",Extension),
    make("SWBG3D07","Polar mesh snapshot regression",PolarSnapshots),make("SWBG3D08","Parallel switch-on RH limit",ParallelShockLimit)};
}
