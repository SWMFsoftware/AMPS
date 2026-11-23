#include "sep.h"
#include "reduced_shock_geometry.h"
#include "util/sep_initialization.h"
#include "util/sep_run_configuration.h"
#include "util/sep_background_runtime.h"
#include "adapters/swcme1d_adapter.h"
#include "transport_common.h"
#include "coefficient_providers.h"
#include "reduced_shock_particles.h"
#include "reduced_shock_boundary.h"
#include "util/sep_source_lifecycle.h"
#include "physics_checkpoint.h"
#include "util/sep_population_production_adapter.h"
#include <cstring>
#include <array>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <sstream>

void install_reduced_background_on_field_lines(double);

// This executable links the ordinary configured AMPS archive. No MPI adapter,
// private pool, substitute geometry or test-only physical background is used.
// The two particles are explicit manufactured initial records, not a new
// production shock injection law. Each is allocated only on its segment owner.
int main(int argc,char** argv) {
  MPI_Init(&argc,&argv);MPI_GLOBAL_COMMUNICATOR=MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD,&PIC::ThisThread);MPI_Comm_size(MPI_COMM_WORLD,&PIC::nTotalThreads);
  namespace FL=PIC::FieldLine;namespace PB=PIC::ParticleBuffer;
  namespace R=SEP::ReducedShock;namespace B=SEP::Background;
  const auto fatal=[&](const std::string& message) {
    if(PIC::ThisThread==0)std::cerr<<"[SINT08] ERROR "<<message<<std::endl;
    MPI_Abort(MPI_COMM_WORLD,2);
  };
  if(argc<3||argc>5)fatal("expected runtime geometry input, core event, optional diagnostic mode and checkpoint path");
  const std::string mode=argc>=4?argv[3]:"geometry";
  const bool checkpoint=mode.find("checkpoint-")==0;
  const bool merge=mode=="source-merge";
  const bool sourceMode=mode.find("source-")==0||checkpoint;
  const std::string checkpointPath=argc==5?argv[4]:"material.checkpoint";
  if(mode!="geometry"&&mode!="connectivity"&&mode!="transport"&&!sourceMode)fatal("unknown diagnostic mode");
  const double epoch=mode=="connectivity"||sourceMode?12000.0:0.0;
  PIC::SimulationTime::SetInitialValue(epoch);
  SEP::Initialization::Configuration initialization;
  auto loaded=SEP::Initialization::LoadFile(argv[1],&initialization);if(!loaded.ok())fatal(loaded.message);
  if(mode=="source-disconnected") {
    // A declared manufactured line lies outside the finite SSE angular cap.
    // Its Parker ambient remains valid; lack of shock connection must yield
    // zero injection instead of borrowing the polar line's source strength.
    const double r=initialization.parkerInitialPointM.z;
    initialization.parkerInitialPointM.x=r*std::sin(1.2);
    initialization.parkerInitialPointM.z=r*std::cos(1.2);
  }
  loaded=SEP::Initialization::Install(initialization);if(!loaded.ok())fatal(loaded.message);
  auto run=SEP::Run::Defaults();run.turbulence.source=SEP::Turbulence::Source::Prescribed;
  run.turbulence.coupling=SEP::Turbulence::CouplingPolicy::Disabled;
  run.turbulence.advectionEnabled=run.turbulence.reflectionEnabled=run.turbulence.cascadeEnabled=run.turbulence.shockInjectionEnabled=false;
  // Freeze every background-affecting global before the first publication.
  // Changing ActiveFlag after InstallBackground must remain a rejected
  // provider handoff; a fixture may not evade that production ownership gate.
  SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag=false;
  SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode=false;
  SEP::AccountAdiabaticCoolingFlag=true;
  if(sourceMode) {
    auto& source=run.reducedParticles;source.enabled=true;
    source.releaseDistanceM=mode=="source-invalid"?1e12:1e8; // explicit verification distances [m]
    // This restart fixture deliberately keeps a nonempty accelerated cohort
    // away from the absorbing front, so an empty-population restart cannot
    // masquerade as particle-history qualification. The declared physical
    // reference distance stays fixed in write, resume and continuous roles.
    if(checkpoint||merge)source.releaseDistanceM=1e9;
    source.minimumEnergyJ=initialization.sourceMinimumEnergyJ;
    source.maximumEnergyJ=initialization.sourceMaximumEnergyJ;
    source.releaseFraction=mode=="source-zero"?0:initialization.sourceNetFraction;
    if(mode=="source-invalid-energy") {
      // Manufactured 1–2 eV protons cannot outrun this front at any pitch.
      // A NET source must reject that complete spectral contract atomically;
      // silently dropping all proposals would not realize its specified rate.
      source.minimumEnergyJ=1.602176634e-19;source.maximumEnergyJ=2*source.minimumEnergyJ;
    }
    source.samplesPerStep=merge?128:checkpoint?16:4;
    run.injection.macroparticlesPerEvent=source.samplesPerStep;
    if(mode=="source-mono")source.verificationMonoenergeticEnergyJ=1.602176634e-14;
    run.coefficients.source=SEP::Transport::Coefficient::SourceMode::Prescribed;
    run.coefficients.pitchAngle=SEP::Transport::Coefficient::PitchAngleKind::Constant;
    run.coefficients.constantDmumuPerS=1e-3;
    run.coefficients.spatial=SEP::Transport::Coefficient::SpatialKind::FromPitchAngle;
    const auto coefficients=SEP::Transport::Coefficient::SetActiveConfiguration(run.coefficients);
    if(!coefficients.ok())fatal(coefficients.message);
    run.mover.value=mode=="source-parker"||mode.find("checkpoint-parker-")==0?SEP::Mover::ProductionMover::Parker:
        mode=="source-fte-mfp"||mode.find("checkpoint-fte-mfp-")==0?SEP::Mover::ProductionMover::FocusedTransportMeanFreePath:
        SEP::Mover::ProductionMover::FocusedTransportDiffusion;
    SEP::Mover::SelectProductionMover(run.mover.value);
  }
  auto installed=SEP::Run::InstallActive(run);if(!installed.ok())fatal(installed.message);
  SEP::SW1DAdapter::ConfigurationRequest request;
  const auto configured=SEP::SW1DAdapter::Configure(request);if(!configured.ok())fatal(configured.detail);
  std::string error;if(!R::Configure(argv[2],&error))fatal(error);
  PB::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;
  auto& allocation=FL::VertexAllocationManager;
  allocation.MagneticField=allocation.PlasmaVelocity=allocation.PlasmaDensity=true;
  allocation.PlasmaTemperature=allocation.PlasmaPressure=allocation.Fluence=true;
  allocation.PreviousVertexData.MagneticField=allocation.PreviousVertexData.PlasmaVelocity=allocation.PreviousVertexData.PlasmaDensity=true;
  allocation.PreviousVertexData.PlasmaTemperature=allocation.PreviousVertexData.PlasmaPressure=true;
  SEP::RequestParticleData();
  FL::Init();PB::Init(checkpoint||merge?1024:64);PIC::MolecularData::SetMass(1.67262192369e-27,0);
  PIC::ParticleWeightTimeStep::GlobalParticleWeight[0]=2;
  const double dt=initialization.timeStepS;
  PIC::ParticleWeightTimeStep::GlobalTimeStep[0]=dt;
  PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=true;
  if(sourceMode)PIC::MolecularData::SetElectricCharge(ElectronCharge,0);
  std::snprintf(PIC::OutputDataFileDirectory,sizeof(PIC::OutputDataFileDirectory),".");
  std::list<SEP::cFieldLine> points;
  auto built=R::Native::BuildInitialLine(&points);if(!built.ok())fatal(built.message);
  SEP::Mesh::InitFieldLineAMPS(&points);
  built=R::Native::RegisterLines(epoch);if(!built.ok())fatal(built.message);
  // Output must name the same reduced authority that actually moves the
  // particles, despite the legacy selector/radius remaining initialized.
  // A non-cadence epoch detects accidental use of the prepared-front cache.
  SEP::ShockDiagnostic diagnostic;const double diagnosticEpoch=epoch+129.5;
  const auto diagnosticStatus=SEP::QueryShockDiagnostic(diagnosticEpoch,&diagnostic);
  if(!diagnosticStatus.ok())fatal(diagnosticStatus.message);
  const auto diagnosticView=R::CapturePhysicalStep(epoch,epoch+dt,R::EventIdentity());
  if(!diagnosticView.ok())fatal(diagnosticView.status.message);
  const auto diagnosticSphere=diagnosticView.value.SphereAt(diagnosticEpoch);
  if(!diagnosticSphere.ok())fatal(diagnosticSphere.status.message);
  if(std::string(diagnostic.modelName)!="reduced-SSE/Parker"||
      std::abs(diagnostic.apexRadiusM/(diagnosticSphere.value.centerM[2]+diagnosticSphere.value.radiusM)-1)>1e-12)
    fatal("output shock metadata disagrees with the actual non-cadence core apex");
  PIC::ParallelFieldLines::StaticDecompositionFieldLineLength(0.005);
  FL::SegmentVolume=SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3;
  install_reduced_background_on_field_lines(epoch);
  auto& line=FL::FieldLinesAll[0];const int original=line.GetTotalSegmentNumber();
  const auto finish=[&](bool passed) {
    line.Delete();FL::nFieldLine=0;PIC::ParallelFieldLines::GenerateThreadSegmentTable();
    MPI_Finalize();return passed?0:1;
  };
  if(sourceMode) {
    namespace T=SEP::Transport;namespace A=T::PICAdapter;
    namespace Source=SEP::ReducedParticles::Native;
    A::SetStandaloneRunPurpose(A::RunPurpose::ControlledDiagnostic);
    SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag=false;
    SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode=false;
    SEP::AccountAdiabaticCoolingFlag=true;
    auto ready=SEP::SourceLifecycle::InitializeNativeFresh();if(!ready.ok())fatal(ready.message);
    const bool resume=checkpoint&&mode.find("-resume")!=std::string::npos;
    const bool topologyResume=checkpoint&&mode.find("-resume8")!=std::string::npos;
    const bool write=checkpoint&&mode.find("-write")!=std::string::npos;
    const bool reject=checkpoint&&mode.find("-reject")!=std::string::npos;
    const bool lateReject=checkpoint&&mode.find("-late-reject")!=std::string::npos;
    std::uint64_t completed=0;
    if(checkpoint) {
      // Three unchanged physical samples form one diagnostic window. Saving
      // after four steps exercises both an already completed window and the
      // next partial window, with alternating native buffer phase and newly
      // inserted zero-age material segments. Resume is a separate MPI process.
      PIC::RequiredSampleLength=3;PIC::CollectingSampleCounter=0;
      if(reject) {
        // Corruption must reject on every owner without publishing even one
        // decoded particle, sample, counter or metadata field. Compare the
        // complete normalized image before/after, rather than only a count.
        std::string before,after;
        ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(0,&before);if(!ready.ok())fatal(ready.message);
        if(lateReject) {
          // This is a VALID, populated writer checkpoint. Fail only on the
          // last MPI owner, after detached geometry/table/source/flux plans
          // and replacement particle allocation exist. No collective belongs
          // inside this rank-local probe: LoadCollective owns the subsequent
          // all-rank readiness barrier, including peers whose staging failed.
          const auto firstVertex=line.GetFirstVertex();
          const auto firstSegment=line.GetFirstSegment();
          const auto particleCount=PB::GetAllPartNum();
          unsigned visited=0;std::uint64_t preparedParticles=0;
          completed=99;
          const auto rejected=SEP::PhysicsCheckpoint::Diagnostics::LoadWithPreparationProbe(
              checkpointPath,&completed,[&](const SEP::PhysicsCheckpoint::Diagnostics::PreparedReplacement& staged) {
                ++visited;preparedParticles=staged.particles;
                if(!staged.material||staged.lines!=1)
                  return T::Status::Error(T::StatusCode::InvalidArgument,"late preparation probe has no detached material owner");
                return PIC::ThisThread==PIC::nTotalThreads-1?
                    T::Status::Error(T::StatusCode::InvalidArgument,"injected late rank-local preparation failure"):
                    T::Status::Ok();
              });
          ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(0,&after);if(!ready.ok())fatal(ready.message);
          int local=!rejected.ok()&&rejected.message.find("injected late rank-local")!=std::string::npos&&
              visited==1&&completed==99&&before==after&&PB::GetAllPartNum()==particleCount&&
              line.GetFirstVertex()==firstVertex&&line.GetFirstSegment()==firstSegment;
          std::uint64_t totalPrepared=0;
          MPI_Allreduce(&preparedParticles,&totalPrepared,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
          local=local&&totalPrepared>0;
          // Production has no failure-injection configuration. Even an
          // explicit caller of the diagnostic API must be rejected before a
          // hook can run. Restore the declared purpose before exact images.
          A::SetStandaloneRunPurpose(A::RunPurpose::Production);visited=0;
          const auto denied=SEP::PhysicsCheckpoint::Diagnostics::LoadWithPreparationProbe(
              checkpointPath,&completed,[&](const SEP::PhysicsCheckpoint::Diagnostics::PreparedReplacement&) {
                ++visited;return T::Status::Ok();
              });
          A::SetStandaloneRunPurpose(A::RunPurpose::ControlledDiagnostic);
          ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(0,&after);if(!ready.ok())fatal(ready.message);
          local=local&&!denied.ok()&&visited==0&&completed==99&&before==after;
          // An immediate ordinary retry exercises table/pool/cache cleanup,
          // not just equality of the old owner's normalized record bytes.
          const auto retried=SEP::PhysicsCheckpoint::LoadCollective(checkpointPath,&completed);
          local=local&&retried.ok()&&completed==4;
          int all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
          if(PIC::ThisThread==0)std::cout<<"MATERIAL_LATE_REJECTION prepared_particles="<<totalPrepared
              <<" failed_rank="<<PIC::nTotalThreads-1<<" retry_step="<<completed<<'\n'
              <<"[MCP12] "<<(all?"PASS":"FAIL")<<" late rank-local restore failure preserves complete state and permits retry\n";
          return finish(all);
        }
        const auto rejected=SEP::PhysicsCheckpoint::LoadCollective(checkpointPath,&completed);
        ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(0,&after);if(!ready.ok())fatal(ready.message);
        int local=!rejected.ok()&&before==after&&completed==0,all=0;
        MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
        if(PIC::ThisThread==0)std::cout<<"[MCP07] "<<(all?"PASS":"FAIL")<<" corrupt material checkpoint rejected without partial publication\n";
        return finish(all);
      }
      if(resume) {
        ready=SEP::PhysicsCheckpoint::LoadCollective(checkpointPath,&completed);
        if(!ready.ok())fatal(ready.message);
        if(completed!=(topologyResume?8:4))fatal("material checkpoint did not restore the requested physical step count");
        if(topologyResume) {
          // The eighth material step of this declared mesh has a different
          // segment count from initialization. Require that change explicitly:
          // an equal-count moving-position test cannot qualify topology restore.
          if(line.GetTotalSegmentNumber()==original)fatal("topology fixture did not restore a changed segment count");
          if(PIC::ThisThread==0)std::cout<<"MATERIAL_TOPOLOGY_RESTORE original_segments="<<original
              <<" restored_segments="<<line.GetTotalSegmentNumber()<<'\n';
        }
      }
    }
    bool valid=true;double initialNumber=0;long merged=0;
    for(int step=static_cast<int>(completed);step<(checkpoint?(write?4:9):2);++step) {
      const double start=epoch+step*dt;
      ready=A::FreezeNativeDriverWindow();if(!ready.ok())fatal(ready.message);
      T::Birth::DriverWindow window;ready=A::GetNativeDriverWindow(&window);if(!ready.ok())fatal(ready.message);
      {
        auto read=B::SnapshotStore::Instance().BeginParticleRead(start);
        const auto history=SEP::SourceLifecycle::CaptureNativeState();const auto count=PB::GetAllPartNum();
        ready=Source::Inject(window);
        if(mode=="source-invalid"||mode=="source-invalid-energy") {
          valid=!ready.ok()&&PB::GetAllPartNum()==count&&
              SEP::SourceLifecycle::CaptureNativeState().counters.size()==history.counters.size();
        }else if(!ready.ok())fatal(ready.message);
        if(merge) {
          // An explicitly dense verification source activates the SAME
          // registered native controller, rather than lowering a production
          // trigger or bypassing its conservation checks. The second pass
          // mixes cooled survivors with fresh births; current and birth energy
          // are then independent invariants. Source accounting below must
          // close without treating zero-weight computational deletion as loss.
          SEP::PopulationControl::RegisterConservativeMerger();
          SEP::PopulationControl::BeginConservativePopulationPass();
          for(auto* segment=line.GetFirstSegment();segment;segment=segment->GetNext())
            if(segment->Thread==PIC::ThisThread)
              PIC::ParticleSplitting::FledLine::MergeController(0,segment->FirstParticleIndex,1,1,1,1.0,8,16,8);
          for(const auto& report:SEP::PopulationControl::ConservativeMergeReports())merged+=report.before-report.after;
        }
        std::vector<long> moving;
        for(auto* segment=line.GetFirstSegment();segment;segment=segment->GetNext())if(segment->Thread==PIC::ThisThread) {
          for(long particle=segment->FirstParticleIndex;particle!=-1;particle=PB::GetNext(particle)) {
            moving.push_back(particle);
            if(step==0) {
              const auto* data=PB::GetParticleDataPointer(particle);T::Birth::Record birth;
              std::memcpy(&birth,data+SEP::Offset::BirthRecord,sizeof(birth));
              valid=valid&&birth.pendingFirstMove&&birth.firstMoveDurationS>0&&birth.firstMoveDurationS<dt;
              initialNumber+=PB::GetIndividualStatWeightCorrection(particle)*PIC::ParticleWeightTimeStep::GlobalParticleWeight[0];
              // Independent polar-axis source count. The core supplies the
              // front state; this oracle writes the dimensional formula out
              // rather than calling the application rate helper under test.
              const double tb=window.endEpochS-birth.firstMoveDurationS;
              const auto view=R::CapturePhysicalStep(start,start+dt,R::EventIdentity());
              if(!view.ok())fatal(view.status.message);
              const auto sphere=view.value.SphereAt(tb);if(!sphere.ok())fatal(sphere.status.message);
              const R::Position apex{{0,0,sphere.value.centerM[2]+sphere.value.radiusM}};
              const auto front=view.value.FrontAt(apex,tb);if(!front.ok())fatal(front.status.message);
              const double expected=run.reducedParticles.releaseFraction*front.value.upstream.protonNumberDensityM3*
                  (front.value.normalSpeedMPerS-front.value.upstream.velocityMPerS[2])*
                  SEP::FieldLine::FluxTubeGeometry::GetMagneticFluxWb(0)/std::abs(front.value.upstream.magneticFieldT[2])*
                  dt/run.reducedParticles.samplesPerStep;
              const double actual=PB::GetIndividualStatWeightCorrection(particle)*PIC::ParticleWeightTimeStep::GlobalParticleWeight[0];
              if(!merge)valid=valid&&std::abs(actual/expected-1)<2e-10;
              SEP::SourceLifecycle::Event sourceEvent;
              sourceEvent.channel.kind=SEP::SourceLifecycle::Kind::Shock;
              sourceEvent.channel.line=sourceEvent.channel.species=0;
              sourceEvent.campaign=A::CampaignRandomSeed;sourceEvent.index=0;
              std::uint64_t identity=0;std::memcpy(&identity,data+SEP::Offset::StableParticleId,sizeof(identity));
              std::uint64_t macro=run.reducedParticles.samplesPerStep;
              for(std::uint64_t j=0;j<run.reducedParticles.samplesPerStep;++j)
                if(SEP::SourceLifecycle::ParticleIdentity(sourceEvent,j)==identity)macro=j;
              valid=valid&&macro<run.reducedParticles.samplesPerStep;
              auto pitch=SEP::SourceLifecycle::Random(sourceEvent,macro,UINT64_C(0x52534d55));
              const double u=pitch.UniformOpen01();
              const double speed=std::hypot(PB::GetVParallel(particle),PB::GetVNormal(particle));
              const auto release=view.value.PrimitiveAt({{0,0,apex[2]+run.reducedParticles.releaseDistanceM}},tb);
              if(!release.ok())fatal(release.status.message);
              const double critical=(front.value.normalSpeedMPerS-release.value.velocityMPerS[2])/speed;
              // The fixture's super-Alfvenic proton energies give a partial
              // outward cone. Its exact affine-flux inverse is independent of
              // the production core CDF table/inversion implementation.
              const bool isotropicProjection=!SEP::Mover::CurrentCapabilities().usesPitchAngleState;
              const double expectedMu=isotropicProjection?2*u-1:critical+(1-critical)*std::sqrt(u);
              valid=valid&&(isotropicProjection||(critical>-1&&critical<1))&&
                  std::abs(PB::GetVParallel(particle)/speed-expectedMu)<1e-10;
              if(macro==0) {
                // Grade the polar performance specialization against the
                // complete advected-polyline quadratic intersection routine.
                // Both endpoints use the accepted INTERNAL interval; the
                // original 300-s physical host duration is never changed.
                A::ParticleContext context;auto status=A::LoadParticle(particle,&context);
                if(!status.ok())fatal(status.message);
                double effective=0;status=A::PrepareParticleStep(&context,dt,&effective);
                if(!status.ok())fatal(status.message);
                SEP::ReducedParticles::Native::BoundaryInterval boundary;
                double internal=effective;
                status=SEP::ReducedParticles::Native::PrepareBoundary(context,0,&internal,&boundary);
                if(!status.ok())fatal(status.message);
                // Repeat the exact query through the worker's immutable-view
                // cache. Cached coordinates must be bitwise identical to the
                // original physical query; the full polyline oracle below
                // still grades them to the unchanged 1-mm absolute bound.
                auto repeated=boundary;double repeatedDuration=internal;
                status=SEP::ReducedParticles::Native::PrepareBoundary(context,0,&repeatedDuration,&repeated);
                if(!status.ok())fatal(status.message);
                valid=valid&&repeatedDuration==internal&&repeated.frontBeginArcM==boundary.frontBeginArcM&&
                    repeated.frontEndArcM==boundary.frontEndArcM;
                const auto mesh=R::Native::CaptureMaterialView();if(!mesh.ok())fatal(mesh.status.message);
                for(int endpoint=0;endpoint<2;++endpoint) {
                  const double queryEpoch=window.endEpochS-(birth.firstMoveDurationS-(endpoint?internal:0));
                  const auto full=R::Native::ConnectionsAt(queryEpoch);
                  if(!full.ok())fatal(full.status.message);
                  valid=valid&&full.value.size()==1;
                  if(full.value.size()==1) {
                    const auto& root=full.value.front();const auto& arc=mesh.value->arc;
                    const double reference=arc[root.segment]+root.fraction*(arc[root.segment+1]-arc[root.segment]);
                    valid=valid&&std::abs(reference-(endpoint?boundary.frontEndArcM:boundary.frontBeginArcM))<1e-3;
                  }
                  // Additional software-equivalence control: retain the old
                  // full-range binary search independently of the optimized
                  // departure-coordinate hint. Every vertex still uses the
                  // actual core flow map and unchanged accuracy controls.
                  // Require bitwise-identical scalar coordinates, in addition
                  // to the independent full-polyline 1-mm physics gate above.
                  const auto& material=*mesh.value;
                  const auto sphere=material.interval.SphereAt(queryEpoch);
                  if(!sphere.ok())fatal(sphere.status.message);
                  const double shock=sphere.value.centerM[2]+sphere.value.radiusM;
                  const auto controls=R::Native::Controls(dt,material.arc[1]);
                  const auto height=[&](std::size_t i) {
                    const auto moved=R::MoveVertex(material.interval,material.vertices[i],
                        material.interval.DepartureS(),queryEpoch,controls);
                    if(!moved.ok())fatal(moved.status.message);
                    return moved.value;
                  };
                  std::size_t left=0,right=material.vertices.size()-1;
                  auto a=height(left),b=height(right);
                  while(right-left>1) {
                    const auto middle=left+(right-left)/2;const auto point=height(middle);
                    if(point[2]<=shock){left=middle;a=point;}else{right=middle;b=point;}
                  }
                  const double exact=material.arc[left]+(shock-a[2])/(b[2]-a[2])*
                      (material.arc[right]-material.arc[left]);
                  valid=valid&&exact==(endpoint?boundary.frontEndArcM:boundary.frontBeginArcM);
                }
              }
            }
          }
          segment->FirstParticleIndex=-1;
        }
        if(mode!="source-invalid"&&mode!="source-invalid-energy")for(long particle:moving) {
          PB::SetNext(-1,particle);PB::SetPrev(-1,particle);
          SEP::Mover::DispatchProductionMover(particle,dt,nullptr);
        }
        // Reproduce the native joined-list contract after direct dispatch.
        // This fixture does not claim to execute the PIC::TimeStep driver;
        // the separately recorded application campaign provides that gate.
        for(auto* segment=line.GetFirstSegment();segment;segment=segment->GetNext())
          segment->FirstParticleIndex=segment->tempFirstParticleIndex.exchange(-1);
      }
      // A mover can finish on a segment owned by another MPI rank. Exactly as
      // in the native driver, join ALL local destination lists first, then
      // migrate their complete particle records to the destination owners.
      // Joining only this rank's owned segments leaves outbound records in a
      // temporary list; skipping migration would instead strand those records
      // outside the next step's owned-segment sweep. Both errors are exposed
      // by the unchanged production CloseNativeDriverWindow invariant.
      PIC::ParallelFieldLines::ExchangeFieldLineParticles();
      FL::Sampling();PIC::SimulationTime::SetInitialValue(start+dt);
      ready=A::CloseNativeDriverWindow();if(!ready.ok())fatal(ready.message);
      if(checkpoint) {
        ++PIC::CollectingSampleCounter;
        if(PIC::CollectingSampleCounter==PIC::RequiredSampleLength) {
          const auto filename="material-output-"+std::to_string(PIC::DataOutputFileNumber)+".dat";
          FL::Output(filename.c_str(),false);
          PIC::LastSampleLength=PIC::CollectingSampleCounter;
          PIC::CollectingSampleCounter=0;++PIC::DataOutputFileNumber;
        }
      }
      if(mode=="source-invalid"||mode=="source-invalid-energy")break;
      long controlledRetirement=-1;long countBeforeRetirement=0;
      if(mode=="source-fte-dmumu"&&step==1&&PIC::nTotalThreads==1) {
        // Independent geometry-owner fixture, AFTER actual source/mover/
        // sampling work: relocate one surviving record into the outgoing tail
        // chord without changing its weight, momentum, birth or RNG history.
        // This manufactured positional control grades deletion of an active
        // source cohort by geometry retirement; it is not a physical particle
        // trajectory or evidence about the source's spatial distribution.
        for(auto* segment=line.GetFirstSegment();segment&&controlledRetirement==-1;segment=segment->GetNext())
          if(segment->FirstParticleIndex!=-1) {
            controlledRetirement=segment->FirstParticleIndex;
            segment->FirstParticleIndex=PB::GetNext(controlledRetirement);
            if(segment->FirstParticleIndex!=-1)PB::SetPrev(-1,segment->FirstParticleIndex);
          }
        valid=valid&&controlledRetirement!=-1;
        if(controlledRetirement!=-1) {
          auto* tail=line.GetLastSegment();
          const double coordinate=line.GetTotalSegmentNumber()-0.5;
          PB::SetFieldLineCoord(coordinate,controlledRetirement);
          double x[3];line.GetCartesian(x,coordinate);PB::SetX(x,controlledRetirement);
          PB::SetPrev(-1,controlledRetirement);PB::SetNext(tail->FirstParticleIndex,controlledRetirement);
          if(tail->FirstParticleIndex!=-1)PB::SetPrev(controlledRetirement,tail->FirstParticleIndex);
          tail->FirstParticleIndex=controlledRetirement;
          countBeforeRetirement=PB::GetAllPartNum();
          std::cout<<"CONTROLLED_GEOMETRY_RELOCATION coordinate="<<coordinate<<" expected_retired_macros=1\n";
        }
      }
      ready=Source::BeforeGeometry();if(!ready.ok())fatal(ready.message);
      const auto moved=R::Native::Advance(start,start+dt);if(!moved.ok())fatal(moved.message);
      ready=Source::AfterGeometry(start+dt);if(!ready.ok())fatal(ready.message);
      if(controlledRetirement!=-1)valid=valid&&PB::GetAllPartNum()==countBeforeRetirement-1;
      valid=valid&&line.GeometryEpochS==start+dt&&PIC::ParticleWeightTimeStep::GlobalTimeStep[0]==dt;
      if(checkpoint) {
        completed=step+1;
        if(!resume&&!write&&completed==8) {
          // The continuous reference never restores. This observational save
          // supplies a SECOND, changed-topology checkpoint for a fresh process
          // that must execute step nine and the next output. Both alternating
          // buffer phases and a different segment count are thereby exercised.
          std::string before,after;
          ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(completed,&before);if(!ready.ok())fatal(ready.message);
          ready=SEP::PhysicsCheckpoint::SaveCollective("material-topology.checkpoint",completed);if(!ready.ok())fatal(ready.message);
          ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(completed,&after);if(!ready.ok())fatal(ready.message);
          valid=valid&&before==after;
        }
        if(write&&completed==4) {
          // Retain the exact replicated-state diagnostic on every rank before
          // attempting publication. An agreement error then has a byte offset
          // and both original inputs for investigation; no failed writer is
          // reclassified as a physics pass or stripped of its replica check.
          std::string replica;ready=SEP::PhysicsCheckpoint::LocalReplicatedPhysicsImage(completed,&replica);
          if(!ready.ok())fatal(ready.message);
          std::ofstream diagnostic("material-replica-rank"+std::to_string(PIC::ThisThread)+".bin",std::ios::binary);
          diagnostic.write(replica.data(),replica.size());diagnostic.close();
          if(diagnostic.fail())fatal("cannot retain replicated checkpoint diagnostic");
          ready=SEP::PhysicsCheckpoint::SaveCollective(checkpointPath,completed);
          if(!ready.ok())fatal(ready.message);
        }
      }
    }
    if(mode=="source-parker") {
      // Grade the real collective writer AFTER the final transport/geometry
      // checks, so its sampling-buffer rotation cannot affect a later fixture
      // step. The production application installs this exact title callback.
      FL::UserDefinedTecplotFileTitle=SEP::TecplotFileTitle;
      FL::Output("native-source-output.dat",false);
      if(PIC::ThisThread==0) {
        std::ifstream output("native-source-output.dat");std::string row;bool title=false,field=false;
        while(std::getline(output,row)) {
          if(row.find("shock model=reduced-SSE/Parker;")!=std::string::npos)title=true;
          std::istringstream values(row);double value=0;
          if(!(values>>value))continue;
          // The documented geometry/count/speed/segment prefix has 9 values;
          // electric-field output is disabled in this controlled profile.
          for(int column=1;column<9;++column)valid=valid&&bool(values>>value);
          const double* b=line.GetFirstVertex()->GetMagneticField();
          for(int component=0;component<3;++component) {
            valid=valid&&bool(values>>value);
            // Independent unit identity: 1 nT = 1e-9 T. Allow only the existing
            // six-digit text-serialization precision, not a physics tolerance.
            valid=valid&&std::abs(value*1e-9-b[component])<=1e-6*std::abs(b[component]);
          }
          field=true;break;
        }
        valid=valid&&title&&field;
      }
    }
    double globalInitialNumber=0;MPI_Allreduce(&initialNumber,&globalInitialNumber,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    if(mode=="source-zero"||mode=="source-disconnected")valid=valid&&PB::GetAllPartNum()==0&&globalInitialNumber==0;
    else if(!resume&&mode!="source-invalid"&&mode!="source-invalid-energy")valid=valid&&globalInitialNumber>0;
    if(checkpoint) {
      std::string image;ready=SEP::PhysicsCheckpoint::LocalPhysicsImage(completed,&image);
      if(!ready.ok())fatal(ready.message);
      std::ofstream output("material-final-rank"+std::to_string(PIC::ThisThread)+".bin",std::ios::binary);
      output.write(image.data(),image.size());output.close();valid=valid&&!output.fail();
      // Export complete owner-neutral records as an independent comparison
      // view. Only allocation links are normalized, exactly as in persistence;
      // IDs, birth energy, all operator streams and MFP hazard remain bytes.
      // The runner combines rank files and sorts IDs, allowing a changed MPI
      // decomposition without weakening the exact particle-history assertion.
      std::ofstream records("material-particles-rank"+std::to_string(PIC::ThisThread)+".bin",std::ios::binary);
      const auto word=[&](std::uint64_t value){for(unsigned byte=0;byte<8;++byte)records.put(static_cast<char>(value>>(8*byte)));};
      word(PB::ParticleDataLength);
      long localCount=0;
      for(auto* segment=line.GetFirstSegment();segment;segment=segment->GetNext())
        for(long pointer=segment->FirstParticleIndex;pointer!=-1;pointer=PB::GetNext(pointer)) {
          std::string record(reinterpret_cast<const char*>(PB::GetParticleDataPointer(pointer)),PB::ParticleDataLength);
          auto* bytes=reinterpret_cast<PB::byte*>(&record[0]);
          PB::SetNext(-1,bytes);PB::SetPrev(-1,bytes);
          std::uint64_t id=0;std::memcpy(&id,bytes+SEP::Offset::StableParticleId,sizeof(id));
          word(id);records.write(record.data(),record.size());++localCount;
        }
      records.close();valid=valid&&!records.fail();
      long totalCount=0;MPI_Allreduce(&localCount,&totalCount,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
      valid=valid&&totalCount>0;
    }
    const char* checkpointId=run.mover.value==SEP::Mover::ProductionMover::Parker?"MCP01":
        run.mover.value==SEP::Mover::ProductionMover::FocusedTransportDiffusion?"MCP02":"MCP03";
    if(merge) {
      long totalMerged=0;MPI_Allreduce(&merged,&totalMerged,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
      valid=valid&&totalMerged>0;
      if(PIC::ThisThread==0)std::cout<<"MATERIAL_MERGE removed_computational_representatives="<<totalMerged<<'\n';
    }
    const char* identity=merge?"MCP08":checkpoint?checkpointId:mode=="source-parker"?"SINT25N":mode=="source-fte-dmumu"?"SINT26N":
        mode=="source-fte-mfp"?"SINT27N":mode=="source-mono"?"SINT28N":mode=="source-zero"?"SINT29N":
        mode=="source-invalid"?"SINT30N":mode=="source-invalid-energy"?"SINT32N":"SINT31N";
    int local=valid,all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(PIC::ThisThread==0)std::cout<<'['<<identity<<"] "<<(all?"PASS":"FAIL")
        <<" native reduced source, birth clock, selected mover and joined energy/number ledger; dt="<<dt<<std::endl;
    return finish(all);
  }
  if(mode=="connectivity") {
    const auto snapshot=B::SnapshotStore::Instance().Current();const auto generation=line.GeometryGeneration;
    const auto now=R::Native::ConnectionsAt(epoch),later=R::Native::ConnectionsAt(epoch+dt);
    bool valid=now.ok()&&later.ok()&&now.value.size()==1&&later.value.size()==1;
    if(valid) {
      const auto& a=now.value.front();const auto& b=later.value.front();
      valid=a.front.acceptedShock&&b.front.acceptedShock&&a.front.epochS==epoch&&b.front.epochS==epoch+dt&&
          a.surfaceResidualM<0.02&&b.surfaceResidualM<0.02&&b.arcLengthM>a.arcLengthM;
      auto* edge=line.GetSegment(static_cast<int>(a.segment));
      const double coordinate=a.segment+a.fraction+0.1;
      const auto estimate=R::Native::CrossingTimeEstimate(coordinate,0);
      // Independent polar-axis cut speed: the front radius evolves in HCI,
      // while this coordinate follows plasma. Only their DIFFERENCE sets the
      // shock approach; U does not enter relative particle displacement.
      const double relative=a.front.normalSpeedMPerS-a.front.upstream.velocityMPerS[2];
      const double expected=0.1*edge->GetLength()/relative;
      valid=valid&&estimate.ok()&&std::abs(estimate.value-expected)<1e-6*expected;
      if(PIC::ThisThread==0)std::cout<<std::setprecision(17)<<"CONNECTIVITY_METRIC epoch_s="<<epoch
          <<" root_s="<<a.segment+a.fraction<<" next_root_s="<<b.segment+b.fraction<<" approach_s="<<estimate.value<<std::endl;
    }
    valid=valid&&B::SnapshotStore::Instance().Current()==snapshot&&line.GeometryGeneration==generation&&line.GeometryEpochS==epoch;
    int local=valid,all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
    if(PIC::ThisThread==0)std::cout<<"[SINT13N] "<<(all?"PASS":"FAIL")<<" native SSE material-time connectivity and relative front speed"<<std::endl;
    return finish(all);
  }
  if(mode=="transport") {
    namespace T=SEP::Transport;namespace A=T::PICAdapter;namespace C=T::Coefficient;
    // These are explicitly manufactured initial conditions on the real core
    // ambient. A slow test ion makes any accidental U*dt unmistakable: correct
    // motion is hundreds of metres, whereas U*300 s is about 10^8 metres.
    // The physical dt remains the input's unchanged 300 s for ALL movers.
    A::SetStandaloneRunPurpose(A::RunPurpose::ControlledDiagnostic);
    SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag=false;
    SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode=false;
    SEP::AccountAdiabaticCoolingFlag=true;
    bool allPassed=true;
    for(const auto& descriptor:SEP::Mover::Registry()) {
      auto typed=C::ActiveConfiguration();typed.source=C::SourceMode::Prescribed;
      typed.pitchAngle=C::PitchAngleKind::Constant;typed.constantDmumuPerS=
          descriptor.mover==SEP::Mover::ProductionMover::Parker?1e-3:0.0;
      typed.spatial=C::SpatialKind::FromPitchAngle;
      typed.meanFreePath=C::MeanFreePathKind::Tenishev2005;
      auto status=C::SetActiveConfiguration(typed);if(!status.ok())fatal(status.message);
      SEP::Mover::SelectProductionMover(descriptor.mover);
      const int edgeIndex=40;const double initialCoordinate=edgeIndex+0.5;
      auto* edge=line.GetSegment(edgeIndex);bool valid=true;
      const auto moverBefore=SEP::Mover::CompletedDispatchCount();
      if(edge->Thread==PIC::ThisThread) {
        const auto particle=PB::GetNewParticle();PB::SetI(0,particle);PB::SetFieldLineId(0,particle);
        PB::SetFieldLineCoord(initialCoordinate,particle);PB::SetVParallel(0.6,particle);PB::SetVNormal(0.8,particle);
        PB::SetIndividualStatWeightCorrection(3,particle);PB::SetNext(-1,particle);PB::SetPrev(-1,particle);
        double initialX[3];line.GetCartesian(initialX,initialCoordinate);PB::SetX(initialX,particle);
        A::InitializeParticleTransportStateWithIdentity(particle,UINT64_C(0x53494e540000)+static_cast<int>(descriptor.mover)+1);
        if(descriptor.mover==SEP::Mover::ProductionMover::FocusedTransportMeanFreePath)
          // A conditioned no-event path, not a replacement waiting-time law.
          // The existing FTEM tests separately grade the exponential clock,
          // residual hazard and wave-frame conservation at actual events.
          *reinterpret_cast<double*>(PB::GetParticleDataPointer(particle)+SEP::Offset::MfpOpticalDepth)=1000;
        auto read=B::SnapshotStore::Instance().BeginParticleRead(epoch);
        A::ParticleContext context;status=A::LoadParticle(particle,&context);if(!status.ok())fatal(status.message);
        double effective=0;status=A::PrepareParticleStep(&context,dt,&effective);if(!status.ok())fatal(status.message);
        A::LocalBackground background;status=A::EvaluateLocalBackground(context,&background);if(!status.ok())fatal(status.message);
        valid=background.coordinateAdvectionMPerS==0&&background.plasmaAdvectionMPerS>1e5&&
            background.generation==read.snapshot().field_line_generation();
        const auto p0=T::MomentumFromSpeed(1.0,context.state.massKg,SpeedOfLight);
        const auto dispatched=SEP::Mover::DispatchProductionMover(particle,dt,nullptr);
        valid=valid&&dispatched==_PARTICLE_MOTION_FINISHED_&&SEP::Mover::CompletedDispatchCount()==moverBefore+1;
        const double movedCoordinate=PB::GetFieldLineCoord(particle);
        const double displacement=(movedCoordinate-initialCoordinate)*edge->GetLength();
        const double speed=std::hypot(PB::GetVParallel(particle),PB::GetVNormal(particle));
        const auto p1=T::MomentumFromSpeed(speed,context.state.massKg,SpeedOfLight);
        valid=valid&&std::isfinite(displacement)&&std::abs(displacement)<10000&&p1.status.ok()&&p1.value<p0.value;
        double derived[3],stored[3];line.GetCartesian(derived,movedCoordinate);PB::GetX(stored,particle);
        for(int d=0;d<3;++d)valid=valid&&derived[d]==stored[d];
        if(descriptor.mover==SEP::Mover::ProductionMover::Parker) {
          const double expected=p0.value*std::exp(-background.velocityDivergencePerS*dt/3);
          valid=valid&&std::abs(p1.value/expected-1)<2e-6;
        } else {
          // Independent RK4 reference for the published gyrotropic ODE at
          // the departure position; the <1 km trajectory changes these core
          // coefficients by <10^-7 here. This independently grades plasma-
          // frame focusing/strain/cooling, without calling a mover kernel to
          // manufacture its own reference. Numerical tolerances are frozen.
          using State=std::array<double,3>;State reference{{0.0,std::log(p0.value),0.6}};
          const auto rate=[&](const State& y) {
            const double mu=y[2],p=std::exp(y[1]);const double speedRef=p/context.state.massKg;
            return State{{speedRef*mu,-0.5*((1-mu*mu)*background.velocityDivergencePerS+
                (3*mu*mu-1)*background.fieldAlignedStrainPerS),
                0.5*(1-mu*mu)*(-speedRef*background.dLnAbsBdsPerM+
                    mu*(background.velocityDivergencePerS-3*background.fieldAlignedStrainPerS))}};
          };
          const auto shifted=[](State a,const State& b,double h){for(int d=0;d<3;++d)a[d]+=h*b[d];return a;};
          const double h=dt/4096;
          for(int i=0;i<4096;++i){auto k1=rate(reference),k2=rate(shifted(reference,k1,h/2)),
              k3=rate(shifted(reference,k2,h/2)),k4=rate(shifted(reference,k3,h));
            for(int d=0;d<3;++d)reference[d]+=h*(k1[d]+2*k2[d]+2*k3[d]+k4[d])/6;}
          valid=valid&&std::abs(displacement-reference[0])<0.1&&
              std::abs(p1.value/std::exp(reference[1])-1)<2e-5&&
              std::abs(PB::GetVParallel(particle)/speed-reference[2])<2e-5;
        }
        std::cout<<std::setprecision(17)<<"TRANSPORT_METRIC mover="<<descriptor.canonicalName<<" rank="<<PIC::ThisThread
            <<" dt_s="<<dt<<" displacement_m="<<displacement<<" momentum_ratio="<<p1.value/p0.value<<std::endl;
        // Direct dispatch attaches only to the standard temporary list. The
        // test owns this one record and removes it before selecting a mover.
        for(auto* segment=line.GetFirstSegment();segment;segment=segment->GetNext())
          if(segment->tempFirstParticleIndex.load()==particle)segment->tempFirstParticleIndex=-1;
        PB::DeleteParticle(particle);
      }
      int local=valid,all=0;MPI_Allreduce(&local,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
      const char* identity=descriptor.mover==SEP::Mover::ProductionMover::Parker?"SINT14N":
          descriptor.mover==SEP::Mover::ProductionMover::FocusedTransportDiffusion?"SINT15N":"SINT16N";
      if(PIC::ThisThread==0)std::cout<<'['<<identity<<"] "<<(all?"PASS":"FAIL")<<" actual native mover on unchanged material core background"<<std::endl;
      allPassed=allPassed&&all;
    }
    return finish(allPassed);
  }
  namespace Tube=SEP::FieldLine::FluxTubeGeometry;
  const double flux=Tube::GetMagneticFluxWb(0);
  const auto volumeConsistent=[&]() {
    auto* segment=line.GetSegment(40);const double b0=segment->GetBegin()->GetMagneticField()[2],
        b1=segment->GetEnd()->GetMagneticField()[2];
    // The maintained native interpolant is componentwise-linear B. On this
    // polar radial tube its exact volume integral is Phi*L*log(B1/B0)/(B1-B0).
    // This independent quadrature oracle grades the ACTUAL discrete measure,
    // without pretending that linear B equals the continuum r^-2 field.
    const double expected=flux*segment->GetLength()*std::log(b1/b0)/(b1-b0);
    const double actual=Tube::SegmentVolumeM3(segment,0);
    const double area=Tube::AreaAtVertexM2(segment->GetBegin(),0);
    return std::abs(actual/expected-1)<2e-8&&std::abs(area*b0/flux-1)<1e-13&&Tube::GetMagneticFluxWb(0)==flux;
  };
  const bool initialVolume=volumeConsistent();
  const auto particle=[&](int segment,double fraction) {
    auto* edge=line.GetSegment(segment);if(edge->Thread!=PIC::ThisThread)return -1L;
    const long int p=PB::GetNewParticle(edge->FirstParticleIndex);PB::SetI(0,p);PB::SetFieldLineCoord(segment+fraction,p);PB::SetFieldLineId(0,p);
    PB::SetIndividualStatWeightCorrection(3,p);PB::SetVParallel(2,p);PB::SetVNormal(3,p);
    double x[3];line.GetCartesian(x,PB::GetFieldLineCoord(p));PB::SetX(x,p);return p;
  };
  const auto kept=particle(0,0.25);particle(original-1,0.75);
  // Non-dyadic coordinates straddle the binary64 exponent boundary 15->16
  // when one head segment is inserted. Unlike the old exact 0.25 fraction,
  // these expose loss of low fractional bits in S+prefix. Cartesian X must be
  // derived from that STORED coordinate; using the pre-addition fraction can
  // leave an inconsistent one-ULP native view and break a correct checkpoint.
  std::vector<long int> precisionParticles;
  for(unsigned i=0;i<32;++i) {
    const auto p=particle(15,(i+0.123456789)/64.0);
    if(p!=-1)precisionParticles.push_back(p);
  }
  bool valid=original==512&&initialVolume;
  // Thread arrivals and MPI allocation must not choose the order of floating
  // additions in a material sample. Reverse the SAME 32 manufactured records,
  // without changing any physical value, and compare all native vertex bytes
  // and material collecting numerators after one actual Sampling() call.
  // Restore every touched exposure, numerator and counter between calls: these
  // are diagnostic trials, not extra physical sampling opportunities. Original
  // list order is restored too, before the geometry/retirement tests continue.
  std::vector<std::string> unsampledVertices;
  for(int v=0;v<=original;++v) {
    auto* vertex=line.GetVertex(v);
    unsampledVertices.emplace_back(vertex->GetAssociatedDataBufferPointer(),vertex->AssociatedDataLength());
  }
  std::vector<FL::Geometry::Exposure> unsampledExposure;
  std::vector<std::vector<double>> unsampledCollecting;
  for(auto* edge=line.GetFirstSegment();edge;edge=edge->GetNext()) {
    unsampledExposure.push_back(edge->SamplingExposure);
    unsampledCollecting.push_back(edge->MaterialCollecting);
  }
  const int unsampledCycle=FL::SampleCycleCounter;
  const auto resetSampling=[&]() {
    for(int v=0;v<=original;++v)
      std::memcpy(line.GetVertex(v)->GetAssociatedDataBufferPointer(),unsampledVertices[v].data(),unsampledVertices[v].size());
    int s=0;
    for(auto* edge=line.GetFirstSegment();edge;edge=edge->GetNext(),++s) {
      edge->SamplingExposure=unsampledExposure[s];edge->MaterialCollecting=unsampledCollecting[s];
    }
    FL::SampleCycleCounter=unsampledCycle;
  };
  const auto sampledImage=[&]() {
    std::string bytes;
    for(int v=0;v<=original;++v) {
      auto* vertex=line.GetVertex(v);
      bytes.append(vertex->GetAssociatedDataBufferPointer(),vertex->AssociatedDataLength());
    }
    for(auto* edge=line.GetFirstSegment();edge;edge=edge->GetNext())
      if(!edge->MaterialCollecting.empty())bytes.append(
          reinterpret_cast<const char*>(edge->MaterialCollecting.data()),edge->MaterialCollecting.size()*sizeof(double));
    return bytes;
  };
  auto* precisionEdge=line.GetSegment(15);std::vector<long int> originalOrder;
  for(long p=precisionEdge->FirstParticleIndex;p!=-1;p=PB::GetNext(p))originalOrder.push_back(p);
  const auto attachOrder=[&](const std::vector<long int>& order) {
    precisionEdge->FirstParticleIndex=order.empty()?-1:order.front();
    for(std::size_t i=0;i<order.size();++i) {
      PB::SetPrev(i?order[i-1]:-1,order[i]);
      PB::SetNext(i+1<order.size()?order[i+1]:-1,order[i]);
    }
  };
  FL::Sampling();const auto forwardSample=sampledImage();resetSampling();
  auto reverseOrder=originalOrder;std::reverse(reverseOrder.begin(),reverseOrder.end());attachOrder(reverseOrder);
  FL::Sampling();const auto reverseSample=sampledImage();resetSampling();attachOrder(originalOrder);
  int localSampleExact=forwardSample==reverseSample,allSampleExact=0;
  MPI_Allreduce(&localSampleExact,&allSampleExact,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  valid=valid&&allSampleExact;
  if(PIC::ThisThread==0)std::cout<<"MATERIAL_SAMPLE_ORDER checked=32 exact="<<allSampleExact<<'\n';
  FL::Sampling();PIC::SimulationTime::SetInitialValue(dt);
  auto advanced=R::Native::Advance(0,dt);if(!advanced.ok())fatal(advanced.message);
  valid=valid&&line.GetTotalSegmentNumber()==original+1&&line.FirstMaterialLabel==-1;
  if(kept!=-1)valid=valid&&PB::GetFieldLineCoord(kept)==1.25&&PB::GetVParallel(kept)==2&&PB::GetVNormal(kept)==3;
  valid=valid&&line.GetFirstSegment()->FirstParticleIndex==-1&&line.GetFirstSegment()->SamplingExposure.collecting==0;
  unsigned localMismatches=0;
  for(auto p:precisionParticles) {
    double actual[3],expected[3];PB::GetX(actual,p);line.GetCartesian(expected,PB::GetFieldLineCoord(p));
    for(unsigned d=0;d<3;++d)if(actual[d]!=expected[d])++localMismatches;
  }
  unsigned mismatches=0;MPI_Allreduce(&localMismatches,&mismatches,1,MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
  valid=valid&&mismatches==0;
  if(PIC::ThisThread==0)std::cout<<"MATERIAL_COORDINATE_CANONICAL checked=32 mismatches="<<mismatches<<'\n';
  FL::Sampling();PIC::SimulationTime::SetInitialValue(2*dt);
  advanced=R::Native::Advance(dt,2*dt);if(!advanced.ok())fatal(advanced.message);
  valid=valid&&line.GetTotalSegmentNumber()==original&&line.GeometryEpochS==2*dt&&line.GeometryGeneration==2&&volumeConsistent();
  valid=valid&&line.GetFirstSegment()->SamplingExposure.collecting==1&&line.GetSegment(2)->SamplingExposure.collecting==2;
  long int local=PB::GetAllPartNum(),global=0;MPI_Allreduce(&local,&global,1,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);
  valid=valid&&global==33&&B::SnapshotStore::Instance().Current()->epoch_seconds()==2*dt;
  int localGood=valid?1:0,good=0;MPI_Allreduce(&localGood,&good,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(PIC::ThisThread==0)std::cout<<"[SINT08] "<<(good?"PASS":"FAIL")<<" actual core-driven PIC material motion, particle loss and sampling; dt="<<dt<<" ranks="<<PIC::nTotalThreads<<std::endl;

  // A rank-local fault in the LAST old vertex reaches the real late trajectory
  // callback after earlier vertices have been staged. Collective readiness
  // must reject on EVERY rank without publishing any geometry or metadata.
  const auto metadata=B::SnapshotStore::Instance().Current();
  const auto first=line.GetFirstVertex()->GetX()[2];
  auto* last=line.GetLastVertex();double saved[3];std::copy(last->GetX(),last->GetX()+3,saved);
  if(PIC::ThisThread==PIC::nTotalThreads-1){double unsupported[3]={0,0,3*149597870700.};last->SetX(unsupported);}
  const auto coreEpoch=R::FrontEpoch();const auto coreGeneration=R::Metadata().generation;
  const auto rejectedInitial=R::Native::InstallBackground(2*dt);
  const auto rejected=R::Native::Advance(2*dt,3*dt);
  bool atomic=!rejectedInitial.ok()&&R::FrontEpoch()==coreEpoch&&R::Metadata().generation==coreGeneration&&!rejected.ok()&&line.GeometryEpochS==2*dt&&line.GeometryGeneration==2&&
      line.GetFirstVertex()->GetX()[2]==first&&B::SnapshotStore::Instance().Current()==metadata;
  last->SetX(saved);
  localGood=atomic?1:0;int allAtomic=0;MPI_Allreduce(&localGood,&allAtomic,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  if(PIC::ThisThread==0)std::cout<<"[SINT09] "<<(allAtomic?"PASS":"FAIL")<<" late rank-local core failure preserves native material state and metadata"<<std::endl;
  line.Delete();FL::nFieldLine=0;PIC::ParallelFieldLines::GenerateThreadSegmentTable();
  MPI_Finalize();return good&&allAtomic?0:1;
}
