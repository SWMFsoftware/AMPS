// Dedicated wave/codec operator verification. The fixed geometry defines an
// operator's coefficient grid; it is never used to integrate a SEP trajectory.
// Full particle/source/material continuation is exercised by CP05C/F in amps.
// Keeping a separate executable prevents production or --all-tests from
// obtaining a frozen-line exception to the material-owner guard.
#include "sep.h"
#include "transport_common.h"
#include "util/sep_cli.h"
#include "util/sep_checkpoint_native_validation.h"
#include "util/sep_background_runtime.h"
#include "util/sep_analytic_radial_background.h"
#include "util/sep_run_configuration.h"
#include "util/sep_source_lifecycle.h"
#include "turbulence_production_adapter.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <list>
#include <sstream>
#include <stdexcept>

int main(int argc,char** argv) {
  SEP::Util::CLI::Options options;
  if(!SEP::Util::CLI::ParseCommandLine(argc,argv,options,std::cout,std::cerr))return 2;
  if(options.printHelp){SEP::Util::CLI::PrintHelp(argv[0],std::cout);return 0;}
  MPI_Init(&argc,&argv);MPI_GLOBAL_COMMUNICATOR=MPI_COMM_WORLD;
  MPI_Comm_rank(MPI_COMM_WORLD,&PIC::ThisThread);MPI_Comm_size(MPI_COMM_WORLD,&PIC::nTotalThreads);
  namespace FL=PIC::FieldLine;namespace PB=PIC::ParticleBuffer;
  namespace W=SEP::AlfvenTurbulence_Kolmogorov;namespace T=SEP::Transport;
  namespace Registry=SEP::Testing;
  try {
    auto descriptors=Registry::OperatorCheckpointDescriptors();
    auto fresh=Registry::OperatorFreshCheckpointDescriptors();
    descriptors.insert(descriptors.end(),fresh.begin(),fresh.end());
    for(auto& descriptor:descriptors) {
      descriptor.group=descriptor.id.find("CP05OF")==0?"checkpoint-operators-fresh-process":"checkpoint-operators";
      descriptor.supportedBuildModes="isolated frozen wave/codec operator; no particle motion";
    }
    Registry::Registry registry(descriptors);
    if(options.listTests){if(PIC::ThisThread==0)registry.PrintList(std::cout);MPI_Finalize();return 0;}
    const auto selected=registry.Select(options.testIds,options.testGroups,options.runAllTests,options.testExpressions);
    if(selected.empty())throw std::runtime_error("select checkpoint-operators or one explicit operator ID");
    if(!options.inputPath.empty()||!options.reducedShockEventPath.empty())
      throw std::runtime_error("operator executable accepts test protocols, not a production initialization");
    Registry::ExecutionContext context;context.inputPath=options.testInputPath;
    context.artifactDirectory=options.testArtifactDirectory;Registry::SetExecutionContext(context);
    T::PICAdapter::SetStandaloneRunPurpose(T::PICAdapter::RunPurpose::ControlledDiagnostic);
    auto run=SEP::Run::Defaults();run.turbulence.source=SEP::Turbulence::Source::SelfConsistentIntegrated;
    run.turbulence.coupling=SEP::Turbulence::CouplingPolicy::Disabled;
    auto installed=SEP::Run::InstallActive(run);if(!installed.ok())throw std::runtime_error(installed.message);
    installed=SEP::Turbulence::SetActiveConfiguration(run.turbulence);
    if(!installed.ok())throw std::runtime_error(installed.message);
    SEP::ShockModelType=SEP::cShockModelType::Analytic1D;
    W::ActiveFlag=true;W::ParticleCouplingMode=false;
    // Register precisely the existing integrated-wave datums before pools
    // allocate buffers. No four-polarization authority is advertised here.
    for(auto* datum:{&W::CellIntegratedWaveEnergy,&W::WaveEnergyDensity,
        &W::IsotropicSEP::S,&W::IsotropicSEP::S_pm,&W::G_plus_streaming,
        &W::G_minus_streaming,&W::gamma_plus_array,&W::gamma_minus_array})
      FL::cFieldLineSegment::AddDatumStored(datum);
    auto& a=FL::VertexAllocationManager;
    a.MagneticField=a.PlasmaVelocity=a.PlasmaDensity=a.PlasmaTemperature=a.PlasmaPressure=true;
    a.PreviousVertexData.MagneticField=a.PreviousVertexData.PlasmaVelocity=
        a.PreviousVertexData.PlasmaDensity=a.PreviousVertexData.PlasmaTemperature=a.PreviousVertexData.PlasmaPressure=true;
    PB::OptionalParticleFieldAllocationManager.MomentumParallelNormal=true;
    SEP::RequestParticleData();FL::Init();PB::Init(2048);
    PIC::SimulationTime::SetInitialValue(0);
    PIC::ParticleWeightTimeStep::GlobalTimeStep[0]=0.01;
    PIC::ParticleWeightTimeStep::GlobalTimeStepInitialized=true;
    PIC::ParticleWeightTimeStep::GlobalParticleWeight[0]=2;
    std::list<SEP::cFieldLine> points;
    // Explicit smooth operator coefficients in SI. B follows r^-2 along this
    // radial ray; U=400 km/s and n=5 cm^-3 at 1 AU use the retained analytic
    // declaration. The wave-only Eulerian operator consumes this fixed grid;
    // the process never calls a particle mover, PIC::TimeStep or amps_time_step.
    for(unsigned i=0;i<=32;++i) {
      SEP::cFieldLine point;point.x[0]=1.495978707e11+i*1e8;
      point.U[0]=400000;point.B[0]=5e-9*std::pow(1.495978707e11/point.x[0],2);
      points.push_back(point);
    }
    SEP::Mesh::InitFieldLineAMPS(&points);
    FL::SegmentVolume=SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3;
    PIC::ParallelFieldLines::StaticDecompositionFieldLineLength(0.005);
    installed=SEP::Background::AnalyticRadial::InitializeFixedNativeHistory();
    if(!installed.ok())throw std::runtime_error(installed.message);
    installed=SEP::SourceLifecycle::InitializeNativeFresh();if(!installed.ok())throw std::runtime_error(installed.message);
    SEP::Background::PublishModelOwnedSnapshot(SEP::Background::Provider::Analytic,0,
        std::numeric_limits<double>::infinity(),"isolated frozen wave/source-record codec operator");
    std::ostringstream suppressed;
    auto summary=registry.Run(selected,PIC::ThisThread==0?std::cout:suppressed);
    int code=summary.ExitCode(),global=0;
    MPI_Allreduce(&code,&global,1,MPI_INT,MPI_MAX,MPI_GLOBAL_COMMUNICATOR);
    int reportFailure=0;
    if(PIC::ThisThread==0) {
      std::string error;
      if(!Registry::WriteJsonSummary(summary,options.testJsonPath,&error)||
          !Registry::WriteJUnitSummary(summary,options.testJunitPath,&error)) {
        std::cerr<<error<<'\n';reportFailure=1;
      }
    }
    MPI_Bcast(&reportFailure,1,MPI_INT,0,MPI_GLOBAL_COMMUNICATOR);
    FL::FieldLinesAll[0].Delete();FL::nFieldLine=0;
    PIC::ParallelFieldLines::GenerateThreadSegmentTable();MPI_Finalize();
    return reportFailure?2:global;
  }catch(const std::exception& error) {
    std::cerr<<"isolated checkpoint operator ERROR: "<<error.what()<<'\n';
    MPI_Abort(MPI_GLOBAL_COMMUNICATOR,2);return 2;
  }
}
