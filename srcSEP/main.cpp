

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <list>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <iostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <time.h>

#include <sys/time.h>
#include <sys/resource.h>

//$Id$


//#include "vt_user.h"
//#include <VT.h>

//the particle class
#include "constants.h"
#include "sep.h"
#include "adapters/swcme1d_adapter.h"
#include "transport_common.h"
#include "turbulence_production_adapter.h"
#include "util/sep_cli.h"
#include "util/sep_background_runtime.h"
#include "util/sep_initialization.h"
#include "util/sep_run_configuration.h"
#include "debug/sep_debug_fieldline_datum.h"

#include "tests.h"

void amps_init();
void amps_init_mesh();
void amps_time_step();

/**
 * Prepare and publish SWCME at the one authoritative PIC clock epoch.
 *
 * The former driver maintained a global elapsed time and a separate static CME
 * launch-time accumulator in addition to PIC::SimulationTime.  Their update
 * order made the SWCME cache lag the particle step by one iteration.  This
 * routine has no clock state of its own: it prepares the physical StepState at
 * the supplied PIC epoch, publishes immutable background metadata for exactly
 * the upcoming global interval, and optionally writes the completed-step
 * diagnostic.
 */
void publish_sw1d_for_particle_step(double epoch_seconds, double dt,
                                    long int iteration) {
  // Build the small time-dependent cache behind the adapter before publishing
  // its metadata.  No SWCME type crosses into this AMPS driver translation
  // unit, exactly as in srcSEP3D's provider adapter.
  const SEP::SW1DAdapter::Status preparation =
      SEP::SW1DAdapter::PrepareState(epoch_seconds);
  if (!preparation.ok()) {
    const SEP::SW1DAdapter::Diagnostics diagnostics =
        SEP::SW1DAdapter::GetDiagnostics();
    std::ostringstream message;
    message << "SWCME preparation failed rank=" << PIC::ThisThread
            << " epoch_s=" << epoch_seconds
            << " previous_source_state_id=" << diagnostics.prepared_state_id
            << " detail=" << preparation.detail;
    exit(__LINE__, __FILE__, message.str().c_str());
  }

  SEP::Background::PublishModelOwnedSnapshot(
      SEP::Background::Provider::Swcme, epoch_seconds, epoch_seconds + dt,
      "standalone SWCME StepState and Parker field-line background");

  // One call does everything: n, V, Br, Bphi, |B|, ∇·V → Tecplot POINT file
  if ((PIC::ThisThread==0)&&((iteration+1)%10==0)) {
    char fname[200];

    const int N = 400;
    static double r[N];
    const double rmin = 1.05*_SUN__RADIUS_, rmax = 2.00*_AU_;     // 0.2–2 AU
    for (int i=0;i<N;++i){
      double t = double(i)/(N-1);
      r[i] = rmin*std::pow(rmax/rmin, t);            // log-spacing (nice for r^-2)
    }

    sprintf(fname,"sw_profile_%li.dat",iteration+1);
    SEP::SW1DAdapter::WriteRadialProfileFromR(
        r, N, fname, epoch_seconds);
  }
}

int main(int argc,char **argv) {
  //      MPI_Init(&argc,&argv);

  // --------------------------------------------------------------------------
  // Parse standalone-driver command-line options before the AMPS/SEP model is
  // initialized.  The CLI is intentionally kept in srcSEP/util so the parsing
  // logic is separated from the physics driver and can be reused or extended
  // without cluttering main.cpp.  If the user requests help, exit immediately
  // before reading input files or allocating AMPS data structures.
  // --------------------------------------------------------------------------
  SEP::Util::CLI::Options cli_options;
  if (!SEP::Util::CLI::ParseCommandLine(argc, argv, cli_options, std::cout, std::cerr)) {
    if (PIC::ThisThread == 0) SEP::Util::CLI::PrintHelp(argv[0], std::cerr);
    return 1;
  }

  if (cli_options.printHelp) {
    if (PIC::ThisThread == 0) SEP::Util::CLI::PrintHelp(argv[0], std::cout);
    return 0;
  }

  // Mover discovery is a dependency-free pre-initialization action, parallel
  // to --list-tests.  The listing comes from the same immutable registry used
  // by parsing and runtime dispatch, so it cannot advertise stale movers.
  if (cli_options.listMovers) {
    if (PIC::ThisThread == 0) SEP::Mover::PrintProductionMovers(std::cout);
    return 0;
  }

  // Listing and selection resolution are intentionally performed before shock
  // configuration, post-compile input, AMPS mesh allocation, field-line setup,
  // turbulence initialization, or any output writer.  Thus malformed test
  // requests cannot accidentally fall through into a long production run.
  const SEP::Testing::Registry& componentTestRegistry = ComponentTestRegistry();
  if (cli_options.listTests) {
    if (PIC::ThisThread == 0) componentTestRegistry.PrintList(std::cout);
    return 0;
  }

  const bool componentTestMode =
      SEP::Util::CLI::IsComponentTestExecutionRequested(cli_options);
  std::vector<const SEP::Testing::Descriptor*> selectedComponentTests;
  if (componentTestMode) {
    try {
      selectedComponentTests = componentTestRegistry.Select(
          cli_options.testIds, cli_options.testGroups,
          cli_options.runAllTests);
    }
    catch (const std::exception& exception) {
      if (PIC::ThisThread == 0) {
        std::cerr << "ERROR: cannot resolve component-test selection: "
                  << exception.what() << '\n';
      }
      return 1;
    }

    if (selectedComponentTests.empty()) {
      if (PIC::ThisThread == 0) {
        std::cerr << "ERROR: component-test selection resolved to no tests.\n";
      }
      return 1;
    }

    if (PIC::ThisThread == 0) {
      std::cout << "SEP component-test mode: selected";
      for (const SEP::Testing::Descriptor* descriptor : selectedComponentTests) {
        std::cout << ' ' << descriptor->id;
      }
      std::cout << '\n';
    }

    // Publish the already-validated case paths once before either the
    // initialization-free or field-line-initialized dispatch point. The
    // registry context is deliberately independent of the CV01 implementation
    // so later end-to-end cases reuse the same production CLI contract.
    SEP::Testing::ExecutionContext testContext;
    testContext.inputPath = cli_options.testInputPath;
    testContext.artifactDirectory = cli_options.testArtifactDirectory;
    SEP::Testing::SetExecutionContext(testContext);

    // Initialization-free callbacks are executed immediately.  This makes
    // component tests of pure formulas genuinely lightweight and proves that
    // test-only execution does not require a mesh merely because other catalog
    // entries do.
    if (RequiredInitializationLevel(selectedComponentTests) ==
        SEP::Testing::InitializationLevel::None) {
      return RunSelectedComponentTests(selectedComponentTests, std::cout,
          cli_options.testJsonPath, cli_options.testJunitPath);
    }
  }


  // Parse and freeze the optional mesh/field-line initialization contract
  // before the historical post-compile parser, SWCME setup, MPI, or AMPS mesh
  // allocation.  Failure is therefore atomic: no partially initialized AMPS
  // state survives a malformed scientific input.  An absent --input leaves
  // HasActive()==false and the legacy mesh path remains unchanged.
  if (!cli_options.inputPath.empty()) {
    SEP::Initialization::Configuration initialization;
    const SEP::Transport::Status loaded =
        SEP::Initialization::LoadFile(cli_options.inputPath, &initialization);
    if (!loaded.ok()) {
      if (PIC::ThisThread == 0)
        std::cerr << "ERROR: initialization input: " << loaded.message << '\n';
      return 1;
    }
    if (!cli_options.initializationOutputDirectory.empty()) {
      const SEP::Transport::Status redirected =
          SEP::Initialization::ApplyOutputDirectoryOverride(
              cli_options.initializationOutputDirectory, &initialization);
      if (!redirected.ok()) {
        if (PIC::ThisThread == 0)
          std::cerr << "ERROR: initialization output directory: "
                    << redirected.message << '\n';
        return 1;
      }
    }
    const SEP::Transport::Status installed =
        SEP::Initialization::Install(initialization);
    if (!installed.ok()) {
      if (PIC::ThisThread == 0)
        std::cerr << "ERROR: initialization install: " << installed.message << '\n';
      return 1;
    }
    if (PIC::ThisThread == 0)
      std::cout << "Initialization input=" << cli_options.inputPath
                << " fingerprint="
                << SEP::Initialization::Fingerprint(initialization) << '\n';

    // Version 2 makes the one-dimensional observer and source sampling count
    // part of the same immutable startup contract as the mesh.  The retained
    // sampler consumes heliocentric radius, while the command line remains a
    // higher-precedence compatibility layer for convergence studies.
    if (initialization.schemaVersion >= 2) {
      SEP::Sampling::SamplingHeliocentricDistanceList.clear();
      SEP::Sampling::SamplingHeliocentricDistanceList.push_back(
          initialization.observerHeliocentricRadiusM);
      if (!cli_options.injectionParticlesProvided)
        cli_options.injectionParticlesPerIteration = static_cast<int>(
            initialization.macroparticlesPerStep);
    }
  }


  //read post-compile input file
  if (PIC::PostCompileInputFileName!="") {
     SEP::Parser::ReadFile(PIC::PostCompileInputFileName);
  }

  // D02 resolves the complete canonical SWCME configuration before the AMPS
  // mesh or species storage exists. The request carries only provider-neutral
  // strings; key ownership, unit conversion, preset expansion, and physics
  // validation remain in src/models/swcme. Coupled hosts can construct the
  // same request directly and never invoke this command-line/PARAM parser.
  SEP::SW1DAdapter::ConfigurationRequest swcmeRequest;
  swcmeRequest.preset=cli_options.slowCmeScenario
      ? SEP::SW1DAdapter::Scenario::Slow
      : SEP::SW1DAdapter::Scenario::Fast;
  if (SEP::Initialization::HasActive() &&
      SEP::Initialization::Active().schemaVersion >= 2) {
    const SEP::Initialization::Configuration& initialization =
        SEP::Initialization::Active();
    for (const SEP::Initialization::SwcmeAssignment& raw :
         initialization.swcmeAssignments) {
      SEP::SW1DAdapter::ParameterAssignment assignment;
      assignment.key = raw.key;
      assignment.value = raw.value;
      assignment.origin = cli_options.inputPath;
      assignment.line = raw.line;
      swcmeRequest.input_assignments.push_back(assignment);
    }
  }
  if (cli_options.cmeScenarioProvided) {
    SEP::SW1DAdapter::ParameterAssignment presetAssignment;
    presetAssignment.key="preset";
    presetAssignment.value=cli_options.slowCmeScenario ? "slow" : "fast";
    presetAssignment.origin="command-line --cme-scenario";
    swcmeRequest.command_line_assignments.push_back(presetAssignment);
  }
  for (std::size_t i=0;i<cli_options.swcmeOverrides.size();++i) {
    const std::string& raw=cli_options.swcmeOverrides[i];
    const std::size_t separator=raw.find('=');
    SEP::SW1DAdapter::ParameterAssignment assignment;
    assignment.key=raw.substr(0,separator);
    assignment.value=raw.substr(separator+1);
    assignment.origin="command-line --swcme-override";
    assignment.line=i+1;
    swcmeRequest.command_line_assignments.push_back(assignment);
  }
  SEP::SW1DAdapter::Status swcmeConfigurationStatus=
      SEP::SW1DAdapter::Configure(swcmeRequest);
  if (!swcmeConfigurationStatus.ok()) {
    if (PIC::ThisThread==0)
      std::cerr<<"ERROR: invalid canonical SWCME configuration: "
               <<swcmeConfigurationStatus.detail<<'\n';
    return 1;
  }
  const SEP::SW1DAdapter::ConfigurationSummary swcmeSummary=
      SEP::SW1DAdapter::GetConfigurationSummary();
  if (SEP::Initialization::HasActive() &&
      SEP::Initialization::Active().schemaVersion >= 2) {
    const SEP::Initialization::Configuration& initialization =
        SEP::Initialization::Active();
    const auto samePhysicalScalar = [](double left, double right) {
      return std::fabs(left - right) <=
          64.0 * std::numeric_limits<double>::epsilon() *
          std::max(std::fabs(left), std::fabs(right));
    };
    const double relativeX = initialization.parkerInitialPointM.x -
        initialization.parkerOriginM.x;
    const double relativeY = initialization.parkerInitialPointM.y -
        initialization.parkerOriginM.y;
    const double relativeZ = initialization.parkerInitialPointM.z -
        initialization.parkerOriginM.z;
    const double sourceRadius = std::sqrt(
        relativeX * relativeX + relativeY * relativeY + relativeZ * relativeZ);
    const double lineSinTheta =
        std::sqrt(relativeX * relativeX + relativeY * relativeY) / sourceRadius;
    if (!samePhysicalScalar(initialization.solarWindSpeedMPerS,
                            swcmeSummary.ambient_wind_speed_m_per_s) ||
        !samePhysicalScalar(initialization.solarRotationRateRadPerS,
                            swcmeSummary.solar_rotation_rate_rad_per_s) ||
        !samePhysicalScalar(initialization.innerRadiusM,
                            swcmeSummary.parker_source_radius_m) ||
        !samePhysicalScalar(lineSinTheta,
                            swcmeSummary.parker_reference_sin_theta) ||
        !std::isfinite(swcmeSummary.parker_radial_field_at_one_au_t) ||
        swcmeSummary.parker_radial_field_at_one_au_t == 0.0) {
      if (PIC::ThisThread == 0)
        std::cerr << "ERROR: initialization Parker wind/rotation/source/"
                     "latitude differs from canonical SWCME, or the resolved "
                     "magnetic normalization is zero\n";
      return 1;
    }
  }

  // Preserve a constant Dmumu value supplied by the post-compile input unless
  // the command line explicitly overrides it.  Keeping the effective value in
  // the options record also makes the startup configuration print truthful.
  if (!cli_options.constantDmumuProvided) {
    cli_options.constantDmumuPerS =
        SEP::Diffusion::ConstPitchAngleDiffusionValue;
    cli_options.coefficients.constantDmumuPerS =
        cli_options.constantDmumuPerS;
  }

  // WP30 freezes the effective run contract after defaults, post-compile input,
  // and CLI overrides have all been resolved.  Downstream code receives only a
  // const view and the fingerprint is emitted before mesh/model initialization.
  SEP::Run::Configuration runConfiguration=SEP::Run::Defaults();
  runConfiguration.mover.value=cli_options.particleMover;
  runConfiguration.mover.source=cli_options.particleMoverProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  runConfiguration.shockModel.value=cli_options.analyticalShock
      ? SEP::Run::ShockModel::Analytical : SEP::Run::ShockModel::Swcme1d;
  if (!cli_options.shockModelProvided)
    runConfiguration.shockModel.value=SEP::ShockModelType==
        SEP::cShockModelType::Analytic1D ? SEP::Run::ShockModel::Analytical
                                        : SEP::Run::ShockModel::Swcme1d;
  runConfiguration.shockModel.source=cli_options.shockModelProvided
      ? SEP::Run::ValueSource::CommandLine
      : (PIC::PostCompileInputFileName!="" ? SEP::Run::ValueSource::InputFile
                                           : SEP::Run::ValueSource::Default);
  runConfiguration.scenario.value=swcmeSummary.preset=="SLOW"
      ? SEP::Run::CmeScenario::Slow : SEP::Run::CmeScenario::Fast;
  runConfiguration.scenario.source=cli_options.cmeScenarioProvided
      ? SEP::Run::ValueSource::CommandLine
      : (PIC::PostCompileInputFileName!="" ? SEP::Run::ValueSource::InputFile
                                           : SEP::Run::ValueSource::Default);
  switch (cli_options.swcmeFailurePolicy) {
    case SEP::Util::CLI::Options::SwcmeFailurePolicy::Strict:
      runConfiguration.swcmeFailurePolicy.value=
          SEP::Run::SwcmeFailurePolicy::Strict;
      break;
    case SEP::Util::CLI::Options::SwcmeFailurePolicy::ClampRadius:
      runConfiguration.swcmeFailurePolicy.value=
          SEP::Run::SwcmeFailurePolicy::ClampRadius;
      break;
    case SEP::Util::CLI::Options::SwcmeFailurePolicy::DiagnosticFallback:
      runConfiguration.swcmeFailurePolicy.value=
          SEP::Run::SwcmeFailurePolicy::DiagnosticFallback;
      break;
  }
  runConfiguration.swcmeFailurePolicy.source=
      cli_options.swcmeFailurePolicyProvided
          ? SEP::Run::ValueSource::CommandLine
          : SEP::Run::ValueSource::Default;
  runConfiguration.swcmeFallbackDensityM3.value=
      cli_options.swcmeFallbackDensityM3;
  runConfiguration.swcmeFallbackDensityM3.source=
      cli_options.swcmeFallbackDensityProvided
          ? SEP::Run::ValueSource::CommandLine
          : SEP::Run::ValueSource::Default;
  runConfiguration.swcmeFallbackSpeedMPerS.value=
      cli_options.swcmeFallbackSpeedMPerS;
  runConfiguration.swcmeFallbackSpeedMPerS.source=
      cli_options.swcmeFallbackSpeedProvided
          ? SEP::Run::ValueSource::CommandLine
          : SEP::Run::ValueSource::Default;
  runConfiguration.swcmeFallbackDivergencePerS.value=
      cli_options.swcmeFallbackDivergencePerS;
  runConfiguration.swcmeFallbackDivergencePerS.source=
      cli_options.swcmeFallbackDivergenceProvided
          ? SEP::Run::ValueSource::CommandLine
          : SEP::Run::ValueSource::Default;
  runConfiguration.swcmeConfigurationFingerprint.value=
      swcmeSummary.fingerprint;
  runConfiguration.swcmeConfigurationFingerprint.source=
      (!cli_options.swcmeOverrides.empty() || cli_options.cmeScenarioProvided)
          ? SEP::Run::ValueSource::CommandLine
          : (PIC::PostCompileInputFileName!=""
              ? SEP::Run::ValueSource::InputFile
              : SEP::Run::ValueSource::Default);
  runConfiguration.totalIterations.value=static_cast<std::uint64_t>(
      cli_options.totalIterations);
  runConfiguration.totalIterations.source=cli_options.totalIterationsProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  runConfiguration.fieldLineSeedAreaM2.value=cli_options.fieldLineSeedAreaM2;
  runConfiguration.fieldLineSeedAreaM2.source=cli_options.fieldLineSeedAreaProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  runConfiguration.shockTurbulenceEfficiency.value=
      cli_options.shockTurbulenceEfficiency;
  runConfiguration.shockTurbulenceEfficiency.source=
      cli_options.shockTurbulenceEfficiencyProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  runConfiguration.shockTurbulencePlusFraction.value=
      cli_options.shockTurbulencePlusFraction;
  runConfiguration.shockTurbulencePlusFraction.source=
      cli_options.shockTurbulencePlusFractionProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  runConfiguration.mergeMinimum.value=cli_options.mergeMinimum;
  runConfiguration.mergeMinimum.source=cli_options.mergeMinimumProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  runConfiguration.mergeMaximum.value=cli_options.mergeMaximum;
  runConfiguration.mergeMaximum.source=cli_options.mergeMaximumProvided
      ? SEP::Run::ValueSource::CommandLine : SEP::Run::ValueSource::Default;
  // WP33 keeps the legacy CLI aliases but freezes the complete population
  // policy in one record.  The PIC merge/split calls below consume this record,
  // so a printed option cannot diverge from the thresholds actually executed.
  runConfiguration.populationControl.minimumParticlesPerCell=
      runConfiguration.mergeMinimum.value;
  runConfiguration.populationControl.maximumParticlesPerCell=
      runConfiguration.mergeMaximum.value;
  runConfiguration.coefficients=cli_options.coefficients;
  runConfiguration.numericalTolerances=cli_options.numericalTolerances;
  runConfiguration.turbulence=cli_options.turbulence;
  runConfiguration.injection.macroparticlesPerEvent=
      static_cast<std::uint64_t>(cli_options.injectionParticlesPerIteration);
  runConfiguration.injection.injectionEfficiency=
      swcmeSummary.source_injection_efficiency;
  const SEP::Transport::Status runStatus=
      SEP::Run::InstallActive(runConfiguration);
  if (!runStatus.ok()) {
    if (PIC::ThisThread==0)
      std::cerr<<"ERROR: invalid final RunConfiguration: "
               <<runStatus.message<<'\n';
    return 1;
  }
  const SEP::Run::FrozenConfiguration& frozenRun=SEP::Run::Active();
  SEP::ShockModelType=frozenRun.get().shockModel.value==
      SEP::Run::ShockModel::Swcme1d ? SEP::cShockModelType::SwCme1d
                                    : SEP::cShockModelType::Analytic1D;
  SEP::FieldLine::FluxTubeGeometry::SetReferenceAreaM2(
      frozenRun.get().fieldLineSeedAreaM2.value);
  if (PIC::ThisThread==0)
    std::cout<<"RunConfiguration fingerprint="<<frozenRun.fingerprint()
             <<" seed-area="<<frozenRun.get().fieldLineSeedAreaM2.value
             <<" m2 source="
             <<SEP::Run::ValueSourceName(
                   frozenRun.get().fieldLineSeedAreaM2.source)<<'\n'
             <<"SWCME configuration fingerprint="<<swcmeSummary.fingerprint
             <<" preset="<<swcmeSummary.preset<<'\n'
             <<swcmeSummary.normalized_manifest;


  //set up shock wave model
  SEP::SW1DAdapter::EnableSheathClamp(true);
  // Bridge the canonical source interval/efficiency into the existing 1-D
  // field-line injector. Energies remain MeV at this public boundary and are
  // converted exactly once to joules inside field_line.cpp.
  SEP::FieldLine::InjectionParameters::emin=swcmeSummary.source_energy_min_MeV;
  SEP::FieldLine::InjectionParameters::emax=swcmeSummary.source_energy_max_MeV;
  SEP::FieldLine::InjectionParameters::InjectionEfficiency=
      swcmeSummary.source_injection_efficiency;
  // Install the already-frozen recovery contract. The fallback sample is
  // validated even in strict mode so a restart fingerprint can never preserve
  // malformed dormant values that become active in a later run.
  SEP::SW1DAdapter::Status adapterStatus=
      SEP::SW1DAdapter::SetDiagnosticFallback(
          SEP::SW1DAdapter::BackgroundSample(
              frozenRun.get().swcmeFallbackDensityM3.value,
              frozenRun.get().swcmeFallbackSpeedMPerS.value,
              frozenRun.get().swcmeFallbackDivergencePerS.value));
  if (!adapterStatus.ok()) {
    if (PIC::ThisThread==0)
      std::cerr<<"ERROR: invalid SWCME fallback: "<<adapterStatus.detail<<'\n';
    return 1;
  }
  SEP::SW1DAdapter::FailurePolicy adapterPolicy=
      SEP::SW1DAdapter::FailurePolicy::Strict;
  if (frozenRun.get().swcmeFailurePolicy.value==
      SEP::Run::SwcmeFailurePolicy::ClampRadius)
    adapterPolicy=SEP::SW1DAdapter::FailurePolicy::ClampRadius;
  else if (frozenRun.get().swcmeFailurePolicy.value==
           SEP::Run::SwcmeFailurePolicy::DiagnosticFallback)
    adapterPolicy=SEP::SW1DAdapter::FailurePolicy::DiagnosticFallback;
  adapterStatus=SEP::SW1DAdapter::SetFailurePolicy(adapterPolicy);
  if (!adapterStatus.ok()) {
    if (PIC::ThisThread==0)
      std::cerr<<"ERROR: invalid SWCME failure policy: "
               <<adapterStatus.detail<<'\n';
    return 1;
  }

  // Prepare the initial SWCME cache before mesh/field-line initialization can
  // query shock geometry.  Metadata publication waits until immediately before
  // the first particle step, after all background-affecting CLI options and
  // field-line/turbulence initialization have completed.
  {
    const double initial_epoch=SEP::Background::SimulationTimeSeconds();
    const SEP::SW1DAdapter::Status preparation=
        SEP::SW1DAdapter::PrepareState(initial_epoch);
    if (!preparation.ok()) {
      if (PIC::ThisThread==0)
        std::cerr << "ERROR: SWCME initial preparation failed rank="
                  << PIC::ThisThread << " epoch_s=" << initial_epoch
                  << " detail=" << preparation.detail << '\n';
      return 1;
    }
  }

  //output parameters of the sshock
  // The production shock diagnostic is unrelated to component-test setup and
  // would create an unrequested shared artifact.  Field-line tests still receive
  // the configured SWCME model, but only a production run writes this file.
  if (!componentTestMode && !cli_options.initializationOnly) {
    SEP::SW1DAdapter::WriteShockVsTime(
        2.0*24.0*3600, 200, "shock_vs_time.dat");
  }

  // --------------------------------------------------------------------------
  // Configure the optional turbulence physics from command-line options before
  // registering AMPS field-line datums.  The selected turbulence model affects
  // which segment datums must be allocated.  In particular, the new
  // wave-number-resolved model needs an additional 2*NK spectral-energy datum.
  // --------------------------------------------------------------------------
  // In `configured` mode the post-compile input remains authoritative.  An
  // explicit CLI provider is applied below.  Do not overwrite either selection
  // here: doing so previously made startup metadata disagree with the source
  // that the user had requested.
  SEP::Util::CLI::ApplyTurbulenceOptions(cli_options);
  if (PIC::ThisThread == 0) {
    SEP::Util::CLI::PrintTurbulenceOptions(cli_options, std::cout);
    SEP::Mover::PrintRuntimeConfiguration(std::cout);
  }

  //setup datum to store the segment's data for the Alfven turbulence model
  if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity);

    // The wave-number-resolved model stores E+(k_j) and E-(k_j) in an
    // additional hidden segment datum.  It is registered only when selected by
    // the CLI to avoid increasing memory in legacy integrated-turbulence runs.
    if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
      PIC::FieldLine::cFieldLineSegment::AddDatumStored(
          &SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy);

      // Hidden diagnostic datum used by the 2-D spectrum writer.  It stores the
      // per-bin wave-energy exchange rates caused by cascade, particle coupling,
      // and reflection during the current main-loop iteration.  It is registered
      // only in wave-number-resolved mode to avoid extra memory in legacy runs.
      PIC::FieldLine::cFieldLineSegment::AddDatumStored(
          &SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergyExchangeRate);
    }

    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S_pm);


    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::gamma_plus_array);
    PIC::FieldLine::cFieldLineSegment::AddDatumStored(&SEP::AlfvenTurbulence_Kolmogorov::gamma_minus_array);

  }

  //set up datum to store distance of a field line vertex to the location of the shock
  PIC::FieldLine::UserDefinedfDataProcessingManager=SEP::FieldLine::CalculateVertexShockDistances;

  amps_init_mesh();
  amps_init();

  // The canonical source record describes one physical species. srcSEP may use
  // several AMPS species indices as numerical populations, but each must match
  // that mass and signed charge. A mixed-species campaign requires the existing
  // explicit SpeciesSource table and is deliberately not inferred from one
  // SWCME source record.
  if (SEP::ShockModelType==SEP::cShockModelType::SwCme1d) {
    const double elementaryChargeC=1.602176634e-19;
    for (int spec=0;spec<PIC::nTotalSpecies;++spec) {
      const double actualMass=PIC::MolecularData::GetMass(spec);
      const double actualCharge=PIC::MolecularData::GetElectricCharge(spec);
      const double expectedCharge=
          swcmeSummary.source_charge_number*elementaryChargeC;
      const double massScale=std::max(
          std::fabs(actualMass),std::fabs(swcmeSummary.source_particle_mass_kg));
      const double chargeScale=std::max(
          std::fabs(actualCharge),std::fabs(expectedCharge));
      if (std::fabs(actualMass-swcmeSummary.source_particle_mass_kg)>
              1.0e-5*std::max(massScale,1.0e-40) ||
          std::fabs(actualCharge-expectedCharge)>
              1.0e-10*std::max(chargeScale,1.0e-30)) {
        std::ostringstream message;
        message<<"SWCME source species does not match AMPS species "<<spec
               <<": expected mass="<<swcmeSummary.source_particle_mass_kg
               <<" kg charge="<<expectedCharge
               <<" C; actual mass="<<actualMass
               <<" kg charge="<<actualCharge<<" C";
        exit(__LINE__,__FILE__,message.str().c_str());
      }
    }
  }

  // This boundary is intentionally after amps_init() and the canonical AMPS
  // species check: mesh construction, block allocation, field-line creation,
  // observer installation, time-step/weight initialization, and initialization
  // Tecplot output have all succeeded.  It is intentionally before turbulence
  // evolution and the first amps_time_step(), so a mesh preview cannot inject
  // or move a particle.  All ranks synchronize before finalizing MPI.
  if (cli_options.initializationOnly) {
    MPI_Barrier(MPI_GLOBAL_COMMUNICATOR);
    if (PIC::ThisThread == 0) {
      const SEP::Initialization::Configuration& initialization =
          SEP::Initialization::Active();
      std::cout << "srcSEP initialization complete; no time steps executed\n"
                << "initialization_mesh="
                << initialization.meshTecplotFile << '\n'
                << "initialization_field_line="
                << initialization.fieldLineTecplotFile << '\n';
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  // Turbulence storage is initialized exactly once below, after the optional
  // test-only exit.  The former ModelInit::Init() call wrote an independent
  // hard-coded profile here and was then overwritten by a second initializer.
  // More importantly, it wrote even when SWMF or a prescribed provider owned
  // the wave state.  Deferring to the source-aware block below gives every
  // field exactly one declared startup owner.

  // Selected component tests own the process after their declared field-line
  // prerequisite is available.  Returning here is the critical test-only
  // boundary: neither the compatibility TestManager path nor the production
  // timestep loop can execute after a --test/--test-group/--all-tests request.
  if (componentTestMode) {
    if (RequiredInitializationLevel(selectedComponentTests) ==
        SEP::Testing::InitializationLevel::FieldLineModel) {
      // Native mover fixtures call the same production coefficient adapters as
      // an ordinary particle step.  Publish one immutable snapshot before the
      // registry starts, but do not enter its read phase here: each fixture must
      // first install and later restore its temporary vertex values.  The
      // fixture itself opens ParticleReadPhase only around actual mover calls.
      const double epochS = SEP::Background::SimulationTimeSeconds();
      if (SEP::Background::ConfiguredProvider() ==
          SEP::Background::Provider::Swcme) {
        publish_sw1d_for_particle_step(
            epochS, PIC::ParticleWeightTimeStep::GlobalTimeStep[0], 0);
      }
      else {
        SEP::Background::PrepareSnapshotForParticleStep();
      }
    }
    return RunSelectedComponentTests(selectedComponentTests, std::cout,
        cli_options.testJsonPath, cli_options.testJunitPath);
  }

  // --------------------------------------------------------------------------
  // Optional development diagnostics.
  //
  // TestManager() performs standalone field-line/model tests and is useful when
  // debugging the SEP/turbulence implementation.  It should not run by default
  // in production simulations because it can add extra diagnostic work/output
  // and may change the intended run flow.  The CLI therefore leaves it OFF
  // unless explicitly requested with one of:
  //   --test-manager on
  //   --testmanager on
  //   --run-test-manager
  // Step 5 makes field-line infrastructure a compile-time precondition, so no
  // runtime fallback can silently skip a requested diagnostic.
  // --------------------------------------------------------------------------
  if (cli_options.runTestManager) {
    TestManager();
  }

  const long int TotalIterations=(_PIC_NIGHTLY_TEST_MODE_==_PIC_MODE_ON_)
      ? static_cast<long int>(PIC::RequiredSampleLength+10)
      : static_cast<long int>(frozenRun.get().totalIterations.value);

  // Initialize wave storage only for a source that is locally owned from the
  // beginning of the run.  Prescribed and SWMF-read-only providers remain
  // authoritative, while SwmfInitialThenEvolveLocal must first copy the
  // imported generation through the production adapter's ImportFieldLine path
  // and record its one-time handoff; treating that mode as an ordinary startup
  // would destroy the very SWMF state it is required to inherit.  The
  // configurable prescribed amplitude is reused as the initial self-consistent
  // deltaB/B; the named 1-AU field constant is only a fallback for malformed or
  // missing local magnetic-field values inside the legacy initializer.
  const SEP::Turbulence::Source turbulenceSource =
      SEP::Turbulence::ActiveConfiguration().source;
  const bool locallyOwnedFromStartup =
      turbulenceSource ==
          SEP::Turbulence::Source::SelfConsistentIntegrated ||
      turbulenceSource ==
          SEP::Turbulence::Source::SelfConsistentSpectral;
  const bool initializeLocalTurbulence =
      SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
      locallyOwnedFromStartup;
  if (initializeLocalTurbulence) {
    const double initialDeltaBOverB =
        SEP::Transport::Coefficient::ActiveConfiguration().
            prescribedDeltaBOverB;

    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,0);
    SEP::AlfvenTurbulence_Kolmogorov::InitializeWaveEnergyFromPhysicalParameters(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
        SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyConstants::TYPICAL_B0_1AU,
        initialDeltaBOverB, initialDeltaBOverB, -2.0, false, true);

    // The initializer writes owned edge segments.  Gather to rank zero and
    // broadcast the complete generation once, before any boundary reservoir is
    // captured, so every rank begins from the same authoritative state.
    PIC::FieldLine::Parallel::MPIGatherDatumStoredAtEdge(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,0);
    PIC::FieldLine::Parallel::MPIBcastDatumStoredAtEdge(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,0);
  }

  // Capture the right-boundary W- initial condition after the turbulence wave
  // energy has been initialized and the edge data have been synchronized.
  //
  // Legacy/integrated model:
  //   store one W- density at the last segment.
  //
  // Wave-number-resolved model:
  //   first expand the integrated E+,E- initial condition into E+(k_j),E-(k_j)
  //   using the Kolmogorov log-bin weights, then store the full W-(k_j) density
  //   at the last segment.  This keeps the right boundary fixed as a
  //   pre-existing spectral turbulence reservoir rather than an artificial
  //   time-growing source.
  if (initializeLocalTurbulence &&
      SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
    SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::InitializeSpectrumFromIntegratedEnergy(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
    SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::ResetRightBoundarySpectrumInitialCondition();
    SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::CaptureRightBoundarySpectrumInitialCondition();
    SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::EnforceRightBoundarySpectrumInitialCondition();
    SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::UpdateIntegratedEnergyFromSpectrum(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
  }
  else if (initializeLocalTurbulence) {
    SEP::AlfvenTurbulence_Kolmogorov::ResetRightBoundaryEminusInitialCondition();
    SEP::AlfvenTurbulence_Kolmogorov::CaptureRightBoundaryEminusInitialCondition(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
    SEP::AlfvenTurbulence_Kolmogorov::EnforceRightBoundaryEminusInitialCondition(
        SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
  }


  //set background plasma density
auto set_background_plasma_density = []() {
    // reference density at 1 AU [m⁻³]
    constexpr double n0 = 5.0e6;
    // 1 astronomical unit in meters
    constexpr double AU = 1.495978707e11;

    // Loop over all field lines
    for (int fldIdx = 0; fldIdx < PIC::FieldLine::nFieldLine; ++fldIdx) {
        auto fieldLine = &PIC::FieldLine::FieldLinesAll[fldIdx];
        if (!fieldLine) continue;

        int nSeg = fieldLine->GetTotalSegmentNumber();
        if (nSeg < 1) continue;

        // Loop over all segments in this field line
        for (int segIdx = 0; segIdx < nSeg; ++segIdx) {
            auto seg = fieldLine->GetSegment(segIdx);
            if (!seg) continue;

            // Process both end‐points (left & right) of the segment
            PIC::FieldLine::cFieldLineVertex* vertices[2] = { seg->GetBegin(), seg->GetEnd() };
            for (auto vtx : vertices) {
                if (!vtx) continue;

                // get pointer to {x,y,z} [m]
                double* X = vtx->GetX();
                // compute radial distance from Sun [m]
                double r = std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);

                // scale density as n0/(r/AU)^2
                double density = n0 / ((r/AU) * (r/AU));

                // store it in the vertex
                vtx->SetPlasmaDensity(density); //SetDatum(density, PIC::FieldLine::DatumAtVertexPlasmaDensity);
            }
        }
    }
};


 // Analytic Parker runs own their plasma profile and may create the r^-2
 // density used by the legacy field-line geometry.  SWCME and SWMF provide
 // authoritative plasma state; overwriting either provider here previously
 // made the snapshot metadata disagree with the arrays consumed by movers.
 if (SEP::Background::ConfiguredProvider() ==
     SEP::Background::Provider::Analytic) {
   set_background_plasma_density();
 }

  // Calculate turbulence wave enregy density from wave energy integrated over the segments of the magnetic tube:
auto CalculateWaveEnergyDensity = [&]() {
    // Focus specifically on SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy
    auto& integrated_energy = SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy;
    auto& energy_density = SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity;

    int local_count = 0;

    for (int fl = 0; fl < PIC::FieldLine::nFieldLine; fl++) {
        auto& field_line = PIC::FieldLine::FieldLinesAll[fl];
        int num_segs = field_line.GetTotalSegmentNumber();

        for (int seg = 0; seg < num_segs; seg++) {
            auto* segment = field_line.GetSegment(seg);
            if (!segment) continue;

            // Only process segments belonging to this thread
            if (segment->Thread == PIC::ThisThread) {
                // Focus on getting CellIntegratedWaveEnergy values
                double* energy_data = segment->GetDatum_ptr(integrated_energy);
                double* density_data = segment->GetDatum_ptr(energy_density);

                // Use SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3 for volume calculation
                double volume = SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3(segment, fl);

                if (energy_data && density_data && volume > 0.0) {
                    // CellIntegratedWaveEnergy stores the conservative variables
                    // E+ and E- integrated over the current magnetic-tube segment.
                    // WaveEnergyDensity is the printable AMPS field-line datum used
                    // in amps.FieldLines.out=*.dat.  Its first two elements must be
                    // the plotted wave-energy densities W+ and W-, while its third
                    // element is the normalized cross helicity sigma_c.  We keep
                    // sigma_c in the same datum as W+ and W- rather than registering
                    // a separate output datum so that the generic field-line writer
                    // prints the three turbulence diagnostics as one adjacent block:
                    //   "W+", "W-", "sigma_c".
                    const double Wplus  = energy_data[0] / volume;
                    const double Wminus = energy_data[1] / volume;

                    density_data[0] = Wplus;
                    density_data[1] = Wminus;

                    // sigma_c is the normalized Elsasser/turbulence imbalance.  It
                    // is bounded by [-1,1] for non-negative W+ and W-.  A zero value
                    // is used when both wave populations vanish to avoid division by
                    // zero and to keep the output finite.
                    const double Wsum = Wplus + Wminus;
                    density_data[2] = (Wsum > 0.0) ? (Wplus - Wminus) / Wsum : 0.0;

                    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
                      validate_numeric(density_data[0],__LINE__,__FILE__);
                      validate_numeric(density_data[1],__LINE__,__FILE__);
                      validate_numeric(density_data[2],__LINE__,__FILE__);
                    }

                    local_count++;
                }
            }
        }
    }

    // MPI operations to gather and broadcast results
    PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(energy_density);

    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,0);
    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,2);

    if (PIC::ThisThread == 0) {
        std::cout << "Compact wave energy density calculation completed with MPI operations" << std::endl;
    }
};


 //calculate the wave energy density
 CalculateWaveEnergyDensity();

    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,0);
    SEP::AlfvenTurbulence_Kolmogorov::TestPrintEPlusValues(SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,2);


    //PIC::FieldLine::Output("fl-edge-test.dat",false);

//hooks for calculating the magnertic tube radius and the volume of the field line segment
PIC::FieldLine::SegmentVolume=SEP::FieldLine::FluxTubeGeometry::SegmentVolumeM3;



  vector<vector<double> > DeltaE_plus, DeltaE_minus;


  // --------------------------------------------------------------------------
  // ResetWaveParticleStreamingAccumulators()
  // --------------------------------------------------------------------------
  // G_plus_streaming and G_minus_streaming are not physical wave-energy state
  // variables.  They are one-time Monte-Carlo accumulators filled by the particle
  // mover during the current AMPS time step and then MPI-summed before the
  // turbulence wave-particle coupling manager is called.  Therefore they must be
  // zeroed on every MPI rank, on every local/ghost copy of every field-line
  // segment, before particles are moved.  If stale values remain on non-owning
  // ranks, MPIAllReduceDatumStoredAtEdge() can repeatedly re-sum old source terms
  // and produce artificially large growth/damping rates.
  //
  // Keep the reset local to main.cpp rather than relying on a coupling manager to
  // clear only owned segments after use.  The particle mover is the producer of
  // these source terms, so the safe place to clear them is immediately before the
  // particle mover is entered through amps_time_step().
  // --------------------------------------------------------------------------
  auto ResetWaveParticleStreamingAccumulators = []() {
    auto ResetDatum = [](PIC::Datum::cDatumStored& Datum) {
      for (int iFieldLine=0; iFieldLine<PIC::FieldLine::nFieldLine; ++iFieldLine) {
        for (PIC::FieldLine::cFieldLineSegment* Segment =
                 PIC::FieldLine::FieldLinesAll[iFieldLine].GetFirstSegment();
             Segment != NULL;
             Segment = Segment->GetNext()) {
          double* data = Segment->GetDatum_ptr(Datum);
          if (data == NULL) continue;

          for (int i=0; i<Datum.length; ++i) data[i]=0.0;
        }
      }
    };

    ResetDatum(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
    ResetDatum(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);
    ResetDatum(SEP::AlfvenTurbulence_Kolmogorov::gamma_plus_array);
    ResetDatum(SEP::AlfvenTurbulence_Kolmogorov::gamma_minus_array);
  };


  //time step
  for (long int niter=0;niter<TotalIterations;niter++) {
    // Freeze the complete background identity for the upcoming particle step.
    // SWCME state, snapshot epoch, validity interval, and the PIC clock are
    // published together.  All later turbulence/shock updates occur after the
    // ParticleReadPhase in amps_time_step() has ended.
    const double global_dt =
        PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
    // SWCME is the only provider whose physical cache must be prepared by the
    // standalone driver before the common application step.  Analytic and SWMF
    // runs are prepared by PrepareSnapshotForParticleStep() in amps_time_step;
    // publishing an SWCME generation for either would violate provider
    // authority and is rejected by PublishModelOwnedSnapshot().
    if (SEP::Background::ConfiguredProvider() ==
        SEP::Background::Provider::Swcme) {
      publish_sw1d_for_particle_step(
          SEP::Background::SimulationTimeSeconds(), global_dt, niter);
    }

    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag &&
        SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode) {
      ResetWaveParticleStreamingAccumulators();
    }

    // ----------------------------------------------------------------------
    // Debug-only validation of field-line datums at the very beginning of the
    // main iteration.
    //
    // This check is intentionally placed before amps_time_step(), shock
    // injection, particle/turbulence coupling, advection, reflection, cascade,
    // and all MPI synchronization performed later in the iteration.  If an AMPS
    // field-line MPI unpack routine later finds a non-finite or overflowed
    // value, these pre-iteration checks help determine whether the bad value
    // already existed in the local field-line state or was created by one of the
    // subsequent operators / MPI reductions in the current iteration.
    //
    // The helper takes the datum as an argument, so additional datums can be
    // checked by adding one more call here.  The calls are protected by the AMPS
    // debugger-mode macro and are therefore inactive in normal production runs.
    // ----------------------------------------------------------------------
    if (_PIC_DEBUGGER_MODE_ == _PIC_DEBUGGER_MODE_ON_) {
      SEP::Debug::ValidateFieldLineDatum(
          SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
          "CellIntegratedWaveEnergy",
          "beginning of main iteration",
          niter,
          __LINE__,
          __FILE__);

      SEP::Debug::ValidateFieldLineDatum(
          SEP::AlfvenTurbulence_Kolmogorov::WaveEnergyDensity,
          "WaveEnergyDensity",
          "beginning of main iteration",
          niter,
          __LINE__,
          __FILE__);

      SEP::Debug::ValidateFieldLineDatum(
          SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,
          "G_plus_streaming",
          "beginning of main iteration",
          niter,
          __LINE__,
          __FILE__);

      SEP::Debug::ValidateFieldLineDatum(
          SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,
          "G_minus_streaming",
          "beginning of main iteration",
          niter,
          __LINE__,
          __FILE__);

      if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
        SEP::Debug::ValidateFieldLineDatum(
            SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy,
            "WaveNumberResolved::SpectralWaveEnergy",
            "beginning of main iteration",
            niter,
            __LINE__,
            __FILE__);

        SEP::Debug::ValidateFieldLineDatum(
            SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergyExchangeRate,
            "WaveNumberResolved::SpectralWaveEnergyExchangeRate",
            "beginning of main iteration",
            niter,
            __LINE__,
            __FILE__);
      }
    }

    static double rsh0=SEP::ParticleSource::ShockWave::Tenishev2005::rShock;

    if (niter==0) {
      switch (SEP::ShockModelType) {
      case SEP::cShockModelType::Analytic1D:
        rsh0=SEP::ParticleSource::ShockWave::Tenishev2005::rShock;
        break;
      case SEP::cShockModelType::SwCme1d:
        rsh0=SEP::SW1DAdapter::ShockRadiusM();
        break;
      }
    }


    amps_time_step();

    // amps_time_step() owns the single post-particle turbulence transaction for
    // both standalone and coupled execution.  Do not advance turbulence again
    // here: doing so previously applied advection, reflection, cascade, and
    // particle-wave exchange twice per global particle step.  Shock-radius
    // history is now captured at that common transaction boundary as well.

    // Retain the former orchestration temporarily as unreachable migration
    // evidence.  Keeping it beside the adapter call makes review against old
    // runs straightforward; no production control path can enter it.
    if (false && SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {
      // Source ownership is independent of the selected particle mover.
      // Prescribed and SWMF-read-only sources may be synchronized and sampled
      // below, but every local mutating operator is gated by this one contract.
      const bool evolve_turbulence_locally = SEP::Turbulence::EvolvesLocally(
          SEP::Turbulence::ActiveConfiguration().source);

      // Function to increment integrated wave energy due to shock passing
      if (evolve_turbulence_locally &&
          SEP::Turbulence::ActiveConfiguration().shockInjectionEnabled &&
          niter!=0) {
         double rsh1;

         switch (SEP::ShockModelType) {
         case SEP::cShockModelType::Analytic1D:
           rsh1=SEP::ParticleSource::ShockWave::Tenishev2005::rShock;
           break;
         case SEP::cShockModelType::SwCme1d:
           rsh1=SEP::SW1DAdapter::ShockRadiusM();
           break;
         }

	 SEP::ParticleSource::ShockWave::ShockTurbulenceEnergyInjection(rsh0, rsh1, PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);

         // The present shock turbulence source is formulated for the legacy
         // branch-integrated E+ and E- datum.  In the wave-number-resolved
         // model, immediately project the updated integrated shock increment
         // back to E±(k_j), preserving the local spectral shape where possible
         // and using a Kolmogorov distribution only when a branch previously had
         // no spectral energy.
         if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
           SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::ProjectIntegratedEnergyToSpectrum(
               SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
         }

	 rsh0=rsh1;
      }

      if (evolve_turbulence_locally &&
          SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
        // Reset per-bin source/sink diagnostics after the shock source has been
        // applied and before particle coupling/reflection/cascade for this
        // iteration.  The requested Tecplot diagnostics describe the exchange
        // rates due to cascade, particle interaction, and reflection only; shock
        // injection and spatial advection are not included in these local-rate
        // columns.
        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::ResetSpectralEnergyExchangeRates();
      }

      // Dispatch policy is described by capabilities, not by comparing raw
      // function addresses.  None of the three production movers evolves wave
      // state directly; movers that accumulate streaming feed this manager.
      if (evolve_turbulence_locally &&
          SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode &&
          !SEP::Mover::CurrentCapabilities().evolvesWaveStateDirectly) {

      // Function to increment integrated wave energy due to shock passing
      //reduce S
      PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
      PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming);

      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatum(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,0,"g+",0);
      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatum(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,0,"g-",0);

      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatumMPI(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,"g+",0);
      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatumMPI(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,"g-",0);


      SEP::AlfvenTurbulence_Kolmogorov::AnalyzeMaxSegmentParticles(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,"G_plus_streaming" ,0);
      SEP::AlfvenTurbulence_Kolmogorov::AnalyzeMaxSegmentParticles(SEP::AlfvenTurbulence_Kolmogorov::G_minus_streaming,"G_minus_streaming" ,0);


//      SEP::AlfvenTurbulence_Kolmogorov::SetDatumAll(0.0,SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming);
//      SEP::AlfvenTurbulence_Kolmogorov::TestPrintDatumMPI(SEP::AlfvenTurbulence_Kolmogorov::G_plus_streaming,"g+",0);

      //couple particles and turbulence
//      SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::UpdateAllSegmentsWaveEnergyWithParticleCoupling(
//		     SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
//		    SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::S,
//		   PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);

      if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
        // New spectral coupling: G±(k_j) modifies the same E±(k_j) bin, and
        // the equal-and-opposite particle energy change is redistributed only
        // to particles whose resonant wave number falls into that bin.
        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::WaveParticleCouplingManager(
            SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
            PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);

        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::EnforceRightBoundarySpectrumInitialCondition();
        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::UpdateIntegratedEnergyFromSpectrum(
            SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
      }
      else {
        // Legacy coupling: growth rates are integrated over k before the two
        // branch-integrated energies E+ and E- are updated.
        SEP::AlfvenTurbulence_Kolmogorov::IsotropicSEP::WaveParticleCouplingManager(
            SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
            PIC::ParticleWeightTimeStep::GlobalTimeStep[0]);

        // The wave-particle coupling operator updates E+ and E- in every segment.
        // Re-apply the fixed right-boundary W- condition immediately afterward so
        // the boundary remains equal to the pre-existing initial turbulence state.
        SEP::AlfvenTurbulence_Kolmogorov::EnforceRightBoundaryEminusInitialCondition(
            SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
      }
      }


      // ----------------------------------------------------------------------
      // Advect turbulence energy along field lines.
      //
      // The explicit wave-energy advection step is limited by an Alfven-speed
      // CFL condition.  The particle global time step can be much larger than
      // the Alfven crossing time of the smallest field-line segment, especially
      // near boundaries.  Advancing the turbulence with the full particle time
      // step can therefore pile energy up at the edge cells.  Subcycle the
      // turbulence advection with the global stable time step returned by
      // GetGlobalMaxStableTimeStep().  The function already includes its own
      // safety factor when estimating the CFL limit.
      //
      // The last segment of each field line is treated as a fixed reservoir of
      // pre-existing inward-propagating turbulence W-.  The initial W- value is
      // captured after the turbulence initial condition is generated and is
      // restored after each turbulence operator.  This prevents the right
      // boundary from either draining away by advection or growing by an
      // artificial source.  The TurbulenceLevelEnd argument is retained only for
      // backward-compatible call signatures and is not used to inject W-.
      // ----------------------------------------------------------------------
      if (evolve_turbulence_locally &&
          SEP::Turbulence::ActiveConfiguration().advectionEnabled) {
        const double dt_total_turbulence = PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
        double dt_done_turbulence = 0.0;

        while (dt_done_turbulence < dt_total_turbulence) {
          double dt_cfl_turbulence = SEP::AlfvenTurbulence_Kolmogorov::GetGlobalMaxStableTimeStep();

          // An invalid CFL limit is a rejected physical update.  Advancing the
          // remaining particle timestep would knowingly violate the wave
          // solver's stability condition and can create an unreported energy
          // pile-up, so fail before mutating another subcycle.
          if (dt_cfl_turbulence <= 0.0 || !std::isfinite(dt_cfl_turbulence)) {
            exit(__LINE__,__FILE__,
                 "No finite positive turbulence-advection CFL timestep");
          }

          const double dt_subcycle = std::min(dt_cfl_turbulence, dt_total_turbulence - dt_done_turbulence);

          if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
            // Spectral advection: every E+(k_j) and E-(k_j) bin is transported
            // independently with the same finite-volume Alfvénic flux geometry.
            // The compact integrated datum is refreshed afterward so the rest of
            // the model and the standard output still see E+ and E-.
            SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::AdvectSpectrumAllFieldLines(
                dt_subcycle,
                0.01,  // inner-boundary total W+ turbulence level, distributed over k
                0.0);  // retained argument; spectral right-boundary W-(k) is fixed

            SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::UpdateIntegratedEnergyFromSpectrum(
                SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
            PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
                SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
          }
          else {
            SEP::AlfvenTurbulence_Kolmogorov::AdvectTurbulenceEnergyAllFieldLines(
                DeltaE_plus, DeltaE_minus,
                SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy,
                dt_subcycle,
                0.01,  // inner-boundary W+ source level
                0.0);  // retained argument; right-boundary W- is fixed to its captured initial value

            // The advection operator updates only locally-owned field-line segments.
            // The next subcycle reads neighbor states to compute finite-volume face
            // fluxes.  Refresh the ghost/edge copies after every subcycle; otherwise
            // MPI-domain interfaces use stale E+/E- values during all later subcycles,
            // producing artificial jumps and incorrect interior evolution of W+/W-.
            PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
                SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);

            // MPI synchronization may overwrite the locally prescribed boundary-cell
            // value on the owner rank.  Restore the fixed right-boundary W- density
            // immediately after synchronization so both diagnostics and subsequent
            // flux calculations use the intended pre-existing turbulence value.
            SEP::AlfvenTurbulence_Kolmogorov::EnforceRightBoundaryEminusInitialCondition(
                SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
          }

          dt_done_turbulence += dt_subcycle;
        }
      }

      //model the effect of wave reflection
      if (evolve_turbulence_locally &&
          SEP::AlfvenTurbulence_Kolmogorov::Reflection::active==true) {
        const double C_reflection =
            SEP::Turbulence::ActiveConfiguration().reflectionCoefficient;

        if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
          // Fully spectral reflection: convert E+(k_j) and E-(k_j) into one
          // another at the same wave-number bin j.  This replaces the older
          // compatibility path that first reflected the branch-integrated E± and
          // then projected the result back onto E±(k).  The integrated datum is
          // refreshed afterward only so the compact W+,W-,sigma_c diagnostics and
          // legacy helper routines continue to see the branch sums.
          SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::ReflectSpectrumAllFieldLines(
              PIC::ParticleWeightTimeStep::GlobalTimeStep[0],C_reflection,0.0,false);

          SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::UpdateIntegratedEnergyFromSpectrum(
              SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
          PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
              SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        }
        else {
          SEP::AlfvenTurbulence_Kolmogorov::Reflection::ReflectTurbulenceEnergyAllFieldLines(
              PIC::ParticleWeightTimeStep::GlobalTimeStep[0],C_reflection,0.0,false);

          // Reflection can convert part of W+ into W- in the boundary segment.
          // Keep the last-segment W- fixed to the captured initial value so the
          // boundary condition remains prescribed rather than dynamically evolved.
          SEP::AlfvenTurbulence_Kolmogorov::EnforceRightBoundaryEminusInitialCondition(
              SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        }
      }

      // Configure cascade
      if (evolve_turbulence_locally &&
          SEP::AlfvenTurbulence_Kolmogorov::Cascade::active==true) {
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetCascadeCoefficient(
            SEP::Turbulence::ActiveConfiguration().cascadeCoefficient);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetDefaultPerpendicularCorrelationLength(
            SEP::Turbulence::ActiveConfiguration().perpendicularCorrelationLengthM);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetDefaultEffectiveArea(1.0);           // V_cell = Δs
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::SetElectronHeatingFraction(
            SEP::Turbulence::ActiveConfiguration().electronHeatingFraction);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableCrossHelicityModulation(false);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableTwoSweepIMEX(false);

        // Advance cascade for all field lines (ΔE arrays accumulate changes)
        // Optional: stronger physics
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableCrossHelicityModulation(true);
        SEP::AlfvenTurbulence_Kolmogorov::Cascade::EnableTwoSweepIMEX(true);

        if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
          // Fully spectral cascade: move energy in log-k space within each
          // segment, from bin j to j+1, using the same model parameters as the
          // integrated cascade configuration above: C_nl=0.8, lambda_perp=1e7 m,
          // cross-helicity modulation ON, and a two-sweep update.  Only the
          // energy that reaches k_max is removed as unresolved dissipation.
          // The compact integrated datum is then refreshed as the sum over k.
          SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::CascadeSpectrumAllFieldLines(
              PIC::ParticleWeightTimeStep::GlobalTimeStep[0],
              SEP::Turbulence::ActiveConfiguration().cascadeCoefficient,
              SEP::Turbulence::ActiveConfiguration().perpendicularCorrelationLengthM,
              true,    // enable per-bin cross-helicity modulation
              true,    // two half-sweeps for a Picard-like update
              false);  // enable_logging

          SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::UpdateIntegratedEnergyFromSpectrum(
              SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
          PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
              SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        }
        else {
          SEP::AlfvenTurbulence_Kolmogorov::Cascade::CascadeTurbulenceEnergyAllFieldLines(
              PIC::ParticleWeightTimeStep::GlobalTimeStep[0],/*enable_logging=*/ false);

          // Nonlinear cascade/damping changes both Elsasser wave populations.
          // Restore only the outer-boundary W- value; all interior segments and
          // the outer-boundary W+ outflow remain governed by the cascade update.
          SEP::AlfvenTurbulence_Kolmogorov::EnforceRightBoundaryEminusInitialCondition(
              SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
        }
      }


      //scatter wave energy
      if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive()) {
        PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
            SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy);
        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::EnforceRightBoundarySpectrumInitialCondition();
        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::UpdateIntegratedEnergyFromSpectrum(
            SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
      }
      else {
        PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);

        // Edge synchronization can refresh data stored on field-line boundaries.
        // Enforce the prescribed right-boundary W- value once more before derived
        // wave-energy-density diagnostics are calculated.
        SEP::AlfvenTurbulence_Kolmogorov::EnforceRightBoundaryEminusInitialCondition(
            SEP::AlfvenTurbulence_Kolmogorov::CellIntegratedWaveEnergy);
      }


      //calculate the wave energy density
      CalculateWaveEnergyDensity();

      // --------------------------------------------------------------------
      // Optional 2-D spectral turbulence diagnostics.
      //
      // This output exists only for the wave-number-resolved model because the
      // legacy model does not store E+(k_j) and E-(k_j).  The diagnostic is
      // intentionally controlled by a CLI cadence rather than written every
      // iteration: one file contains all field lines, all segments, and all
      // 128 wave-number bins, so output every time step can become large.
      //
      // The call is placed here, after all wave-energy source/transport terms
      // and after CalculateWaveEnergyDensity(), so the compact W+,W-,sigma_c
      // field-line output and the spectral W±(s,k),sigma_c(s,k) output refer
      // to the same turbulence state.  The shock model was advanced near the
      // beginning of this iteration, so the output title includes the current
      // shock location at this simulation time.
      // --------------------------------------------------------------------
      if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive() &&
          cli_options.spectralOutputInterval > 0 &&
          ((niter+1) % cli_options.spectralOutputInterval == 0)) {
        PIC::FieldLine::Parallel::MPIAllGatherDatumStoredAtEdge(
            SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::SpectralWaveEnergy);

        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::OutputSpectrumTecplot2D(
            niter+1,SEP::Background::SimulationTimeSeconds());
      }

/*
      if ((niter+1)%2==0)  {
         char fname[300];
	 sprintf(fname,"fl-%i-%e.dat",PIC::ThisThread,rsh0/_AU_);

	 //PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(&PIC::FieldLine::DatumAtVertexParticleWeight);
         //PIC::FieldLine::Parallel::MPIAllReduceDatumStoredAtVertex(&PIC::FieldLine::DatumAtVertexParticleCosPitchAngle);

	 PIC::FieldLine::Output(fname,false);

	 //PIC::FieldLine::Parallel::SetDatumStoredAtVertex(0.0,&PIC::FieldLine::DatumAtVertexParticleWeight);
	 //PIC::FieldLine::Parallel::SetDatumStoredAtVertex(0.0,&PIC::FieldLine::DatumAtVertexParticleCosPitchAngle);
      }
*/

    }

    // The adapter already exported these derived values.  Recomputing through
    // the historical output helper is harmless and preserves byte-compatible
    // diagnostic formatting while the old block remains for migration review.
    if (SEP::AlfvenTurbulence_Kolmogorov::ActiveFlag) {
      CalculateWaveEnergyDensity();
      if (SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::IsActive() &&
          cli_options.spectralOutputInterval > 0 &&
          ((niter + 1) % cli_options.spectralOutputInterval == 0)) {
        SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::
            OutputSpectrumTecplot2D(
                niter + 1, SEP::Background::SimulationTimeSeconds());
      }
    }


    //PIC::ParticleSplitting::Split::SplitWithVelocityShift_FL(10,200);
    //
    //PIC::ParticleSplitting::FledLine::WeightedParticleMerging(20,20,20,500,800);
    //PIC::ParticleSplitting::FledLine::WeightedParticleSplitting(20,20,20,500,800);
  }


  char fname[400];

  sprintf(fname,"%s/test_SEP.dat",PIC::OutputDataFileDirectory);
  PIC::RunTimeSystemState::GetMeanParticleMicroscopicParameters(fname);

  {
    const SEP::SW1DAdapter::Diagnostics localDiagnostics=
        SEP::SW1DAdapter::GetDiagnostics();
    const unsigned long long localStateId=
        static_cast<unsigned long long>(localDiagnostics.prepared_state_id);
    unsigned long long minimumStateId=0,maximumStateId=0;
    double minimumEpoch=0.0,maximumEpoch=0.0;
    MPI_Allreduce(&localStateId,&minimumStateId,1,MPI_UNSIGNED_LONG_LONG,
                  MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&localStateId,&maximumStateId,1,MPI_UNSIGNED_LONG_LONG,
                  MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(&localDiagnostics.prepared_epoch_seconds,&minimumEpoch,1,
                  MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    MPI_Allreduce(&localDiagnostics.prepared_epoch_seconds,&maximumEpoch,1,
                  MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if (minimumStateId!=maximumStateId || minimumEpoch!=maximumEpoch) {
      std::ostringstream message;
      message<<"SWCME MPI state consensus failed: state_id_range=["
             <<minimumStateId<<','<<maximumStateId<<"] epoch_range=["
             <<minimumEpoch<<','<<maximumEpoch<<']';
      exit(__LINE__,__FILE__,message.str().c_str());
    }

    // Query counts depend on each rank's particle ownership, so sum rather than
    // require equality. State identity/epoch above must agree exactly. The D03
    // native gate parses this explicit consensus marker and archives the global
    // recovery totals for each decomposition.
    const unsigned long long localCounters[4]={
        static_cast<unsigned long long>(localDiagnostics.successful_queries),
        static_cast<unsigned long long>(localDiagnostics.failed_queries),
        static_cast<unsigned long long>(localDiagnostics.radius_clamps),
        static_cast<unsigned long long>(localDiagnostics.diagnostic_fallbacks)};
    unsigned long long globalCounters[4]={0,0,0,0};
    MPI_Reduce(localCounters,globalCounters,4,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,
               MPI_COMM_WORLD);
    const unsigned long long localCompletedDispatches=
        static_cast<unsigned long long>(SEP::Mover::CompletedDispatchCount());
    unsigned long long globalCompletedDispatches=0;
    MPI_Reduce(&localCompletedDispatches,&globalCompletedDispatches,1,
               MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    if (PIC::ThisThread==0) {
    cout << "SWCME background summary: policy="
         << SEP::SW1DAdapter::FailurePolicyName(
                SEP::SW1DAdapter::GetFailurePolicy())
         << " successful_queries=" << globalCounters[0]
         << " failed_queries=" << globalCounters[1]
         << " radius_clamps=" << globalCounters[2]
         << " diagnostic_fallbacks=" << globalCounters[3]
         << " final_source_state_id=" << maximumStateId
         << " completed_particle_dispatches=" << globalCompletedDispatches
         << " mpi_consensus=pass"
         << endl;
    cout << "End of the run:" << PIC::nTotalSpecies << endl;
    }
  }

  MPI_Finalize();
  return EXIT_SUCCESS;
}
