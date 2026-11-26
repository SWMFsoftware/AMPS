
#include "sep.h"
#include "adapters/swcme1d_adapter.h"
#include "adapters/reduced_shock_background_adapter.h"
#include "transport_common.h"
#include "util/sep_population_production_adapter.h"

int SEP::Offset::MeanFreePath=-1;
int SEP::Offset::TransportSchema=-1;
int SEP::Offset::StableParticleId=-1;
int SEP::Offset::MfpOpticalDepth=-1;
int SEP::Offset::MfpEventIndex=-1;
int SEP::Offset::OperatorRandomStates=-1;
int SEP::Offset::SplitEventCounter=-1;
int SEP::Offset::BirthRecord=-1;
int SEP::Offset::ReducedSourceRecord=-1;

//selector of the shock wave model
SEP::cShockModelType SEP::ShockModelType=SEP::cShockModelType::Analytic1D;

//in the case the model is run as a part of the SWMF, FreezeTimeSimulationMHD  is the sumulation time starting which the control of the
//model run is not returned to the SWMF and the sumulation continues with AMPS only and "freezed" MHD solar wind
double SEP::FreezeSolarWindModelTime=-1.0;

//composition table of the GCR composition
cCompositionGroupTable *SEP::CompositionGroupTable=NULL;
int *SEP::CompositionGroupTableIndex=NULL;
int SEP::nCompositionGroups=0;

cInternalSphericalData* SEP::InnerBoundary=NULL;

//the type of the equations that is soleved
int SEP::ModelEquation=SEP::ModelEquationFTE;

//IMF used in the calcualtions: either Parker spiral or background magnetif field
int SEP::ModeIMF=SEP::ModeIMF_background;


//types of the differentiation of the pitch angle diffusion coeffcient
int SEP::Diffusion::PitchAngleDifferentialMode=SEP::Diffusion::PitchAngleDifferentialModeAnalytical;

//parameters of the scattering model
int SEP::Scattering::Tenishev2005AIAA::status=SEP::Scattering::Tenishev2005AIAA::_disabled;
double SEP::Scattering::Tenishev2005AIAA::alpha=1.0/3.0;
double SEP::Scattering::Tenishev2005AIAA::beta=2.0/3.0;
double SEP::Scattering::Tenishev2005AIAA::lambda0=0.4*_AU_;

//the limit to switch from solving FTE to the Parker Equation when the D_{\mu\mu} is to high
double SEP::TimeStepRatioSwitch_FTE2PE=-1.0;

//min/max particle number limit during a run
int SEP::MinParticleLimit=10,SEP::MaxParticleLimit=20;

//the model for solar wind density
int SEP::ParticleSource::ShockWaveSphere::SolarWindDensityMode=SEP::ParticleSource::ShockWaveSphere::SolarWindDensityMode_analytic;


SEP::Transport::Status SEP::QueryShockDiagnostic(double epochS,ShockDiagnostic* result) {
  using namespace Transport;
  if(!result||!std::isfinite(epochS))return Status::Error(StatusCode::InvalidArgument,
      "shock output metadata requires a finite physical epoch and an output record");
  ShockDiagnostic candidate;
  if(ReducedShock::Enabled()) {
    // The reduced provider takes precedence over the retained legacy selector.
    // Its disabled Tenishev radius can contain an initialization sentinel;
    // reporting that value beside real core transport would misidentify the
    // physical front. Query the CORE leading ray at this output's exact epoch,
    // not its separate prepared-cache cadence. No state is published here.
    const auto view=ReducedShock::CapturePhysicalStep(epochS,epochS,ReducedShock::EventIdentity());
    if(!view.ok())return Status::Error(StatusCode::UnsupportedConfiguration,view.status.message);
    const auto sphere=view.value.SphereAt(epochS);
    if(!sphere.ok())return Status::Error(StatusCode::UnsupportedConfiguration,sphere.status.message);
    const auto apex=view.value.LeadingPointOnRay(sphere.value.direction,epochS);
    if(!apex.ok())return Status::Error(StatusCode::UnsupportedConfiguration,apex.status.message);
    candidate.apexRadiusM=0;
    for(double x:apex.value)candidate.apexRadiusM=std::hypot(candidate.apexRadiusM,x);
    candidate.modelName=ReducedShock::HeliocentricSphereFront()?
        "reduced-sphere/Parker":"reduced-SSE/Parker";
    if(!std::isfinite(candidate.apexRadiusM)||candidate.apexRadiusM<=0)
      return Status::Error(StatusCode::InvalidParticleState,"core output apex distance is nonfinite/nonpositive");
  }else if(ShockModelType==cShockModelType::Analytic1D) {
    candidate.apexRadiusM=ParticleSource::ShockWave::Tenishev2005::rShock;candidate.modelName="analytic-1D";
  }else if(ShockModelType==cShockModelType::SwCme1d) {
    candidate.apexRadiusM=SW1DAdapter::ShockRadiusM();candidate.modelName="SW-CME-1D";
  }
  *result=candidate;return Status::Ok();
}

//title that will be printed in Tecplot output file (simulation time and shock location)
void SEP::TecplotFileTitle(char* title) {
  // This function is installed below as
  //   PIC::FieldLine::UserDefinedTecplotFileTitle
  // and is called directly from PIC::FieldLine::Output() at the moment when a
  // file such as amps.FieldLines.out=<counter>.dat is being created.  Placing
  // the SEP-specific title extension here is therefore preferable to scanning
  // the output directory after every iteration: the title is generated once,
  // for the file currently being written, with no post-processing and no race
  // with the generic AMPS field-line writer.

  ShockDiagnostic shock;
  const auto metadata=QueryShockDiagnostic(SEP::Background::SimulationTimeSeconds(),&shock);
  if(!metadata.ok())exit(__LINE__,__FILE__,metadata.message.c_str());
  const double rShock=shock.apexRadiusM;const char* ShockModelName=shock.modelName;

  // Keep the original information, time=<...>, because existing Tecplot macros
  // and post-processing scripts may rely on it.  The additional R_sh fields are
  // redundant on purpose: meters are the internal SI units, AU is convenient for
  // heliospheric plots, and solar radii are convenient near the Sun.
  if (rShock>0.0) {
    sprintf(title,
            "time=%e; shock model=%s; R_sh=%e m = %e AU = %e R_s",
            SEP::Background::SimulationTimeSeconds(),ShockModelName,
            rShock,rShock/_AU_,rShock/_SUN__RADIUS_);
  }
  else {
    // Fallback for initialization phases when a valid shock radius has not yet
    // been published.  This preserves the previous title format rather than
    // putting a nonphysical negative radius into the Tecplot file.
    sprintf(title,"time=%e",SEP::Background::SimulationTimeSeconds());
  }
}


void SEP::Init() {
  //title that will be printed inn Tecplot output file (simuation time)
  PIC::FieldLine::UserDefinedTecplotFileTitle=TecplotFileTitle;

  //composition of the GCRs
  nCompositionGroups=1;
  CompositionGroupTable=new cCompositionGroupTable[nCompositionGroups];
  CompositionGroupTableIndex=new int[PIC::nTotalSpecies];

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) CompositionGroupTableIndex[spec]=0; //all simulated model species are hydrogen
  CompositionGroupTable[0].FistGroupSpeciesNumber=0;
  CompositionGroupTable[0].nModelSpeciesGroup=PIC::nTotalSpecies;

  CompositionGroupTable[0].minVelocity=Relativistic::E2Speed(SEP::BoundingBoxInjection::minEnergy,PIC::MolecularData::GetMass(0));
  CompositionGroupTable[0].maxVelocity=Relativistic::E2Speed(SEP::BoundingBoxInjection::maxEnergy,PIC::MolecularData::GetMass(0));

  CompositionGroupTable[0].GroupVelocityStep=(CompositionGroupTable[0].maxVelocity-CompositionGroupTable[0].minVelocity)/CompositionGroupTable[0].nModelSpeciesGroup;

  if (PIC::ThisThread==0) {
    cout << "$PREFIX: Composition Group Velocity and Energy Characteristics:\nspec\tmin Velocity [m/s]\tmax Velovity[m/s]\t min Energy[eV]\tmax Energy[eV]" << endl;

    for (int s=0;s<PIC::nTotalSpecies;s++) {
      double minV,maxV,minE,maxE,mass;

      mass=PIC::MolecularData::GetMass(s);

      minV=::SEP::CompositionGroupTable[0].GetMinVelocity(s);
      maxV=::SEP::CompositionGroupTable[0].GetMaxVelocity(s);

      //convert velocity into energy and distribute energy of a new particles
      minE=Relativistic::Speed2E(minV,mass);
      maxE=Relativistic::Speed2E(maxV,mass);

      cout << s << "\t" << minV << "\t" << maxV << "\t" << minE*J2eV << "\t" <<  maxE*J2eV << endl;
    }
  }

  //init source models of SEP and GCR
  if (_PIC_EARTH_GCR__MODE_==_PIC_MODE_ON_) BoundingBoxInjection::GCR::Init();
}


void SEP::RequestParticleData() {
  long int offset;

  // These are scalar field-line transport state, not a second Cartesian
  // position/momentum representation.  AMPS checkpoints include requested
  // particle storage, so stochastic event continuity survives restart.
  PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double),alignof(double));
  Offset::MeanFreePath=offset;
  PIC::ParticleBuffer::RequestDataStorage(
      offset,sizeof(std::uint64_t),alignof(std::uint64_t));
  Offset::TransportSchema=offset;
  PIC::ParticleBuffer::RequestDataStorage(
      offset,sizeof(std::uint64_t),alignof(std::uint64_t));
  Offset::StableParticleId=offset;
  PIC::ParticleBuffer::RequestDataStorage(offset,sizeof(double),alignof(double));
  Offset::MfpOpticalDepth=offset;
  PIC::ParticleBuffer::RequestDataStorage(
      offset,sizeof(std::uint64_t),alignof(std::uint64_t));
  Offset::MfpEventIndex=offset;
  PIC::ParticleBuffer::RequestDataStorage(
      offset,3*sizeof(SEP::Transport::OperatorRandomState),
      alignof(SEP::Transport::OperatorRandomState));
  Offset::OperatorRandomStates=offset;
  PIC::ParticleBuffer::RequestDataStorage(
      offset,sizeof(std::uint64_t),alignof(std::uint64_t));
  Offset::SplitEventCounter=offset;
  PIC::ParticleBuffer::RequestDataStorage(
      offset,sizeof(SEP::Transport::Birth::Record),
      alignof(SEP::Transport::Birth::Record));
  Offset::BirthRecord=offset;
  // Two eight-byte scalar slots: schema and physical birth kinetic energy [J].
  // The latter remains unchanged by cooling/scattering and exact splitting.
  PIC::ParticleBuffer::RequestDataStorage(offset,16,alignof(double));
  Offset::ReducedSourceRecord=offset;
  SEP::Transport::PICAdapter::RegisterPopulationControlHooks();
  SEP::PopulationControl::RegisterConservativeMerger();
}
