
#include "pic.h"
#include "sep.h"
#include "adapters/swcme1d_adapter.h"
#include "transport_common.h"
#include "amps2swmf.h"

// field_line.cpp is copied to build/main and compiled by AMPS's generic
// object rule.  That rule always exposes the AMPS root, but it does not inherit
// include flags appended by the earlier srcSEP submake.  Resolve these shared
// model headers through the same canonical-path adapter used by public srcSEP
// headers so the source works in both source/srcSEP and copied build/main.
#include "util/sep_common_header_path.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_species_source.h)
#include SRCSEP_SEP_COMMON_HEADER(sep_injection_spectrum.h)

#include <cmath>
#include <cstring>
#include <string>


int SEP::FieldLine::InjectionParameters::nParticlesPerIteration=300;
double SEP::FieldLine::InjectionParameters::PowerIndex=4.0;
double SEP::FieldLine::InjectionParameters::emin=0.1,SEP::FieldLine::InjectionParameters::emax=500;
double SEP::FieldLine::InjectionParameters::InjectionEfficiency=3.4E-4; //Sokolov-2004-AJ

double SEP::FieldLine::InjectionParameters::ConstEnergyInjectionValue=0.0;
double SEP::FieldLine::InjectionParameters::ConstSpeedInjectionValue=0.0;
double SEP::FieldLine::InjectionParameters::ConstMuInjectionValue=0.5;

#if _SEP_FIELD_LINE_INJECTION_ == _SEP_FIELD_LINE_INJECTION__SHOCK_
int SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectShockLocations;
#else
int SEP::FieldLine::InjectionParameters::InjectLocation=SEP::FieldLine::InjectionParameters::_InjectBegginingFL;
#endif



int SEP::FieldLine::InjectionParameters::InjectionMomentumModel=SEP::FieldLine::InjectionParameters::_tenishev2005aiaa;
int SEP::FieldLine::InjectionParameters::UseAnalyticShockModel=SEP::FieldLine::InjectionParameters::AnalyticShockModel_Tenishev2005;

namespace {

std::uint64_t InjectionEventKey() {
  const double epoch=SEP::Background::SimulationTimeSeconds();
  std::uint64_t bits=0;
  static_assert(sizeof(bits)==sizeof(epoch),"double event key must be 64 bits");
  std::memcpy(&bits,&epoch,sizeof(bits));
  return bits;
}

SEP::Transport::KeyedRandomStream SourceRandom(
    int fieldLine,int species,std::uint64_t macro,
    SEP::Injection::RandomPurpose purpose) {
  SEP::Injection::RandomKey key;
  key.campaign=0; // WP30 replaces this with the frozen campaign seed.
  key.event=InjectionEventKey();
  key.fieldLine=static_cast<std::uint64_t>(fieldLine);
  key.species=static_cast<std::uint64_t>(species);
  key.macroparticle=macro;
  key.purpose=purpose;
  return SEP::Injection::MakeRandomStream(key);
}

void SampleIsotropicMomentum(double magnitude,
                             SEP::Transport::RandomStream& random,
                             double p[3]) {
  const double mu=2.0*random.UniformOpen01()-1.0;
  const double phi=2.0*Pi*random.UniformOpen01();
  const double perpendicular=std::sqrt(std::max(0.0,1.0-mu*mu));
  p[0]=magnitude*perpendicular*std::cos(phi);
  p[1]=magnitude*perpendicular*std::sin(phi);
  p[2]=magnitude*mu;
}

void SamplePitchAngleMomentum(double magnitude,double mu,const double l[3],
                              SEP::Transport::RandomStream& random,double p[3]) {
  // Construct a deterministic perpendicular frame, then randomize only the
  // physical gyrophase with its purpose-tagged stream.  This removes the hidden
  // global RNG consumed by Vector3D::GetRandomNormFrame.
  double reference[3]={0.0,0.0,1.0};
  if (std::fabs(l[2])>0.9) { reference[0]=1.0;reference[2]=0.0; }
  double e0[3]={l[1]*reference[2]-l[2]*reference[1],
                l[2]*reference[0]-l[0]*reference[2],
                l[0]*reference[1]-l[1]*reference[0]};
  const double norm=Vector3D::Length(e0);
  for (int d=0;d<3;++d) e0[d]/=norm;
  const double e1[3]={l[1]*e0[2]-l[2]*e0[1],
                      l[2]*e0[0]-l[0]*e0[2],
                      l[0]*e0[1]-l[1]*e0[0]};
  const double phi=2.0*Pi*random.UniformOpen01();
  const double perpendicular=std::sqrt(std::max(0.0,1.0-mu*mu));
  for (int d=0;d<3;++d)
    p[d]=magnitude*(mu*l[d]+perpendicular*(std::cos(phi)*e0[d]+
                                           std::sin(phi)*e1[d]));
}

SEP::Transport::SpeciesSource::Configuration ResolveSpeciesSource(int spec) {
  namespace SS = SEP::Transport::SpeciesSource;
  SS::Configuration source;
  if (SS::HasActiveConfiguration()) {
    const SEP::Transport::Status status =
        SS::FindActiveConfiguration(spec, &source);
    if (!status.ok()) exit(__LINE__, __FILE__, status.message.c_str());
    return source;
  }

  // Preserve legacy input behavior when no explicit WP18 source table was
  // installed, while still returning a complete typed species record to the
  // injection calculations below.  The fallback is intentionally visible and
  // cannot be used after an explicit table is installed: a missing configured
  // species then fails before it can inherit proton-like parameters.
  source.species.modelSpecies = spec;
  source.species.name = "PIC-species-" + std::to_string(spec);
  source.species.signedChargeC =
      PIC::MolecularData::GetElectricCharge(spec);
  source.species.restMassKg = PIC::MolecularData::GetMass(spec);
  const double protonMassKg = 1.67262192369e-27;
  const double massRatio = source.species.restMassKg / protonMassKg;
  source.species.nucleonCount =
      massRatio >= 0.5 ? std::floor(massRatio + 0.5) : 0.0;
  source.abundanceFraction = 1.0;
  source.injectionEfficiency =
      SEP::FieldLine::InjectionParameters::InjectionEfficiency;
  source.spectralIndex = SEP::FieldLine::InjectionParameters::PowerIndex;
  source.energyConvention = SS::EnergyConvention::TotalKinetic;
  return source;
}

double ResolveTotalEnergyJ(
    const SEP::Transport::SpeciesSource::Configuration& source,
    double configuredEnergyMeV) {
  const SEP::Transport::ScalarResult total =
      SEP::Transport::SpeciesSource::TotalKineticEnergyJ(
          source, SEP::Units::EnergyFromMeV(configuredEnergyMeV).Value());
  if (!total.status.ok())
    exit(__LINE__, __FILE__, total.status.message.c_str());
  return total.value;
}

double ExplicitSourceScale(
    const SEP::Transport::SpeciesSource::Configuration& source) {
  if (!SEP::Transport::SpeciesSource::HasActiveConfiguration()) return 1.0;
  const double legacyEfficiency =
      SEP::FieldLine::InjectionParameters::InjectionEfficiency;
  if (!(legacyEfficiency > 0.0))
    exit(__LINE__, __FILE__,
         "explicit species source cannot scale a non-positive legacy efficiency");
  // Analytic legacy rates already include the historical global efficiency.
  // Replace that factor with the species-local efficiency, then partition the
  // event by normalized abundance.  SWMF count construction below uses the
  // local efficiency directly and therefore does not call this helper.
  return source.abundanceFraction *
      source.injectionEfficiency / legacyEfficiency;
}

}  // namespace



long int SEP::FieldLine::InjectParticleFieldLineBeginning(int spec,int iFieldLine) {
  namespace FL = PIC::FieldLine;

  long int newParticle;
  PIC::ParticleBuffer::byte *newParticleData;
  int nInjectedParticles=0;
  int npart;
  double l[3],pAbs,p[3],ParticleWeightCorrectionFactor=1.0;

  npart=100;
  const SEP::Transport::SpeciesSource::Configuration speciesSource =
      ResolveSpeciesSource(spec);
  ParticleWeightCorrectionFactor=ExplicitSourceScale(speciesSource);
  pAbs=Relativistic::Energy2Momentum(
      ResolveTotalEnergyJ(speciesSource,100.0),
      PIC::MolecularData::GetMass(spec));

  FL::FieldLinesAll[iFieldLine].GetSegment(0)->GetDir(l);

  for (int i=0;i<npart;i++) {
    //generate a particle
    SEP::Transport::KeyedRandomStream random=SourceRandom(
        iFieldLine,spec,static_cast<std::uint64_t>(i),
        SEP::Injection::RandomPurpose::PitchAngle);
    SampleIsotropicMomentum(pAbs,random,p);

    if (Vector3D::DotProduct(p,l)<0.0) for (int idim=0;idim<3;idim++) p[idim]=-p[idim];

    if ((newParticle=PIC::FieldLine::InjectParticle_default(spec,p,ParticleWeightCorrectionFactor,iFieldLine,0))!=-1) {
      SEP::Transport::PICAdapter::InitializeParticleTransportState(
          newParticle, UINT64_C(1),
          (static_cast<std::uint64_t>(iFieldLine) << 32) |
              static_cast<std::uint64_t>(i),
          p);
      nInjectedParticles++;

    }
  }

  return nInjectedParticles;
}

long int InjectSolarWindIons(int spec,int iFieldLine) {
  namespace FL = PIC::FieldLine;

  double InjectionArea,n_sw,t_sw,v_sw[3];
  auto Segment=FL::FieldLinesAll[iFieldLine].GetFirstSegment();
  FL::cFieldLineVertex* FirstVertex=Segment->GetBegin();

  //determine the parameters of of the solar wind at the beginning of the field line
  FirstVertex->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw);
  FirstVertex->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw);
  FirstVertex->GetPlasmaVelocity(v_sw);

  InjectionArea=SEP::FieldLine::FluxTubeGeometry::AreaAtVertexM2(
      FirstVertex,iFieldLine);

  //inject model partiles
  return PIC::FieldLine::InjectMaxwellianLineBeginning(spec,n_sw,t_sw,v_sw,InjectionArea,iFieldLine,200);
}


long int SEP::FieldLine::InjectParticlesSingleFieldLine(int spec,int iFieldLine) {
  namespace FL = PIC::FieldLine;

  int iShockFieldLine,npart;
  double xInjection[3]={0.0,0.0,0.0},S=0.0,anpart,p[3],ParticleWeightCorrectionFactor;
  int nInjectedParticles=0;
  bool injection_coordinate_was_set=false;

  // Resolve the species contract once for the complete injection event.  Its
  // charge/mass metadata, abundance, efficiency, spectral index, and energy
  // convention therefore cannot change between count and momentum sampling.
  const SEP::Transport::SpeciesSource::Configuration speciesSource =
      ResolveSpeciesSource(spec);


  //determine the filed line to inject particles
  iShockFieldLine=0;

  if (InjectionParameters::InjectLocation==InjectionParameters::_InjectInputFileAMPS) {
    #ifdef _SEP_SHOCK_LOCATION_COUPLER_TABLE_
    #if _SEP_SHOCK_LOCATION_COUPLER_TABLE_ == _PIC_MODE_ON_
    if ((iShockFieldLine=AMPS2SWMF::ShockData[iFieldLine].iSegmentShock)==-1) return 0;
    #endif
    #endif
  }
  else {
    switch (InjectionParameters::InjectLocation) {
    case InjectionParameters::_InjectShockLocations:
      #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
      if (AMPS2SWMF::ShockData==NULL) {
        exit(__LINE__,__FILE__,"Error: the shock location table is not allocated");
      }
      else {
        if ((iShockFieldLine=AMPS2SWMF::ShockData[iFieldLine].iSegmentShock)==-1) return 0;
      }
      #else
      switch (InjectionParameters::UseAnalyticShockModel) {
      case InjectionParameters::AnalyticShockModel_Tenishev2005:
	iShockFieldLine=SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionLocation(iFieldLine,S,xInjection);
	injection_coordinate_was_set=(iShockFieldLine>=0);
	break;
      case InjectionParameters::AnalyticShockModel_none:
        iShockFieldLine=0;
	break;
      default:
	exit(__LINE__,__FILE__,"Error: the option is unknown");
      }
      #endif

      break;
    case  InjectionParameters::_InjectBegginingFL:
      iShockFieldLine=0;
      break;
    }
  }

  // Resolve a physical injection coordinate for every configuration branch.
  // The analytic shock already supplies its precise local intersection.  SWMF,
  // file-driven, and no-shock branches identify only a segment, so they use its
  // midpoint; beginning injection uses its beginning.  This prevents an
  // uninitialized field-line coordinate from being stored in a new particle.
  FL::cFieldLineSegment* Segment=FL::FieldLinesAll[iFieldLine].GetSegment(iShockFieldLine);

  if (Segment==NULL) return 0;
  if (Segment->Thread!=PIC::ThisThread) return 0;

  if (!injection_coordinate_was_set) {
    const double local_fraction =
        (InjectionParameters::InjectLocation==
         InjectionParameters::_InjectBegginingFL) ? 0.0 : 0.5;
    S=iShockFieldLine+local_fraction;
    Segment->GetCartesian(xInjection,local_fraction);
  }

  //determine the volume swept by the shock wave during the time step
  double xBegin[3],xEnd[3],xMiddle[3];

  Segment->GetBegin()->GetX(xBegin);
  Segment->GetEnd()->GetX(xEnd);

  for (int idim=0;idim<3;idim++) xMiddle[idim]=0.5*(xBegin[idim]+xEnd[idim]);

  //velocity of the shock wave
  double swept_volume_m3;
  cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node=PIC::Mesh::Search::FindBlock(xMiddle);
  if (node==NULL || node->block==NULL) return 0;

  const double local_fraction=S-iShockFieldLine;

  #if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  if (AMPS2SWMF::ShockData[iFieldLine].ShockSpeed>AMPS2SWMF::MinShockSpeed) {
    swept_volume_m3=SEP::FieldLine::FluxTubeGeometry::SweptVolumeM3(
        Segment,iFieldLine,local_fraction,
        AMPS2SWMF::ShockData[iFieldLine].ShockSpeed,
        node->block->GetLocalTimeStep(spec));
  }
  else {
    if (AMPS2SWMF::MinShockSpeed==0.0) exit(__LINE__,__FILE__,"Error: AMPS2SWMF::MinShockSpeed is not set");

    swept_volume_m3=SEP::FieldLine::FluxTubeGeometry::SweptVolumeM3(
        Segment,iFieldLine,local_fraction,AMPS2SWMF::MinShockSpeed,
        node->block->GetLocalTimeStep(spec));
  }
  #else
    switch (InjectionParameters::UseAnalyticShockModel) {
    case InjectionParameters::AnalyticShockModel_Tenishev2005:
      swept_volume_m3=SEP::ParticleSource::ShockWave::Tenishev2005::GetShockSpeed();
      break;
    case InjectionParameters::AnalyticShockModel_none:
      swept_volume_m3=1.0;
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }

    double LocalTimeStep=-1;

    switch( _SIMULATION_TIME_STEP_MODE_) {
    case _SPECIES_DEPENDENT_LOCAL_TIME_STEP_:
      LocalTimeStep=node->block->GetLocalTimeStep(spec);
      break;
    case  _SPECIES_DEPENDENT_GLOBAL_TIME_STEP_:
      LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[spec];
      break;
    case  _SINGLE_GLOBAL_TIME_STEP_:
      LocalTimeStep=PIC::ParticleWeightTimeStep::GlobalTimeStep[0];
      break;
    default:
      exit(__LINE__,__FILE__,"not implemented");
    }


    swept_volume_m3=SEP::FieldLine::FluxTubeGeometry::SweptVolumeM3(
        Segment,iFieldLine,local_fraction,swept_volume_m3,LocalTimeStep);
  #endif


  //determine the number of particles to inject
  double t_sw_begin,t_sw_end; //=Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaTemperature);
  double n_sw_begin,n_sw_end; //=Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaDensity);
  Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_begin);
  Segment->GetBegin()->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_begin);

  Segment->GetEnd()->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t_sw_end);
  Segment->GetEnd()->GetDatum(FL::DatumAtVertexPlasmaDensity,&n_sw_end);


#if _PIC_COUPLER_MODE_ == _PIC_COUPLER_MODE__SWMF_
  n_sw_end=AMPS2SWMF::ShockData[iFieldLine].DownStreamDensity;

  const double injected_physical_particles=
      SEP::FieldLine::FluxTubeGeometryCore::InjectedPhysicalParticleCount(
          SEP::Units::NumberDensityPerM3(n_sw_end),
          SEP::Units::VolumeM3(swept_volume_m3),
          speciesSource.injectionEfficiency);
  // abundanceFraction partitions the total source across explicitly
  // configured species.  Legacy mode returns one here and is bit-compatible
  // with the former single-parameter source.
  anpart=speciesSource.abundanceFraction*injected_physical_particles;
  anpart/=node->block->GetLocalParticleWeight(spec);
#else
  n_sw_end=1.0;

  switch (InjectionParameters::UseAnalyticShockModel) {
  case InjectionParameters::AnalyticShockModel_Tenishev2005:
    anpart=ExplicitSourceScale(speciesSource)*swept_volume_m3*
        SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate()/
        node->block->GetLocalParticleWeight(spec);
    cout << "Shock locaiton=" << Vector3D::Length(xInjection)/_AU_ << "[AU], Source Rate=" << SEP::ParticleSource::ShockWave::Tenishev2005::GetInjectionRate() << endl << flush;
    break;
  case InjectionParameters::AnalyticShockModel_none:
    anpart=InjectionParameters::nParticlesPerIteration;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }
#endif

  double GlobalWeightCorrectionFactor=1.0;

  if (anpart==0.0) return 0.0;
  else if (anpart<InjectionParameters::nParticlesPerIteration) {
    GlobalWeightCorrectionFactor=anpart/InjectionParameters::nParticlesPerIteration;
    anpart=InjectionParameters::nParticlesPerIteration;
  }
  else if (anpart>10.0*InjectionParameters::nParticlesPerIteration) {
    GlobalWeightCorrectionFactor=anpart/(10.0*InjectionParameters::nParticlesPerIteration);
    anpart=10*InjectionParameters::nParticlesPerIteration;
  }

  //in case particle are injected at the beginning of the field line, the actual plasma density is not used -> set the particle weight == 1
  if (InjectionParameters::InjectLocation==InjectionParameters::_InjectBegginingFL) {
    GlobalWeightCorrectionFactor=1.0;
  }

  npart=(int)anpart;
  SEP::Transport::KeyedRandomStream countRandom=SourceRandom(
      iFieldLine,spec,0,SEP::Injection::RandomPurpose::EventCount);
  if (anpart-npart>countRandom.UniformOpen01()) npart++;

  auto GetMomentum_Tenishev2005AIAA = [&] (double *pAbsTable,double *WeightCorrectionTable,int nParticles) -> bool {
    // Convert MeV or MeV/nucleon to total joules before the relativistic
    // momentum transform.  The convention is part of the per-species source
    // record and is therefore unambiguous for alpha particles and heavy ions.
    double emin=ResolveTotalEnergyJ(speciesSource,InjectionParameters::emin);
    double emax=ResolveTotalEnergyJ(speciesSource,InjectionParameters::emax);

    double s;

    switch (_PIC_COUPLER_MODE_) {
    case _PIC_COUPLER_MODE__SWMF_:
      s=AMPS2SWMF::ShockData[iFieldLine].CompressionRatio;
      break;
    default:
      s=SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio();    //InjectionParameters::PowerIndex;

      switch (SEP::ShockModelType) {
      case SEP::cShockModelType::Analytic1D:
        s=SEP::ParticleSource::ShockWave::Tenishev2005::GetCompressionRatio();
        break;
      case SEP::cShockModelType::SwCme1d:
        s=SEP::SW1DAdapter::CompressionRatio();
        break;
      default:
        exit(__LINE__,__FILE__,"Error: the case is not known");
      }


    }

    if (s>SEP::ParticleSource::ShockWave::MaxLimitCompressionRatio) s=SEP::ParticleSource::ShockWave::MaxLimitCompressionRatio;

    if (s==1.0) return false;

    double q=3.0*s/(s-1.0);
    double pAbs,pmin,pmax,speed,pvect[3];
    double mass=PIC::MolecularData::GetMass(spec);

    if (q<1.0) q=1.0;

    pmin=Relativistic::Energy2Momentum(emin,mass);
    pmax=Relativistic::Energy2Momentum(emax,mass);

    speed=Relativistic::E2Speed(emin,PIC::MolecularData::GetMass(spec));
    pmin=Relativistic::Speed2Momentum(speed,mass);

    speed=Relativistic::E2Speed(emax,PIC::MolecularData::GetMass(spec));
    pmax=Relativistic::Speed2Momentum(speed,mass);

    // Diffusive-shock acceleration supplies phase-space f(p) proportional to
    // p^-q.  The number spectrum is dN/dp proportional to p^(2-q), hence the
    // explicit density index q-2 below.  Sampling its normalized inverse CDF
    // eliminates the ambiguous historical log-p proposal weight.
    SEP::Injection::Spectrum spectrum;
    spectrum.measure=SEP::Injection::Measure::Momentum;
    spectrum.minimum=pmin;
    spectrum.maximum=pmax;
    spectrum.powerIndex=q-2.0;
    for (int i=0;i<nParticles;i++) {
      SEP::Transport::KeyedRandomStream random=SourceRandom(
          iFieldLine,spec,static_cast<std::uint64_t>(i),
          SEP::Injection::RandomPurpose::Spectrum);
      const SEP::Transport::ScalarResult sample=SEP::Injection::InverseCdf(
          spectrum,random.UniformOpen01());
      if (!sample.status.ok()) exit(__LINE__,__FILE__,sample.status.message.c_str());
      pAbsTable[i]=sample.value;
      WeightCorrectionTable[i]=1.0;
    }

    return true;
  };

  auto GetMomentum_Sokolov2004AJ = [&] (double *pAbsTable,double *WeightCorrectionTable,int nParticles) {
    // Configuration energies are MeV; relativistic momentum helpers require
    // joules.  The typed conversion makes this unit boundary explicit.
    const double energy_min_J=ResolveTotalEnergyJ(
        speciesSource,SEP::FieldLine::InjectionParameters::emin);
    const double energy_max_J=ResolveTotalEnergyJ(
        speciesSource,SEP::FieldLine::InjectionParameters::emax);
    double p_injection_min=Relativistic::Energy2Momentum(
        energy_min_J,PIC::MolecularData::GetMass(spec));
    double p_injection_max=Relativistic::Energy2Momentum(
        energy_max_J,PIC::MolecularData::GetMass(spec));

    if (speciesSource.spectralIndex<=1.0) exit(__LINE__,__FILE__,"species source spectral index is out of range: must be greater than one");
    SEP::Injection::Spectrum spectrum;
    spectrum.measure=SEP::Injection::Measure::Momentum;
    spectrum.minimum=p_injection_min;
    spectrum.maximum=p_injection_max;
    spectrum.powerIndex=speciesSource.spectralIndex;
    for (int i=0;i<nParticles;i++) {
      SEP::Transport::KeyedRandomStream random=SourceRandom(
          iFieldLine,spec,static_cast<std::uint64_t>(i),
          SEP::Injection::RandomPurpose::Spectrum);
      const SEP::Transport::ScalarResult sample=SEP::Injection::InverseCdf(
          spectrum,random.UniformOpen01());
      if (!sample.status.ok()) exit(__LINE__,__FILE__,sample.status.message.c_str());
      pAbsTable[i]=sample.value;
      WeightCorrectionTable[i]=1.0;
    }
  };

  double *pAbsTable=new double [npart];
  double *WeightCorrectionTable=new double [npart];
  double p_const;
  bool shock_injects_particles=true;

  switch (InjectionParameters::InjectionMomentumModel) {
  case InjectionParameters::_tenishev2005aiaa:
    shock_injects_particles=GetMomentum_Tenishev2005AIAA(pAbsTable,WeightCorrectionTable,npart);
    break;
  case InjectionParameters::_sokolov2004aj:
    GetMomentum_Sokolov2004AJ(pAbsTable,WeightCorrectionTable,npart);
    break;
  case InjectionParameters::_const_speed:
    p_const=Relativistic::Speed2Momentum(SEP::FieldLine::InjectionParameters::ConstSpeedInjectionValue,PIC::MolecularData::GetMass(spec));

    for (int i=0;i<npart;i++) pAbsTable[i]=p_const,WeightCorrectionTable[i]=1.0;
    break;
  case InjectionParameters::_const_energy:
     // ConstEnergyInjectionValue is part of the legacy input contract and is
     // expressed in MeV, not joules.
     p_const=Relativistic::Energy2Momentum(
         ResolveTotalEnergyJ(
             speciesSource,
             SEP::FieldLine::InjectionParameters::ConstEnergyInjectionValue),
         PIC::MolecularData::GetMass(spec));

    for (int i=0;i<npart;i++) pAbsTable[i]=p_const,WeightCorrectionTable[i]=1.0;
    break;
  default:
    exit(__LINE__,__FILE__,"Error: the option is unknown");
  }

  if (shock_injects_particles==true) for (int i=0;i<npart;i++) {
    if ((InjectionParameters::InjectionMomentumModel==InjectionParameters::_const_speed)||(InjectionParameters::InjectionMomentumModel==InjectionParameters::_const_energy)) {
      double l[3];
      double mu;

      Segment->GetDir(l);
      mu=SEP::FieldLine::InjectionParameters::ConstMuInjectionValue;
      SEP::Transport::KeyedRandomStream gyrophaseRandom=SourceRandom(
          iFieldLine,spec,static_cast<std::uint64_t>(i),
          SEP::Injection::RandomPurpose::Gyrophase);
      SamplePitchAngleMomentum(pAbsTable[i],mu,l,gyrophaseRandom,p);
    }
    else {
      SEP::Transport::KeyedRandomStream directionRandom=SourceRandom(
          iFieldLine,spec,static_cast<std::uint64_t>(i),
          SEP::Injection::RandomPurpose::PitchAngle);
      SampleIsotropicMomentum(pAbsTable[i],directionRandom,p);
    }

    long int newParticle;

    if ((newParticle=PIC::FieldLine::InjectParticle_default(spec,p,GlobalWeightCorrectionFactor*WeightCorrectionTable[i],iFieldLine,iShockFieldLine))!=-1) {
      SEP::Transport::PICAdapter::InitializeParticleTransportState(
          newParticle, UINT64_C(2),
          (static_cast<std::uint64_t>(iFieldLine) << 32) |
              static_cast<std::uint64_t>(i),
          p);
      nInjectedParticles++;

      //Set the local coordinte to the shock location
      PIC::ParticleBuffer::SetFieldLineCoord(S,newParticle);

      //set the initiali distance of the particle from the assigned magnetic field line
    }
  }

  delete [] pAbsTable;
  delete [] WeightCorrectionTable;

  return nInjectedParticles;
}

long int SEP::FieldLine::InjectParticles() {
  long int res=0;

  for (int spec=0;spec<PIC::nTotalSpecies;spec++) for (int iFieldLine=0;iFieldLine<PIC::FieldLine::nFieldLine;iFieldLine++) {
    switch (SEP::FieldLine::InjectionParameters::InjectLocation) {
    case SEP::FieldLine::InjectionParameters::_InjectShockLocations:
      res+=InjectParticlesSingleFieldLine(spec,iFieldLine);
      break;

    case SEP::FieldLine::InjectionParameters::_InjectBegginingFL:
      if (InjectionParameters::InjectionMomentumModel==SEP::FieldLine::InjectionParameters::_background_sw_temperature) {
        res+=InjectSolarWindIons(spec,iFieldLine);
      }
      else {
        res+=InjectParticleFieldLineBeginning(spec,iFieldLine);
      }
      break;
    default:
      exit(__LINE__,__FILE__,"Error: the option is unknown");
    }
  }

  return res;
}

//=============================================================================
    // DeleteAllParticles
    //=============================================================================
    // Purpose: Delete all particles from all magnetic field lines and segments
    //
    // Description:
    //   Traverses the entire field line structure and removes all particles from
    //   every segment of every field line. This function provides a clean way to
    //   reset the particle population, useful for restarting simulations, clearing
    //   initialization states, or preparing for new injection scenarios.
    //
    // Physics:
    //   - Removes all computational particles while preserving field line geometry
    //   - Maintains field line structure and background plasma data
    //   - Resets particle lists to empty state for fresh initialization
    //
    // Algorithm:
    //   1. Loop through all field lines (0 to nFieldLine-1)
    //   2. For each field line, traverse all segments
    //   3. For each segment, delete all attached particles using PIC framework
    //   4. Reset segment particle list pointers to -1 (empty)
    //   5. Accumulate total count of deleted particles
    //
    // Parameters:
    //   None
    //
    // Returns:
    //   Total number of particles deleted across all field lines
    //
    // Usage:
    //   // Clear all particles before reinitialization:
    //   long int deletedCount = SEP::SolarWind::DeleteAllParticles();
    //
    //   // Then reinitialize with new parameters:
    //   SEP::SolarWind::InitializeSolarWindPopulation(spec, nParticles);
    //
    // Notes:
    //   - Only processes segments assigned to current MPI thread
    //   - Uses PIC::ParticleBuffer::DeleteParticle() for proper memory management
    //   - Preserves field line geometry and background plasma data
    //   - Thread-safe operation respects MPI domain decomposition
    //   - Resets FirstParticleIndex to -1 for each segment
    //
    // Warning:
    //   This function permanently removes all particles. Make sure this is
    //   the intended behavior before calling, especially in production runs.
    //=============================================================================

long int SEP::FieldLine::DeleteAllParticles() {
    namespace FL = PIC::FieldLine;
    namespace PB = PIC::ParticleBuffer;

    long int totalDeletedParticles = 0;

    // Check if field line mode is active and particles are attached to segments
    if ((_PIC_PARTICLE_LIST_ATTACHING_ != _PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_) ||
        (_PIC_FIELD_LINE_MODE_ != _PIC_MODE_ON_)) {
        // Field line particle management not active - return 0
        return 0;
    }

    // Loop through all field lines
    for (int iFieldLine = 0; iFieldLine < FL::nFieldLine; iFieldLine++) {
        FL::cFieldLineSegment* Segment = FL::FieldLinesAll[iFieldLine].GetFirstSegment();

        // Loop through all segments in current field line
        while (Segment != nullptr) {
            // Only process segments assigned to this thread
            if (Segment->Thread == PIC::ThisThread) {
                long int ptr = Segment->FirstParticleIndex;
                long int ptr_next;

                // Delete all particles in this segment
                while (ptr != -1) {
                    ptr_next = PB::GetNext(ptr);
                    PB::DeleteParticle(ptr);
                    totalDeletedParticles++;
                    ptr = ptr_next;
                }

                // Reset segment particle list pointer
                Segment->FirstParticleIndex = -1;
            }

            // Move to next segment
            Segment = Segment->GetNext();
        }
    }

    return totalDeletedParticles;
}
