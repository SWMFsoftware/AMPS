#include "thermal_population.h"
#include "sep.h"
#include "transport_common.h"
#include "util/sep_background_runtime.h"
#include "util/sep_source_lifecycle.h"
#include "util/sep_maxwell_juttner_volume.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>
#include <vector>

namespace SEP { namespace ParticleSource { namespace {
namespace FL=PIC::FieldLine;
namespace V=Sources::MaxwellJuttnerVolume;
Transport::Status Invalid(const char* text){return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,text);}
struct State {
  // Fluid-rest thermodynamic inputs and the actual received Cartesian bulk
  // velocity. gamma is the local simulation-volume Lorentz normalization.
  double density=0,temperature=0,gamma=1;
  std::array<double,3> bulk{{0,0,0}};
};
State Interpolate(FL::cFieldLineSegment* segment,double fraction) {
  // Use the same piecewise-linear endpoint state as transport. In particular,
  // neither a radial temperature law nor a fixed 400 km/s drift replaces the
  // SWMF n/T/U data selected by the user for this initial population.
  State state;double n[2],t[2],u[2][3];
  auto* begin=segment->GetBegin();auto* end=segment->GetEnd();
  begin->GetDatum(FL::DatumAtVertexPlasmaDensity,&n[0]);end->GetDatum(FL::DatumAtVertexPlasmaDensity,&n[1]);
  begin->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t[0]);end->GetDatum(FL::DatumAtVertexPlasmaTemperature,&t[1]);
  begin->GetPlasmaVelocity(u[0]);end->GetPlasmaVelocity(u[1]);
  state.density=(1-fraction)*n[0]+fraction*n[1];state.temperature=(1-fraction)*t[0]+fraction*t[1];
  for(int d=0;d<3;++d)state.bulk[d]=(1-fraction)*u[0][d]+fraction*u[1][d];
  const auto status=V::BulkGamma(state.bulk,&state.gamma);
  if(!status.ok())throw std::logic_error(status.message);
  if(!std::isfinite(state.density)||state.density<0||!std::isfinite(state.temperature)||state.temperature<=0)
    throw std::logic_error("imported thermal volume requires finite nonnegative SWMF density and positive SWMF temperature");
  return state;
}
double DensityPerCoordinate(FL::cFieldLineSegment* segment,int line,double fraction) {
  const State state=Interpolate(segment,fraction);
  // Flux-tube area is the existing authoritative Phi/|B| measure. Native
  // spatial coordinates/areas describe the simulation-frame volume, whose
  // equilibrium number density is gamma_U*n0. Uniform coordinate proposals
  // are corrected by this full Jacobian, not an assumed uniform tube area.
  return state.density*state.gamma*FieldLine::FluxTubeGeometry::AreaAtSegmentFractionM2(segment,line,fraction)*segment->GetLength();
}
} // namespace

Transport::Status PopulateImportedThermalInitialCondition() {
  using namespace Transport;
  if(Background::ConfiguredProvider()!=Background::Provider::Swmf||
     PICAdapter::NativeDriverWindowActive()||PIC::ParticleBuffer::GetAllPartNum()!=0)
    return Invalid("imported thermal prepopulation requires a declared fresh empty population before source-window freeze");
  // CON_bline explicitly exports rho/cProtonMass. The requested single H+
  // campaign therefore has a declared proton density; a compiled electron or
  // heavy-ion list does not authorize inventing ambient species abundances.
  if(PIC::nTotalSpecies!=1)return Status::Error(StatusCode::UnsupportedConfiguration,
      "imported multi-species thermal prepopulation requires explicit SWMF composition/abundances");
  std::string symbol=PIC::MolecularData::GetChemSymbol(0);
  std::transform(symbol.begin(),symbol.end(),symbol.begin(),[](unsigned char c){return static_cast<char>(std::toupper(c));});
  if(symbol!="H_PLUS"&&symbol!="H+")return Status::Error(StatusCode::UnsupportedConfiguration,
      "the declared SWMF proton-density initial condition requires a compiled H+ species or explicit composition mapping");
  const int species=0;const double mass=PIC::MolecularData::GetMass(species);
  const double nominalWeight=PIC::ParticleWeightTimeStep::GlobalParticleWeight[species];
  if(!std::isfinite(mass)||mass<=0||!std::isfinite(nominalWeight)||nominalWeight<=0)
    return Invalid("thermal prepopulation requires initialized positive species mass and global statistical weight");
  if(_INDIVIDUAL_PARTICLE_WEIGHT_MODE_!=_INDIVIDUAL_PARTICLE_WEIGHT_ON_)
    return Status::Error(StatusCode::UnsupportedConfiguration,
        "spatially varying imported thermal volume requires individual statistical weight corrections");
  Birth::Record birth;auto status=Birth::InitialCondition(Background::SimulationTimeSeconds(),&birth);
  if(!status.ok())return status;

  // Validate the complete native source input before reserving an event,
  // sampling momentum or allocating records. The midpoint/end evaluations
  // also validate the geometry measure used by count proposals and weights.
  std::vector<std::vector<double> > proposals(FL::nFieldLine);
  try {
    for(int line=0;line<FL::nFieldLine;++line) {
      const int count=FL::FieldLinesAll[line].GetTotalSegmentNumber();proposals[line].resize(count);
      for(int segment=0;segment<count;++segment) {
        auto* s=FL::FieldLinesAll[line].GetSegment(segment);
        const double integral=(DensityPerCoordinate(s,line,0)+4*DensityPerCoordinate(s,line,0.5)+DensityPerCoordinate(s,line,1))/6;
        const double expected=std::min(200.0,integral/nominalWeight);
        if(!std::isfinite(expected)||expected<0)throw std::logic_error("thermal population proposal weight is invalid");
        proposals[line][segment]=expected;
      }
    }
  }catch(const std::exception& error){return Status::Error(StatusCode::InvalidParticleState,error.what());}

  for(int line=0;line<FL::nFieldLine;++line) {
    // Reserve one persisted event per physical line on every PIC rank before
    // owner filtering. Macro keys use initial segment/sample indices, never
    // allocation pointers, sampled momenta, successful counts or rank ids.
    const auto event=SourceLifecycle::BeginNativeEvent(SourceLifecycle::Kind::SolarWindPopulation,species,line);
    for(std::size_t segment=0;segment<proposals[line].size();++segment) {
      auto* s=FL::FieldLinesAll[line].GetSegment(static_cast<int>(segment));
      if(s->Thread!=PIC::ThisThread)continue;
      const double expected=proposals[line][segment];if(expected==0)continue;
      const std::uint64_t segmentKey=static_cast<std::uint64_t>(segment)<<32;
      auto countRandom=SourceLifecycle::Random(event,segmentKey|UINT64_C(0xffffffff),UINT64_C(0x564f4c434e54));
      unsigned count=static_cast<unsigned>(expected);
      if(countRandom.UniformOpen01()<expected-count)++count;
      for(unsigned macro=0;macro<count;++macro) {
        const auto key=segmentKey|macro;
        auto random=SourceLifecycle::Random(event,key,UINT64_C(0x564f4c554d45));
        const double fraction=random.UniformOpen01();
        State state;
        try{state=Interpolate(s,fraction);}catch(const std::exception& error){return Status::Error(StatusCode::InvalidParticleState,error.what());}
        V::Input input;input.properTemperatureK=state.temperature;input.restMassKg=mass;input.bulkVelocityMPerS=state.bulk;
        V::Sample sampled;status=V::Draw(input,&random,&sampled);if(!status.ok())return status;
        // Simpson's value above chooses computational effort only. The actual
        // importance weight uses n(s)*gamma(s)*A(s)*L / E[N], so count rounding
        // and spatial interpolation do not renormalize the physical density.
        // All momenta are local-plasma values; full 3-D sampling followed by
        // gyro projection gives the volume marginal, not a boundary-flux law.
        double correction=0;
        try{correction=DensityPerCoordinate(s,line,fraction)/(expected*nominalWeight);}
        catch(const std::exception& error){return Status::Error(StatusCode::InvalidParticleState,error.what());}
        if(!std::isfinite(correction)||correction<=0)return Invalid("thermal volume statistical correction is nonpositive or unrepresentable");
        double momentum[3];for(int d=0;d<3;++d)momentum[d]=sampled.restMomentumKgMPerS[d];
        const long particle=FL::InjectParticle_default(species,momentum,correction,line,
            static_cast<int>(segment),segment+fraction);
        if(particle<0)return Invalid("thermal volume source failed to allocate its fresh particle");
        status=PICAdapter::InitializeParticleTransportStateWithIdentity(particle,
            SourceLifecycle::ParticleIdentity(event,key),birth);
        if(!status.ok())return status;
      }
    }
  }
  return Status::Ok();
}
} }
