#include "CometRuntimeInput.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "CometBrightness.h"
#include "CometChemistry.h"
#include "CometDiagnostics.h"
#include "CometModel.h"
#include "CometPhysics.h"
#include "CometRadiativeCooling.h"
#include "RadialMesh.h"
#include "RadialSimulation.h"
#include "pic.h"
#include "src/general/input/amps_input.h"

bool Comet::RuntimeInput::Required=false;

namespace {
namespace In=AMPS::Input;

std::string LoadedInputPath;
std::string LoadedCaseName;

bool Fail(const std::string &message,std::string *error) {
  if (error!=NULL) *error=message;
  return false;
}

bool FailAt(const In::Assignment &item,const std::string &message,
    std::string *error) {
  return Fail(In::ErrorAt(item.where,message,item.text).message,error);
}

std::string Trim(const std::string &text) {
  std::size_t first=0,last=text.size();
  while (first<last && std::isspace(static_cast<unsigned char>(text[first])))
    first++;
  while (last>first &&
      std::isspace(static_cast<unsigned char>(text[last-1]))) last--;
  return text.substr(first,last-first);
}

std::string Upper(std::string text) {
  for (char &c:text)
    c=static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return text;
}

std::vector<std::string> SplitPipe(const std::string &text) {
  std::vector<std::string> result;
  std::size_t begin=0;
  for (;;) {
    const std::size_t delimiter=text.find('|',begin);
    result.push_back(Trim(text.substr(begin,delimiter==std::string::npos ?
        std::string::npos : delimiter-begin)));
    if (delimiter==std::string::npos) break;
    begin=delimiter+1;
  }
  return result;
}

bool ParseDouble(const In::Assignment &item,double *value,std::string *error) {
  char *end=NULL;
  errno=0;
  const double parsed=std::strtod(item.value.c_str(),&end);
  if (errno!=0 || end==item.value.c_str() || *end!='\0' ||
      !std::isfinite(parsed))
    return FailAt(item,"expected a finite number",error);
  *value=parsed;
  return true;
}

bool ParseLong(const In::Assignment &item,long int *value,std::string *error) {
  char *end=NULL;
  errno=0;
  const long parsed=std::strtol(item.value.c_str(),&end,10);
  if (errno!=0 || end==item.value.c_str() || *end!='\0')
    return FailAt(item,"expected an integer",error);
  *value=parsed;
  return true;
}

bool ParseBool(const In::Assignment &item,bool *value,std::string *error) {
  const std::string normalized=Upper(item.value);
  if (normalized=="TRUE" || normalized=="ON") {
    *value=true;
    return true;
  }
  if (normalized=="FALSE" || normalized=="OFF") {
    *value=false;
    return true;
  }
  return FailAt(item,"expected true/false (on/off is also accepted)",error);
}

bool ParsePipeDoubles(const In::Assignment &item,std::size_t count,
    std::vector<double> *values,std::string *error) {
  const std::vector<std::string> fields=SplitPipe(item.value);
  if (fields.size()!=count)
    return FailAt(item,"expected "+std::to_string(count)+
        " pipe-separated numerical fields",error);
  values->clear();
  for (const std::string &field:fields) {
    In::Assignment scalar=item;
    scalar.value=field;
    double value=0.0;
    if (!ParseDouble(scalar,&value,error)) return false;
    values->push_back(value);
  }
  return true;
}

int FindSpecies(const std::string &symbol) {
  char compiled[_MAX_STRING_LENGTH_PIC_];
  for (int species=0;species<PIC::nTotalSpecies;species++) {
    PIC::MolecularData::GetChemSymbol(compiled,species);
    if (symbol==compiled) return species;
  }
  return -1;
}

// The reader makes unknown keys fatal. This is important for scientific
// receipts: a typo must not leave a compiled default active while appearing
// in the archived input. Every successful Take() records ownership.
class cAssignments {
 public:
  explicit cAssignments(const In::Block &block) {
    for (const In::Assignment &item:block.assignments) Items[item.key]=&item;
  }

  const In::Assignment* Take(const std::string &key,std::string *error) {
    const auto found=Items.find(key);
    if (found==Items.end()) {
      Fail("runtime comet input is missing required key '"+key+"'",error);
      return NULL;
    }
    Consumed.insert(key);
    return found->second;
  }

  bool Finish(std::string *error) const {
    for (const auto &entry:Items)
      if (Consumed.count(entry.first)==0)
        return FailAt(*entry.second,"unknown comet runtime key '"+
            entry.first+"'",error);
    return true;
  }

 private:
  std::map<std::string,const In::Assignment*> Items;
  std::set<std::string> Consumed;
};

bool ReadPositiveDouble(cAssignments &input,const std::string &key,
    double *value,std::string *error,bool allowZero=false) {
  const In::Assignment *item=input.Take(key,error);
  if (item==NULL || !ParseDouble(*item,value,error)) return false;
  if (allowZero ? *value<0.0 : !(*value>0.0))
    return FailAt(*item,allowZero ? "value must be non-negative" :
        "value must be positive",error);
  return true;
}

bool ReadPositiveLong(cAssignments &input,const std::string &key,
    long int *value,std::string *error,bool allowZero=false) {
  const In::Assignment *item=input.Take(key,error);
  if (item==NULL || !ParseLong(*item,value,error)) return false;
  if (allowZero ? *value<0 : !(*value>0))
    return FailAt(*item,allowZero ? "value must be a non-negative integer" :
        "value must be a positive integer",error);
  return true;
}

bool ConfigureEmission(cAssignments &input,const std::string &key,int line,
    std::string *error) {
  const In::Assignment *item=input.Take(key,error);
  if (item==NULL) return false;
  const std::vector<std::string> field=SplitPipe(item->value);
  if (field.size()!=11)
    return FailAt(*item,"brightness emission expects SPECIES | TRANSITION | "
        "G_FACTOR_OR_MODEL | G_PROVENANCE | THIN | ISOTROPIC_4PI | "
        "NONE_OR_LINE_CENTER | CROSS_SECTION_M2 | TAU_THRESHOLD | "
        "TAU_PROVENANCE | EFFECTIVE_WAVELENGTH_NM",error);

  const std::string speciesName=Upper(field[0]);
  const int species=FindSpecies(speciesName);
  if (species<0) return FailAt(*item,"unknown compiled emission species '"+
      speciesName+"'",error);
  if (Upper(field[4])!="THIN" || Upper(field[5])!="ISOTROPIC_4PI")
    return FailAt(*item,"only THIN | ISOTROPIC_4PI brightness is implemented",error);

  int model=Comet::Brightness::GFactor::Fixed;
  double fixedG=0.0;
  const std::string modelName=Upper(field[2]);
  if (modelName=="H_LYA_IKL2020") {
    if (speciesName!="H")
      return FailAt(*item,"H_LYA_IKL2020 requires species H",error);
    model=Comet::Brightness::GFactor::HydrogenLymanAlphaIKL2020;
  }
  else if (modelName=="OH_AX00_SA1988_LOWELL_UNQUENCHED") {
    if (speciesName!="OH")
      return FailAt(*item,
          "OH_AX00_SA1988_LOWELL_UNQUENCHED requires species OH",error);
    model=Comet::Brightness::GFactor::HydroxylAX00SchleicherAhearn1988Unquenched;
  }
  else {
    In::Assignment scalar=*item;
    scalar.value=field[2];
    if (!ParseDouble(scalar,&fixedG,error) || !(fixedG>0.0))
      return FailAt(*item,"unknown g-factor model or non-positive fixed g factor",error);
  }

  double crossSection=0.0,threshold=0.0,wavelengthNm=0.0;
  In::Assignment scalar=*item;
  scalar.value=field[7];
  if (!ParseDouble(scalar,&crossSection,error)) return false;
  scalar.value=field[8];
  if (!ParseDouble(scalar,&threshold,error)) return false;
  scalar.value=field[10];
  if (!ParseDouble(scalar,&wavelengthNm,error)) return false;

  const std::string tau=Upper(field[6]);
  int tauMode=Comet::Brightness::OpticalDepthDiagnostic::None;
  if (tau=="NONE") {
    if (crossSection!=0.0 || threshold!=0.0 || Upper(field[9])!="NONE")
      return FailAt(*item,"tau NONE requires zero values and provenance NONE",error);
  }
  else if (tau=="LINE_CENTER") {
    if (!(crossSection>0.0) || !(threshold>0.0) || Upper(field[9])=="NONE")
      return FailAt(*item,"LINE_CENTER requires positive cross section and "
          "threshold plus non-NONE provenance",error);
    tauMode=Comet::Brightness::OpticalDepthDiagnostic::LineCenter;
  }
  else return FailAt(*item,"tau mode must be NONE or LINE_CENTER",error);

  if (model==Comet::Brightness::GFactor::HydroxylAX00SchleicherAhearn1988Unquenched) {
    if (!(wavelengthNm>0.0))
      return FailAt(*item,"the OH energy g factor requires a positive wavelength",error);
  }
  else if (wavelengthNm!=0.0)
    return FailAt(*item,"this g factor requires zero effective wavelength",error);

  Comet::Brightness::EmissionSpecies[line]=species;
  Comet::Brightness::PhotonGFactor[line]=fixedG;
  Comet::Brightness::GFactorModel[line]=model;
  Comet::Brightness::EffectivePhotonWavelength[line]=wavelengthNm*1.0e-9;
  Comet::Brightness::TauDiagnosticMode[line]=tauMode;
  Comet::Brightness::LineCenterCrossSection[line]=crossSection;
  Comet::Brightness::TauFlagThreshold[line]=threshold;
  return Comet::Brightness::SetEmissionMetadata(line,field[1].c_str(),
      field[3].c_str(),field[9].c_str(),error);
}

bool ApplyCometBlock(const In::Block &block,const std::string &selectedCase,
    std::string *error) {
  if (!block.blocks.empty())
    return Fail(In::ErrorAt(block.blocks.front().begin,
        "the resolved comet application block may contain assignments only").message,error);
  cAssignments input(block);

  const In::Assignment *item=input.Take("schema",error);
  if (item==NULL || item->value!="1")
    return item==NULL ? false : FailAt(*item,"supported comet runtime schema is 1",error);
  item=input.Take("geometry",error);
  if (item==NULL || Upper(item->value)!="SPHERICAL")
    return item==NULL ? false : FailAt(*item,
        "this DIM1 executable requires geometry=spherical",error);
  Comet::Model::GeometryMode=Comet::Model::Geometry::Spherical;

  item=input.Take("case.index",error);
  if (item==NULL || item->value!=selectedCase)
    return item==NULL ? false : FailAt(*item,
        "case.index must match the explicit -case selector",error);
  // These receipt fields do not alter the solver, but consuming and requiring
  // them prevents observational provenance from being detached from a case.
  if (input.Take("case.campaign_phase",error)==NULL ||
      input.Take("case.observation_start_utc",error)==NULL ||
      input.Take("case.observation_end_utc",error)==NULL ||
      input.Take("case.exposure_mean_f107p_sfu",error)==NULL) return false;

  if (!ReadPositiveDouble(input,"heliocentric_distance_au",
      &Comet::Model::HeliocentricDistanceAU,error)) return false;
  Comet::Chemistry::HeliocentricDistanceAU=Comet::Model::HeliocentricDistanceAU;

  long int seed=0;
  if (!ReadPositiveLong(input,"random_seed",&seed,error)) return false;
  Comet::Model::RandomSeed=static_cast<std::uint64_t>(seed);
  Comet::Model::RandomSeedConfigured=true;
  Comet::Chemistry::ResidualPhaseSeed=Comet::Model::RandomSeed;

  for (int species=0;species<PIC::nTotalSpecies;species++) {
    Comet::Model::SurfaceProductionRate[species]=0.0;
    Comet::Model::SurfaceDistributionMode[species]=
        Comet::Model::SurfaceDistribution::Disabled;
    Comet::Model::SurfaceTemperature[species]=0.0;
    Comet::Model::SurfaceBulkSpeed[species]=0.0;
    Comet::Model::ExtendedProductionRate[species]=0.0;
    Comet::Model::ExtendedDistributionMode[species]=
        Comet::Model::ExtendedDistribution::Disabled;
  }
  const char *parents[]={"H2O","CO","CO2"};
  for (const char *symbol:parents) {
    const int species=FindSpecies(symbol);
    if (species<0) return Fail(std::string("runtime I3 input requires compiled species ")+symbol,error);
    const std::string prefix=std::string("source.")+symbol+".";
    item=input.Take(prefix+"distribution",error);
    if (item==NULL || Upper(item->value)!="UNIFORM")
      return item==NULL ? false : FailAt(*item,"I3 direct sources require uniform distribution",error);
    Comet::Model::SurfaceDistributionMode[species]=
        Comet::Model::SurfaceDistribution::Uniform;
    if (!ReadPositiveDouble(input,prefix+"rate_s-1",
        &Comet::Model::SurfaceProductionRate[species],error) ||
        !ReadPositiveDouble(input,prefix+"temperature_k",
        &Comet::Model::SurfaceTemperature[species],error) ||
        !ReadPositiveDouble(input,prefix+"bulk_speed_m_s",
        &Comet::Model::SurfaceBulkSpeed[species],error,true)) return false;
  }

  if (!ReadPositiveDouble(input,"collision.temperature_floor_k",
      &Comet::Physics::CollisionTemperatureFloor,error) ||
      !ReadPositiveDouble(input,"collision.minimum_sample_particles",
      &Comet::Physics::CollisionTemperatureMinimumSampleParticleNumber,error,true))
    return false;

  bool enabled=false;
  item=input.Take("chemistry.enabled",error);
  if (item==NULL || !ParseBool(*item,&enabled,error)) return false;
  Comet::Chemistry::Enabled=enabled;
  if (!enabled) return FailAt(*item,"the compiled I3 network requires chemistry.enabled=true",error);
  item=input.Take("chemistry.residual_sampling",error);
  if (item==NULL) return false;
  const std::string residual=Upper(item->value);
  if (residual=="SYSTEMATIC")
    Comet::Chemistry::ResidualSamplingMode=
        Comet::Chemistry::ResidualSampling::RandomPhaseSystematic;
  else if (residual=="BERNOULLI")
    Comet::Chemistry::ResidualSamplingMode=
        Comet::Chemistry::ResidualSampling::IndependentBernoulli;
  else return FailAt(*item,"chemistry.residual_sampling must be systematic or bernoulli",error);
  for (std::size_t reaction=0;
       reaction<Comet::Chemistry::GetCompiledReactionCount();reaction++) {
    const char *id=Comet::Chemistry::GetCompiledReactionId(reaction);
    if (id==NULL) return Fail("compiled chemistry contains an unnamed reaction",error);
    double rate=0.0;
    if (!ReadPositiveDouble(input,std::string("chemistry.rate.")+id+
        "_1au_s-1",&rate,error)) return false;
    std::string chemistryError;
    if (!Comet::Chemistry::ConfigureReactionRate(id,rate,&chemistryError))
      return Fail(chemistryError,error);
  }

  item=input.Take("diagnostics.enabled",error);
  if (item==NULL || !ParseBool(*item,&Comet::Diagnostics::Enabled,error)) return false;

  std::fill(Comet::RadiativeCooling::SpeciesModel,
      Comet::RadiativeCooling::SpeciesModel+PIC::nTotalSpecies,0);
  item=input.Take("cooling.H2O.model",error);
  if (item==NULL) return false;
  const int water=FindSpecies("H2O");
  if (water<0) return Fail("runtime I3 input requires compiled species H2O",error);
  if (Upper(item->value)!="GOMBOSI1986_H2O")
    return FailAt(*item,"supported H2O cooling model is GOMBOSI1986_H2O",error);
  Comet::RadiativeCooling::SpeciesModel[water]=static_cast<int>(
      CometaryPhysics::RadiativeCooling::Model::Gombosi1986Water);

  if (!ReadPositiveDouble(input,"nucleus.radius_m",
      &Comet::RadialModel::NucleusRadius,error) ||
      !ReadPositiveDouble(input,"mesh.outer_radius_m",
      &Comet::RadialModel::DomainOuterRadius,error)) return false;
  Comet::RadialModel::DomainRadialExtent=-1.0;
  long int radialCells=0;
  if (!ReadPositiveLong(input,"mesh.radial_cell_count",&radialCells,error) ||
      radialCells>std::numeric_limits<int>::max()) return false;
  Comet::RadialModel::RequestedCellNumber=static_cast<int>(radialCells);
  if (!ReadPositiveDouble(input,"mesh.surface_resolution_m",
      &Comet::RadialModel::SurfaceResolution,error)) return false;

  long int referenceInjection=0;
  if (!ReadPositiveLong(input,"simulation.reference_injection_particles",
      &referenceInjection,error)) return false;
  PIC::ParticleWeightTimeStep::maxReferenceInjectedParticleNumber=
      static_cast<double>(referenceInjection);
  Comet::RadialSimulation::Enabled=true;
  if (!ReadPositiveDouble(input,"simulation.courant_number",
      &Comet::RadialSimulation::CourantNumber,error) ||
      !ReadPositiveDouble(input,"simulation.inner_boundary_temperature_k",
      &Comet::RadialSimulation::InnerBoundaryTemperature,error) ||
      !ReadPositiveLong(input,"simulation.sampling_length",
      &Comet::RadialSimulation::SamplingLength,error) ||
      !ReadPositiveLong(input,"simulation.weight_adaptation_interval",
      &Comet::RadialSimulation::ParticleWeightAdaptationInterval,error) ||
      !ReadPositiveLong(input,"simulation.load_balance_interval",
      &Comet::RadialSimulation::LoadBalanceInterval,error) ||
      !ReadPositiveLong(input,"simulation.max_iterations",
      &Comet::RadialSimulation::MaxIterations,error)) return false;

  item=input.Take("simulation.weight_adaptation_enabled",error);
  if (item==NULL || !ParseBool(*item,
      &Comet::RadialSimulation::ParticleWeightAdaptationEnabled,error)) return false;
  if (!ReadPositiveLong(input,"simulation.particles_per_cell.minimum",
      &Comet::RadialSimulation::MinimumParticlesPerCell,error) ||
      !ReadPositiveLong(input,"simulation.particles_per_cell.maximum",
      &Comet::RadialSimulation::MaximumParticlesPerCell,error)) return false;
  long int chemicalMinimum=0,chemicalMaximum=0;
  if (!ReadPositiveLong(input,"simulation.chemical_only_particles_per_cell.minimum",
      &chemicalMinimum,error) ||
      !ReadPositiveLong(input,"simulation.chemical_only_particles_per_cell.maximum",
      &chemicalMaximum,error)) return false;
  bool lockDirect=false;
  item=input.Take("simulation.lock_direct_source_weights",error);
  if (item==NULL || !ParseBool(*item,&lockDirect,error)) return false;

  if (!ReadPositiveDouble(input,"simulation.emergency_buffer_utilization",
      &Comet::RadialSimulation::EmergencyParticleBufferUtilization,error) ||
      !ReadPositiveDouble(input,"simulation.allocation_low_watermark",
      &Comet::RadialSimulation::AllocationLowWatermark,error) ||
      !ReadPositiveDouble(input,"simulation.allocation_reserve_fraction",
      &Comet::RadialSimulation::AllocationReserveFraction,error) ||
      !ReadPositiveLong(input,"simulation.max_face_transfer_multiplicity",
      &Comet::RadialSimulation::MaximumFaceTransferMultiplicity,error) ||
      !ReadPositiveLong(input,"simulation.max_photochemical_products_per_reaction",
      &Comet::RadialSimulation::MaximumPhotochemicalProductsPerReaction,error) ||
      !ReadPositiveDouble(input,"simulation.convergence_tolerance",
      &Comet::RadialSimulation::ConvergenceTolerance,error)) return false;
  long int convergenceWindows=0;
  if (!ReadPositiveLong(input,"simulation.convergence_windows",
      &convergenceWindows,error) || convergenceWindows>std::numeric_limits<int>::max())
    return false;
  Comet::RadialSimulation::RequiredConvergedWindows=
      static_cast<int>(convergenceWindows);
  Comet::RadialSimulation::IndividualCorrectionWindowsEnabled=false;
  Comet::RadialSimulation::MaximumIndividualCorrectionWindowSplits=-1;
  Comet::RadialSimulation::Stage9ScienceDiagnosticsEnabled=false;

  for (int species=0;species<PIC::nTotalSpecies;species++) {
    char symbol[_MAX_STRING_LENGTH_PIC_];
    PIC::MolecularData::GetChemSymbol(symbol,species);
    if (!ReadPositiveDouble(input,std::string("transport.characteristic_speed.")+
        symbol+"_m_s",&Comet::RadialSimulation::CharacteristicSpeed[species],error))
      return false;
    const bool direct=Comet::Model::SurfaceProductionRate[species]>0.0;
    Comet::RadialSimulation::SpeciesMinimumParticlesPerCell[species]=
        direct ? Comet::RadialSimulation::MinimumParticlesPerCell : chemicalMinimum;
    Comet::RadialSimulation::SpeciesMaximumParticlesPerCell[species]=
        direct ? Comet::RadialSimulation::MaximumParticlesPerCell : chemicalMaximum;
    Comet::RadialSimulation::ParticleWeightLocked[species]=direct && lockDirect ? 1 : 0;
    Comet::RadialSimulation::PopulationPolicySource[species]=direct ?
        Comet::RadialSimulation::PolicySource::DirectSource :
        Comet::RadialSimulation::PolicySource::ChemicalOnly;
    Comet::RadialSimulation::PopulationClass[species]=direct ?
        Comet::RadialSimulation::PolicyClass::DirectSource :
        Comet::RadialSimulation::PolicyClass::ChemicalOnly;
    Comet::RadialSimulation::IndividualCorrectionWindowLower[species]=0.0;
    Comet::RadialSimulation::IndividualCorrectionWindowTarget[species]=0.0;
    Comet::RadialSimulation::IndividualCorrectionWindowUpper[species]=0.0;
    Comet::RadialSimulation::IndividualCorrectionWindowPolicySource[species]=0;
  }

  if (!Comet::Brightness::ResetRuntimeConfiguration(error)) return false;
  item=input.Take("brightness.enabled",error);
  if (item==NULL || !ParseBool(*item,&Comet::Brightness::Enabled,error)) return false;
  if (!Comet::Brightness::Enabled)
    return FailAt(*item,"the I3 runtime model requires brightness.enabled=true",error);
  item=input.Take("brightness.impact_parameter_grid",error);
  if (item==NULL) return false;
  std::vector<std::string> grid=SplitPipe(item->value);
  if (grid.size()!=4)
    return FailAt(*item,"brightness grid expects LOG|LINEAR | MIN_M | MAX_M | COUNT",error);
  const std::string gridMode=Upper(grid[0]);
  if (gridMode=="LOG")
    Comet::Brightness::ImpactParameterGridMode=Comet::Brightness::GridMode::Logarithmic;
  else if (gridMode=="LINEAR")
    Comet::Brightness::ImpactParameterGridMode=Comet::Brightness::GridMode::Linear;
  else return FailAt(*item,"brightness grid mode must be LOG or LINEAR",error);
  In::Assignment scalar=*item;
  scalar.value=grid[1];
  if (!ParseDouble(scalar,&Comet::Brightness::MinimumImpactParameter,error)) return false;
  scalar.value=grid[2];
  if (!ParseDouble(scalar,&Comet::Brightness::MaximumImpactParameter,error)) return false;
  scalar.value=grid[3];
  if (!ParseLong(scalar,&Comet::Brightness::ImpactParameterCount,error)) return false;

  item=input.Take("brightness.observation_geometry",error);
  std::vector<double> geometry;
  if (item==NULL || !ParsePipeDoubles(*item,9,&geometry,error)) return false;
  const double radial[3]={geometry[0],geometry[1],geometry[2]};
  const double los[3]={geometry[3],geometry[4],geometry[5]};
  const double impact[3]={geometry[6],geometry[7],geometry[8]};
  if (!Comet::Brightness::SetObservationGeometry(radial,los,impact,error)) return false;
  double radialVelocityKmS=0.0;
  item=input.Take("brightness.heliocentric_radial_velocity_km_s",error);
  if (item==NULL || !ParseDouble(*item,&radialVelocityKmS,error)) return false;
  double cometVelocity[3];
  for (int component=0;component<3;component++)
    cometVelocity[component]=radialVelocityKmS*1.0e3*radial[component];
  if (!Comet::Brightness::SetCometHeliocentricVelocity(cometVelocity,error))
    return false;

  double solarFluxCm2=0.0;
  if (!ReadPositiveDouble(input,
      "brightness.solar_lya_flux_1au_photon_cm-2_s-1",&solarFluxCm2,error))
    return false;
  Comet::Brightness::SolarLymanAlphaFlux1AU=solarFluxCm2*1.0e4;
  item=input.Take("brightness.velocity_domain_policy",error);
  if (item==NULL) return false;
  const std::string velocityPolicy=Upper(item->value);
  if (velocityPolicy=="ABORT")
    Comet::Brightness::OutOfRangeVelocityPolicy=
        Comet::Brightness::VelocityDomainPolicy::Abort;
  else if (velocityPolicy=="EXCLUDE_AND_REPORT")
    Comet::Brightness::OutOfRangeVelocityPolicy=
        Comet::Brightness::VelocityDomainPolicy::ExcludeAndReport;
  else return FailAt(*item,
      "brightness.velocity_domain_policy must be ABORT or EXCLUDE_AND_REPORT",error);

  if (!ConfigureEmission(input,"brightness.emission.OH",0,error) ||
      !ConfigureEmission(input,"brightness.emission.H",1,error)) return false;
  Comet::Brightness::EmissionLineCount=2;

  bool dust=false;
  item=input.Take("dust.enabled",error);
  if (item==NULL || !ParseBool(*item,&dust,error)) return false;
  if (dust) return FailAt(*item,"dust is not compiled into the I3 executable",error);

  return input.Finish(error);
}

std::uint64_t HashCanonical(const std::string &text) {
  std::uint64_t hash=UINT64_C(1469598103934665603);
  for (unsigned char byte:text) {
    hash^=byte;
    hash*=UINT64_C(1099511628211);
  }
  return hash;
}
}  // namespace

bool Comet::RuntimeInput::LoadAndApply(const std::string &path,
    const std::string &caseName,std::string *error) {
  In::Document document;
  In::Status status=In::ParseFile(path,&document);
  std::string localError=status.ok() ? std::string() : status.message;
  const In::Block *application=NULL;
  In::Block resolved;
  if (status.ok()) {
    status=In::RequireNoAssignments(document.root);
    if (status.ok()) status=In::RequireOnlyBlocks(document.root,{"application"},
        "the AMPS runtime-input root");
    if (status.ok()) {
      if (document.root.blocks.size()!=1 ||
          document.root.blocks.front().kind!="application" ||
          document.root.blocks.front().name!="comet")
        status=In::Status::Error(
            "runtime comet input requires exactly '#begin application: comet'");
      else application=&document.root.blocks.front();
    }
    if (status.ok()) status=In::ResolveCase(*application,caseName,&resolved);
    if (!status.ok()) localError=status.message;
  }

  int localOK=status.ok() ? 1 : 0,globalOK=0;
  MPI_Allreduce(&localOK,&globalOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!globalOK)
    return Fail(localOK ? "another MPI rank could not parse/select the runtime input" :
        localError,error);

  const std::string canonical=In::ToCanonicalText(resolved);
  const unsigned long long localHash=HashCanonical(canonical);
  unsigned long long minimumHash=0,maximumHash=0;
  MPI_Allreduce(&localHash,&minimumHash,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,
      MPI_GLOBAL_COMMUNICATOR);
  MPI_Allreduce(&localHash,&maximumHash,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,
      MPI_GLOBAL_COMMUNICATOR);
  if (minimumHash!=maximumHash)
    return Fail("MPI ranks resolved different runtime comet inputs",error);

  std::string applyError;
  localOK=ApplyCometBlock(resolved,caseName,&applyError) ? 1 : 0;
  MPI_Allreduce(&localOK,&globalOK,1,MPI_INT,MPI_MIN,MPI_GLOBAL_COMMUNICATOR);
  if (!globalOK)
    return Fail(localOK ? "another MPI rank rejected the runtime comet input" :
        applyError,error);

  LoadedInputPath=path;
  LoadedCaseName=caseName;
  return true;
}

const std::string& Comet::RuntimeInput::LoadedPath() { return LoadedInputPath; }
const std::string& Comet::RuntimeInput::SelectedCase() { return LoadedCaseName; }
