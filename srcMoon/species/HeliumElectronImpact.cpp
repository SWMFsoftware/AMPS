#include "HeliumElectronImpact.h"
#include "PhysicalReactionEvents.h"
#include <filesystem>

namespace Moon { namespace Species { namespace HeliumElectronImpact {
namespace D=ReactionData;
namespace K=ReducedKinematics;
namespace Math=::AMPS::ChemicalDatabase::ReducedModels;
namespace {
constexpr auto He=ElectronImpact::NeutralSpecies::He;
constexpr double PilotThermalEnergyEv=10;
// D05 frozen original used by U16. Retain the fit's own energy threshold
// and validate the caller's declaration before reading its coefficients.
constexpr char D05Sha256[]="470ea3078c64159eece6a469314dde2356e2d6ccaf618f4eaecb46649d61a5cd";
constexpr char EnergyReference[]="P02:separated-neutral-ground-isotope-atoms-and-free-electron-reference-zero";

bool Convert(B::Real value,double* output){
  const double projected=static_cast<double>(value);
  if(!std::isfinite(projected)||(value!=0&&projected==0))return false;
  *output=projected;
  return true;
}
}

Status Provider::Prepare(const PhysicalStateProfile::Registry& profile,
    const std::string& nativeD05Path,const std::string& nativeD05Sha256,
    int neutralCompiledSpecies,int ionCompiledSpecies,Provider* output){
  if(!output||neutralCompiledSpecies<0||ionCompiledSpecies<0||neutralCompiledSpecies==ionCompiledSpecies)
    return {Code::InvalidInput,"distinct explicit compiled indices for helium and its ion are required"};
  if(!profile.loaded()||!std::filesystem::is_regular_file(nativeD05Path))
    return {Code::Unavailable,"authenticated helium state profile and original D05 rate table required"};
  std::string error;
  if(nativeD05Sha256!=D05Sha256||!D::VerifyFileSha256(nativeD05Path,nativeD05Sha256,&error))
    return {Code::CorruptSource,"frozen D05 rate table identity mismatch: "+error};
  Provider candidate;
  if(!candidate.rates_.LoadNativeTable(nativeD05Path,&error)||
     !candidate.rates_.GetCoefficient(He,&candidate.source_,&error))
    return {Code::CorruptSource,error};
  PhysicalStateProfile::Constant elementaryCharge;
  if(!profile.GetConstant("elementary charge","C",&elementaryCharge,&error)||elementaryCharge.value!=1.602176634e-19L)
    return {Code::CorruptSource,"exact SI elementary charge unavailable from the state profile"};
  candidate.electronVoltJ_=elementaryCharge.value;
  candidate.profileIdentity_=profile.identity();
  candidate.d05Identity_="D05/Voronov:cfit.dat:"+nativeD05Sha256+":row"+std::to_string(candidate.source_.sourceLine);

  constexpr const char* ids[]={"4He:ground","4He+:ground","electron:free"};
  constexpr std::int64_t charges[]={0,1,-1};
  const int compiled[]={neutralCompiledSpecies,ionCompiledSpecies,-1};
  candidate.catalog_.resize(3);
  for(unsigned i=0;i<3;++i){
    K::ResolvedState state;
    K::MassValue physical;
    if(!profile.GetState(ids[i],&state,&error)||!profile.GetPhysicalRestMass(ids[i],&physical,&error))
      return {Code::Unavailable,"required isotope/state/physical mass unavailable: "+std::string(ids[i])+"; "+error};
    const auto& definition=state.definition();
    const auto& nuclei=definition.nuclei;
    const bool expectedNuclei=i==2?nuclei.empty():
        nuclei.size()==1&&nuclei.front().identity=="4He"&&nuclei.front().count==1;
    if(!expectedNuclei||definition.elementaryCharge!=charges[i]||
       definition.excitation.joules!=0||
       (i==1?definition.formation.joules<=0:definition.formation.joules!=0))
      return {Code::CorruptSource,"helium/electron isotope, charge or ground-reference energy mismatch"};
    double internal=0;
    const auto status=PhysicalReactionEvents::ReferenceCatalogEntry(state,compiled[i],&candidate.catalog_[i],&internal);
    if(!status.ok())return {Code::CorruptSource,status.message};
    auto& bathState=candidate.states_[i];
    bathState.identity=definition.identity;
    bathState.state=definition.state;
    bathState.sourceIdentity=definition.formation.sourceIdentity;
    bathState.energyReferenceIdentity=EnergyReference;
    bathState.kineticMassKg=state.referenceMassKg();
    bathState.physicalRestMassKg=physical.kg;
    bathState.stateEnergyJ=definition.formation.joules;
    bathState.chargeElementary=charges[i];
    for(const auto& n:nuclei)bathState.nuclei.push_back({n.identity,n.count});
  }
  candidate.reaction_={"D05-prescribed-4He-electron-impact-ground-ion",
      {{0,1,R::Population::Tracked},{2,1,R::Population::Bath}},
      {{1,1,R::Population::Tracked},{2,2,R::Population::Bath}}};
  if(!R::ValidateReaction(candidate.catalog_,candidate.reaction_,&error))return {Code::CorruptSource,error};
  *output=std::move(candidate);
  return {};
}

B::MaxwellianInput Provider::Incoming(const Bath& bath,const Evaluation& evaluation)const {
  B::MaxwellianInput input;
  input.derivativeSourceIdentity=d05Identity_;
  input.temperatureDomainIdentity="prescribed-stationary-helium-10eV-on-D05-native-domain";
  input.thermalEnergyJ=bath.thermalEnergyEv*electronVoltJ_; // k_B T in J
  input.minimumThermalEnergyJ=source_.temperatureMinEv*electronVoltJ_;
  input.maximumThermalEnergyJ=source_.temperatureMaxEv*electronVoltJ_;
  input.logarithmicRateDerivative=evaluation.logSlope;
  input.thresholdJ=source_.ionizationEnergyEv*electronVoltJ_;
  input.bulkVelocityParentRestMS=bath.bulkVelocityMS;
  // These conditions define this narrow reference. Evaluate requires
  // measured drift and target motion to use a separate, reviewed provider.
  input.isotropicUnshiftedParentRest=true;
  input.crossSectionTemperatureIndependent=true;
  input.collisionWeightRotationallyInvariant=true;
  return input;
}

Status Provider::Evaluate(const Bath& bath,Evaluation* output) const {
  if(!output||!std::isfinite(bath.thermalEnergyEv)||!std::isfinite(bath.densityM3)||
      bath.thermalEnergyEv<=0||bath.densityM3<0||bath.frameIdentity.empty()||
      !Math::Finite(bath.bulkVelocityMS)||!Math::Finite(bath.targetVelocityMS))
    return {Code::InvalidInput,"finite positive temperature, nonnegative density and explicit inertial bath required"};
  if(!rates_.loaded())return {Code::Unavailable,"helium pilot source is not prepared"};
  if(bath.thermalEnergyEv!=PilotThermalEnergyEv||bath.bulkVelocityMS!=B::Vector3{}||bath.targetVelocityMS!=B::Vector3{})
    return {Code::OutsideDomain,"helium pilot requires 10 eV and stationary bath and target in the declared frame"};
  Evaluation value;
  std::string error;
  if(!rates_.RateCoefficientM3S(He,bath.thermalEnergyEv,&value.coefficientM3S,&error)||
     !rates_.IonizationRateS(He,bath.thermalEnergyEv,bath.densityM3,&value.rateS,&error)||
     !rates_.LogTemperatureDerivative(He,bath.thermalEnergyEv,&value.logSlope,&error))
    return {Code::NumericalFailure,error};
  value.fitThresholdEv=source_.ionizationEnergyEv;
  value.physicalThresholdJ=states_[1].stateEnergyJ;
  const auto status=B::MaxwellianCollisionMeanEnergy(Incoming(bath,value),&value.meanIncomingEnergyJ);
  if(!status.ok())return status;
  *output=value;
  return {};
}

Status Provider::Compose(const Bath& bath,const Admission& admission,const R::ParentState& parent,
    const R::TraceInterval& interval,double eventEpochS,double ionReferenceWeight,
    const R::OpticalDepthClock& ionClock,Composition* output) const {
  if(!output)return {Code::InvalidInput,"helium composition output required"};
  Evaluation evaluation;
  auto status=Evaluate(bath,&evaluation);
  if(!status.ok())return status;
  const auto& evidence=admission.evidence;
  if(!admission.planner.qualified||admission.crosswalkPath.empty()||admission.modelErrorPath.empty()||
     !std::isfinite(admission.maximumAbsoluteLogSlopeError)||admission.maximumAbsoluteLogSlopeError<0||
     !std::isfinite(admission.thresholdCrosswalkMaximumDifferenceJ)||admission.thresholdCrosswalkMaximumDifferenceJ<0||
     !std::filesystem::is_regular_file(admission.crosswalkPath)||!std::filesystem::is_regular_file(admission.modelErrorPath))
    return {Code::Unavailable,"reviewed helium threshold/state, conditional-moment and heavy-recoil evidence required"};
  std::string error;
  if(evidence.sourceProfileIdentity!=profileIdentity_||evidence.sourceProfileSha256!=profileIdentity_||
     !D::VerifyFileSha256(admission.crosswalkPath,evidence.stateCrosswalkSha256,&error)||
     !D::VerifyFileSha256(admission.modelErrorPath,evidence.modelErrorEvidenceSha256,&error))
    return {Code::CorruptSource,"helium admission does not bind the prepared state profile and evidence: "+error};
  if(evaluation.rateS==0||!interval.valid()||parent.species!=0||parent.velocityMS!=R::Vector3{}||
     parent.internalEnergyJ!=0||evidence.frameIdentity!=bath.frameIdentity||
     admission.planner.eventFrame!=bath.frameIdentity)
    return {Code::OutsideDomain,"actual accepted nonzero-rate event on the stationary helium reference required"};

  B::ElectronImpactInput request;
  request.evidence=evidence;
  request.numerical=admission.numerical;
  request.event.positionM={parent.positionM[0],parent.positionM[1],parent.positionM[2]};
  request.event.physicalStartS=interval.startS();
  request.event.elapsedS=interval.elapsedS();
  request.event.remainingS=interval.remainingS();
  request.event.physicalEventS=eventEpochS;
  request.event.representedNumber=static_cast<B::Real>(parent.referenceWeight)*parent.individualCorrection;
  request.event.frameIdentity=bath.frameIdentity;
  request.event.timeScaleIdentity=evidence.timeScaleIdentity;
  request.parent=states_[0];request.ion=states_[1];request.electron=states_[2];
  request.incoming=Incoming(bath,evaluation);
  request.sourceThresholdCrosswalkMaximumDifferenceJ=admission.thresholdCrosswalkMaximumDifferenceJ;
  request.elementaryChargeC=electronVoltJ_;
  // Check both uncertainty ends before composing the central prescription.
  // The lower end must have eligible incident energy. The upper end makes
  // the omitted-recoil bound conservative for the whole supported interval.
  const auto central=request.incoming.logarithmicRateDerivative;
  request.incoming.logarithmicRateDerivative=central-admission.maximumAbsoluteLogSlopeError;
  B::Real minimum=0;
  status=B::MaxwellianCollisionMeanEnergy(request.incoming,&minimum);
  if(!status.ok())return status;
  if(minimum<evaluation.physicalThresholdJ)
    return {Code::OutsideDomain,"admitted moment uncertainty reaches below physical helium threshold"};
  request.incoming.logarithmicRateDerivative=central+admission.maximumAbsoluteLogSlopeError;
  B::Result upper;
  status=B::ElectronImpactInheritedVelocity(request,&upper);
  if(!status.ok())return status;
  request.incoming.logarithmicRateDerivative=central;
  Composition value;
  status=B::ElectronImpactInheritedVelocity(request,&value.reservoir);
  if(!status.ok())return status;
  value.reservoir.error=upper.error;

  value.policy=admission.planner;
  if(!value.policy.products.empty()||value.policy.bath.supplied)
    return {Code::InvalidInput,"caller policy must supply only authorization and arithmetic bounds"};
  // EventPlan computes tracked changes from the source catalog. The bath
  // moments below come independently from the two outgoing and one incoming
  // electrons, in the same frame and per one physical reaction.
  auto& exchange=value.policy.bath;
  exchange.supplied=true;
  if(!Convert(value.reservoir.bathMassChangeKg,&exchange.modelMassChangeKg)||
     !Convert(value.reservoir.bathEnergyChangeJ,&exchange.energyChangeJ))
    return {Code::NumericalFailure,"electron bath mass/energy exceeds the planner representation"};
  for(unsigned i=0;i<3;++i)
    if(!Convert(value.reservoir.bathMomentumChangeKgMS[i],&exchange.momentumChangeKgMS[i]))
      return {Code::NumericalFailure,"electron bath momentum exceeds the planner representation"};
  R::ProductKinematics ion;
  ion.term=0;
  ion.referenceWeight=ionReferenceWeight;
  ion.velocityMS=parent.velocityMS;
  ion.nextClock=ionClock;
  if(!Convert(states_[1].stateEnergyJ,&ion.internalEnergyJ))
    return {Code::NumericalFailure,"helium ion state energy exceeds the planner representation"};
  value.policy.products.push_back(ion);
  if(!R::MakeEventPlan(catalog_,reaction_,parent,interval,R::PlanScope::Transport,&value.policy,&value.plan,&error))
    return {Code::NumericalFailure,error};
  *output=std::move(value);
  return {};
}
}}}
