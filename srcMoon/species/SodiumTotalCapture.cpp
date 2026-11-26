#include "SodiumTotalCapture.h"
#if __has_include("../../models/chemical_database/DataIO.h")
#include "../../models/chemical_database/DataIO.h"
#else
#include "../../src/models/chemical_database/DataIO.h"
#endif
#include <algorithm>
#include <cmath>
#include <limits>
#include <locale>
#include <sstream>
#include <utility>

namespace Moon { namespace Species { namespace SodiumTotalCapture {
const char ModelId[]="hplus-na3s-total-capture-forward-1to5kev-u-v1";
namespace {
using Vector=K::Vector3;
namespace CX=ChargeExchange;
// Exact bytes of the reviewed prepare_total_capture.py outputs (2026-11-25).
const char ShingalDigest[]="462760edaa5e13fea4b489a1cb81498167e7042d4101b9a44207142928c8d89a";
const char LevelsDigest[]="1d9f99878f08a25940d5f1fbc72e953a69b913e703f6a800c4565167fdd12199";
constexpr std::size_t ShingalBytes=262,LevelsBytes=528;
bool Positive(Real value){return std::isfinite(value)&&value>0;}
Real Norm(const Vector& v){return std::hypot(v[0],v[1],v[2]);}
template<class Bound> bool BoundValid(const Bound& b){
  return std::isfinite(b.absolute)&&std::isfinite(b.relative)&&b.absolute>=0&&b.relative>=0&&b.relative<1;
}
DB::Status Fail(DB::Code code,const std::string& message){return {code,message};}
// Locale-independent conversion of one complete native decimal token.
bool Decimal(const std::string& token,Real* value){
  std::istringstream stream(token);stream.imbue(std::locale::classic());
  Real parsed=0;
  if(!(stream>>parsed)||!stream.eof()||!std::isfinite(parsed))return false;
  *value=parsed;return true;
}
// Statistical weight 2J+1 from the native half-integer token "k/2".
bool WeightFromJ(const std::string& token,long* weight){
  const auto slash=token.find('/');
  if(slash==std::string::npos||token.substr(slash)!="/2")return false;
  try{const long twiceJ=std::stol(token.substr(0,slash));if(twiceJ<1)return false;*weight=twiceJ+1;}
  catch(...){return false;}
  return true;
}
}

DB::Status Provider::Prepare(const Definition& d,Provider* output){
  if(!output)return Fail(DB::Code::InvalidInput,"output storage is required");
  const auto& n=d.numerical;
  if(d.frame.empty()||d.constantsIdentity.empty()||!Positive(d.atomicMassUnitKg)||d.electronVoltJ!=1.602176634e-19L||
     !Positive(n.maximumClassicalSpeedMS)||n.maximumClassicalSpeedMS>=K::SpeedOfLightMS||
     !std::isfinite(n.unitVectorTolerance)||n.unitVectorTolerance<0||n.unitVectorTolerance>=1||
     !BoundValid(n.referenceMassKg)||!BoundValid(n.momentumKgMS)||!BoundValid(n.energyJ)||!BoundValid(n.thresholdJ)||
     !BoundValid(d.ledgerBounds.modelMassKg)||!BoundValid(d.ledgerBounds.momentumKgMS)||!BoundValid(d.ledgerBounds.energyJ))
    return Fail(DB::Code::InvalidInput,"explicit SI constants, common frame, classical and arithmetic controls required");
  if(d.sodiumSpecies<0||d.sodiumIonSpecies<0||d.hydrogenSpecies<0||
     d.sodiumSpecies==d.sodiumIonSpecies||d.sodiumSpecies==d.hydrogenSpecies||d.sodiumIonSpecies==d.hydrogenSpecies)
    return Fail(DB::Code::InvalidInput,"Na, Na+ and H require three distinct compiled species slots");
  Provider candidate;candidate.definition_=d;std::string error;
  // Morgan total: the shared registry authenticates both v1 publications.
  if(!candidate.morgan_.LoadPublishedReferences(CX::FrozenV1References(d.morganPublicationDirectory),&error))
    return Fail(DB::Code::CorruptSource,"Morgan 1985 total-capture publication: "+error);
  CX::ChannelMetadata morgan;
  if(!candidate.morgan_.GetMetadata(CX::ChannelId::ProtonNaMorgan1985LogLog,&morgan,&error))
    return Fail(DB::Code::Unavailable,error);
  if(morgan.minimumEnergyEvPerU>MinimumEnergyEvPerU||morgan.maximumEnergyEvPerU<MaximumEnergyEvPerU)
    return Fail(DB::Code::CorruptSource,"approved 1-5 keV/u envelope lies outside the Morgan source domain");
  // Shingal Table 4: E [keV/u], sigma_1s, 2s, 2p, 3s, 3p, 3d [A^2].
  DB::ReactionData::Csv table;
  if(!DB::ReactionData::ReadPinnedCsv(d.shingalTablePath,ShingalDigest,ShingalBytes,
      {"energy_kev_u","sigma_1s_A2","sigma_2s_A2","sigma_2p_A2","sigma_3s_A2","sigma_3p_A2","sigma_3d_A2"},
      &table,&error))return Fail(DB::Code::CorruptSource,"Shingal Table 4 projection: "+error);
  for(const auto& row:table){
    Real cells[7];
    for(std::size_t i=0;i<7;++i)
      if(!Decimal(row[i],&cells[i])||!Positive(cells[i]))
        return Fail(DB::Code::CorruptSource,"Shingal Table 4 contains a nonpositive or invalid value");
    Row value;value.energyEvPerU=1000*cells[0];
    // Shell partials: n=1 (1s), n=2 (2s+2p), n=3 (3s+3p+3d); A^2 = 1e-20 m^2.
    value.sigmaM2={{1e-20L*cells[1],1e-20L*(cells[2]+cells[3]),1e-20L*(cells[4]+cells[5]+cells[6])}};
    if(!candidate.shingal_.empty()&&value.energyEvPerU<=candidate.shingal_.back().energyEvPerU)
      return Fail(DB::Code::CorruptSource,"Shingal Table 4 energies are not strictly increasing");
    candidate.shingal_.push_back(value);
  }
  if(candidate.shingal_.empty()||candidate.shingal_.front().energyEvPerU>MinimumEnergyEvPerU||
     candidate.shingal_.back().energyEvPerU<MaximumEnergyEvPerU)
    return Fail(DB::Code::CorruptSource,"Shingal rows do not bracket the approved 1-5 keV/u envelope");
  // NIST ASD H I levels: shell excitation = sum(g E) / sum(g) over nlj levels.
  DB::ReactionData::Csv levels;
  if(!DB::ReactionData::ReadPinnedCsv(d.hydrogenLevelsPath,LevelsDigest,LevelsBytes,
      {"n","l","configuration","term","J","g","level_eV","uncertainty_eV","reference"},&levels,&error))
    return Fail(DB::Code::CorruptSource,"NIST H I level projection: "+error);
  std::array<Real,Shells> weighted{},weights{};
  for(const auto& row:levels){
    long g=0,parsed=0;Real level=0;
    const int shell=row[0].size()==1?row[0][0]-'0':0;
    if(shell<1||shell>Shells||!WeightFromJ(row[4],&g)||!Decimal(row[6],&level)||level<0)
      return Fail(DB::Code::CorruptSource,"invalid NIST H I level row");
    try{parsed=std::stol(row[5]);}catch(...){parsed=0;}
    if(parsed!=g)return Fail(DB::Code::CorruptSource,"NIST statistical weight differs from 2J+1");
    weighted[shell-1]+=g*level;weights[shell-1]+=g;
  }
  for(int s=0;s<Shells;++s){
    if(!(weights[s]>0))return Fail(DB::Code::CorruptSource,"a hydrogen shell has no NIST level");
    candidate.excitationJ_[s]=(weighted[s]/weights[s])*d.electronVoltJ;
  }
  if(candidate.excitationJ_[0]!=0||!(candidate.excitationJ_[1]>0)||!(candidate.excitationJ_[2]>candidate.excitationJ_[1]))
    return Fail(DB::Code::CorruptSource,"hydrogen ground level must be zero and shells increasing");
  // States: identical convention and checks as the 4 keV/u ground-capture case.
  const char* names[]={"Na","H+","Na+","H"};
  const char* states[]={"3s","bare-proton","2p6","1s"};
  const char* nuclei[]={"23Na","1H","23Na","1H"};
  const int charges[]={0,1,1,0},slots[]={d.sodiumSpecies,-1,d.sodiumIonSpecies,d.hydrogenSpecies};
  for(std::size_t i=0;i<4;++i){
    const auto& state=d.states[i];
    if(!state.valid())return Fail(DB::Code::Unavailable,"missing source-resolved capture state");
    const auto& s=state.definition();
    if(s.identity!=names[i]||s.state!=states[i]||s.nuclei.size()!=1||s.nuclei[0].identity!=nuclei[i]||
       s.nuclei[0].count!=1||s.elementaryCharge!=charges[i]||s.excitation.joules!=0||
       (charges[i]? !Positive(s.formation.joules):s.formation.joules!=0)||
       state.referenceInventoryIdentity()!=d.states[0].referenceInventoryIdentity())
      return Fail(DB::Code::InvalidInput,"wrong isotope/state or incompatible separated-atom energy/mass convention");
    candidate.catalog_.push_back({s.identity,s.state,slots[i],s.nuclei,s.elementaryCharge,static_cast<double>(state.referenceMassKg())});
    candidate.internal_[i]=static_cast<double>(s.formation.joules);
  }
  if(!K::FormationQ({d.states[0],d.states[1]},{d.states[2],d.states[3]},&candidate.formation_,&error))
    return Fail(DB::Code::InvalidInput,error);
  candidate.reaction_={ModelId,{{0,1,R::Population::Tracked},{1,1,R::Population::Bath}},
      {{2,1,R::Population::Tracked},{3,1,R::Population::Tracked}}};
  if(!R::ValidateReaction(candidate.catalog_,candidate.reaction_,&error))return Fail(DB::Code::InvalidInput,error);
  *output=std::move(candidate);return {};
}

DB::Status Provider::Evaluate(const Beam& beam,const R::Vector3& sodium,Evaluation* output) const{
  if(!output||catalog_.empty()||!std::isfinite(beam.protonDensityM3)||beam.protonDensityM3<0)
    return Fail(DB::Code::InvalidInput,"prepared model, output and finite nonnegative proton density [m^-3] required");
  Vector relative;
  for(unsigned k=0;k<3;++k)relative[k]=beam.protonVelocityMS[k]-sodium[k];
  const Real speed=Norm(relative);
  if(!Positive(speed)||!std::isfinite(Norm(beam.protonVelocityMS))||
      Norm(beam.protonVelocityMS)>definition_.numerical.maximumClassicalSpeedMS||
      std::hypot(sodium[0],sodium[1],sodium[2])>definition_.numerical.maximumClassicalSpeedMS)
    return Fail(DB::Code::OutsideDomain,"invalid velocity or beyond the declared classical speed domain");
  // Energy per atomic mass unit in the target rest frame [eV/u].
  const Real energy=(0.5L*definition_.atomicMassUnitKg/definition_.electronVoltJ)*speed*speed;
  if(!std::isfinite(energy)||energy<MinimumEnergyEvPerU||energy>MaximumEnergyEvPerU)
    return Fail(DB::Code::OutsideDomain,"total capture is admitted only for 1000-5000 eV/u in the target rest frame");
  Evaluation value;value.energyEvPerU=energy;value.relativeSpeedMS=speed;
  double sigma=0;std::string error;
  if(!morgan_.CrossSectionM2(CX::ChannelId::ProtonNaMorgan1985LogLog,static_cast<double>(energy),&sigma,&error))
    return Fail(DB::Code::OutsideDomain,error);
  value.totalCrossSectionM2=sigma;
  // Shell fractions: log-log interpolation of each shell partial between the
  // adjacent Shingal rows, then normalization. A row energy uses that row.
  const auto upper=std::lower_bound(shingal_.begin(),shingal_.end(),energy,
      [](const Row& row,Real e){return row.energyEvPerU<e;});
  if(upper==shingal_.end()||(upper==shingal_.begin()&&upper->energyEvPerU!=energy))
    return Fail(DB::Code::OutsideDomain,"energy outside the Shingal rows; no extrapolation");
  std::array<Real,Shells> partial{};
  if(upper->energyEvPerU==energy)partial=upper->sigmaM2;
  else{
    const Row& lo=*(upper-1);const Row& hi=*upper;
    const Real t=std::log(energy/lo.energyEvPerU)/std::log(hi.energyEvPerU/lo.energyEvPerU);
    for(int s=0;s<Shells;++s)
      partial[s]=std::exp(std::log(lo.sigmaM2[s])+t*(std::log(hi.sigmaM2[s])-std::log(lo.sigmaM2[s])));
  }
  const Real sum=partial[0]+partial[1]+partial[2];
  if(!Positive(sum))return Fail(DB::Code::NumericalFailure,"shell partials are not representable");
  for(int s=0;s<Shells;++s)value.shellFraction[s]=partial[s]/sum;
  value.rateS=beam.protonDensityM3*value.totalCrossSectionM2*speed;
  if(!std::isfinite(value.rateS)||(beam.protonDensityM3>0&&value.rateS==0))
    return Fail(DB::Code::NumericalFailure,"beam event rate overflow/underflow");
  *output=value;return {};
}

DB::Status Provider::Compose(const Beam& beam,const R::ParentState& parent,const R::TraceInterval& interval,
    Real uniform,const std::array<double,2>& weights,const std::array<R::OpticalDepthClock,2>& clocks,
    Composition* output) const{
  if(!output||catalog_.empty()||parent.species!=0||parent.internalEnergyJ!=internal_[0])
    return Fail(DB::Code::InvalidInput,"capture requires a prepared ground-state 23Na parent and output");
  if(!(uniform>0&&uniform<1))
    return Fail(DB::Code::InvalidInput,"shell selection requires a strictly open uniform in (0,1)");
  Composition c;std::string error;auto status=Evaluate(beam,parent.velocityMS,&c.evaluation);
  if(!status.ok())return status;
  if(!Positive(c.evaluation.rateS))return Fail(DB::Code::OutsideDomain,"no collision can occur in a zero-density beam");
  const auto& f=c.evaluation.shellFraction;
  c.shell=uniform<f[0]?1:(uniform<f[0]+f[1]?2:3);
  const Real excitation=excitationJ_[c.shell-1];
  Vector incident;
  for(unsigned k=0;k<3;++k)incident[k]=(beam.protonVelocityMS[k]-parent.velocityMS[k])/c.evaluation.relativeSpeedMS;
  const auto& s=definition_.states;
  const Real mNa=s[0].referenceMassKg(),mp=s[1].referenceMassKg();
  const Real mIon=s[2].referenceMassKg(),mH=s[3].referenceMassKg();
  const Real mu=mp*(mNa/(mp+mNa));
  K::TwoBodyInput solve;solve.products={{s[2],s[3]}};solve.controls=definition_.numerical;
  solve.context={definition_.frame,ModelId};
  // Prompt H(n): Q(n) = IE(H) - IE(Na) - E_exc(n). The relative motion is
  // separated from uniform COM transport, as in the ground-capture case.
  c.shellQJ=formation_.qJ-excitation;
  solve.availableCOMKineticJ=.5L*mu*c.evaluation.relativeSpeedMS*c.evaluation.relativeSpeedMS+c.shellQJ;
  for(unsigned k=0;k<3;++k){
    solve.totalMomentumKgMS[k]=mp*beam.protonVelocityMS[k]+mNa*parent.velocityMS[k];
    // Zero-angle exchange: product 1 (H) moves along +incident in the COM.
    solve.relativeDirection[k]=-incident[k];
  }
  if(!K::GeneralTwoBodyFinalState(solve,&c.products,&error))return Fail(DB::Code::NumericalFailure,error);
  c.products.formationBalance=formation_;c.products.stateQJ=c.shellQJ;c.products.modelIdentity=ModelId;
  c.products.sourceIdentities={ShingalDigest,LevelsDigest,
      CX::FrozenV1References("").sodiumPublicationSha256};
  // Omitted-physics bounds for the approved model scope.
  const Real muFinal=mIon*(mH/(mIon+mH));
  c.finalComMomentumKgMS=std::sqrt(2*muFinal*solve.availableCOMKineticJ);
  const Real theta=MaximumComAngleDeg*(std::acos(-1.0L)/180);
  c.omittedImpulseBoundKgMS=2*c.finalComMomentumKgMS*std::sin(theta/2);
  c.sodiumIonRecoilBoundJ=c.omittedImpulseBoundKgMS*c.omittedImpulseBoundKgMS/(2*mIon);
  c.radiatedEnergyJ=excitation;
  c.photonRecoilSpeedBoundMS=excitation/(mH*K::SpeedOfLightMS);
  auto& policy=c.policy;policy.qualified=true; // scoped authority: approved 2026-11-25 I24 model
  policy.identity=ModelId;
  policy.evidenceIdentity=std::string("Morgan1985 total + Shingal1986 shell split + NIST ASD centroids; "
      "zero-angle COM, theta_max 0.4 deg; ")+ShingalDigest+";"+LevelsDigest;
  policy.eventFrame=definition_.frame;policy.maximumClassicalSpeedMS=definition_.numerical.maximumClassicalSpeedMS;
  policy.bounds=definition_.ledgerBounds;
  // Untracked environment per event: it loses the incoming proton and gains
  // the radiated excitation energy. Both moments come from the prescribed
  // projectile and the selected shell, never from a residual fit.
  policy.bath.supplied=true;policy.bath.modelMassChangeKg=static_cast<double>(-mp);
  const Real vp=Norm(beam.protonVelocityMS);
  policy.bath.energyChangeJ=static_cast<double>(-(.5L*mp*vp*vp+s[1].definition().formation.joules)+excitation);
  for(unsigned k=0;k<3;++k)policy.bath.momentumChangeKgMS[k]=static_cast<double>(-mp*beam.protonVelocityMS[k]);
  for(unsigned i=0;i<2;++i){
    R::ProductKinematics p;p.term=i;p.referenceWeight=weights[i];p.nextClock=clocks[i];p.internalEnergyJ=internal_[i+2];
    for(unsigned k=0;k<3;++k)p.velocityMS[k]=static_cast<double>(c.products.products[i].velocityMS[k]);
    policy.products.push_back(p);
  }
  if(!R::MakeEventPlan(catalog_,reaction_,parent,interval,R::PlanScope::Transport,&policy,&c.plan,&error))
    return Fail(DB::Code::InvalidInput,error);
  *output=std::move(c);return {};
}
}}}
