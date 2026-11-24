#include "SodiumGroundCapture.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <locale>
#include <sstream>
#include <utility>

namespace Moon { namespace Species { namespace SodiumGroundCapture {
const char ModelId[]="hplus-na3s-ground-capture-4kev-u-isotropic-v1";
namespace {
using Vector=K::Vector3;
bool Positive(Real value){return std::isfinite(value)&&value>0;}
Real Norm(const Vector& v){return std::hypot(v[0],v[1],v[2]);}
Vector Cross(const Vector& a,const Vector& b){
  return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
template<class Bound> bool BoundValid(const Bound& b){
  return std::isfinite(b.absolute)&&std::isfinite(b.relative)&&b.absolute>=0&&b.relative>=0&&b.relative<1;
}
DB::Status Fail(DB::Code code,const std::string& message){return {code,message};}
const char PointDigest[]="59a2d1d7e5a9342dcbd8ba1a8251c1864650ae9dc42624dca42e8befde69029a";
}

// The case consumes a published partial for the PROMPT ground products.
// Never query the D119990 total outer-shell rate or assume inheritance of
// Poppe's bulk pickup-ion birth velocity. Neither has this state's meaning.
DB::Status Provider::Prepare(const DB::Catalogue& database,const std::string& path,
    const Definition& d,Provider* output){
  if(!output||!database.loaded())return Fail(DB::Code::InvalidInput,"an open catalogue and output are required");
  if(database.identity().releaseSha256!="26876e4c09e251220674a65bd43efa16317d1c354a7550b9a06667254637a86c"||
      database.identity().exportSha256!="585edadbf8e4df35be9c3db1a33583074932df14cea8811e44832535a96135e9")
    return Fail(DB::Code::Unavailable,"select the frozen r14 source-bound scenario for sodium ground capture");
  std::shared_ptr<const DB::ReactionData::Csv> selections;
  auto status=database.TableRows("amps-model-source-bindings-v1","table-access",&selections);
  if(!status.ok())return status;
  if(std::count_if(selections->begin(),selections->end(),[](const DB::ReactionData::Row& r){
    return r.size()==2&&r[0]==ModelId;})!=1)return Fail(DB::Code::CorruptSource,"missing or ambiguous r14 model binding");
  DB::ReactionData::Csv points;std::string error;
  if(!DB::ReactionData::ReadPinnedCsv(path,PointDigest,290,
      {"model","source","energy","energy_unit","energy_frame","cross_section","area_unit","angular_model","recoil_model"},
      &points,&error))return Fail(DB::Code::CorruptSource,error);
  if(points.size()!=1||points[0][0]!=ModelId||points[0][3]!="keV/u"||points[0][6]!="cm^2")
    return Fail(DB::Code::CorruptSource,"r14 projection has wrong model, energy axis or area unit");
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
  // r14 pins the complete lexical record. Convert using the native decimal
  // tokens after authentication. The keV/u and cm² scale factors are exact.
  Provider candidate;candidate.definition_=d;
  std::istringstream energy(points[0][2]),area(points[0][5]);
  energy.imbue(std::locale::classic());area.imbue(std::locale::classic());
  Real nativeEnergy=0,nativeArea=0;
  if(!(energy>>nativeEnergy)||!energy.eof()||!(area>>nativeArea)||!area.eof()||
     !Positive(nativeEnergy)||!Positive(nativeArea))
    return Fail(DB::Code::CorruptSource,"ground-capture point contains an invalid native energy or cross section");
  candidate.energyEvPerU_=1000*nativeEnergy;
  candidate.crossSectionM2_=1e-4L*nativeArea;
  const char* names[]={"Na","H+","Na+","H"};
  const char* levels[]={"3s","bare-proton","2p6","1s"};
  const char* nuclei[]={"23Na","1H","23Na","1H"};
  const int charges[]={0,1,1,0},slots[]={d.sodiumSpecies,-1,d.sodiumIonSpecies,d.hydrogenSpecies};
  for(std::size_t i=0;i<4;++i){
    const auto& state=d.states[i];
    if(!state.valid())return Fail(DB::Code::Unavailable,"missing source-resolved capture state");
    const auto& s=state.definition();
    if(s.identity!=names[i]||s.state!=levels[i]||s.nuclei.size()!=1||s.nuclei[0].identity!=nuclei[i]||
       s.nuclei[0].count!=1||s.elementaryCharge!=charges[i]||s.excitation.joules!=0||
       (charges[i]? !Positive(s.formation.joules):s.formation.joules!=0)||
       state.referenceInventoryIdentity()!=d.states[0].referenceInventoryIdentity())
      return Fail(DB::Code::InvalidInput,"wrong isotope/state or incompatible separated-atom energy/mass convention");
    candidate.catalog_.push_back({s.identity,s.state,slots[i],s.nuclei,s.elementaryCharge,static_cast<double>(state.referenceMassKg())});
    candidate.internal_[i]=static_cast<double>(s.formation.joules);
  }
  // The additive mass convention preserves integer material, while Q follows
  // from the ion formation energies. Adding an ion mass defect to Q again
  // would count the same electronic binding energy twice.
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
  const Real energy=(0.5L*definition_.atomicMassUnitKg/definition_.electronVoltJ)*speed*speed;
  // Relative energy is specified per atomic mass unit, independent of the
  // projectile isotope mass used below for true momentum. Allow only binary64
  // coordinate roundoff on this single node, never a finite physical bin.
  if(!std::isfinite(energy)||std::abs(energy-energyEvPerU_)>32*std::numeric_limits<double>::epsilon()*energyEvPerU_)
    return Fail(DB::Code::OutsideDomain,"H(1s) sodium theory partial is defined only at 4.0 keV/u in the target rest frame");
  const Real rate=beam.protonDensityM3*crossSectionM2_*speed;
  if(!std::isfinite(rate)||(beam.protonDensityM3>0&&rate==0))
    return Fail(DB::Code::NumericalFailure,"beam event rate overflow/underflow");
  *output={energy,speed,crossSectionM2_,rate};return {};
}

DB::Status Provider::Compose(const Beam& beam,const R::ParentState& parent,const R::TraceInterval& interval,
    Real polar,Real azimuth,const std::array<double,2>& weights,const std::array<R::OpticalDepthClock,2>& clocks,
    Composition* output) const{
  if(!output||catalog_.empty()||parent.species!=0||parent.internalEnergyJ!=internal_[0])
    return Fail(DB::Code::InvalidInput,"capture requires a prepared ground-state 23Na parent and output");
  Composition c;std::string error;auto status=Evaluate(beam,parent.velocityMS,&c.evaluation);
  if(!status.ok())return status;
  if(!Positive(c.evaluation.rateS))return Fail(DB::Code::OutsideDomain,"no collision can occur in a zero-density beam");
  Vector incident,local;
  for(unsigned k=0;k<3;++k)incident[k]=(beam.protonVelocityMS[k]-parent.velocityMS[k])/c.evaluation.relativeSpeedMS;
  if(!K::IsotropicDirectionFromOpenUniforms(polar,azimuth,&local,&error))return Fail(DB::Code::InvalidInput,error);
  // Define the polar axis by the incident relative velocity and both azimuthal
  // axes deterministically. Choosing the least-aligned Cartesian axis avoids
  // ill-conditioning for beams parallel to a coordinate axis. theta describes
  // the outgoing H projectile; the first tracked product, Na+, recoils oppositely.
  Vector axis{};unsigned least=0;
  for(unsigned k=1;k<3;++k)if(std::abs(incident[k])<std::abs(incident[least]))least=k;
  axis[least]=1;
  Vector tangent=Cross(axis,incident);const Real length=Norm(tangent);
  for(auto& x:tangent)x/=length;
  const Vector bitangent=Cross(incident,tangent);
  const auto& s=definition_.states;
  const Real mNa=s[0].referenceMassKg(),mp=s[1].referenceMassKg();
  const Real mu=mp*(mNa/(mp+mNa));
  K::TwoBodyInput solve;solve.products={{s[2],s[3]}};solve.controls=definition_.numerical;
  solve.context={definition_.frame,ModelId};
  // Separate incoming relative motion from uniform transport of the COM.
  // This avoids subtracting large bulk kinetic energies in a moving frame.
  solve.availableCOMKineticJ=.5L*mu*c.evaluation.relativeSpeedMS*c.evaluation.relativeSpeedMS+formation_.qJ;
  for(unsigned k=0;k<3;++k){
    solve.totalMomentumKgMS[k]=mp*beam.protonVelocityMS[k]+mNa*parent.velocityMS[k];
    solve.relativeDirection[k]=-(local[0]*tangent[k]+local[1]*bitangent[k]+local[2]*incident[k]);
  }
  if(!K::GeneralTwoBodyFinalState(solve,&c.products,&error))return Fail(DB::Code::NumericalFailure,error);
  c.products.formationBalance=formation_;c.products.stateQJ=formation_.qJ;c.products.modelIdentity=ModelId;
  c.products.sourceIdentities.push_back(PointDigest);
  auto& policy=c.policy;policy.qualified=true; // scoped authority for this modeled development event
  policy.identity=ModelId;policy.evidenceIdentity=std::string("r14 theory + approved isotropic sensitivity; ")+PointDigest;
  policy.eventFrame=definition_.frame;policy.maximumClassicalSpeedMS=definition_.numerical.maximumClassicalSpeedMS;
  policy.bounds=definition_.ledgerBounds;
  // The bath loses exactly the incoming proton. Derive its moments from the
  // prescribed projectile itself; a residual-fitted bath would hide errors in
  // tracked product mass, velocity or internal excitation.
  policy.bath.supplied=true;policy.bath.modelMassChangeKg=static_cast<double>(-mp);
  const Real vp=Norm(beam.protonVelocityMS);
  policy.bath.energyChangeJ=static_cast<double>(-(.5L*mp*vp*vp+s[1].definition().formation.joules));
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
