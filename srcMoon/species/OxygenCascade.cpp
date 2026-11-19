#include "OxygenCascade.h"
#include "../../src/models/chemical_database/ReducedModels.h"
#include <algorithm>
#include <cmath>
#include <locale>
#include <new>
#include <set>
#include <sstream>

namespace Moon { namespace Species { namespace OxygenCascade {
namespace DB=AMPS::ChemicalDatabase;
namespace RM=DB::ReducedModels;
namespace {
// These are source identities, not fitted calibration constants. Pinning the
// original files AND their maintained lexical projections prevents a partial
// exit list from masquerading as a complete level lifetime.
const std::string levelSha="8b0bfcbae54365ed668efa35f337e7b9ed6a60a55d0b4fbfe94693e22660c192";
const std::string linkSha="bf4e17705e69cbaebbabd07f624e6034c5d96817b10aeb51ee8a77c96a21ed4f";
const DB::ReactionData::Row levelColumns={"level_index","configuration_native","label_native","spin_multiplicity_native","L_native","J_native","observed_energy_cm_inverse_native","observed_energy_status","theoretical_energy_cm_inverse_native","theoretical_energy_status","uncertainty_status","native_line","native_row_sha256"};
const DB::ReactionData::Row linkColumns={"lower_index","upper_index","signed_wavelength_angstrom_native","wavelength_status","gf_native","decay_rate_s_inverse_native","rate_kind","uncertainty_status","label_native","native_line","native_row_sha256"};
Status error(DB::Code code,const std::string& message){return {code,message};}
bool finite(Real x){return std::isfinite(x);}
bool number(const std::string& text,Real* value){
  std::istringstream s(text);s.imbue(std::locale::classic());Real v=0;s>>v;
  if(!s||!finite(v))return false;
  s>>std::ws;
  if(!s.eof())return false;
  *value=v;return true;
}
bool integer(const std::string& text,unsigned* value){
  Real x=0;if(!number(text,&x)||x<1||x>100||std::floor(x)!=x)return false;
  *value=static_cast<unsigned>(x);return true;
}
bool valid(const Bound& b){return finite(b.absolute)&&finite(b.relative)&&b.absolute>=0&&b.relative>=0;}
bool within(Real a,Real b,const Bound& bound){
  const Real allowed=bound.absolute+bound.relative*std::max(std::fabs(a),std::fabs(b));
  return finite(a)&&finite(b)&&finite(allowed)&&std::fabs(a-b)<=allowed;
}
Status pinnedSource(const Catalogue& catalogue,const std::string& id,const std::string& sha){
  DB::SourceRecord r;auto s=catalogue.GetSource(id,&r);if(!s.ok())return s;
  if(r.sha256!=sha)return error(DB::Code::CorruptSource,"selected CHIANTI original source identity changed: "+id);
  std::string path;s=catalogue.SourcePath(id,"source-access",&path);if(!s.ok())return s;
  std::string why;if(!DB::ReactionData::VerifyFileSha256(path,sha,&why))return error(DB::Code::CorruptSource,why);
  return {};
}
}
Status Network::LoadFrozenSelectedChiantI(const Catalogue& catalogue){
  if(!catalogue.loaded())return error(DB::Code::Unavailable,"an authenticated shared Catalogue is required");
  try {
    auto s=pinnedSource(catalogue,"chianti-oi:levels-native","230bbdcd95c269fc0eea04a1fd37f534c1a42979d2e29c1b346e49dd71c688da");if(!s.ok())return s;
    s=pinnedSource(catalogue,"chianti-oi:radiative-native","f24efc2bd0fa5120baca6dce06ffb7831ae2a59e17cc186c233a43b5145866a8");if(!s.ok())return s;
    DB::TableRecord lm,tm;s=catalogue.GetTable("chianti-oi-levels-five",&lm);if(!s.ok())return s;
    s=catalogue.GetTable("chianti-oi-radiative-lowest-five",&tm);if(!s.ok())return s;
    if(lm.sha256!=levelSha||tm.sha256!=linkSha)return error(DB::Code::CorruptSource,"frozen CHIANTI selection changed; no weak exit may be removed");
    if(lm.schemaId!="amps-chianti-oi-levels-five-lexical-v1"||tm.schemaId!="amps-chianti-oi-radiative-lowest-five-lexical-v1"||lm.columns!=levelColumns||tm.columns!=linkColumns)
      return error(DB::Code::UnsupportedSchema,"selected CHIANTI column/schema mismatch");
    std::shared_ptr<const DB::ReactionData::Csv> lr,tr;
    s=catalogue.TableRows(lm.id,"table-access",&lr);if(!s.ok())return s;
    s=catalogue.TableRows(tm.id,"table-access",&tr);if(!s.ok())return s;
    if(lr->size()!=5||tr->size()!=9)return error(DB::Code::CorruptSource,"complete selected five-level/nine-link closure required");
    Network next;next.identity_=catalogue.identity();
    for(std::size_t n=0;n<lr->size();++n){
      const auto& r=(*lr)[n];Level l;
      if(r.size()!=13||!integer(r[0],&l.index)||l.index!=n+1||r[7]!="REPORTED"||r[9]!="REPORTED"||r[10]!="NOT_SUPPLIED_IN_FILE"||!number(r[6],&l.observedCmInverse)||!number(r[8],&l.theoreticalCmInverse)||l.observedCmInverse<0||l.theoreticalCmInverse<0||!DB::ReactionData::DigestSyntax(r[12]))
        return error(DB::Code::CorruptSource,"invalid/missing selected level metadata; no theoretical fallback");
      l.configuration=r[1];l.multiplicity=r[3];l.L=r[4];l.J=r[5];l.observedNative=r[6];l.theoreticalNative=r[8];l.rowSha256=r[12];next.levels_.push_back(l);
    }
    // Native indices encode the selected downward ordering. We retain native
    // rate/wavelength text as well as SI-valued rates; no negative wavelength
    // is silently turned into observed data, and no missing A is replaced.
    std::set<std::pair<unsigned,unsigned>> pairs;
    for(const auto& r:*tr){
      Link l;
      if(r.size()!=11||!integer(r[0],&l.lower)||!integer(r[1],&l.upper)||!integer(r[9],&l.nativeLine)||l.lower>=l.upper||l.upper>5||!pairs.insert({l.lower,l.upper}).second||!number(r[2],&l.signedVacuumWavelengthAngstrom)||l.signedVacuumWavelengthAngstrom==0||!number(r[5],&l.rateS)||!(l.rateS>0)||r[6]!="SINGLE_PHOTON_RADIATIVE"||r[7]!="NOT_SUPPLIED_IN_FILE"||!DB::ReactionData::DigestSyntax(r[10]))
        return error(DB::Code::CorruptSource,"invalid/missing selected radiative link; zero/missing A is not invented");
      if(!(next.levels_[l.upper-1].observedCmInverse>next.levels_[l.lower-1].observedCmInverse)||!(next.levels_[l.upper-1].theoreticalCmInverse>next.levels_[l.lower-1].theoreticalCmInverse))
        return error(DB::Code::CorruptSource,"non-downward source energies");
      l.rateNative=r[5];l.wavelengthNative=r[2];l.rowSha256=r[10];l.uncertaintyStatus=r[7];next.links_.push_back(l);
    }
    const std::set<std::pair<unsigned,unsigned>> expected={{1,2},{1,3},{2,3},{1,4},{2,4},{3,4},{1,5},{2,5},{4,5}};
    if(pairs!=expected)return error(DB::Code::CorruptSource,"selected radiative graph is not downward complete");
    // This is the same coherent A/lifetime/branching implementation exposed
    // by the shared production database. There is no second rate sum law.
    for(unsigned upper=1;upper<=5;++upper){
      auto& law=next.laws_[upper-1];std::vector<Real> rates;
      for(std::size_t i=0;i<next.links_.size();++i)if(next.links_[i].upper==upper){law.links.push_back(i);rates.push_back(next.links_[i].rateS);}
      if(rates.empty()){law.stableWithinSelectedNetwork=true;continue;}
      RM::RadiativeSet set;s=RM::RadiativeBranches(rates,&set);if(!s.ok())return s;
      law.totalRateS=set.totalRateS;law.lifetimeS=set.lifetimeS;law.branchProbabilities=set.branches;
      for(std::size_t i=0;i<law.links.size();++i)next.links_[law.links[i]].branchProbability=set.branches[i];
    }
    next.loaded_=true;*this=std::move(next);return {};
  }catch(const std::bad_alloc&){return error(DB::Code::ResourceLimit,"O I network allocation failed; reduce concurrent resources");}
}
Status Network::GetDecayLaw(unsigned level,DecayLaw* out)const{
  if(!out||level<1||level>5)return error(DB::Code::InvalidInput,"explicit native level1..5/output required");
  if(!loaded_)return error(DB::Code::Unavailable,"O I network is not loaded");
  try {
    DecayLaw temporary=laws_[level-1];*out=std::move(temporary);return {};
  }catch(const std::bad_alloc&){return error(DB::Code::ResourceLimit,"decay-law allocation failed");}
}
Status Network::SelectDecay(unsigned level,Real u,std::size_t* out)const{
  if(!out||!finite(u)||u<0||u>=1)return error(DB::Code::InvalidInput,"decay CDF requires uniform in[0,1)");
  DecayLaw law;auto s=GetDecayLaw(level,&law);if(!s.ok())return s;
  if(law.stableWithinSelectedNetwork)return error(DB::Code::OutsideDomain,"stable selected level has no decay channel");
  // Native link order defines the conditional CDF. The final positive branch
  // owns the representable upper tail; no weak channel is removed or fitted.
  Real cumulative=0;
  for(std::size_t i=0;i<law.links.size();++i){cumulative+=law.branchProbabilities[i];if(u<cumulative||i+1==law.links.size()){*out=law.links[i];return {};}}
  return error(DB::Code::NumericalFailure,"conditional radiative CDF could not be selected");
}
Status Network::WaitingTime(unsigned level,Real u,Real* out)const{
  if(!out||!finite(u)||u<=0||u>=1)return error(DB::Code::InvalidInput,"decay wait requires an open uniform");
  DecayLaw law;auto s=GetDecayLaw(level,&law);if(!s.ok())return s;
  if(!law.lifetimeS)return error(DB::Code::OutsideDomain,"stable selected level has no finite wait");
  // Exponential survival is the constant-A continuous-time Markov law, not
  // a lifetime cutoff. Open uniforms avoid log(0) and zero-duration events.
  const Real t=-std::log(u)*(*law.lifetimeS);
  if(!finite(t)||!(t>0))return error(DB::Code::NumericalFailure,"radiative wait is not representable");
  *out=t;return {};
}
Status Network::GetStateGap(std::size_t index,EnergyBasis basis,const PhotonConstants& c,StateGap* out)const{
  if(!out||index>=9||!finite(c.planckJS)||!finite(c.speedOfLightMS)||c.planckJS<=0||c.speedOfLightMS<=0||c.planckSource.empty()||c.lightSource.empty()||(basis!=EnergyBasis::Observed&&basis!=EnergyBasis::Theoretical))return error(DB::Code::InvalidInput,"explicit energy basis, positive SI h/c and source identities required");
  if(!loaded_)return error(DB::Code::Unavailable,"O I network is not loaded");
  try {
  const auto& l=links_[index];const auto& lo=levels_[l.lower-1];const auto& hi=levels_[l.upper-1];
  const bool observed=basis==EnergyBasis::Observed;
  const Real low=observed?lo.observedCmInverse:lo.theoreticalCmInverse,high=observed?hi.observedCmInverse:hi.theoreticalCmInverse;
  const Real conversion=c.planckJS*c.speedOfLightMS*100; // exact cm^-1 to m^-1 ratio
  StateGap v;v.lower=l.lower;v.upper=l.upper;v.basis=basis;v.lowerEnergyJ=low*conversion;v.upperEnergyJ=high*conversion;v.releasedStateEnergyJ=(high-low)*conversion;v.signedNativeVacuumWavelengthAngstrom=l.signedVacuumWavelengthAngstrom;
  if(!finite(conversion)||!(conversion>0)||!finite(v.lowerEnergyJ)||!finite(v.upperEnergyJ)||!finite(v.releasedStateEnergyJ)||!(v.releasedStateEnergyJ>0))return error(DB::Code::NumericalFailure,"source state-gap conversion is not representable");
  v.lowerNativeCmInverse=observed?lo.observedNative:lo.theoreticalNative;v.upperNativeCmInverse=observed?hi.observedNative:hi.theoreticalNative;v.planckJS=c.planckJS;v.speedOfLightMS=c.speedOfLightMS;v.levelTableSha256=levelSha;v.linkTableSha256=linkSha;v.planckSource=c.planckSource;v.lightSource=c.lightSource;*out=std::move(v);return {};
  }catch(const std::bad_alloc&){return error(DB::Code::ResourceLimit,"state-gap metadata allocation failed");}
}
Status Network::AdvanceRadiativeOnly(unsigned initial,Real duration,const std::vector<UniformPair>& draws,History* out)const{
  if(!out||!finite(duration)||duration<0||initial<1||initial>5)return error(DB::Code::InvalidInput,"finite nonnegative interval and initial level1..5 required");
  if(!loaded_)return error(DB::Code::Unavailable,"O I network is not loaded");
  try {
  History h;h.initialLevel=initial;h.finalLevel=initial;h.durationS=duration;Real elapsed=0;
  // A strict downward graph bounds this loop by four accepted events. The
  // first wait extending beyond the declared interval is censored; its branch
  // uniform is not consumed. No unrepresented tiny event is silently skipped.
  while(h.finalLevel!=1&&elapsed<duration){
    if(h.waitingDrawsUsed>=draws.size())return error(DB::Code::InvalidInput,"radiative history lacks explicit waiting/branch uniforms");
    const auto& u=draws[h.waitingDrawsUsed++];Real delay=0;auto s=WaitingTime(h.finalLevel,u.waiting,&delay);if(!s.ok())return s;
    if(delay>duration-elapsed)break;
    const Real next=elapsed+delay;if(!(next>elapsed)||!finite(next))return error(DB::Code::NumericalFailure,"radiative event offset cannot be represented");
    std::size_t link=0;s=SelectDecay(h.finalLevel,u.branch,&link);if(!s.ok())return s;
    ++h.branchDrawsUsed;++h.photonCounts[link];h.events.push_back({next,link});h.finalLevel=links_[link].lower;elapsed=next;
  }
  *out=std::move(h);return {};
  }catch(const std::bad_alloc&){return error(DB::Code::ResourceLimit,"radiative history allocation failed");}
}
Status SolveRestFramePhotonRecoil(const RecoilInput& in,RecoilLedger* out){
  const Real mi=in.physicalParentMassKg,mf=in.physicalDaughterMassKg,c=in.constants.speedOfLightMS,d=in.stateGap.releasedStateEnergyJ;
  if(!out||!finite(mi)||!finite(mf)||!finite(c)||!finite(d)||mi<=mf||mf<=0||c<=0||d<=0||!finite(in.constants.planckJS)||in.constants.planckJS<=0||in.constants.planckSource.empty()||in.constants.lightSource.empty()||in.constants.planckJS!=in.stateGap.planckJS||c!=in.stateGap.speedOfLightMS||in.constants.planckSource!=in.stateGap.planckSource||in.constants.lightSource!=in.stateGap.lightSource||in.parentMassSource.empty()||in.daughterMassSource.empty()||in.isotopeIdentity.empty()||in.frameIdentity.empty()||!valid(in.massEnergyConsistencyJ)||!valid(in.energyJ)||!valid(in.momentumKgMS)||!finite(in.unitDirectionTolerance)||in.unitDirectionTolerance<0||in.unitDirectionTolerance>=1||!RM::Finite(in.photonDirection))return error(DB::Code::InvalidInput,"rest-emitter recoil needs explicit physical masses/state/units/sources/arithmetic bounds");
  try {
  // DeltaE is an internal-state/rest-energy difference, counted exactly once.
  // Photon plus relativistic atom recoil must close BOTH the supplied level
  // gap and independent physical mass gap within caller-frozen arithmetic
  // bounds. Such bounds are not permission to invent a physical mass profile.
  const Real rest=(mi*c)*c,massGap=((mi-mf)*c)*c;
  if(!finite(rest)||!finite(massGap)||!(rest>0)||!(massGap>0)||d>=rest)return error(DB::Code::OutsideDomain,"rest mass/state release is outside the photon-emission domain");
  if(!within(massGap,d,in.massEnergyConsistencyJ))return error(DB::Code::OutsideDomain,"physical mass difference is inconsistent with supplied state release");
  const Real norm=std::hypot(in.photonDirection[0],in.photonDirection[1],in.photonDirection[2]);
  if(!finite(norm)||norm==0||std::fabs(norm-1)>in.unitDirectionTolerance)return error(DB::Code::InvalidInput,"nonzero unit photon direction required within explicit bound");
  RecoilLedger v;v.stateGap=in.stateGap;v.constants=in.constants;v.physicalParentMassKg=mi;v.physicalDaughterMassKg=mf;v.stateEnergyReleaseJ=d;v.massEnergyReleaseJ=massGap;
  // Rest-frame two-body conservation: E_gamma=(Mi^2-Mf^2)c^2/(2Mi).
  // The DeltaE form avoids subtracting nearly equal squared atomic masses.
  // Ratio-first arithmetic avoids overflow of 2*Mi*c^2.
  v.photonEnergyJ=d*(1-(d/rest)/2);v.photonMomentumMagnitudeKgMS=v.photonEnergyJ/c;
  const Real gammaMass=std::hypot(mf,v.photonMomentumMagnitudeKgMS/c),denominator=gammaMass+mf;
  // Rationalized relativistic KE avoids gammaMass-minus-Mf cancellation.
  v.atomicKineticEnergyJ=(v.photonMomentumMagnitudeKgMS/denominator)*v.photonMomentumMagnitudeKgMS;
  v.energyResidualJ=RM::Sum3(v.photonEnergyJ,v.atomicKineticEnergyJ,-d);
  v.massEnergyResidualJ=RM::Sum3(v.photonEnergyJ,v.atomicKineticEnergyJ,-massGap);
  const Real energyAllowed=in.energyJ.absolute+in.energyJ.relative*std::max(std::fabs(d),std::fabs(massGap));
  if(!finite(energyAllowed)||!finite(v.energyResidualJ)||!finite(v.massEnergyResidualJ)||std::fabs(v.energyResidualJ)>energyAllowed||std::fabs(v.massEnergyResidualJ)>energyAllowed)return error(DB::Code::NumericalFailure,"recoil ledger exceeds explicit state/rest-mass energy bounds");
  if(!finite(v.photonEnergyJ)||!(v.photonEnergyJ>0)||!finite(v.photonMomentumMagnitudeKgMS)||!(v.photonMomentumMagnitudeKgMS>0)||!finite(gammaMass)||!finite(denominator)||!finite(v.atomicKineticEnergyJ)||!(v.atomicKineticEnergyJ>0))return error(DB::Code::NumericalFailure,"finite-emitter photon/recoil energy is not represented within explicit bounds");
  for(unsigned i=0;i<3;++i){
    const Real unit=in.photonDirection[i]/norm;
    v.photonMomentumKgMS[i]=v.photonMomentumMagnitudeKgMS*unit;v.atomicMomentumKgMS[i]=-v.photonMomentumKgMS[i];v.atomicVelocityMS[i]=v.atomicMomentumKgMS[i]/gammaMass;
    if(!finite(v.atomicVelocityMS[i])||(unit!=0&&v.photonMomentumKgMS[i]==0)||!within(v.photonMomentumKgMS[i]+v.atomicMomentumKgMS[i],0,in.momentumKgMS))return error(DB::Code::NumericalFailure,"photon/recoil vector is outside explicit representation bounds");
  }
  // The represented Cartesian vector is part of the actual ledger. A scalar
  // E/c check alone cannot certify it: normalizing a non-axis direction and
  // multiplying its components introduces another rounding operation. Check
  // its norm, massless dispersion and recoil energy against the SAME explicit
  // caller bounds. Zero bounds therefore refuse a non-exact representation;
  // they never acquire an implicit epsilon.
  const Real representedP=std::hypot(v.photonMomentumKgMS[0],v.photonMomentumKgMS[1],v.photonMomentumKgMS[2]);
  const Real representedGammaMass=std::hypot(mf,representedP/c);
  const Real representedKE=(representedP/(representedGammaMass+mf))*representedP;
  const Real photonDispersionResidual=RM::Sum3(v.photonEnergyJ,-representedP*c,0);
  if(!finite(representedP)||!finite(representedGammaMass)||!finite(representedKE)||!finite(photonDispersionResidual)||!within(representedP,v.photonMomentumMagnitudeKgMS,in.momentumKgMS)||std::fabs(photonDispersionResidual)>energyAllowed||std::fabs(RM::Sum3(representedKE,-v.atomicKineticEnergyJ,0))>energyAllowed)
    return error(DB::Code::NumericalFailure,"represented photon vector or atomic recoil is outside explicit dispersion/energy bounds");
  v.atomicKineticEnergyJ=representedKE;
  v.energyResidualJ=RM::Sum3(v.photonEnergyJ,representedKE,-d);
  v.massEnergyResidualJ=RM::Sum3(v.photonEnergyJ,representedKE,-massGap);
  if(std::fabs(v.energyResidualJ)>energyAllowed||std::fabs(v.massEnergyResidualJ)>energyAllowed)
    return error(DB::Code::NumericalFailure,"represented vector energy ledger does not close within explicit bounds");
  for(unsigned i=0;i<3;++i){
    v.atomicVelocityMS[i]=v.atomicMomentumKgMS[i]/representedGammaMass;
    if(!finite(v.atomicVelocityMS[i])||!within(v.atomicVelocityMS[i]*representedGammaMass,v.atomicMomentumKgMS[i],in.momentumKgMS))
      return error(DB::Code::NumericalFailure,"represented atomic velocity does not reconstruct recoil momentum within explicit bounds");
  }
  if(!(std::hypot(v.atomicVelocityMS[0],v.atomicVelocityMS[1],v.atomicVelocityMS[2])<c))
    return error(DB::Code::NumericalFailure,"represented recoiling atom is not subluminal");
  v.isotopeIdentity=in.isotopeIdentity;v.frameIdentity=in.frameIdentity;v.parentMassSource=in.parentMassSource;v.daughterMassSource=in.daughterMassSource;*out=std::move(v);return {};
  }catch(const std::bad_alloc&){return error(DB::Code::ResourceLimit,"recoil metadata allocation failed");}
}
}}}
