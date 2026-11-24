// Production primitive verification only. All masses/labels/velocities below
// are explicitly analytical mechanical controls, not lunar species, measured
// cross sections, daughter distributions, observations or kinetic authority.
#include "srcMoon/species/ReactionExecution.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <type_traits>

namespace R=Moon::Species::ReactionExecution;
static_assert(std::is_trivially_copyable<R::TraceInterval>::value,"local interval must support PIC memcpy storage");
static_assert(std::is_trivially_copyable<R::OpticalDepthClock>::value,"clock must support PIC memcpy storage");
static_assert(std::is_nothrow_move_assignable<R::EventPlan>::value,"plan publication must not partially modify an output on allocation failure");
int checks=0,failures=0;
void Check(bool condition,const char* name) {
  ++checks;if (!condition) ++failures;
  std::cout<<(condition ? "PASS ":"FAIL ")<<name<<'\n';
}
bool Near(double a,double b,double relative=1e-13) {
  return std::abs(a-b)<=relative*std::max(1.0,std::max(std::abs(a),std::abs(b)));
}
R::OpticalDepthClock Fresh(double u) {
  R::OpticalDepthClock clock;std::string e;
  if (!R::OpticalDepthClock::FromOpenUniform(u,&clock,&e)) throw std::runtime_error(e);
  return clock;
}
R::TraceInterval Interval(double start,double duration,double elapsed=0) {
  R::TraceInterval value;std::string e;
  if (!R::TraceInterval::Create(start,duration,elapsed,&value,&e)) throw std::runtime_error(e);
  return value;
}
struct Fixture {
  // Exact integer composition and convenient model kg values make independent
  // identities transparent; they do not claim an atomic/molecular mass model.
  std::vector<R::SpeciesIdentity> catalog{
    {"control-dimer","ground",0,{{"X",1},{"Y",1}},0,2},
    {"control-X","ground",1,{{"X",1}},0,1},
    {"control-Y","ground",2,{{"Y",1}},0,1}};
  R::Reaction reaction{"mechanical-control",{{0,1,R::Population::Tracked}},
      {{1,1,R::Population::Tracked},{2,1,R::Population::Tracked}}};
  R::ParentState parent;
  R::KinematicPolicy policy;
  Fixture() {
    parent.species=0;parent.referenceWeight=2.5;parent.individualCorrection=.5;
    parent.positionM={{1,2,3}};
    policy.qualified=true;policy.identity="analytical-mechanical-control";
    policy.evidenceIdentity="CONTROL-NOT-LUNAR-KINETICS";policy.eventFrame="declared-Cartesian-control";
    policy.maximumClassicalSpeedMS=10;
    policy.bounds.modelMassKg={1e-14,1e-14};policy.bounds.momentumKgMS={1e-14,1e-14};policy.bounds.energyJ={1e-14,1e-14};
    policy.bath.supplied=true;policy.bath.energyChangeJ=-1;
    R::ProductKinematics a;a.term=0;a.referenceWeight=.25;a.velocityMS={{1,0,0}};a.nextClock=Fresh(.2);
    R::ProductKinematics b;b.term=1;b.referenceWeight=4;b.velocityMS={{-1,0,0}};b.nextClock=Fresh(.8);
    policy.products={a,b};
  }
  R::EventPlan Plan(R::PlanScope scope=R::PlanScope::Transport) {
    R::EventPlan plan;std::string e;
    if (!R::MakeEventPlan(catalog,reaction,parent,Interval(1000,10,3),scope,
        scope==R::PlanScope::RateOnly ? nullptr:&policy,&plan,&e)) throw std::runtime_error(e);
    return plan;
  }
};
// This adapter implements the same reserve/stage/atomic-publish contract that
// a PIC backend must subsequently verify with actual slots/lists. Only the
// adapter is a control fixture; all physics/plan/clock/weight code is production.
class Store:public R::TransactionBackend {
 public:
  std::size_t capacity=10,reserved=0;bool alive=true,stageFail=false,commitFail=false,throwReserve=false;
  int commits=0,rollbacks=0,reservations=0;
  std::vector<R::ProductPacket> published,staged;
  std::unique_ptr<R::EventLedger> ledger,stagedLedger;
  R::ReservationStatus Reserve(const R::EventPlan& p,std::string* error) override {
    ++reservations;reserved=std::min(capacity,p.products().size());
    if (throwReserve) throw std::runtime_error("control partial-reserve exception");
    if (reserved<p.products().size()) {*error="control capacity deferral";return R::ReservationStatus::Deferred;}
    return R::ReservationStatus::Ready;
  }
  bool Stage(const R::EventPlan& p,std::string* error) override {
    staged=p.products();stagedLedger.reset(new R::EventLedger(p.ledger()));
    if (stageFail) {*error="control stage refusal";return false;}return true;
  }
  bool Commit(const R::EventPlan&,std::string* error) override {
    if (commitFail) {*error="control commit refusal before publish";return false;}
    published.swap(staged);ledger.swap(stagedLedger);alive=false;++commits;reserved=0;
    return true; // no allocations/fallible operation after physical publication
  }
  void Rollback() noexcept override {reserved=0;staged.clear();stagedLedger.reset();++rollbacks;}
};

int main() {
  std::string error;Fixture f;
  Check(R::ValidateReaction(f.catalog,f.reaction,&error),"integer topology and upstream structural validation");
  auto bad=f.reaction;bad.products[0].multiplicity=2;
  Check(!R::ValidateReaction(f.catalog,bad,&error),"wrong nuclei refused");
  auto catalog=f.catalog;catalog[1].elementaryCharge=1;
  Check(!R::ValidateReaction(catalog,f.reaction,&error),"wrong integer charge refused");
  catalog=f.catalog;catalog[2].identity=catalog[1].identity;catalog[2].state=catalog[1].state;
  Check(!R::ValidateCatalog(catalog,&error),"duplicate exact species state refused");
  catalog=f.catalog;catalog[1].modelMassKg=0;
  Check(!R::ValidateCatalog(catalog,&error),"zero model mass refused");
  catalog=f.catalog;catalog[0].nuclei[0].count=std::numeric_limits<std::int64_t>::max();
  bad=f.reaction;bad.reactants[0].multiplicity=2;
  Check(!R::ValidateReaction(catalog,bad,&error),"nonunit tracked parent refused");
  catalog=f.catalog;catalog[1].nuclei[0].count=std::numeric_limits<std::int64_t>::max();bad=f.reaction;bad.products[0].multiplicity=3;
  Check(!R::ValidateReaction(catalog,bad,&error),"integer composition multiplication overflow refused");
  auto interval=Interval(1000,10,3),unchanged=interval;
  Check(!interval.Advance(8,&unchanged,&error) && unchanged.elapsedS()==3,"invalid local advance leaves output unchanged");
  R::TraceInterval next;
  Check(interval.Advance(7,&next,&error) && next.startS()==1000 && next.durationS()==10 && next.elapsedS()==10,"local interval endpoint retained");
  R::TraceInterval copied;std::memcpy(&copied,&interval,sizeof(copied));
  Check(copied.startS()==1000 && copied.remainingS()==7,"interval byte-storage roundtrip");
  auto clock=Fresh(.5),saved=clock;
  Check(!R::OpticalDepthClock::FromOpenUniform(0,&saved,&error) && saved.remaining()==clock.remaining(),"open zero draw refused transactionally");
  Check(!R::OpticalDepthClock::FromOpenUniform(1,&saved,&error),"open one draw refused");
  Check(!R::OpticalDepthClock::FromOpenUniform(std::numeric_limits<double>::quiet_NaN(),&saved,&error),"nonfinite draw refused");
  R::OpticalDepthClock copyClock;std::memcpy(&copyClock,&clock,sizeof(clock));
  Check(copyClock.fresh() && copyClock.remaining()==clock.remaining(),"clock byte-storage roundtrip");
  // Decimal ln(2) from the alternating/atanh log series, independent of
  // native libm. A binary64 log alone loses about 2.3e-17 at this input.
  // Compare both retained components, including after an ordinary PIC-style
  // byte copy, rather than comparing their rounded binary64 projection.
  const auto extendedMode=R::OpticalDepthPrecision::CompensatedExtendedLog;
  R::OpticalDepthClock extended;
  if (std::numeric_limits<long double>::radix==2 &&
      std::numeric_limits<long double>::digits>std::numeric_limits<double>::digits &&
      std::numeric_limits<long double>::digits<=2*std::numeric_limits<double>::digits) {
    double high=0,low=0;
    const long double ln2=0.693147180559945309417232121458176568075500134360255254120680L;
    Check(R::OpticalDepthClock::FromOpenUniform(.5,&extended,&error,extendedMode) &&
        extended.PendingDepthComponents(&high,&low,&error) && low!=0 &&
        std::abs(extended.remainingExtended()-ln2)<std::abs(static_cast<long double>(clock.remaining())-ln2),
        "compensated source log resolves independently known ln(2) more closely");
    R::OpticalDepthClock migrated;
    std::memcpy(&migrated,&extended,sizeof(extended));
    Check(migrated.fresh() && migrated.SamePendingDepth(extended) && !migrated.SamePendingDepth(clock),
        "extended clock preserves both physical depth components across an opaque copy");
  }
  else {
    Check(!R::OpticalDepthClock::FromOpenUniform(.5,&extended,&error,extendedMode) && !extended.valid(),
        "unsupported host cannot silently claim a higher precision optical clock");
  }
  const auto beforeInvalid=extended;
  Check(!R::OpticalDepthClock::FromOpenUniform(0,&extended,&error,extendedMode) &&
      (!beforeInvalid.valid() || extended.SamePendingDepth(beforeInvalid)),
      "invalid explicit physical source draw leaves extended clock unchanged");
  Check(!R::OpticalDepthClock::FromOpenUniform(.5,&extended,&error,
      static_cast<R::OpticalDepthPrecision>(-1)), "unknown optical-clock precision is refused");
  double delay=-1;
  const long double exactDelay=static_cast<long double>(clock.remaining())/3;
  Check(R::ConstantHazardEventDelay({{0,1},{1,2}},clock,&delay,&error) &&
      static_cast<long double>(delay)>=exactDelay &&
      static_cast<long double>(std::nextafter(delay,0.0))<exactDelay,
      "pure constant-hazard delay is first nonpremature binary64 time");
  const auto previewClock=clock;
  Check(R::ConstantHazardEventDelay({{0,0}},clock,&delay,&error) && std::isinf(delay) &&
      clock.fresh() && clock.remaining()==previewClock.remaining(),"zero-rate pure preview changes no clock");
  Check(!R::ConstantHazardEventDelay({{0,1e-320}},clock,&delay,&error),"unrepresentable positive delay refused");
  R::HazardClockAdvance reconstructed,clockSentinel;
  const auto endpointClock=Fresh(std::exp(-1.0));
  Check(R::AdvanceConstantHazardClock({{0,.25},{1,.75}},1.0L,endpointClock,.2,&reconstructed,&error) &&
      reconstructed.event && reconstructed.reaction==0 && reconstructed.clock.consumed(),
      "origin reconstruction consumes exact caller endpoint independently of rounded segment partitions");
  const auto futureClock=Fresh(std::nextafter(std::exp(-1.0),0.0));
  Check(futureClock.remaining()==std::nextafter(1.0,std::numeric_limits<double>::infinity()) &&
      R::AdvanceConstantHazardClock({{0,1}},1.0L,futureClock,NAN,&reconstructed,&error) &&
      !reconstructed.event && reconstructed.clock.remaining()>0,
      "nextafter-future optical endpoint remains pending without a tolerance-based event");
  const long double preciseTime=.3L+.2L;
  Check(R::AdvanceConstantHazardClock({{0,1}},preciseTime,endpointClock,NAN,&reconstructed,&error) &&
      !reconstructed.event && std::abs(static_cast<long double>(reconstructed.clock.remaining())-(1-preciseTime))<1e-18L,
      "total extended-precision elapsed reconstructs a pending origin clock without substep accumulation");
  Check(R::AdvanceConstantHazardClock({{0,0}},0,endpointClock,NAN,&reconstructed,&error) &&
      !reconstructed.event && reconstructed.clock.remaining()==endpointClock.remaining(),
      "zero elapsed and zero hazard inspect no channel draw");
  clockSentinel=reconstructed;
  Check(!R::AdvanceConstantHazardClock({{0,1}},-1,endpointClock,.5,&reconstructed,&error) &&
      reconstructed.clock.remaining()==clockSentinel.clock.remaining(),
      "invalid anchored elapsed leaves output unchanged");
  R::HazardClockAdvance consumed;
  if(!R::AdvanceConstantHazardClock({{0,1}},1,endpointClock,.5,&consumed,&error))throw std::runtime_error(error);
  Check(!R::AdvanceConstantHazardClock({{0,1}},0,consumed.clock,.5,&reconstructed,&error),
      "consumed clock cannot become an origin clock");
  R::HazardAdvance outcome;
  Check(R::AdvanceConstantHazards({{0,.25},{1,.75}},Interval(7,2),clock,.2,&outcome,&error) && outcome.event && outcome.reaction==0 && Near(outcome.eventOffsetS,std::log(2.0)) && outcome.clock.consumed(),"analytic event time and first channel");
  Check(outcome.interval.startS()==7 && Near(outcome.interval.remainingS(),2-std::log(2.0)),"event preserves original local interval");
  Check(R::AdvanceConstantHazards({{0,.25},{1,.75}},Interval(7,2),clock,.8,&outcome,&error) && outcome.reaction==1,"analytic competing second channel");
  Check(R::AdvanceConstantHazards({{0,0},{1,0}},Interval(7,2),clock,std::numeric_limits<double>::quiet_NaN(),&outcome,&error) && !outcome.event && outcome.interval.remainingS()==0 && outcome.clock.remaining()==clock.remaining(),"zero hazard consumes no random decision and no optical depth");
  Check(R::AdvanceConstantHazards({{0,.1}},Interval(7,2),clock,0,&outcome,&error) && !outcome.event && Near(outcome.clock.remaining(),std::log(2.0)-.2) && !outcome.clock.fresh(),"no-event optical continuation");
  Check(!R::ConsumeExposure(clock,clock.remaining(),&saved,&error),"varying-rate crossing must be located by actual trajectory");
  Check(R::ConsumeExposure(clock,.1,&saved,&error) && Near(saved.remaining(),std::log(2.0)-.1) && !saved.fresh(),"caller-integrated no-event exposure");
  auto tiny=clock;
  for (int i=0;i<10000;++i) {
    R::OpticalDepthClock updated;
    if (!R::ConsumeExposure(tiny,1e-18,&updated,&error)) throw std::runtime_error(error);
    tiny=updated;
  }
  Check(std::abs(tiny.remaining()-(clock.remaining()-1e-14))<2e-16,"sub-ULP segment exposures accumulate without erased hazard");
  Check(!R::AdvanceConstantHazards({{0,1},{0,1}},Interval(0,1),clock,.5,&outcome,&error),"duplicate hazards refused");
  Check(!R::AdvanceConstantHazards({{0,-1}},Interval(0,1),clock,.5,&outcome,&error),"negative hazard refused");
  double probability=-1;
  Check(R::EventProbability({{0,1e-12}},1e-6,&probability,&error) && std::abs(probability/1e-18-1)<1e-14,"small optical probability uses stable expm1");
  Check(R::EventProbability({{0,1e300}},1e300,&probability,&error) && probability==1,"large finite optical exposure saturates physically");
  Check(R::AdvanceConstantHazards({{0,1}},Interval(0,clock.remaining()),clock,.5,&outcome,&error) && outcome.event && outcome.interval.remainingS()==0,"event exactly at local endpoint");
  auto plan=f.Plan();const auto& ledger=plan.ledger();
  Check(ledger.integerConservationVerified && ledger.reactantChargeElementary==ledger.productChargeElementary && ledger.reactantNuclei.size()==2,"exact integer material ledger retained");
  Check(plan.products().size()==2 && Near(ledger.parentNumber,1.25) && Near(ledger.productNumber,2.5),"physical represented numbers and multiplicities");
  Check(Near(plan.products()[0].individualCorrection,5) && Near(plan.products()[1].individualCorrection,.3125),"different reference weights without timestep factor");
  auto fractionalParent=f.parent;fractionalParent.referenceWeight=.1;fractionalParent.individualCorrection=1;
  auto fractionalPolicy=f.policy;
  fractionalPolicy.products[0].referenceWeight=.7;fractionalPolicy.products[1].referenceWeight=.7;
  R::EventPlan fractionalPlan;
  Check(R::MakeEventPlan(f.catalog,f.reaction,fractionalParent,interval,R::PlanScope::Transport,&fractionalPolicy,&fractionalPlan,&error) && fractionalPlan.products()[0].targetRepresentedNumber==.1 && fractionalPlan.products()[0].representedNumber==.7*(.1/.7) && fractionalPlan.ledger().productNumber==fractionalPlan.products()[0].representedNumber+fractionalPlan.products()[1].representedNumber,"ledger records actual rounded reference-weight reconstruction rather than ideal quotient");
  Check(plan.products()[0].interval.startS()==1000 && plan.products()[0].interval.elapsedS()==3 && plan.products()[0].interval.remainingS()==7 && plan.products()[0].positionM==f.parent.positionM,"products inherit physical event position and local continuation");
  Check(plan.products()[0].clock.fresh() && plan.products()[1].clock.fresh() && plan.products()[0].clock.remaining()!=plan.products()[1].clock.remaining(),"fresh daughter optical depths retained");
  Check(Near(ledger.trackedEnergyChangeJ,1.25) && Near(ledger.bathEnergyChangeJ,-1.25) && Near(ledger.parentModelMassKg,ledger.productModelMassKg),"explicit bath energy and model mass ledger");
  R::EventPlan localPlan;
  Check(R::MakeEventPlan(f.catalog,f.reaction,f.parent,Interval(3,.125,.03125),R::PlanScope::Transport,&f.policy,&localPlan,&error) && localPlan.products()[0].interval.durationS()==.125 && localPlan.products()[0].interval.remainingS()==.09375,"arbitrary per-particle local tracing interval is preserved");
  auto repeatedCatalog=f.catalog;repeatedCatalog[0].nuclei={{"X",2}};
  auto repeatedReaction=f.reaction;repeatedReaction.products={{1,2,R::Population::Tracked}};
  auto repeatedPolicy=f.policy;repeatedPolicy.products[1].term=0;repeatedPolicy.products[1].occurrence=1;
  R::EventPlan repeatedPlan;
  Check(R::MakeEventPlan(repeatedCatalog,repeatedReaction,f.parent,interval,R::PlanScope::Transport,&repeatedPolicy,&repeatedPlan,&error) && repeatedPlan.products().size()==2 && repeatedPlan.products()[0].species==1 && repeatedPlan.products()[1].species==1 && repeatedPlan.products()[0].velocityMS[0]!=repeatedPlan.products()[1].velocityMS[0],"repeated products retain separate correlated carriers and occurrence identity");
  repeatedPolicy.products.pop_back();
  Check(!R::MakeEventPlan(repeatedCatalog,repeatedReaction,f.parent,interval,R::PlanScope::Transport,&repeatedPolicy,&repeatedPlan,&error),"double-weight single carrier cannot replace two specified daughter velocities");
  auto states=f.catalog;states.push_back({"control-X","excited",1,{{"X",1}},0,1});
  auto stateReaction=f.reaction;stateReaction.products[0].species=3;
  R::EventPlan statePlan;
  Check(R::MakeEventPlan(states,stateReaction,f.parent,interval,R::PlanScope::Transport,&f.policy,&statePlan,&error) && statePlan.products()[0].species==3 && statePlan.products()[0].compiledSpecies==1 && statePlan.products()[0].stateIdentity=="excited","distinct internal state survives shared compiled PIC species index");
  states[3].state="CHANGED-AFTER-PLAN";
  Check(statePlan.products()[0].stateIdentity=="excited" && statePlan.parentStateIdentity()=="ground","pending plan owns resolved state identities independently of later catalog edits");
  auto charged=f.catalog;
  charged.push_back({"control-electron","free",-1,{},-1,.001});
  charged.push_back({"control-dimer-ion","ground",3,{{"X",1},{"Y",1}},1,1.999});
  R::Reaction ionization{"control-bath-ionization",{{0,1,R::Population::Tracked},{3,1,R::Population::Bath}},{{4,1,R::Population::Tracked},{3,2,R::Population::Bath}}};
  auto ionPolicy=f.policy;ionPolicy.bath.energyChangeJ=0;ionPolicy.bath.modelMassChangeKg=.001;
  ionPolicy.products.resize(1);ionPolicy.products[0].velocityMS={{0,0,0}};
  R::EventPlan ionPlan;
  Check(R::MakeEventPlan(charged,ionization,f.parent,interval,R::PlanScope::Transport,&ionPolicy,&ionPlan,&error) && ionPlan.ledger().reactantChargeElementary==-1 && ionPlan.ledger().productChargeElementary==-1 && ionPlan.ledger().trackedChargeChangeElementary==1 && ionPlan.ledger().bathChargeChangeElementary==-1,"full binary electron bath stoichiometry and exact charge close");
  Check(ionPlan.ledger().bathTerms.size()==2 && ionPlan.ledger().bathTerms[0].reactant && Near(ionPlan.ledger().bathTerms[0].representedNumber,1.25) && Near(ionPlan.ledger().bathTerms[1].representedNumber,2.5),"bath incoming and outgoing multiplicities are explicit");
  ionPolicy.bath.modelMassChangeKg=0;
  Check(!R::MakeEventPlan(charged,ionization,f.parent,interval,R::PlanScope::Transport,&ionPolicy,&ionPlan,&error),"bath model mass cannot be omitted or silently repaired");
  R::EventPlan sentinel=plan;
  Check(!R::MakeEventPlan(f.catalog,f.reaction,f.parent,interval,R::PlanScope::Transport,nullptr,&sentinel,&error) && sentinel.products().size()==2,"missing kinematic policy refuses transport without changing output");
  auto policy=f.policy;policy.products.pop_back();
  Check(!R::MakeEventPlan(f.catalog,f.reaction,f.parent,interval,R::PlanScope::Transport,&policy,&sentinel,&error),"missing product occurrence refused");
  policy=f.policy;policy.products[0].nextClock=saved;
  Check(!R::MakeEventPlan(f.catalog,f.reaction,f.parent,interval,R::PlanScope::Transport,&policy,&sentinel,&error),"continued parent clock cannot masquerade as fresh daughter clock");
  policy=f.policy;policy.bath.energyChangeJ=0;
  Check(!R::MakeEventPlan(f.catalog,f.reaction,f.parent,interval,R::PlanScope::Transport,&policy,&sentinel,&error),"missing bath energy does not get fitted from residual");
  policy=f.policy;policy.products[0].velocityMS[0]=11;
  Check(!R::MakeEventPlan(f.catalog,f.reaction,f.parent,interval,R::PlanScope::Transport,&policy,&sentinel,&error),"caller-qualified classical domain enforced");
  auto parent=f.parent;parent.referenceWeight=std::numeric_limits<double>::max();parent.individualCorrection=2;
  Check(!R::MakeEventPlan(f.catalog,f.reaction,parent,interval,R::PlanScope::Transport,&f.policy,&sentinel,&error),"represented inventory overflow refused");
  auto rate=f.Plan(R::PlanScope::RateOnly);Store store;R::PendingEvent pendingRate(rate);
  Check(rate.valid() && rate.products().empty() && Near(rate.ledger().productNumber,2.5) && pendingRate.TryCommit(&store,&error)==R::AttemptStatus::Error && store.reservations==0,"rate-only topology cannot commit particles");
  R::PendingEvent pending(plan);store.capacity=1;
  const double daughterClock=pending.plan().products()[0].clock.remaining();
  Check(pending.TryCommit(&store,&error)==R::AttemptStatus::Deferred && store.alive && store.published.empty() && store.reserved==0 && store.rollbacks==1 && !pending.committed(),"capacity deferral rolls back partial reservation");
  store.capacity=2;store.stageFail=true;
  Check(pending.TryCommit(&store,&error)==R::AttemptStatus::Error && store.alive && store.published.empty() && !store.stagedLedger && store.staged.empty(),"stage failure rolls back without parent or ledger mutation");
  store.stageFail=false;store.commitFail=true;
  Check(pending.TryCommit(&store,&error)==R::AttemptStatus::Error && store.alive && store.commits==0,"atomic commit refusal retains physical state");
  store.commitFail=false;store.throwReserve=true;
  Check(pending.TryCommit(&store,&error)==R::AttemptStatus::Error && store.reserved==0 && store.alive,"backend exception rolls back partial capacity");
  store.throwReserve=false;
  Check(pending.TryCommit(&store,&error)==R::AttemptStatus::Committed && !store.alive && store.published.size()==2 && store.commits==1 && pending.committed() && pending.plan().products()[0].clock.remaining()==daughterClock,"same immutable plan succeeds after retries without redraw");
  const int reservations=store.reservations;
  Check(pending.TryCommit(&store,&error)==R::AttemptStatus::Error && store.reservations==reservations && store.commits==1,"double commit refused before backend mutation");
  // Fixed-probability sampler control, not a physical channel table. Its
  // predeclared 99% simultaneous Hoeffding interval supplements exact CDF
  // boundary checks and does not claim full U17/I24 qualification.
  std::mt19937_64 rng(731993);std::uniform_real_distribution<double> uniform(0,1);
  const int n=100000;int first=0;
  for (int i=0;i<n;++i) {
    double u;do {u=uniform(rng);} while (!(u>0 && u<1));
    if (!R::AdvanceConstantHazards({{0,.25},{1,.75}},Interval(0,2),clock,u,&outcome,&error)) throw std::runtime_error(error);
    if (outcome.reaction==0) ++first;
  }
  const double bound=std::sqrt(.5*n*std::log(4/.01));
  Check(std::abs(first-.25*n)<=bound,"100000 competing-channel draws satisfy predeclared 99-percent control bound");
  std::cout<<"SUMMARY checks="<<checks<<" failed="<<failures<<" first_channel="<<first<<" seed=731993\n";
  return failures ? 1:0;
}
