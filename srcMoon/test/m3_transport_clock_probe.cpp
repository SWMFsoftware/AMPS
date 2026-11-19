// Analytical TIME/PROTOCOL controls only. No observation, lunar reaction rate,
// field, particle source distribution or sampling correction is fabricated.
#include "MoonTransportClock.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
namespace T=Moon::TransportTime;
int Checks=0,Failures=0;
void Check(bool good,const char* name){++Checks;if(!good)++Failures;std::cout<<(good?"PASS ":"FAIL ")<<name<<'\n';}
int main() {
  static_assert(std::is_trivially_copyable<T::Step>::value,"step tokens have no process-local pointers");
  static_assert(std::is_trivially_copyable<T::Clock>::value,"candidate clock publication is nonallocating");
  std::string error;T::Clock clock;T::Step step,sentinel;sentinel.sequence=77;
  Check(!clock.Begin(1,&sentinel,&error) && sentinel.sequence==77,"uninitialized Begin preserves caller output");
  Check(!clock.Initialize(NAN,&error) && !clock.initialized(),"nonfinite origin refused without initialization");
  Check(clock.Initialize(1000,&error) && clock.elapsedS()==0 && clock.completedSteps()==0,"explicit origin initializes empty clock");
  Check(!clock.Initialize(0,&error) && clock.originEtS()==1000,"one-shot origin cannot silently reinitialize");
  Check(!clock.Begin(0,&step,&error) && !clock.active(),"zero outer step refused");
  Check(!clock.Begin(-1,&step,&error) && !clock.active(),"negative outer step refused");
  Check(!clock.Begin(INFINITY,&step,&error) && !clock.active(),"infinite outer step refused");
  Check(clock.Begin(1,&step,&error) && step.sequence==1 && step.relativeBeginS==0 && step.relativeEndS==1 && step.beginEtS==1000 && step.endEtS==1001,
      "immutable first physical step has exact relative and projected endpoints");
  Check(!clock.Begin(1,&sentinel,&error) && sentinel.sequence==77,"nested Begin refused without token mutation");
  T::BirthTime birth,previous;previous.offsetS=77;birth=previous;
  Check(clock.BirthFromOffset(0,&birth,&error) && birth.relativeS==0 && birth.epochEtS==1000 && birth.remainingS==1,"zero offset is beginning birth with full interval");
  Check(clock.BirthFromOffset(.25,&birth,&error) && birth.relativeS==.25L && birth.epochEtS==1000.25 && birth.remainingS==.75,"source birth uses relative offset before absolute projection");
  Check(clock.BirthFromOffset(1,&birth,&error) && birth.relativeS==1 && birth.epochEtS==1001 && birth.remainingS==0,"endpoint source has zero unused seconds");
  birth=previous;Check(!clock.BirthFromOffset(-.1,&birth,&error) && birth.offsetS==77,"negative source offset refused transactionally");
  Check(!clock.BirthFromOffset(std::nextafter(1.0,INFINITY),&birth,&error),"strict future source offset cannot be clamped to endpoint");
  Check(!clock.BirthFromOffset(NAN,&birth,&error),"nonfinite departure offset refused");
  Check(clock.BirthFromRemaining(.75,&birth,&error) && birth.relativeS==.25L && birth.epochEtS==1000.25 && birth.offsetS==.25,"remaining-interval source has same declared birth");
  Check(clock.BirthFromRemaining(0,&birth,&error) && birth.relativeS==1 && birth.epochEtS==1001,"zero remaining source is exact endpoint");
  Check(clock.BirthFromRemaining(1,&birth,&error) && birth.relativeS==0 && birth.epochEtS==1000,"full remaining source is exact beginning");
  Check(!clock.BirthFromRemaining(-1,&birth,&error) && !clock.BirthFromRemaining(2,&birth,&error),"negative or oversized source remaining interval refused");
  Check(!clock.BirthFromRemaining(1e-30,&birth,&error),"unrepresentable positive source residual does not become endpoint birth");
  auto forged=step;forged.sequence++;
  Check(!clock.Commit(forged,&error) && clock.active() && clock.elapsedS()==0,"forged sequence cannot publish elapsed time");
  forged=step;forged.durationS=.5;
  Check(!clock.AbortUnaccepted(forged,&error) && clock.active(),"forged duration cannot abort another active step");
  Check(clock.AbortUnaccepted(step,&error) && !clock.active() && clock.completedSteps()==0 && clock.elapsedS()==0,"PRE-WORK abort retains original accepted chronology");
  Check(!clock.AbortUnaccepted(step,&error),"duplicate abort refused");
  Check(clock.Begin(1,&step,&error) && step.sequence==1 && clock.Commit(step,&error) && clock.elapsedS()==1 && clock.completedSteps()==1,"fresh accepted step publishes exactly once");
  Check(!clock.Commit(step,&error) && clock.elapsedS()==1,"duplicate commit cannot advance physical time");
  Check(!clock.BirthFromOffset(0,&birth,&error),"source outside active physical step refused");
  T::Step second;Check(clock.Begin(.5,&second,&error) && second.sequence==2 && second.beginEtS==1001 && second.endEtS==1001.5,"next step starts at accepted transport endpoint instead of frozen origin");
  Check(!clock.Commit(step,&error) && clock.active(),"stale preceding token cannot commit next interval");
  Check(clock.BirthFromOffset(.25,&birth,&error) && birth.epochEtS==1001.25 && clock.Commit(second,&error),"repeated source birth advances while external origin remains fixed");
  Check(clock.originEtS()==1000 && clock.elapsedS()==1.5L,"external origin is immutable across transport commits");
  T::Clock accumulated;accumulated.Initialize(3e8,&error);double naive=3e8;long double maximumRelativeError=0;
  bool all=true;constexpr int Count=100000;constexpr double Dt=.1;
  for(int i=1;i<=Count;++i) {
    T::Step next;
    all=all && accumulated.Begin(Dt,&next,&error) && accumulated.Commit(next,&error);
    const long double reference=static_cast<long double>(i)*static_cast<long double>(Dt);
    maximumRelativeError=std::max(maximumRelativeError,std::abs(accumulated.elapsedS()-reference));
    naive+=Dt;
  }
  const long double reference=static_cast<long double>(Count)*Dt;
  const long double arithmeticBound=4*std::numeric_limits<long double>::epsilon()*reference;
  Check(all && maximumRelativeError<=arithmeticBound,"100000 intervals match independent binary-dt product within long-double arithmetic bound");
  double projected=0;accumulated.Project(accumulated.elapsedS(),&projected,&error);
  const double exactProjected=static_cast<double>(3e8L+reference);
  Check(projected==exactProjected && std::abs(naive-exactProjected)>std::abs(projected-exactProjected),"relative compensated clock avoids repeated large-ET addition drift");
  T::Clock huge;huge.Initialize(1e30,&error);
  Check(!huge.Begin(.1,&step,&error) && !huge.active() && huge.elapsedS()==0,"outer step below absolute epoch representation is explicit refusal");
  T::Clock negative;negative.Initialize(-2,&error);
  Check(negative.Begin(.25,&step,&error) && step.endEtS==-1.75 && negative.Commit(step,&error),"negative explicit TDB epoch is valid and chronologically increasing");
  const double one=1,next1=std::nextafter(one,INFINITY),next2=std::nextafter(next1,INFINITY),next3=std::nextafter(next2,INFINITY);
  Check(T::EpochProjectionAgrees(one,one) && T::EpochProjectionAgrees(next2,one),"derived two-addition epoch projection envelope accepted");
  Check(!T::EpochProjectionAgrees(next3,one) && !T::EpochProjectionAgrees(NAN,one),"larger or nonfinite endpoint discrepancy cannot silently align");
  std::cout<<std::setprecision(21)<<"METRIC relative_error_s="<<maximumRelativeError<<" arithmetic_bound_s="<<arithmeticBound
    <<" naive_absolute_drift_s="<<std::abs(naive-exactProjected)<<'\n';
  std::cout<<"SUMMARY checks="<<Checks<<" failed="<<Failures<<'\n';return Failures?1:0;
}
