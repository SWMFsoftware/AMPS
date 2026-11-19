#ifndef AMPS_CHEMICAL_DATABASE_REDUCED_MODELS_H
#define AMPS_CHEMICAL_DATABASE_REDUCED_MODELS_H

#include "Runtime.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <new>
#include <optional>
#include <vector>

namespace AMPS { namespace ChemicalDatabase { namespace ReducedModels {

// These explicit-parameter mathematical services select no species, material,
// source profile, units conversion, RNG or application. Every momentum/energy
// in one call belongs to the SAME caller-declared inertial or local event frame.
// Q is a supplied net available energy in that frame's declared state reference;
// do not count formation/rest energy in Q twice. Runtime Status means numerical
// evaluation only, never scientific readiness or physical source authorization.
using Vector3 = std::array<long double,3>;
inline Status Refuse(Code c,const char* message) { return Status{c,message}; }
inline bool Finite(long double x) { return std::isfinite(x); }
inline bool Finite(const Vector3& x) { return Finite(x[0])&&Finite(x[1])&&Finite(x[2]); }
inline long double Sum3(long double a,long double b,long double c) {
  // Neumaier compensation retains a small separately supplied Q while two
  // large kinetic/bulk contributions cancel. No physical tolerance is added.
  long double sum=a,correction=0;
  for(long double value: {b,c}) {
    const long double next=sum+value;
    correction+=(std::fabs(sum)>=std::fabs(value))?(sum-next)+value:(value-next)+sum;
    sum=next;
  }
  return sum+correction;
}
inline long double Sum4(long double a,long double b,long double c,long double d) {
  long double sum=a,correction=0;
  for(long double value: {b,c,d}) {
    const long double next=sum+value;
    correction+=(std::fabs(sum)>=std::fabs(value))?(sum-next)+value:(value-next)+sum;
    sum=next;
  }
  return sum+correction;
}
inline long double Kinetic(const Vector3& p,long double mass) {
  // Divide before multiplying to avoid an unnecessary p^2 intermediate.
  return Sum3((p[0]/mass)*p[0]/2,(p[1]/mass)*p[1]/2,(p[2]/mass)*p[2]/2);
}

struct IonMass {
  long double referenceMassKg,bindingMassKg,physicalMassKg;
};
inline Status IonPhysicalMass(long double neutralMassKg,long double electronMassKg,
    long double ionizationEnergyJ,long double speedOfLightMS,IonMass* output) {
  if(!output||!Finite(neutralMassKg)||!Finite(electronMassKg)||!Finite(ionizationEnergyJ)||
      !Finite(speedOfLightMS)||electronMassKg<=0||neutralMassKg<=electronMassKg||
      ionizationEnergyJ<0||speedOfLightMS<=0)
    return Refuse(Code::InvalidInput,"explicit positive masses/c and nonnegative ionization energy required");
  IonMass value{neutralMassKg-electronMassKg,
                (ionizationEnergyJ/speedOfLightMS)/speedOfLightMS,0};
  value.physicalMassKg=value.referenceMassKg+value.bindingMassKg;
  if(!Finite(value.bindingMassKg)||!Finite(value.physicalMassKg)||value.physicalMassKg<=0||
      (ionizationEnergyJ>0&&value.bindingMassKg==0))
    return Refuse(Code::NumericalFailure,"ion mass/binding contribution is outside numeric representation");
  *output=value; return Status{};
}

struct RadiativeSet {
  long double totalRateS,lifetimeS;
  std::vector<long double> branches;
};
inline Status RadiativeBranches(const std::vector<long double>& ratesS,RadiativeSet* output) {
  if(!output||ratesS.empty()) return Refuse(Code::InvalidInput,"explicit nonempty coherent A-set/output required");
  long double sum=0,correction=0;
  for(long double a:ratesS) {
    if(!Finite(a)||a<0) return Refuse(Code::InvalidInput,"A-values must be finite nonnegative s^-1; missing is not zero");
    const long double next=sum+a;
    if(!Finite(next)) return Refuse(Code::NumericalFailure,"A-set sum overflow");
    correction+=(std::fabs(sum)>=std::fabs(a))?(sum-next)+a:(a-next)+sum;
    sum=next;
  }
  sum+=correction;
  if(!(sum>0)) return Refuse(Code::OutsideDomain,"zero total A: no finite radiative lifetime");
  RadiativeSet value{sum,1/sum,{}};
  if(!Finite(sum)||!Finite(value.lifetimeS)||value.lifetimeS==0)
    return Refuse(Code::NumericalFailure,"radiative sum/lifetime is not representable");
  try {
    value.branches.reserve(ratesS.size());
    for(long double a:ratesS) {
      const long double branch=a/sum;
      if(a>0&&branch==0) return Refuse(Code::NumericalFailure,"positive radiative branch underflow; no channel deleted");
      value.branches.push_back(branch);
    }
  } catch(const std::bad_alloc&) {
    return Refuse(Code::ResourceLimit,"radiative branch allocation failed; reduce supplied state set or increase memory");
  }
  *output=std::move(value); return Status{};
}

inline Status OpalSecondaryEnergy(long double incidentEnergyJ,long double thresholdJ,
    long double widthJ,long double quantile,long double* output) {
  // Opal lower-energy marginal only, wmax=(E-I)/2. This inverse contains no
  // angular law, correlated partner energy or full ion-recoil event closure.
  if(!output||!Finite(incidentEnergyJ)||!Finite(thresholdJ)||!Finite(widthJ)||!Finite(quantile)||
      thresholdJ<0||incidentEnergyJ<=thresholdJ||widthJ<=0||quantile<0||quantile>1)
    return Refuse(Code::InvalidInput,"Opal requires E>I>=0, W>0 and an explicit quantile in [0,1]");
  const long double maximum=(incidentEnergyJ-thresholdJ)/2;
  if(!(maximum>0)) return Refuse(Code::NumericalFailure,"Opal support is not representable");
  if(quantile==0) {*output=0; return Status{};}
  if(quantile==1) {*output=maximum; return Status{};}
  const long double angle=std::atan2(maximum,widthJ);
  // Complementary-angle evaluation avoids loss of a tiny distance from the
  // upper endpoint when atan(maximum/W) rounds near pi/2. No epsilon extends
  // the physical support or changes a supplied quantile.
  // For maximum/W<=1 the direct angle is at most pi/4 and well conditioned
  // even near quantile1; using its complement would instead lose tiny energies
  // through cancellation near pi/2. Choose solely by numerical conditioning.
  const long double value=(maximum<=widthJ||quantile<0.5L) ? widthJ*std::tan(quantile*angle) :
      widthJ/std::tan(std::atan2(widthJ,maximum)+(1-quantile)*angle);
  if(!Finite(value)||value<=0||value>maximum)
    return Refuse(Code::NumericalFailure,"Opal quantile cannot be represented inside its declared support");
  *output=value; return Status{};
}

struct ConservationBounds {
  long double momentumAbsoluteKgMS,momentumRelative,energyAbsoluteJ,energyRelative,unitVectorTolerance;
};
inline bool Valid(const ConservationBounds& b) {
  return Finite(b.momentumAbsoluteKgMS)&&Finite(b.momentumRelative)&&Finite(b.energyAbsoluteJ)&&
    Finite(b.energyRelative)&&Finite(b.unitVectorTolerance)&&b.momentumAbsoluteKgMS>=0&&b.momentumRelative>=0&&
    b.energyAbsoluteJ>=0&&b.energyRelative>=0&&b.unitVectorTolerance>=0&&b.unitVectorTolerance<1;
}
struct TwoBodyState {
  std::array<Vector3,2> momentumKgMS,velocityMS;
  std::array<long double,2> kineticEnergyJ;
  long double centreOfMassKineticEnergyJ,totalKineticEnergyJ;
};
inline Status TwoBodyProducts(const Vector3& totalMomentumKgMS,long double incomingEnergyJ,
    long double netQJ,long double firstMassKg,long double secondMassKg,const Vector3& direction,
    const ConservationBounds& bounds,TwoBodyState* output) {
  if(!output||!Finite(totalMomentumKgMS)||!Finite(direction)||!Finite(incomingEnergyJ)||!Finite(netQJ)||
      !Finite(firstMassKg)||!Finite(secondMassKg)||firstMassKg<=0||secondMassKg<=0||incomingEnergyJ<0||!Valid(bounds))
    return Refuse(Code::InvalidInput,"two-body explicit SI inputs/arithmetic bounds are invalid");
  const long double norm=std::hypot(direction[0],direction[1],direction[2]);
  if(!Finite(norm)||norm==0||std::fabs(norm-1)>bounds.unitVectorTolerance)
    return Refuse(Code::InvalidInput,"two-body direction must be nonzero and unit within explicit caller bound");
  const long double totalMass=firstMassKg+secondMassKg;
  if(!Finite(totalMass)) return Refuse(Code::NumericalFailure,"two-body total mass overflow");
  if(totalMass==firstMassKg||totalMass==secondMassKg)
    return Refuse(Code::NumericalFailure,"two-body mass sum loses a positive mass; use greater arithmetic precision");
  const long double bulk=Kinetic(totalMomentumKgMS,totalMass);
  const long double com=Sum3(incomingEnergyJ,netQJ,-bulk);
  if(!Finite(bulk)||!Finite(com)) return Refuse(Code::NumericalFailure,"two-body available energy overflow");
  if(com<0) return Refuse(Code::OutsideDomain,"two-body state has negative available COM kinetic energy; no tolerance clamp");
  if(com==0) {
    // A rounded zero bulk/Q difference must not erase a representable recoil.
    // Certify the special zero-COM case with error-free product/sum residuals;
    // otherwise refuse instead of inventing a mass-ratio or energy tolerance.
    const long double z=totalMass-firstMassKg;
    const long double massTail=(firstMassKg-(totalMass-z))+(secondMassKg-z);
    const long double target=incomingEnergyJ+netQJ;
    const long double qz=target-incomingEnergyJ;
    const long double energyTail=(incomingEnergyJ-(target-qz))+(netQJ-qz);
    // Below this derived arithmetic floor an exact product residual can itself
    // underflow, so fma(...,-roundedProduct)==0 no longer certifies exactness.
    // This is a floating-representation limit, never a physical energy cutoff.
    const long double productFloor=std::scalbn(std::numeric_limits<long double>::min(),
                                              std::numeric_limits<long double>::digits);
    long double squareSum=0;bool exact=massTail==0&&energyTail==0&&Finite(2*totalMass);
    for(long double p:totalMomentumKgMS) {
      const long double square=p*p,next=squareSum+square;
      const long double sz=next-squareSum;
      const long double squareTail=(squareSum-(next-sz))+(square-sz);
      // A product and its fma residual can BOTH underflow to zero. Such a
      // nonzero momentum is not proof of zero kinetic energy.
      exact=exact&&!(p!=0&&std::fabs(square)<productFloor)&&Finite(square)&&Finite(next)&&
        std::fma(p,p,-square)==0&&squareTail==0;
      squareSum=next;
    }
    const long double energyProduct=target*(2*totalMass);
    exact=exact&&Finite(target)&&Finite(energyProduct)&&
      !(target!=0&&std::fabs(energyProduct)<productFloor)&&
      std::fma(target,2*totalMass,-squareSum)==0;
    if(!exact) return Refuse(Code::NumericalFailure,"zero COM energy cannot be certified after bulk/Q cancellation; increase precision");
  }
  const long double reduced=(firstMassKg/totalMass)*secondMassKg;
  if(!(reduced>0)) return Refuse(Code::NumericalFailure,"positive two-body reduced mass underflow");
  const long double relative=std::sqrt(reduced)*std::sqrt(com)*std::sqrt(2.0L);
  if(!Finite(relative)||(com>0&&relative==0)) return Refuse(Code::NumericalFailure,"two-body relative momentum is not representable");
  TwoBodyState value{}; value.centreOfMassKineticEnergyJ=com;
  const std::array<long double,2> masses{{firstMassKg,secondMassKg}};
  for(int i=0;i<3;++i) {
    const long double recoil=relative*(direction[i]/norm);
    value.momentumKgMS[0][i]=(firstMassKg/totalMass)*totalMomentumKgMS[i]+recoil;
    value.momentumKgMS[1][i]=(secondMassKg/totalMass)*totalMomentumKgMS[i]-recoil;
    for(int j=0;j<2;++j) value.velocityMS[j][i]=value.momentumKgMS[j][i]/masses[j];
  }
  long double momentumError=0,momentumScale=0;
  for(int i=0;i<3;++i) {
    if(!Finite(value.momentumKgMS[0][i])||!Finite(value.momentumKgMS[1][i])||
       !Finite(value.velocityMS[0][i])||!Finite(value.velocityMS[1][i]))
      return Refuse(Code::NumericalFailure,"two-body product vector is not representable");
    momentumError=std::max(momentumError,std::fabs(Sum3(value.momentumKgMS[0][i],value.momentumKgMS[1][i],-totalMomentumKgMS[i])));
    momentumScale=std::max({momentumScale,std::fabs(totalMomentumKgMS[i]),std::fabs(value.momentumKgMS[0][i]),std::fabs(value.momentumKgMS[1][i])});
  }
  for(int j=0;j<2;++j) value.kineticEnergyJ[j]=Kinetic(value.momentumKgMS[j],masses[j]);
  value.totalKineticEnergyJ=value.kineticEnergyJ[0]+value.kineticEnergyJ[1];
  const long double target=incomingEnergyJ+netQJ;
  const long double energyError=std::fabs(Sum4(value.kineticEnergyJ[0],value.kineticEnergyJ[1],-incomingEnergyJ,-netQJ));
  const long double allowedMomentum=bounds.momentumAbsoluteKgMS+bounds.momentumRelative*momentumScale;
  const long double allowedEnergy=bounds.energyAbsoluteJ+bounds.energyRelative*std::max(std::fabs(target),std::fabs(value.totalKineticEnergyJ));
  if(!Finite(target)||!Finite(value.totalKineticEnergyJ)||!Finite(energyError)||
      !Finite(allowedMomentum)||!Finite(allowedEnergy)||momentumError>allowedMomentum||energyError>allowedEnergy)
    return Refuse(Code::NumericalFailure,"two-body output violates explicit arithmetic conservation bounds");
  *output=value; return Status{};
}

struct InheritedState {
  Vector3 heavyVelocityMS,heavyMomentumKgMS,signedReservoirMomentumKgMS;
  long double heavyKineticEnergyJ,signedReservoirEnergyJ,signedReservoirChargeC;
};
inline Status InheritParentVelocity(long double parentMassKg,long double heavyMassKg,
    long double parentChargeC,long double heavyChargeC,long double incomingBathChargeC,
    const Vector3& parentVelocityMS,const Vector3& incomingBathMomentumKgMS,long double incomingBathEnergyJ,
    long double netQJ,InheritedState* output) {
  if(!output||!Finite(parentMassKg)||!Finite(heavyMassKg)||parentMassKg<=0||heavyMassKg<=0||
      !Finite(parentChargeC)||!Finite(heavyChargeC)||!Finite(incomingBathChargeC)||
      !Finite(parentVelocityMS)||!Finite(incomingBathMomentumKgMS)||!Finite(incomingBathEnergyJ)||!Finite(netQJ))
    return Refuse(Code::InvalidInput,"inherit-velocity requires explicit positive masses and finite same-frame inputs");
  InheritedState value{};value.heavyVelocityMS=parentVelocityMS;
  Vector3 parentMomentum{};
  for(int i=0;i<3;++i) {
    parentMomentum[i]=parentMassKg*parentVelocityMS[i];
    value.heavyMomentumKgMS[i]=heavyMassKg*parentVelocityMS[i];
    value.signedReservoirMomentumKgMS[i]=Sum3(parentMomentum[i],incomingBathMomentumKgMS[i],-value.heavyMomentumKgMS[i]);
  }
  value.heavyKineticEnergyJ=Kinetic(value.heavyMomentumKgMS,heavyMassKg);
  const long double initial=Kinetic(parentMomentum,parentMassKg);
  value.signedReservoirEnergyJ=Sum4(initial,incomingBathEnergyJ,netQJ,-value.heavyKineticEnergyJ);
  value.signedReservoirChargeC=Sum3(parentChargeC,incomingBathChargeC,-heavyChargeC);
  if(!Finite(value.heavyMomentumKgMS)||!Finite(value.signedReservoirMomentumKgMS)||
     !Finite(value.heavyKineticEnergyJ)||!Finite(value.signedReservoirEnergyJ)||!Finite(value.signedReservoirChargeC)||!Finite(initial))
    return Refuse(Code::NumericalFailure,"inherit-velocity signed reservoir arithmetic overflow");
  *output=value;return Status{};
}

struct ThermalState {
  long double rateS,remainingPopulation,releasedPopulation;
  std::optional<long double> integratedHazard,logIntegratedHazard;
  bool rateUnderflowed,exposureUnderflowed,releasedUnderflowed;
};
inline Status ThermalRelease(long double attemptFrequencyS,long double bindingEnergyJ,long double boltzmannJK,
    long double temperatureK,long double initialPopulation,long double durationS,ThermalState* output) {
  if(!output||!Finite(attemptFrequencyS)||!Finite(bindingEnergyJ)||!Finite(boltzmannJK)||!Finite(temperatureK)||
      !Finite(initialPopulation)||!Finite(durationS)||attemptFrequencyS<0||bindingEnergyJ<0||boltzmannJK<=0||
      temperatureK<=0||initialPopulation<0||durationS<0)
    return Refuse(Code::InvalidInput,"thermal release requires explicit finite SI inputs, T/kB>0 and nonnegative nu/E/N/dt");
  ThermalState value{0,initialPopulation,0,0,std::nullopt,false,false,false};
  if(attemptFrequencyS==0) {*output=value;return Status{};}
  const long double thermalEnergy=boltzmannJK*temperatureK;
  long double exponent=0;
  if(bindingEnergyJ>0) {
    exponent=(Finite(thermalEnergy)&&thermalEnergy>0)?bindingEnergyJ/thermalEnergy:
       std::exp(std::log(bindingEnergyJ)-std::log(boltzmannJK)-std::log(temperatureK));
  }
  if(!Finite(exponent)) return Refuse(Code::NumericalFailure,"activation exponent outside numeric representation");
  // Do not evaluate nu*exp(-exponent): the exponential factor can underflow
  // even when multiplication by a large explicit nu would give a finite rate.
  // The exact zero-barrier branch also avoids exp(log(nu)) rounding overflow
  // when nu is at the largest representable value.
  value.rateS=exponent==0?attemptFrequencyS:std::exp(Sum3(std::log(attemptFrequencyS),-exponent,0));
  value.rateUnderflowed=value.rateS==0;
  if(durationS==0) {*output=value;return Status{};}
  const long double logHazard=Sum3(std::log(attemptFrequencyS),std::log(durationS),-exponent);
  if(!Finite(logHazard)) return Refuse(Code::NumericalFailure,"integrated activation exposure is not representable");
  value.logIntegratedHazard=logHazard;
  const long double hazard=std::exp(logHazard);
  if(!Finite(hazard)) {
    // A hazard above the floating range exhausts every representable finite
    // population. Nullopt honestly records the unrepresentable hazard itself.
    value.integratedHazard=std::nullopt;value.remainingPopulation=0;value.releasedPopulation=initialPopulation;
  } else if(hazard==0) {
    value.integratedHazard=0;value.exposureUnderflowed=true;
    // exp(logHazard) can underflow while N0*exposure is representable. In this
    // regime 1-exp(-x)=x to much better than machine relative precision; retain
    // that analytic limit in log space, rather than silently losing a stock.
    value.releasedPopulation=initialPopulation==0?0:std::exp(std::log(initialPopulation)+logHazard);
    value.remainingPopulation=initialPopulation-value.releasedPopulation;
  } else {
    value.integratedHazard=hazard;
    value.remainingPopulation=initialPopulation*std::exp(-hazard);
    value.releasedPopulation=initialPopulation*(-std::expm1(-hazard));
  }
  value.releasedUnderflowed=initialPopulation>0&&value.releasedPopulation==0;
  if(!Finite(value.rateS)||!Finite(value.remainingPopulation)||!Finite(value.releasedPopulation)||
      value.remainingPopulation<0||value.releasedPopulation<0||value.releasedPopulation>initialPopulation)
    return Refuse(Code::NumericalFailure,"thermal finite-inventory arithmetic is not representable");
  *output=value;return Status{};
}

} } }
#endif
