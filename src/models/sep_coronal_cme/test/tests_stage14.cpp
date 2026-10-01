#include "test_framework.h"
#include "sep_coronal_cme/research_extensions.h"
#include <cmath>

namespace SCCMTest {
using namespace SEP::CoronalCME;
using namespace SEP::CoronalCME::Research;
namespace {
void close(double a,double b,double tol,const std::string& message) {
  Require(std::abs(a-b)<=tol*std::max(1.0,std::max(std::abs(a),std::abs(b))),message);
}
DynamicGeometryState state(double t) {
  DynamicGeometryState s; s.centerM={10+2*t,3*t,0}; s.centerVelocityMPerS={2,3,0};
  s.rotationAxis={0,0,1}; s.angleRad=0.2*t; s.angularRateRadPerS=0.2;
  s.axesM={3+t,2+0.2*t,1+0.1*t}; s.axisRatesMPerS={1,0.2,0.1}; s.epochS=t; s.frameId="manufactured-inertial"; return s;
}
}
void RegisterStage14(Registry* tests) {
  (*tests)["MFP3D08"]=[] {
    ForeshockProxy p; p.identity={6,"foreshock-distance-proxy","compact-c2-v1","coefficients",1,7,true};
    p.referenceDistanceM=2; p.supportWidthM=10; p.minimumFactor=0.2;
    auto first=EvaluateForeshockProxy(p,2,7,100); Require(first.ok(),"proxy start"); close(first.value.factor,0.2,1e-14,"minimum factor");
    for (int j=0;j<=100;j++) {
      const double d=2+0.1*j;
      // A plane front x=10: volume geometry obtains its signed distance from
      // position/normal, while a line crossing obtains it from arclength.
      const Vec3 physical{10+d,3,-2},front{10,3,-2};
      const double volumeDistance=Dot(physical-front,{1,0,0});
      const double lineDistance=Norm(physical-front);
      auto a=EvaluateForeshockProxy(p,volumeDistance,7,100),b=EvaluateForeshockProxy(p,lineDistance,7,100);
      Require(a.ok()&&b.ok()&&a.value.factor>0&&a.value.factor<=1,"positive bounded reduction");
      Require(a.value.meanFreePathM==b.value.meanFreePathM,"3-D/line physical point parity");
    }
    close(EvaluateForeshockProxy(p,12,7,100).value.meanFreePathM,100,1e-14,"far recovery");
    Require(!EvaluateForeshockProxy(p,1,7,100).ok(),"no retrospective unresolved-layer scattering");
    Require(!EvaluateForeshockProxy(p,2,8,100).ok(),"stale front generation rejected");
    p.minimumFactor=1.1; Require(!EvaluateForeshockProxy(p,2,7,100).ok(),"above-one reduction rejected");
    p.minimumFactor=0; Require(!EvaluateForeshockProxy(p,2,7,100).ok(),"zero reduction rejected");
    p.identity.enabled=false; Require(EvaluateForeshockProxy(p,1,7,100).value.meanFreePathM==100,"disabled bitwise recovery");
  };
  (*tests)["ELL3D11"]=[] {
    const double t=0.5,h=1e-5; auto center=DynamicEllipsoid::Build(state(t)),plus=DynamicEllipsoid::Build(state(t+h)),minus=DynamicEllipsoid::Build(state(t-h));
    Require(center.ok()&&plus.ok()&&minus.ok(),"history constructed");
    const double theta=1.1,phi=0.7; Vec3 x=center.value.Point(theta,phi);
    const Vec3 numerical=(plus.value.Point(theta,phi)-minus.value.Point(theta,phi))/(2*h);
    Require(Norm(numerical-center.value.SurfaceVelocity(theta,phi))<1e-8,"independent position derivative includes center deflection and rotation");
    auto evaluation=center.value.Evaluate(x); Require(evaluation.ok(),"dynamic normal");
    const double ft=(plus.value.Evaluate(x).value.implicitValue-minus.value.Evaluate(x).value.implicitValue)/(2*h);
    const double dx=1e-5;
    const double gradient=(center.value.Evaluate(x+dx*evaluation.value.outwardNormal).value.implicitValue-
                           center.value.Evaluate(x-dx*evaluation.value.outwardNormal).value.implicitValue)/(2*dx);
    close(evaluation.value.normalSpeedMPerS,-ft/gradient,1e-8,"Q_dot normal-speed reference");
    const Vec3 e1=center.value.Rotate({1,0,0}),e2=center.value.Rotate({0,1,0}),e3=center.value.Rotate({0,0,1});
    close(Dot(e1,Cross(e2,e3)),1,1e-14,"proper rotation determinant"); close(Dot(e1,e2),0,1e-14,"orthogonal attitude");
    auto stationary=state(0); stationary.centerVelocityMPerS={}; stationary.angularRateRadPerS=0; stationary.axisRatesMPerS={0,0,0};
    auto fixed=DynamicEllipsoid::Build(stationary); Require(Norm(fixed.value.SurfaceVelocity(theta,phi))==0,"fixed history exact reduction");
    // Noncommuting rotations have a changing inertial angular-velocity axis.
    // This independently differentiated matrix history exercises the general
    // SO(3) authority, not the fixed-axis Rodrigues reduction.
    auto general=[](double time){auto s=state(time);s.fullAttitude=true;
      const double z=0.2*time,x=0.3*time,cz=std::cos(z),sz=std::sin(z),cx=std::cos(x),sx=std::sin(x);
      s.bodyToInertial={cz,-sz*cx,sz*sx,sz,cz*cx,-cz*sx,0,sx,cx};
      s.inertialAngularVelocityPerS={0.3*cz,0.3*sz,0.2};return s;};
    auto g=DynamicEllipsoid::Build(general(t)),gp=DynamicEllipsoid::Build(general(t+h)),gm=DynamicEllipsoid::Build(general(t-h));
    Require(g.ok()&&gp.ok()&&gm.ok(),"general attitude valid");
    Require(Norm((gp.value.Point(theta,phi)-gm.value.Point(theta,phi))/(2*h)-g.value.SurfaceVelocity(theta,phi))<1e-8,"noncommuting attitude derivative");
    auto reflected=general(t);reflected.bodyToInertial[0]*=-1;
    Require(!DynamicEllipsoid::Build(reflected).ok(),"improper/nonorthogonal attitude rejected");
  };
  (*tests)["SRC3D21"]=[] {
    ReferenceStratum a={"proton",1,2,1,4,0,10},b=a; b.momentumLowerSi=2; b.momentumUpperSi=3;
    auto geometry=[](double offset){return SEP::Core::Result<double>::Success(2*(1+offset));};
    auto x=SolveReferenceFamilyMember(a,2,10,1e-10,[](double){return 0.5;},geometry);
    auto y=SolveReferenceFamilyMember(b,2,10,1e-10,[](double){return 1.0;},geometry);
    Require(x.ok()&&y.ok(),"positive continuous root family"); close(x.value.offsetM,4,1e-9,"constant-coefficient independent root"); close(y.value.offsetM,2,1e-9,"momentum-dependent placement");
    Require(ValidateReferenceFamily({x.value,y.value}).ok(),"joint family contiguous support");
    y.value.stratum.momentumLowerSi=2.1; Require(!ValidateReferenceFamily({x.value,y.value}).ok(),"gap rejected");
    auto badGeometry=[](double){return SEP::Core::Result<double>::Failure(SEP::Core::StatusCode::OutOfDomain,"fold/overlap/mask/clearance");};
    Require(!SolveReferenceFamilyMember(a,2,10,1e-8,[](double){return 1.0;},badGeometry).ok(),"uncertified geometry rejects allocation");
    Require(!SolveReferenceFamilyMember(a,2,1,1e-8,[](double){return 1.0;},geometry).ok(),"root outside support rejected");
    Require(!SolveReferenceFamilyMember(a,2,10,1e-8,[](double){return -1.0;},geometry).ok(),"sign certification");
  };
  (*tests)["SRC3D22"]=[] {
    constexpr double c=299792458.0; const double mass=1; Vec3 p={2*c,0,0};
    RenewalBranch absorbed; absorbed.outcome=RenewalOutcome::Absorbed; absorbed.probability=0.3; absorbed.outgoingMomentumKgMPerS=p;
    RenewalBranch down=absorbed; down.outcome=RenewalOutcome::Downstream; down.probability=0.2;
    RenewalBranch renewed=absorbed; renewed.outcome=RenewalOutcome::ReReleased; renewed.probability=0.5; renewed.residenceTimeS=3;
    renewed.outgoingMomentumKgMPerS=2*p; renewed.declaredShockImpulseKgMPerS=p;
    renewed.declaredShockWorkJ=std::hypot(mass*c*c,Norm(2*p)*c)-std::hypot(mass*c*c,Norm(p)*c);
    auto result=EvaluateRenewal(10,mass,p,{absorbed,down,renewed},true); Require(result.ok(),"normalized conditional transition");
    close(result.value.incomingNumber,result.value.absorbedNumber+result.value.downstreamNumber+result.value.reReleasedNumber,1e-14,"number closure");
    close(result.value.outgoingEnergyJ-result.value.incomingEnergyJ,result.value.shockWorkJ,1e-14,"energy/shock-work closure");
    Require(Norm(result.value.outgoingMomentum-result.value.incomingMomentum-result.value.shockImpulse)<1e-5,"momentum ledger closure");
    renewed.residenceTimeS=0; Require(!EvaluateRenewal(10,mass,p,{absorbed,down,renewed},true).ok(),"zero-delay re-emission rejected");
    Require(!EvaluateRenewal(10,mass,p,{absorbed,down},true).ok(),"unnormalized kernel rejected");
    Require(!EvaluateRenewal(10,mass,p,{absorbed},false).ok(),"missing work authority rejected");
  };
  (*tests)["SRC3D20"]=[] {
    ImpulsiveLedger ledger; ImpulsiveBirth b={1,0,10,3,ImpulsiveFrontRelation::NoFrontYet,"ImpulsiveCoronalRelease"};
    Require(ledger.Birth(b).ok(),"pre-front impulsive origin"); Require(!ledger.Birth(b).ok(),"duplicate physical birth rejected");
    Require(ledger.AbsorbFrontEncounter(1,7).ok(),"front overtake absorbed"); Require(!ledger.AbsorbFrontEncounter(1,7).ok(),"exactly once encounter");
    Require(ledger.BornNumber()==10&&ledger.FrontEncounterNumber()==10&&ledger.ShockFirstPassageNumber()==0,"origin/interaction/shock ledgers disjoint");
    Require(ledger.BornKineticEnergyJ()==30&&ledger.FrontEncounterKineticEnergyJ()==30,"disjoint physical kinetic-energy census");
    b.immutableBirthId=2; b.frontRelation=ImpulsiveFrontRelation::Downstream; Require(!ledger.Birth(b).ok(),"unsupported downstream birth rejected");
    b.frontRelation=ImpulsiveFrontRelation::Upstream; Require(!ledger.Birth(b).ok(),"front generation required");
  };
}
}
