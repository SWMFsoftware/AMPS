#include "sep_corona_swcme/sheath_model.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>

namespace SEP { namespace CoronaSwcme {
namespace {

using CoronalCME::Cross;
using CoronalCME::Dot;
using CoronalCME::FixedOrientationEllipsoid;
using CoronalCME::MhdPrimitiveState;
using CoronalCME::Norm;
using CoronalCME::Vec3;

constexpr double kPi=3.141592653589793238462643383279502884;
constexpr double kMu0=1.25663706212e-6;

struct BirthState {
  FixedOrientationEllipsoid front;
  Vec3 position;
  Vec3 surfaceVelocity;
  Vec3 normal;
  Vec3 deficit;
  MhdPrimitiveState upstream;
  CoronalCME::MhdShockSolution jump;
  double normalSpeed=0.0;
  double w1=0.0,w2=0.0;
};

bool Finite(Vec3 value) {
  return std::isfinite(value.x)&&std::isfinite(value.y)&&std::isfinite(value.z);
}

Vec3 PointPolarDerivative(const FixedOrientationEllipsoid& shape,
    double polar,double azimuth) {
  const auto& k=shape.Kinematics();const auto& b=shape.Basis();
  return -k.radialSemiAxisM.value*std::sin(polar)*b.radial+
      k.firstLateralSemiAxisM.value*std::cos(polar)*std::cos(azimuth)*b.firstLateral+
      k.secondLateralSemiAxisM.value*std::cos(polar)*std::sin(azimuth)*b.secondLateral;
}

Vec3 PointAzimuthDerivative(const FixedOrientationEllipsoid& shape,
    double polar,double azimuth) {
  const auto& k=shape.Kinematics();const auto& b=shape.Basis();
  return -k.firstLateralSemiAxisM.value*std::sin(polar)*std::sin(azimuth)*b.firstLateral+
      k.secondLateralSemiAxisM.value*std::sin(polar)*std::cos(azimuth)*b.secondLateral;
}

double Determinant(const std::array<Vec3,3>& columns) {
  return Dot(columns[0],Cross(columns[1],columns[2]));
}

Vec3 ApplyMap(const std::array<Vec3,3>& current,
    const std::array<Vec3,3>& birth,Vec3 value) {
  const double determinant=Determinant(birth);
  const double c0=Dot(value,Cross(birth[1],birth[2]))/determinant;
  const double c1=Dot(value,Cross(birth[2],birth[0]))/determinant;
  const double c2=Dot(value,Cross(birth[0],birth[1]))/determinant;
  return c0*current[0]+c1*current[1]+c2*current[2];
}

std::array<Vec3,3> CartesianMapColumns(const std::array<Vec3,3>& current,
    const std::array<Vec3,3>& birth) {
  return {ApplyMap(current,birth,{1,0,0}),ApplyMap(current,birth,{0,1,0}),
      ApplyMap(current,birth,{0,0,1})};
}

Core::Result<BirthState> EvaluateBirth(const EventConfiguration& event,
    const AmbientModel& ambient,double polar,double azimuth,double time) {
  using Return=Core::Result<BirthState>;
  const auto evolution=event.At(time);
  if(!evolution.ok())return Return::Failure(evolution.status.code,evolution.status.message);
  const auto shape=FixedOrientationEllipsoid::FromCenter(event.basis,
      evolution.value.ellipsoid,event.support.solarRadiusM);
  if(!shape.ok())return Return::Failure(shape.status.code,shape.status.message);
  BirthState result;result.front=shape.value;
  result.position=shape.value.Point(polar,azimuth);
  result.surfaceVelocity=shape.value.SurfaceVelocity(polar,azimuth);
  const auto surface=shape.value.Evaluate(result.position);
  if(!surface.ok())return Return::Failure(surface.status.code,surface.status.message);
  result.normal=surface.value.outwardNormal;result.normalSpeed=surface.value.normalSpeedMPerS;
  const auto upstream=ambient.Evaluate(result.position,time);
  if(!upstream.ok())return Return::Failure(upstream.status.code,upstream.status.message);
  result.upstream={upstream.value.plasma.massDensityKgM3,
      upstream.value.plasma.pressurePa,upstream.value.velocityMPerS,
      upstream.value.magneticFieldT};
  const auto characteristic=CoronalCME::EvaluateMhdCharacteristics(result.upstream,
      result.normal,event.composition.gammaAdiabatic);
  if(!characteristic.ok())return Return::Failure(characteristic.status.code,
      characteristic.status.message);
  result.w1=result.normalSpeed-Dot(result.upstream.velocityMPerS,result.normal);
  if(!(result.w1>characteristic.value.fastSpeedMPerS))return Return::Failure(
      Core::StatusCode::UnsupportedCapability,"material label is on a sub-fast front patch");
  const auto jump=CoronalCME::SolveObliqueFastShock(result.upstream,result.normal,
      result.normalSpeed,event.composition.gammaAdiabatic);
  if(!jump.ok())return Return::Failure(jump.status.code,jump.status.message);
  result.jump=jump.value;
  result.w2=result.normalSpeed-Dot(result.jump.downstream.velocityMPerS,result.normal);
  if(!(result.w2>0&&std::isfinite(result.w2)))return Return::Failure(
      Core::StatusCode::InvalidState,"RH downstream does not enter the sheath");
  result.deficit=result.jump.downstream.velocityMPerS-result.surfaceVelocity;
  return Return::Success(std::move(result));
}

template<class Sampler>
Core::Result<Vec3> OneSidedDerivative(double x,double step,double lower,double upper,
    const Vec3& center,const Sampler& sampler) {
  using Return=Core::Result<Vec3>;
  Core::Result<Vec3> minus=Return::Failure(Core::StatusCode::OutOfDomain,"minus unavailable");
  Core::Result<Vec3> plus=Return::Failure(Core::StatusCode::OutOfDomain,"plus unavailable");
  Core::Result<Vec3> minus2=Return::Failure(Core::StatusCode::OutOfDomain,"minus2 unavailable");
  Core::Result<Vec3> plus2=Return::Failure(Core::StatusCode::OutOfDomain,"plus2 unavailable");
  if(x-step>=lower)minus=sampler(x-step);
  if(x+step<=upper)plus=sampler(x+step);
  if(x-2*step>=lower)minus2=sampler(x-2*step);
  if(x+2*step<=upper)plus2=sampler(x+2*step);
  // A fourth-order centered stencil materially reduces derivative noise in F
  // while retaining one-sided stencils at coverage/fast-support boundaries.
  if(minus2.ok()&&minus.ok()&&plus.ok()&&plus2.ok())return Return::Success(
      (minus2.value-8*minus.value+8*plus.value-plus2.value)/(12*step));
  if(minus.ok()&&plus.ok())return Return::Success((plus.value-minus.value)/(2*step));
  if(plus.ok()&&plus2.ok())return Return::Success(
      (-3*center+4*plus.value-plus2.value)/(2*step));
  if(minus.ok()&&minus2.ok())return Return::Success(
      (3*center-4*minus.value+minus2.value)/(2*step));
  if(plus.ok())return Return::Success((plus.value-center)/step);
  if(minus.ok())return Return::Success((center-minus.value)/step);
  return Return::Failure(Core::StatusCode::UnsupportedCapability,
      "shock-fed derivative has no same-fast-branch stencil");
}

double Relative(double a,double b) {
  return std::abs(a-b)/std::max({std::abs(a),std::abs(b),1e-300});
}

// Diagnostic closure, not the selected production law.  Post-shock material
// retains the exact RH relative drift at age zero, then relaxes smoothly
// toward a nonzero frozen fraction.  L(0)=0 and L'(0)=1 preserve the shock
// position and velocity limits; L'>=kappa>0 avoids the singular old-cohort
// volume produced by a fully arrested drift.  Nothing in mass conservation or
// the Cauchy induction identity selects kappa or T.  Their implied acceleration
// is the large, measured sustaining-force residual of this empirical map, so
// changing them to reduce a test residual would be model fitting rather than
// numerical convergence.  The reviewed production replacement instead makes
// a single ejecta contact the piston boundary and solves radial momentum and
// energy on finite ambient mass cells; see docs/BG3D4_PISTON_CLOSURE.md.
double DriftTime(const EventConfiguration& event,double ageS) {
  const double kappa=event.regional.sheathDriftAsymptoteFraction;
  const double relaxation=event.regional.sheathDriftRelaxationTimeS;
  return kappa*ageS+(1-kappa)*relaxation*(-std::expm1(-ageS/relaxation));
}

double DriftRate(const EventConfiguration& event,double ageS) {
  const double kappa=event.regional.sheathDriftAsymptoteFraction;
  return kappa+(1-kappa)*std::exp(-ageS/
      event.regional.sheathDriftRelaxationTimeS);
}

Core::Result<Vec3> EvaluateMapPosition(const EventConfiguration& event,
    const AmbientModel& ambient,const SheathMaterialLabel& label,double time) {
  // This intentionally evaluates the RH deficit at the query epoch.  That is
  // why the diagnostic cannot evolve old material after the prescribed front
  // becomes sub-fast: no admissible current downstream state exists.  It would
  // be incorrect to fall back to the birth deficit silently because that is a
  // different evolution law.  The production piston model separates ongoing
  // material evolution from shock admission and therefore does not use this
  // helper.
  const auto flow=EvaluateBirth(event,ambient,label.polarRad,label.azimuthRad,time);
  if(!flow.ok())return Core::Result<Vec3>::Failure(
      flow.status.code,flow.status.message);
  return Core::Result<Vec3>::Success(flow.value.position+
      DriftTime(event,time-label.crossingTimeS)*flow.value.deficit);
}

struct PositionClassification {
  Vec3 position;
  double solarSigned=0.0;   // positive outside the physical Sun
  double outerSigned=0.0;   // positive inside the shock front
};

Core::Result<PositionClassification> ClassifyPosition(
    const EventConfiguration& event,const AmbientModel& ambient,
    const SheathMaterialLabel& label,double time) {
  using Return=Core::Result<PositionClassification>;
  const auto flow=EvaluateBirth(event,ambient,label.polarRad,label.azimuthRad,time);
  if(!flow.ok())return Return::Failure(flow.status.code,flow.status.message);
  const auto position=EvaluateMapPosition(event,ambient,label,time);
  if(!position.ok())return Return::Failure(position.status.code,position.status.message);
  PositionClassification result;
  result.position=position.value;
  const auto frontValue=flow.value.front.Evaluate(result.position);
  if(!frontValue.ok())return Return::Failure(
      Core::StatusCode::InvalidState,"cannot evaluate a sheath boundary event function");
  result.solarSigned=Norm(result.position)-event.support.solarRadiusM;
  result.outerSigned=-frontValue.value.implicitValue;
  return Return::Success(result);
}

struct LocatedDisposition {
  SheathCellDisposition disposition=SheathCellDisposition::Retained;
  double time=0.0;
};

Core::Result<LocatedDisposition> LocateFirstDisposition(
    const EventConfiguration& event,const AmbientModel& ambient,
    const SheathMaterialLabel& label,double epoch) {
  using Return=Core::Result<LocatedDisposition>;
  // These are geometric event functions, independent of F/J.  Locating an
  // exit before evaluating deformation prevents an already-removed parcel
  // from being extrapolated into a later nonphysical ballistic fold.
  const int intervals=128;
  const double begin=label.crossingTimeS;
  const double epsilon=std::max(event.support.rootToleranceS,
      1e-9*std::max(1.0,epoch-begin));
  double previousTime=std::min(epoch,begin+epsilon);
  auto previous=ClassifyPosition(event,ambient,label,previousTime);
  if(!previous.ok())return Return::Failure(previous.status.code,previous.status.message);
  auto scalar=[](const PositionClassification& state,
      SheathCellDisposition kind) {
    if(kind==SheathCellDisposition::SolarExit)return state.solarSigned;
    return state.outerSigned;
  };
  for(int i=1;i<=intervals;++i) {
    const double time=begin+(epoch-begin)*i/intervals;
    if(time<=previousTime)continue;
    const auto current=ClassifyPosition(event,ambient,label,time);
    if(!current.ok())return Return::Failure(current.status.code,current.status.message);
    for(const auto kind:{SheathCellDisposition::SolarExit,
        SheathCellDisposition::OuterExit}) {
      if(scalar(previous.value,kind)>0&&scalar(current.value,kind)<=0) {
        double lower=previousTime,upper=time;
        for(int iteration=0;iteration<80&&upper-lower>event.support.rootToleranceS;
            ++iteration) {
          const double middle=0.5*(lower+upper);
          const auto state=ClassifyPosition(event,ambient,label,middle);
          if(!state.ok())return Return::Failure(state.status.code,state.status.message);
          if(scalar(state.value,kind)>0)lower=middle;else upper=middle;
        }
        return Return::Success({kind,0.5*(lower+upper)});
      }
    }
    previousTime=time;previous=current;
  }
  return Return::Success({SheathCellDisposition::Retained,epoch});
}

} // namespace

Core::Result<std::shared_ptr<ShockFedSheathModel>> ShockFedSheathModel::Create(
    std::shared_ptr<const EventConfiguration> event,
    std::shared_ptr<const AmbientModel> ambient) {
  using Return=Core::Result<std::shared_ptr<ShockFedSheathModel>>;
  if(!event||!ambient)return Return::Failure(Core::StatusCode::InvalidConfiguration,
      "shock-fed sheath requires event and ambient authorities");
  if(event->physicsFingerprint!=ambient->Identity())
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "sheath and ambient authorities have different event identities");
  // Keep the selector closed.  The future per-ray-lagrangian-piston value must
  // construct a different type with a finite ambient inventory, a checksummed
  // ray quadrature and a contact/ejecta driver.  Accepting it here would falsely
  // advertise the empirical front-driven map under a production identity.
  if(event->sheathModel!="rh-relaxing-material-map-v1")return Return::Failure(
      Core::StatusCode::UnsupportedCapability,"unsupported analytical sheath model");
  auto model=std::shared_ptr<ShockFedSheathModel>(new ShockFedSheathModel);
  model->event_=std::move(event);model->ambient_=std::move(ambient);
  return Return::Success(std::move(model));
}

Core::Result<SheathMappedState> ShockFedSheathModel::Evaluate(
    const SheathMaterialLabel& label,double epochS) const {
  using Return=Core::Result<SheathMappedState>;
  if(!(std::isfinite(label.polarRad)&&std::isfinite(label.azimuthRad)&&
      std::isfinite(label.crossingTimeS)&&std::isfinite(epochS)&&
      label.polarRad>0&&label.polarRad<kPi&&label.azimuthRad>=0&&
      label.azimuthRad<2*kPi&&label.crossingTimeS>=event_->regional.sheathAdmissionStartS&&
      epochS>=label.crossingTimeS&&epochS<=event_->support.endS))return Return::Failure(
      Core::StatusCode::OutOfDomain,"invalid sheath material label or epoch coverage");
  const auto birth=EvaluateBirth(*event_,*ambient_,label.polarRad,label.azimuthRad,
      label.crossingTimeS);
  if(!birth.ok())return Return::Failure(birth.status.code,birth.status.message);
  const auto nowKinematics=event_->At(epochS);
  if(!nowKinematics.ok())return Return::Failure(nowKinematics.status.code,
      nowKinematics.status.message);
  const auto nowShape=FixedOrientationEllipsoid::FromCenter(event_->basis,
      nowKinematics.value.ellipsoid,event_->support.solarRadiusM);
  if(!nowShape.ok())return Return::Failure(nowShape.status.code,nowShape.status.message);
  const auto currentFlow=epochS==label.crossingTimeS?birth:
      EvaluateBirth(*event_,*ambient_,label.polarRad,label.azimuthRad,epochS);
  if(!currentFlow.ok())return Return::Failure(currentFlow.status.code,
      "current sheath-driving state unavailable: "+currentFlow.status.message);
  const double age=epochS-label.crossingTimeS;
  const double driftTime=DriftTime(*event_,age);
  const double driftRate=DriftRate(*event_,age);
  SheathMappedState result;result.label=label;result.eventIdentity=event_->physicsFingerprint;
  result.upstreamRelativeNormalSpeedMPerS=birth.value.w1;
  result.downstreamRelativeNormalSpeedMPerS=birth.value.w2;
  result.shockResidual=birth.value.jump.residuals.maximum;
  result.birthAreaDensityM2PerRad2=Norm(Cross(
      PointPolarDerivative(birth.value.front,label.polarRad,label.azimuthRad),
      PointAzimuthDerivative(birth.value.front,label.polarRad,label.azimuthRad)));
  const auto mappedPosition=EvaluateMapPosition(*event_,*ambient_,label,epochS);
  if(!mappedPosition.ok())return Return::Failure(mappedPosition.status.code,
      mappedPosition.status.message);
  result.positionM=mappedPosition.value;
  // This is the centre basis of a finite physical reference prism
  // X0(theta,phi,s)=X_sh(theta,phi,tau)-s*deficit(theta,phi,tau).
  // Its determinant is dA_sh/dtheta/dphi times the downstream shock-frame
  // crossing speed w2.  It supplies a real 3-D reference volume; the angular
  // and admission-time labels are not treated as Cartesian identity axes.
  const std::array<Vec3,3> birthBasis={
      PointPolarDerivative(birth.value.front,label.polarRad,label.azimuthRad),
      PointAzimuthDerivative(birth.value.front,label.polarRad,label.azimuthRad),
      -1.0*birth.value.deficit};
  const double birthDet=Determinant(birthBasis);
  const double birthScale=std::max(1.0,
      Norm(birthBasis[0])*Norm(birthBasis[1])*Norm(birthBasis[2]));
  if(!(std::isfinite(birthDet)&&birthDet>1e-14*birthScale))return Return::Failure(
      Core::StatusCode::InvalidState,
      "shock-fed cohort reference prism is singular or reversed");
  result.referenceDerivativeColumns=birthBasis;
  result.referenceVolumeDensityM3PerRad2S=birthDet;
  if(age<=event_->support.rootToleranceS*1e-6) {
    result.jacobian=1;result.deformationColumns={Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}};
    result.labelDerivativeColumns=birthBasis;
    result.currentVolumeDensityM3PerRad2S=birthDet;
    result.primitive=birth.value.jump.downstream;
    return Return::Success(std::move(result));
  }

  const double angularStep=5e-4;
  auto polarSample=[&](double value) {
    const auto sample=EvaluateBirth(*event_,*ambient_,value,label.azimuthRad,
        epochS);
    return sample.ok()?Core::Result<Vec3>::Success(sample.value.deficit):
        Core::Result<Vec3>::Failure(sample.status.code,sample.status.message);
  };
  const auto dPolar=OneSidedDerivative(label.polarRad,angularStep,1e-8,kPi-1e-8,
      currentFlow.value.deficit,polarSample);
  auto azimuthSample=[&](double value) {
    value=std::fmod(value+2*kPi,2*kPi);
    const auto sample=EvaluateBirth(*event_,*ambient_,label.polarRad,value,
        epochS);
    return sample.ok()?Core::Result<Vec3>::Success(sample.value.deficit):
        Core::Result<Vec3>::Failure(sample.status.code,sample.status.message);
  };
  const auto azMinus=azimuthSample(label.azimuthRad-angularStep);
  const auto azPlus=azimuthSample(label.azimuthRad+angularStep);
  const auto azMinus2=azimuthSample(label.azimuthRad-2*angularStep);
  const auto azPlus2=azimuthSample(label.azimuthRad+2*angularStep);
  if(!dPolar.ok()||(!azMinus.ok()&&!azPlus.ok()))return Return::Failure(
      Core::StatusCode::UnsupportedCapability,
      "shock-fed angular derivative crosses unsupported fast boundary");
  const Vec3 dAzimuth=azMinus2.ok()&&azMinus.ok()&&azPlus.ok()&&azPlus2.ok()?
      (azMinus2.value-8*azMinus.value+8*azPlus.value-azPlus2.value)/(12*angularStep):
      azMinus.ok()&&azPlus.ok()?(azPlus.value-azMinus.value)/(2*angularStep):
      (azPlus.ok()?(azPlus.value-currentFlow.value.deficit)/angularStep:
          (currentFlow.value.deficit-azMinus.value)/angularStep);
  // The RH deficit is obtained through several nonlinear characteristic and
  // jump operations.  A 0.5-s fourth-order stencil avoids subtractive loss at
  // heliospheric times while remaining far below the 100-s handoff scale.
  // Where both centered stencils fit, Richardson extrapolation of h and h/2
  // removes their leading fourth-order term.  Near coverage endpoints the
  // one-sided h/2 result is retained rather than assuming centered support.
  const double timeStep=std::max(0.5,500*event_->support.rootToleranceS);
  auto positionSample=[&](double value) {
    return EvaluateMapPosition(*event_,*ambient_,label,value);
  };
  const double timeLower=std::max(event_->support.startS,label.crossingTimeS);
  const auto velocityCoarse=OneSidedDerivative(epochS,timeStep,timeLower,
      event_->support.endS,result.positionM,positionSample);
  const auto velocityFine=OneSidedDerivative(epochS,0.5*timeStep,timeLower,
      event_->support.endS,result.positionM,positionSample);
  if(!velocityCoarse.ok()||!velocityFine.ok())return Return::Failure(
      !velocityCoarse.ok()?velocityCoarse.status.code:velocityFine.status.code,
      !velocityCoarse.ok()?velocityCoarse.status.message:velocityFine.status.message);
  result.primitive.velocityMPerS=velocityFine.value;
  if(epochS-2*timeStep>=timeLower&&
      epochS+2*timeStep<=event_->support.endS)
    result.primitive.velocityMPerS=(16*velocityFine.value-velocityCoarse.value)/15;

  const std::array<Vec3,3> currentBasis={
      PointPolarDerivative(nowShape.value,label.polarRad,label.azimuthRad)+
          driftTime*dPolar.value,
      PointAzimuthDerivative(nowShape.value,label.polarRad,label.azimuthRad)+
          driftTime*dAzimuth,
      -driftRate*currentFlow.value.deficit};
  const double currentDet=Determinant(currentBasis);
  if(!std::isfinite(currentDet))return Return::Failure(
      Core::StatusCode::NumericalFailure,"shock-fed current volume metric is nonfinite");
  result.jacobian=currentDet/birthDet;
  if(!(std::isfinite(result.jacobian)&&result.jacobian>=event_->regional.minimumJacobian))
    return Return::Failure(Core::StatusCode::InvalidState,
        "shock-fed relaxing map folded or crossed its minimum Jacobian");
  result.deformationColumns=CartesianMapColumns(currentBasis,birthBasis);
  result.labelDerivativeColumns=currentBasis;
  result.currentVolumeDensityM3PerRad2S=currentDet;
  result.primitive.massDensityKgM3=birth.value.jump.downstream.massDensityKgM3/result.jacobian;
  result.primitive.pressurePa=birth.value.jump.downstream.pressurePa*
      std::pow(result.jacobian,-event_->composition.gammaAdiabatic);
  result.primitive.magneticFieldT=ApplyMap(currentBasis,birthBasis,
      birth.value.jump.downstream.magneticFieldT)/result.jacobian;
  if(!(result.primitive.massDensityKgM3>0&&result.primitive.pressurePa>0&&
      Finite(result.positionM)&&Finite(result.primitive.velocityMPerS)&&
      Finite(result.primitive.magneticFieldT)))return Return::Failure(
      Core::StatusCode::NumericalFailure,"shock-fed material state is nonfinite");
  return Return::Success(std::move(result));
}

Core::Result<SheathContactState> ShockFedSheathModel::EvaluateContact(
    double polarRad,double azimuthRad,double epochS) const {
  using Return=Core::Result<SheathContactState>;
  const double start=std::max(event_->support.startS,
      event_->regional.sheathAdmissionStartS);
  if(!(std::isfinite(polarRad)&&std::isfinite(azimuthRad)&&
      std::isfinite(epochS)&&epochS>=start&&epochS<=event_->support.endS))
    return Return::Failure(Core::StatusCode::OutOfDomain,
        "material contact query is outside event support");
  // Diagnostic-map rear boundary only.  One global oldest cohort is selected
  // so inventory cannot reset at each local first-fast time.  It is material
  // under this map, but it is not the production contact authority: BG3D-3's
  // fixed-fraction surface disagrees with it and the replacement model uses
  // the reconstructed ejecta surface directly.  Failure at the global start
  // is therefore an unsupported diagnostic patch, never permission to choose
  // a later cohort and erase the preceding physical inventory.
  SheathMaterialLabel label{polarRad,azimuthRad,start,0};
  const auto mapped=Evaluate(label,epochS);
  if(!mapped.ok())return Return::Failure(mapped.status.code,mapped.status.message);
  Vec3 areaVector=Cross(mapped.value.labelDerivativeColumns[0],
      mapped.value.labelDerivativeColumns[1]);
  const double area=Norm(areaVector);
  if(!(std::isfinite(area)&&area>0))return Return::Failure(
      Core::StatusCode::InvalidState,"material contact has a singular surface metric");
  Vec3 normal=areaVector/area;
  const auto eventState=event_->At(epochS);
  if(!eventState.ok())return Return::Failure(eventState.status.code,
      eventState.status.message);
  const Vec3 center=eventState.value.ellipsoid.centerDistanceM.value*
      event_->basis.radial;
  if(Dot(normal,mapped.value.positionM-center)<0)
    normal=-1.0*normal;
  SheathContactState result;
  result.label=label;result.positionM=mapped.value.positionM;
  result.outwardNormal=normal;result.velocityMPerS=mapped.value.primitive.velocityMPerS;
  result.primitive=mapped.value.primitive;result.areaDensityM2PerRad2=area;
  result.normalSpeedMPerS=Dot(result.velocityMPerS,normal);
  // The contact is the image of a fixed material-label surface, so its
  // parameterization velocity is exactly U.  Keeping the dimensional value
  // explicit allows independent finite-difference boundary-velocity tests.
  result.relativeMassFluxKgM2S=result.primitive.massDensityKgM3*
      Dot(result.primitive.velocityMPerS-result.velocityMPerS,normal);
  result.normalMagneticFieldT=Dot(result.primitive.magneticFieldT,normal);
  const double totalPressure=result.primitive.pressurePa+
      Dot(result.primitive.magneticFieldT,result.primitive.magneticFieldT)/(2*kMu0);
  result.tractionPa=totalPressure*normal-
      (result.normalMagneticFieldT/kMu0)*result.primitive.magneticFieldT;
  result.eventIdentity=event_->physicsFingerprint;
  return Return::Success(std::move(result));
}

Core::Result<SheathMappedState> ShockFedSheathModel::EvaluateAtPosition(
    Vec3 target,double epochS,double toleranceM,int maximumIterations) const {
  using Return=Core::Result<SheathMappedState>;
  if(!current_||current_->epochS!=epochS||!Finite(target)||
      !(std::isfinite(toleranceM)&&toleranceM>0)||maximumIterations<1)
    return Return::Failure(Core::StatusCode::InvalidState,
        "spatial sheath query requires a matching committed inventory and tolerance");
  const SheathAdmissionCell* seed=nullptr;double best=std::numeric_limits<double>::infinity();
  for(const auto& cell:current_->cells)if(
      cell.disposition==SheathCellDisposition::Retained) {
    const double distance=Norm(cell.state.positionM-target);
    if(distance<best){best=distance;seed=&cell;}
  }
  if(seed==nullptr)return Return::Failure(Core::StatusCode::OutOfDomain,
      "committed sheath inventory contains no retained material support");
  SheathMaterialLabel label=seed->label;
  for(int iteration=0;iteration<maximumIterations;++iteration) {
    const auto state=Evaluate(label,epochS);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    const Vec3 residual=state.value.positionM-target;
    const double norm=Norm(residual);
    if(norm<=toleranceM)return state;
    const auto& g=state.value.labelDerivativeColumns;
    const double determinant=Determinant(g);
    const double scale=std::max(1.0,Norm(g[0])*Norm(g[1])*Norm(g[2]));
    if(!(std::isfinite(determinant)&&std::abs(determinant)>1e-14*scale))
      return Return::Failure(Core::StatusCode::NumericalFailure,
          "spatial sheath inverse encountered a singular material Jacobian");
    const Vec3 delta={Dot(residual,Cross(g[1],g[2]))/determinant,
        Dot(residual,Cross(g[2],g[0]))/determinant,
        Dot(residual,Cross(g[0],g[1]))/determinant};
    bool accepted=false;
    for(double damping=1;damping>=1.0/128;damping*=0.5) {
      SheathMaterialLabel trial=label;
      trial.polarRad=std::max(1e-8,std::min(kPi-1e-8,
          label.polarRad-damping*delta.x));
      trial.azimuthRad=std::fmod(label.azimuthRad-damping*delta.y+2*kPi,2*kPi);
      trial.crossingTimeS=std::max(event_->regional.sheathAdmissionStartS,
          std::min(epochS,label.crossingTimeS-damping*delta.z));
      const auto candidate=Evaluate(trial,epochS);
      if(candidate.ok()&&Norm(candidate.value.positionM-target)<norm) {
        label=trial;accepted=true;break;
      }
    }
    if(!accepted)return Return::Failure(Core::StatusCode::NumericalFailure,
        "spatial sheath inverse cannot decrease its Cartesian residual");
  }
  return Return::Failure(Core::StatusCode::NumericalFailure,
      "spatial sheath inverse exceeded its iteration budget");
}

Core::Result<SheathMappedState> ShockFedSheathModel::QueryCommitted(
    const SheathMaterialLabel& label) const {
  using Return=Core::Result<SheathMappedState>;
  if(!current_)return Return::Failure(Core::StatusCode::InvalidState,
      "material query requires a committed sheath inventory");
  for(const auto& cell:current_->cells) {
    if(cell.label.shockPatchLineage==label.shockPatchLineage&&
        cell.label.crossingTimeS==label.crossingTimeS&&
        cell.label.polarRad==label.polarRad&&
        cell.label.azimuthRad==label.azimuthRad)
      return Return::Success(cell.state);
  }
  return Return::Failure(Core::StatusCode::OutOfDomain,
      "material label is absent from the committed sheath inventory");
}

Core::Result<std::shared_ptr<const SheathInventory>>
ShockFedSheathModel::PrepareInventory(double epochS,
    std::uint64_t backgroundGeneration,int polarCells,int azimuthCells,
    int admissionTimeCells) {
  using Return=Core::Result<std::shared_ptr<const SheathInventory>>;
  const double begin=std::max(event_->support.startS,
      event_->regional.sheathAdmissionStartS);
  if(!(backgroundGeneration>0&&polarCells>=2&&azimuthCells>=3&&
      admissionTimeCells>=1&&epochS>=begin&&epochS<=event_->support.endS))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "invalid sheath inventory epoch, generation, or resolution");
  auto candidate=std::shared_ptr<SheathInventory>(new SheathInventory);
  candidate->epochS=epochS;candidate->backgroundGeneration=backgroundGeneration;
  candidate->eventIdentity=event_->physicsFingerprint;
  candidate->minimumJacobian=std::numeric_limits<double>::infinity();
  candidate->zeroVolumeStartup=epochS==begin;
  if(epochS==begin) {
    // Diagnostic admission limit only.  The old map has no finite material
    // volume at its first admission time, so attempting to invert A0 or adding
    // a density/J floor would manufacture mass.  This branch is not a valid
    // production CME startup.  The replacement piston model instead starts
    // with a finite ambient column adjacent to a velocity-compatible contact;
    // the disturbed sheath and shock then form dynamically.
    candidate->minimumJacobian=1.0;
    current_=candidate;
    return Return::Success(std::move(candidate));
  }
  const double dt=(epochS-begin)/admissionTimeCells;
  const double dPolar=kPi/polarCells;
  const double dAzimuth=2*kPi/azimuthCells;
  std::set<std::tuple<std::uint64_t,int>> identities;
  for(int it=0;it<admissionTimeCells;++it) {
    const double crossing=begin+(it+0.5)*dt;
    const auto birthEvent=event_->At(crossing);
    if(!birthEvent.ok())return Return::Failure(birthEvent.status.code,birthEvent.status.message);
    const auto birthShape=FixedOrientationEllipsoid::FromCenter(event_->basis,
        birthEvent.value.ellipsoid,event_->support.solarRadiusM);
    if(!birthShape.ok())return Return::Failure(birthShape.status.code,birthShape.status.message);
    const auto patches=birthShape.value.Tessellate(polarCells,azimuthCells);
    if(!patches.ok())return Return::Failure(patches.status.code,patches.status.message);
    for(const auto& patch:patches.value) {
      const auto birth=EvaluateBirth(*event_,*ambient_,patch.polarParameterRad,
          patch.azimuthParameterRad,crossing);
      if(!birth.ok()) {
        if(birth.status.code==Core::StatusCode::UnsupportedCapability) {
          candidate->subfastAreaTimeM2S+=patch.areaM2*dt;continue;
        }
        return Return::Failure(birth.status.code,birth.status.message);
      }
      const auto globalContact=EvaluateBirth(*event_,*ambient_,
          patch.polarParameterRad,patch.azimuthParameterRad,begin);
      if(!globalContact.ok())return Return::Failure(
          Core::StatusCode::UnsupportedCapability,
          "fast sheath patch lacks a compatible globally initialized material contact");
      candidate->fastAreaTimeM2S+=patch.areaM2*dt;
      SheathMaterialLabel label{patch.polarParameterRad,patch.azimuthParameterRad,
          crossing,patch.physicalId};
      const auto located=LocateFirstDisposition(*event_,*ambient_,label,epochS);
      if(!located.ok())return Return::Failure(located.status.code,
          "sheath patch "+std::to_string(patch.physicalId)+": "+located.status.message);
      const auto state=Evaluate(label,located.value.time);
      if(!state.ok())return Return::Failure(state.status.code,
          "sheath patch "+std::to_string(patch.physicalId)+": "+state.status.message);
      if(!identities.insert({patch.physicalId,it}).second)return Return::Failure(
          Core::StatusCode::DataIntegrityFailure,"duplicate sheath admission lineage/time cell");
      // Midpoint integration uses the exact curved surface metric at the
      // admitted label.  The corresponding current volume uses det(A), not an
      // unweighted normal column.  Refinement, rather than a density floor or
      // rescaled shell thickness, controls the quadrature error.
      const double birthAreaM2=state.value.birthAreaDensityM2PerRad2*
          dPolar*dAzimuth;
      const double upstreamMass=birth.value.upstream.massDensityKgM3*birth.value.w1*
          birthAreaM2*dt;
      const double downstreamMass=birth.value.jump.downstream.massDensityKgM3*
          birth.value.w2*birthAreaM2*dt;
      candidate->maximumAdmissionMassResidual=std::max(
          candidate->maximumAdmissionMassResidual,Relative(upstreamMass,downstreamMass));
      SheathAdmissionCell cell;cell.label=label;cell.intervalS=dt;
      cell.evaluationTimeS=located.value.time;
      cell.birthAreaM2=birthAreaM2;cell.admittedMassKg=downstreamMass;
      cell.currentVolumeM3=state.value.currentVolumeDensityM3PerRad2S*
          dPolar*dAzimuth*dt;
      cell.state=state.value;
      cell.disposition=located.value.disposition;
      if(cell.disposition==SheathCellDisposition::SolarExit) {
        candidate->solarExitMassKg+=cell.admittedMassKg;
      } else if(cell.disposition==SheathCellDisposition::OuterExit) {
        candidate->outerExitMassKg+=cell.admittedMassKg;
      } else {
        const auto contact=EvaluateContact(label.polarRad,label.azimuthRad,
            located.value.time);
        if(!contact.ok())return Return::Failure(contact.status.code,
            "retained sheath cell has no material rear boundary: "+contact.status.message);
        const double ordered=Dot(state.value.positionM-contact.value.positionM,
            contact.value.outwardNormal);
        if(!(ordered>=-1e-10*std::max(1.0,Norm(state.value.positionM))))
          return Return::Failure(Core::StatusCode::InvalidState,
              "shock-fed map crossed its oldest-cohort material contact");
        candidate->retainedMassKg+=cell.admittedMassKg;
      }
      const double reconstructed=state.value.primitive.massDensityKgM3*
          cell.currentVolumeM3;
      candidate->maximumInventoryMassAbsoluteResidualKg=std::max(
          candidate->maximumInventoryMassAbsoluteResidualKg,
          std::abs(cell.admittedMassKg-reconstructed));
      candidate->maximumInventoryMassResidual=std::max(
          candidate->maximumInventoryMassResidual,
          Relative(cell.admittedMassKg,reconstructed));
      candidate->minimumJacobian=std::min(candidate->minimumJacobian,
          state.value.jacobian);
      candidate->admittedMassKg+=cell.admittedMassKg;
      candidate->cells.push_back(std::move(cell));
    }
  }
  const double closed=candidate->initialMassKg+candidate->retainedMassKg+
      candidate->solarExitMassKg+candidate->outerExitMassKg;
  candidate->maximumInventoryMassResidual=std::max(
      candidate->maximumInventoryMassResidual,Relative(candidate->admittedMassKg,closed));
  if(candidate->cells.empty())candidate->minimumJacobian=1.0;
  if(!(candidate->minimumJacobian>=event_->regional.minimumJacobian))return Return::Failure(
      Core::StatusCode::InvalidState,"sheath inventory contains a nonpositive-J parcel");
  current_=candidate;
  return Return::Success(std::move(candidate));
}

const char* Name(SheathCellDisposition value) noexcept {
  switch(value) {
    case SheathCellDisposition::Retained:return "retained";
    case SheathCellDisposition::SolarExit:return "solar-exit";
    case SheathCellDisposition::OuterExit:return "outer-exit";
  }
  return "unknown";
}

} } // namespace SEP::CoronaSwcme
