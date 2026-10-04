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
  if(plus.ok())return Return::Success((plus.value-center)/step);
  if(minus.ok())return Return::Success((center-minus.value)/step);
  return Return::Failure(Core::StatusCode::UnsupportedCapability,
      "shock-fed derivative has no same-fast-branch stencil");
}

double Relative(double a,double b) {
  return std::abs(a-b)/std::max({std::abs(a),std::abs(b),1e-300});
}

struct PositionClassification {
  Vec3 position;
  double contactSigned=0.0; // positive in the sheath, zero at contact
  double solarSigned=0.0;   // positive outside the physical Sun
  double outerSigned=0.0;   // positive inside the shock front
};

Core::Result<PositionClassification> ClassifyPosition(
    const EventConfiguration& event,const BirthState& birth,
    const SheathMaterialLabel& label,double time) {
  using Return=Core::Result<PositionClassification>;
  const auto evolution=event.At(time);
  if(!evolution.ok())return Return::Failure(evolution.status.code,evolution.status.message);
  const auto front=FixedOrientationEllipsoid::FromCenter(event.basis,
      evolution.value.ellipsoid,event.support.solarRadiusM);
  const auto contactKinematics=ContactKinematics(event,evolution.value.ellipsoid);
  const auto contact=contactKinematics.ok()?FixedOrientationEllipsoid::FromCenter(
      event.basis,contactKinematics.value,event.support.solarRadiusM):
      Core::Result<FixedOrientationEllipsoid>::Failure(
          contactKinematics.status.code,contactKinematics.status.message);
  if(!front.ok()||!contact.ok())return Return::Failure(Core::StatusCode::InvalidState,
      !front.ok()?front.status.message:contact.status.message);
  PositionClassification result;
  result.position=front.value.Point(label.polarRad,label.azimuthRad)+
      (time-label.crossingTimeS)*birth.deficit;
  const auto frontValue=front.value.Evaluate(result.position);
  const auto contactValue=contact.value.Evaluate(result.position);
  if(!frontValue.ok()||!contactValue.ok())return Return::Failure(
      Core::StatusCode::InvalidState,"cannot evaluate a sheath boundary event function");
  result.contactSigned=contactValue.value.implicitValue;
  result.solarSigned=Norm(result.position)-event.support.solarRadiusM;
  result.outerSigned=-frontValue.value.implicitValue;
  return Return::Success(result);
}

struct LocatedDisposition {
  SheathCellDisposition disposition=SheathCellDisposition::Retained;
  double time=0.0;
};

Core::Result<LocatedDisposition> LocateFirstDisposition(
    const EventConfiguration& event,const BirthState& birth,
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
  auto previous=ClassifyPosition(event,birth,label,previousTime);
  if(!previous.ok())return Return::Failure(previous.status.code,previous.status.message);
  auto scalar=[](const PositionClassification& state,
      SheathCellDisposition kind) {
    if(kind==SheathCellDisposition::ContactExit)return state.contactSigned;
    if(kind==SheathCellDisposition::SolarExit)return state.solarSigned;
    return state.outerSigned;
  };
  for(int i=1;i<=intervals;++i) {
    const double time=begin+(epoch-begin)*i/intervals;
    if(time<=previousTime)continue;
    const auto current=ClassifyPosition(event,birth,label,time);
    if(!current.ok())return Return::Failure(current.status.code,current.status.message);
    for(const auto kind:{SheathCellDisposition::ContactExit,
        SheathCellDisposition::SolarExit,SheathCellDisposition::OuterExit}) {
      if(scalar(previous.value,kind)>0&&scalar(current.value,kind)<=0) {
        double lower=previousTime,upper=time;
        for(int iteration=0;iteration<80&&upper-lower>event.support.rootToleranceS;
            ++iteration) {
          const double middle=0.5*(lower+upper);
          const auto state=ClassifyPosition(event,birth,label,middle);
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
  if(event->physicsFingerprint!=ambient->Event().physicsFingerprint)
    return Return::Failure(Core::StatusCode::DataIntegrityFailure,
        "sheath and ambient authorities have different event identities");
  if(event->sheathModel!="rh-ballistic-material-map-v1")return Return::Failure(
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
  const double age=epochS-label.crossingTimeS;
  SheathMappedState result;result.label=label;result.eventIdentity=event_->physicsFingerprint;
  result.upstreamRelativeNormalSpeedMPerS=birth.value.w1;
  result.downstreamRelativeNormalSpeedMPerS=birth.value.w2;
  result.shockResidual=birth.value.jump.residuals.maximum;
  result.birthAreaDensityM2PerRad2=Norm(Cross(
      PointPolarDerivative(birth.value.front,label.polarRad,label.azimuthRad),
      PointAzimuthDerivative(birth.value.front,label.polarRad,label.azimuthRad)));
  result.positionM=nowShape.value.Point(label.polarRad,label.azimuthRad)+
      age*birth.value.deficit;
  result.primitive.velocityMPerS=nowShape.value.SurfaceVelocity(
      label.polarRad,label.azimuthRad)+birth.value.deficit;
  if(age<=event_->support.rootToleranceS*1e-6) {
    result.jacobian=1;result.deformationColumns={Vec3{1,0,0},Vec3{0,1,0},Vec3{0,0,1}};
    result.labelDerivativeColumns={
        PointPolarDerivative(birth.value.front,label.polarRad,label.azimuthRad),
        PointAzimuthDerivative(birth.value.front,label.polarRad,label.azimuthRad),
        -1.0*birth.value.deficit};
    result.primitive=birth.value.jump.downstream;
    return Return::Success(std::move(result));
  }

  const double angularStep=5e-4;
  auto polarSample=[&](double value) {
    const auto sample=EvaluateBirth(*event_,*ambient_,value,label.azimuthRad,
        label.crossingTimeS);
    return sample.ok()?Core::Result<Vec3>::Success(sample.value.deficit):
        Core::Result<Vec3>::Failure(sample.status.code,sample.status.message);
  };
  const auto dPolar=OneSidedDerivative(label.polarRad,angularStep,1e-8,kPi-1e-8,
      birth.value.deficit,polarSample);
  auto azimuthSample=[&](double value) {
    value=std::fmod(value+2*kPi,2*kPi);
    const auto sample=EvaluateBirth(*event_,*ambient_,label.polarRad,value,
        label.crossingTimeS);
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
      (azPlus.ok()?(azPlus.value-birth.value.deficit)/angularStep:
          (birth.value.deficit-azMinus.value)/angularStep);
  const double timeStep=std::max(0.5,500*event_->support.rootToleranceS);
  auto timeSample=[&](double value) {
    const auto sample=EvaluateBirth(*event_,*ambient_,label.polarRad,
        label.azimuthRad,value);
    return sample.ok()?Core::Result<Vec3>::Success(sample.value.deficit):
        Core::Result<Vec3>::Failure(sample.status.code,sample.status.message);
  };
  const auto dTime=OneSidedDerivative(label.crossingTimeS,timeStep,
      std::max(event_->support.startS,event_->regional.sheathAdmissionStartS),
      event_->support.endS,birth.value.deficit,timeSample);
  if(!dTime.ok())return Return::Failure(dTime.status.code,dTime.status.message);

  const std::array<Vec3,3> birthBasis={
      PointPolarDerivative(birth.value.front,label.polarRad,label.azimuthRad),
      PointAzimuthDerivative(birth.value.front,label.polarRad,label.azimuthRad),
      -1.0*birth.value.deficit};
  const std::array<Vec3,3> currentBasis={
      PointPolarDerivative(nowShape.value,label.polarRad,label.azimuthRad)+age*dPolar.value,
      PointAzimuthDerivative(nowShape.value,label.polarRad,label.azimuthRad)+age*dAzimuth,
      -1.0*birth.value.deficit+age*dTime.value};
  const double birthDet=Determinant(birthBasis),currentDet=Determinant(currentBasis);
  const double scale=std::max(1.0,Norm(birthBasis[0])*Norm(birthBasis[1])*Norm(birthBasis[2]));
  if(!(std::isfinite(birthDet)&&std::isfinite(currentDet)&&birthDet>1e-14*scale))
    return Return::Failure(Core::StatusCode::InvalidState,
        "shock-fed birth map has singular or reversed crossing volume");
  result.jacobian=currentDet/birthDet;
  if(!(std::isfinite(result.jacobian)&&result.jacobian>=event_->regional.minimumJacobian))
    return Return::Failure(Core::StatusCode::InvalidState,
        "shock-fed ballistic map folded or crossed its minimum Jacobian");
  result.deformationColumns=CartesianMapColumns(currentBasis,birthBasis);
  result.labelDerivativeColumns=currentBasis;
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

Core::Result<std::shared_ptr<const SheathInventory>>
ShockFedSheathModel::PrepareInventory(double epochS,
    std::uint64_t backgroundGeneration,int polarCells,int azimuthCells,
    int admissionTimeCells) {
  using Return=Core::Result<std::shared_ptr<const SheathInventory>>;
  const double begin=std::max(event_->support.startS,
      event_->regional.sheathAdmissionStartS);
  if(!(backgroundGeneration>0&&polarCells>=2&&azimuthCells>=3&&
      admissionTimeCells>=1&&epochS>begin&&epochS<=event_->support.endS))
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "invalid sheath inventory epoch, generation, or resolution");
  auto candidate=std::shared_ptr<SheathInventory>(new SheathInventory);
  candidate->epochS=epochS;candidate->backgroundGeneration=backgroundGeneration;
  candidate->eventIdentity=event_->physicsFingerprint;
  candidate->minimumJacobian=std::numeric_limits<double>::infinity();
  const double dt=(epochS-begin)/admissionTimeCells;
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
      candidate->fastAreaTimeM2S+=patch.areaM2*dt;
      SheathMaterialLabel label{patch.polarParameterRad,patch.azimuthParameterRad,
          crossing,patch.physicalId};
      const auto located=LocateFirstDisposition(*event_,birth.value,label,epochS);
      if(!located.ok())return Return::Failure(located.status.code,
          "sheath patch "+std::to_string(patch.physicalId)+": "+located.status.message);
      const auto state=Evaluate(label,located.value.time);
      if(!state.ok())return Return::Failure(state.status.code,
          "sheath patch "+std::to_string(patch.physicalId)+": "+state.status.message);
      if(!identities.insert({patch.physicalId,it}).second)return Return::Failure(
          Core::StatusCode::DataIntegrityFailure,"duplicate sheath admission lineage/time cell");
      const double upstreamMass=birth.value.upstream.massDensityKgM3*birth.value.w1*
          patch.areaM2*dt;
      const double downstreamMass=birth.value.jump.downstream.massDensityKgM3*
          birth.value.w2*patch.areaM2*dt;
      candidate->maximumAdmissionMassResidual=std::max(
          candidate->maximumAdmissionMassResidual,Relative(upstreamMass,downstreamMass));
      SheathAdmissionCell cell;cell.label=label;cell.intervalS=dt;
      cell.evaluationTimeS=located.value.time;
      cell.birthAreaM2=patch.areaM2;cell.admittedMassKg=downstreamMass;
      cell.currentVolumeM3=state.value.jacobian*birth.value.w2*patch.areaM2*dt;
      cell.state=state.value;
      cell.disposition=located.value.disposition;
      if(cell.disposition==SheathCellDisposition::SolarExit) {
        candidate->solarExitMassKg+=cell.admittedMassKg;
      } else if(cell.disposition==SheathCellDisposition::OuterExit) {
        candidate->outerExitMassKg+=cell.admittedMassKg;
      } else if(cell.disposition==SheathCellDisposition::ContactExit) {
        candidate->contactExitMassKg+=cell.admittedMassKg;
      } else {
        candidate->retainedMassKg+=cell.admittedMassKg;
      }
      const double reconstructed=state.value.primitive.massDensityKgM3*
          cell.currentVolumeM3;
      candidate->maximumInventoryMassResidual=std::max(
          candidate->maximumInventoryMassResidual,
          Relative(cell.admittedMassKg,reconstructed));
      candidate->minimumJacobian=std::min(candidate->minimumJacobian,
          state.value.jacobian);
      candidate->admittedMassKg+=cell.admittedMassKg;
      candidate->cells.push_back(std::move(cell));
    }
  }
  const double closed=candidate->retainedMassKg+candidate->contactExitMassKg+
      candidate->solarExitMassKg+candidate->outerExitMassKg;
  candidate->maximumInventoryMassResidual=std::max(
      candidate->maximumInventoryMassResidual,Relative(candidate->admittedMassKg,closed));
  if(candidate->cells.empty()||!(candidate->minimumJacobian>=
      event_->regional.minimumJacobian))return Return::Failure(
      Core::StatusCode::InvalidState,"sheath inventory has no admissible positive-J parcels");
  current_=candidate;
  return Return::Success(std::move(candidate));
}

const char* Name(SheathCellDisposition value) noexcept {
  switch(value) {
    case SheathCellDisposition::Retained:return "retained";
    case SheathCellDisposition::ContactExit:return "contact-exit";
    case SheathCellDisposition::SolarExit:return "solar-exit";
    case SheathCellDisposition::OuterExit:return "outer-exit";
  }
  return "unknown";
}

} } // namespace SEP::CoronaSwcme
