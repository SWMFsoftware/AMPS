#include "sep_coronal_cme/discontinuity_transport.h"
#include "sep_coronal_cme/constants.h"
#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronalCME {
namespace {
bool Finite(Vec3 v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
bool UnitVector(Vec3 v) { return Finite(v) && std::abs(Norm(v)-1.0) < 1.0e-12; }
bool Valid(const MovingPlane& p) {
  return !p.identity.stableId.empty() && UnitVector(p.normal) && Finite(p.originM) &&
      std::isfinite(p.epochS) && std::isfinite(p.normalSpeedMPerS);
}
double Distance(const MovingPlane& p, Vec3 x, double time) {
  return Dot(x-p.originM,p.normal)-p.normalSpeedMPerS*(time-p.epochS);
}
double Relative(double a,double b) {
  return std::abs(a-b)/std::max({std::abs(a),std::abs(b),1.0e-300});
}
double EnergyDensity(const MhdPrimitiveState& s,double gamma) {
  return s.pressurePa/(gamma-1.0)+0.5*s.massDensityKgM3*Dot(s.velocityMPerS,s.velocityMPerS)+
      Dot(s.magneticFieldT,s.magneticFieldT)/(2.0*Constants::kVacuumPermeabilityHPerM);
}
double EnergyFlux(const MhdPrimitiveState& s,Vec3 n,double gamma) {
  const double mu=Constants::kVacuumPermeabilityHPerM;
  return (EnergyDensity(s,gamma)+s.pressurePa+Dot(s.magneticFieldT,s.magneticFieldT)/(2.0*mu))*
      Dot(s.velocityMPerS,n)-Dot(s.velocityMPerS,s.magneticFieldT)*Dot(s.magneticFieldT,n)/mu;
}
double Gamma(Vec3 p,double mass) {
  return std::hypot(1.0,Norm(p)/(mass*Constants::kSpeedOfLightMPerS));
}
Vec3 Velocity(Vec3 p,double mass) { return p/(Gamma(p,mass)*mass); }
// Drift-rotate-drift relativistic Boris map. For E=0 this rotation preserves
// |p| to roundoff, and the two drifts make the map time reversible. The field
// is evaluated at the intermediate position, not at a cell-averaged sector.
OrbitState Boris(const FiniteHcsSheet& sheet,OrbitState s,double mass,double charge,double dt,double grid) {
  const Vec3 mid=s.positionM+0.5*dt*Velocity(s.momentumKgMPerS,mass);
  Vec3 b=sheet.Evaluate(mid).value.plasma.magneticFieldT;
  if (grid>0.0) {
    const auto& p=sheet.Parameters();
    const double d=Dot(mid-p.sheet.originM,p.sheet.normal),lo=grid*std::floor(d/grid);
    const double fraction=(d-lo)/grid;
    const Vec3 left=sheet.Evaluate(mid+(lo-d)*p.sheet.normal).value.plasma.magneticFieldT;
    const Vec3 right=sheet.Evaluate(mid+(lo+grid-d)*p.sheet.normal).value.plasma.magneticFieldT;
    b=p.fieldMagnitudeT*Unit((1.0-fraction)*left+fraction*right);
  }
  const Vec3 t=charge*dt/(2.0*mass*Gamma(s.momentumKgMPerS,mass))*b;
  const Vec3 rotation=2.0*t/(1.0+Dot(t,t));
  const Vec3 prime=s.momentumKgMPerS+Cross(s.momentumKgMPerS,t);
  s.momentumKgMPerS=s.momentumKgMPerS+Cross(prime,rotation);
  s.positionM=mid+0.5*dt*Velocity(s.momentumKgMPerS,mass);
  s.timeS+=dt;
  return s;
}
OrbitState UniformBoris(OrbitState s,const MhdPrimitiveState& plasma,double mass,double charge,double dt) {
  const Vec3 electric=-1.0*Cross(plasma.velocityMPerS,plasma.magneticFieldT);
  const Vec3 mid=s.positionM+0.5*dt*Velocity(s.momentumKgMPerS,mass);
  const Vec3 minus=s.momentumKgMPerS+0.5*charge*dt*electric;
  const Vec3 t=charge*dt/(2.0*mass*Gamma(minus,mass))*plasma.magneticFieldT;
  const Vec3 prime=minus+Cross(minus,t);
  const Vec3 rotated=minus+Cross(prime,2.0*t/(1.0+Dot(t,t)));
  s.momentumKgMPerS=rotated+0.5*charge*dt*electric;
  s.positionM=mid+0.5*dt*Velocity(s.momentumKgMPerS,mass);
  s.timeS+=dt;
  return s;
}

int EventPriority(TransportSurfaceKind kind) {
  // An unresolved exclusion is processed before any operator that could pass
  // the particle through it; a qualified magnetic rotation precedes the jump.
  if (kind==TransportSurfaceKind::PfssScsTransition || kind==TransportSurfaceKind::Separatrix) return 0;
  if (kind==TransportSurfaceKind::Shock) return 2;
  return 1;
}
}

Core::Result<SegmentEvent> LocatePlaneCrossing(const MovingPlane& plane,
    Vec3 a,Vec3 b,double ta,double tb) {
  if (!Valid(plane) || !Finite(a) || !Finite(b) || !std::isfinite(ta) || !std::isfinite(tb) || ta==tb)
    return Core::Result<SegmentEvent>::Failure(Core::StatusCode::InvalidState,"invalid crossing segment/surface");
  SegmentEvent result;
  const double da=Distance(plane,a,ta),db=Distance(plane,b,tb);
  if (da==0.0 || (da>0.0 && db>0.0) || (da<0.0 && db<0.0))
    return Core::Result<SegmentEvent>::Success(result);
  const double fraction=da/(da-db);
  if (!(fraction>0.0 && fraction<=1.0)) return Core::Result<SegmentEvent>::Success(result);
  result.crossed=true;
  result.event={plane.identity,ta+fraction*(tb-ta),fraction,
      a+fraction*(b-a),da<0.0?InterfaceSide::Minus:InterfaceSide::Plus,
      da<0.0?InterfaceSide::Plus:InterfaceSide::Minus};
  return Core::Result<SegmentEvent>::Success(result);
}

std::vector<CrossingEvent> OrderCrossingEvents(std::vector<CrossingEvent> events,double tolerance) {
  std::sort(events.begin(),events.end(),[](const CrossingEvent& a,const CrossingEvent& b){
    return std::tie(a.timeS,a.surface.stableId,a.surface.generation)<
           std::tie(b.timeS,b.surface.stableId,b.surface.generation);
  });
  // Form clusters from their earliest time; a tolerance comparator directly
  // inside std::sort would violate transitivity and yield undefined ordering.
  for (std::size_t begin=0;begin<events.size();) {
    std::size_t end=begin+1;
    while (end<events.size() && events[end].timeS-events[begin].timeS<=std::max(0.0,tolerance)) ++end;
    std::sort(events.begin()+begin,events.begin()+end,[](const CrossingEvent& a,const CrossingEvent& b){
      return std::make_tuple(EventPriority(a.surface.kind),a.surface.stableId,a.surface.generation)<
             std::make_tuple(EventPriority(b.surface.kind),b.surface.stableId,b.surface.generation);
    });
    begin=end;
  }
  return events;
}

Core::Result<std::shared_ptr<const FiniteHcsSheet>> FiniteHcsSheet::Create(const FiniteHcsParameters& p) {
  using R=Core::Result<std::shared_ptr<const FiniteHcsSheet>>;
  if (!Valid(p.sheet) || p.sheet.identity.kind!=TransportSurfaceKind::FiniteHcs ||
      p.sheet.normalSpeedMPerS!=0.0 || !UnitVector(p.tangent) ||
      std::abs(Dot(p.tangent,p.sheet.normal))>1.0e-12 ||
      !(p.halfThicknessM>0.0 && p.fieldMagnitudeT>0.0 && p.massDensityKgM3>0.0 && p.pressurePa>0.0) ||
      !std::isfinite(p.halfThicknessM+p.fieldMagnitudeT+p.massDensityKgM3+p.pressurePa) ||
      !(p.outwardWaveEnergyJPerM3>=0.0 && p.inwardWaveEnergyJPerM3>=0.0) ||
      !std::isfinite(p.outwardWaveEnergyJPerM3+p.inwardWaveEnergyJPerM3))
    return R::Failure(Core::StatusCode::InvalidConfiguration,"finite exterior HCS requires a stationary orthonormal force-free profile and positive physical thickness/state");
  auto result=std::shared_ptr<FiniteHcsSheet>(new FiniteHcsSheet);
  result->parameters_=p;
  return R::Success(result);
}
Core::Result<FiniteHcsState> FiniteHcsSheet::Evaluate(Vec3 x) const {
  if (!Finite(x)) return Core::Result<FiniteHcsState>::Failure(Core::StatusCode::InvalidState,"nonfinite HCS coordinate");
  const auto& p=parameters_;
  FiniteHcsState s;
  s.signedDistanceM=Distance(p.sheet,x,p.sheet.epochS);
  const double z=s.signedDistanceM/p.halfThicknessM;
  // sech computed without cosh overflow, including the thin-sheet limit.
  const double e=std::exp(-std::abs(z));
  const double sech=2.0*e/(1.0+e*e);
  s.plasma={p.massDensityKgM3,p.pressurePa,{},p.fieldMagnitudeT*(
      std::tanh(z)*p.tangent+sech*Cross(p.sheet.normal,p.tangent))};
  s.waves=BuildDirectionalWaveState(p.outwardWaveEnergyJPerM3,p.inwardWaveEnergyJPerM3,
      z<0.0?MagneticSector::Negative:MagneticSector::Positive).value;
  s.signedWaveLabelsValid=z!=0.0;
  if (!s.signedWaveLabelsValid) {
    s.waves.parallelJPerM3=s.waves.antiparallelJPerM3=s.waves.fieldCrossHelicity=0.0;
  }
  return Core::Result<FiniteHcsState>::Success(s);
}
Core::Result<HcsOrbitResult> FiniteHcsSheet::Advance(const OrbitState& initial,
    double mass,double charge,double duration,const OrbitControls& c) const {
  using R=Core::Result<HcsOrbitResult>;
  if (!Finite(initial.positionM) || !Finite(initial.momentumKgMPerS) || !std::isfinite(initial.timeS) ||
      !(mass>0.0) || !std::isfinite(mass+charge+duration) ||
      !(c.maximumGyroAngleRad>0.0 && c.maximumGyroAngleRad<1.0 &&
        c.maximumThicknessFraction>0.0 && c.maximumThicknessFraction<1.0 && c.eventTimeToleranceS>0.0) ||
      !std::isfinite(c.maximumGyroAngleRad+c.maximumThicknessFraction+c.eventTimeToleranceS+c.fieldGridSpacingM) ||
      c.fieldGridSpacingM<0.0 || c.fieldGridSpacingM>parameters_.halfThicknessM || c.maximumSubsteps==0)
    return R::Failure(Core::StatusCode::InvalidConfiguration,"invalid finite-HCS full-orbit state/controls");
  HcsOrbitResult result; result.state=initial;
  if (duration==0.0) return R::Success(result);
  const double speed=Norm(Velocity(initial.momentumKgMPerS,mass));
  const double omega=std::abs(charge)*parameters_.fieldMagnitudeT/(mass*Gamma(initial.momentumKgMPerS,mass));
  const double bound=std::min(omega>0.0?c.maximumGyroAngleRad/omega:std::numeric_limits<double>::infinity(),
      speed>0.0?c.maximumThicknessFraction*parameters_.halfThicknessM/speed:std::numeric_limits<double>::infinity());
  const double needed=std::max(1.0,std::ceil(std::abs(duration)/bound));
  if (!std::isfinite(needed) || needed>static_cast<double>(c.maximumSubsteps))
    return R::Failure(Core::StatusCode::NumericalFailure,"HCS orbit exceeds declared substep budget");
  const auto count=static_cast<std::uint64_t>(needed);
  const double dt=duration/static_cast<double>(count);
  for (std::uint64_t k=0;k<count;++k) {
    const OrbitState start=result.state;
    OrbitState end=Boris(*this,start,mass,charge,dt,c.fieldGridSpacingM);
    const double a=Distance(parameters_.sheet,start.positionM,start.timeS);
    const double b=Distance(parameters_.sheet,end.positionM,end.timeS);
    if (a!=0.0 && ((a<0.0 && b>=0.0)||(a>0.0 && b<=0.0))) {
      // Bisect the actual Boris map, rather than an interpolated mesh field.
      // Then split at the event and recompute the remainder from that state.
      double lo=0.0,hi=1.0;
      for (int iteration=0;iteration<100 && (hi-lo)*std::abs(dt)>c.eventTimeToleranceS;++iteration) {
        const double mid=0.5*(lo+hi);
        const auto trial=Boris(*this,start,mass,charge,mid*dt,c.fieldGridSpacingM);
        const double d=Distance(parameters_.sheet,trial.positionM,trial.timeS);
        if ((a<0.0 && d<0.0)||(a>0.0 && d>0.0)) lo=mid; else hi=mid;
      }
      const double f=0.5*(lo+hi);
      auto at=Boris(*this,start,mass,charge,f*dt,c.fieldGridSpacingM);
      // Project only the event-location roundoff (<=v*event tolerance), not a
      // physical step. The momentum and its energy are never projected.
      at.positionM=at.positionM-Distance(parameters_.sheet,at.positionM,at.timeS)*parameters_.sheet.normal;
      result.events.push_back({parameters_.sheet.identity,at.timeS,f,at.positionM,
          a<0.0?InterfaceSide::Minus:InterfaceSide::Plus,a<0.0?InterfaceSide::Plus:InterfaceSide::Minus});
      end=Boris(*this,at,mass,charge,(1.0-f)*dt,c.fieldGridSpacingM);
    }
    result.state=end; ++result.substeps;
  }
  result.state.timeS=initial.timeS+duration;
  return R::Success(result);
}

Core::Status CrossingLedger::Restore(const std::set<Key>& keys) {
  for (const auto& key:keys)
    if (std::get<0>(key).empty() || std::get<1>(key).empty())
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,"empty crossing lineage/surface identity");
  keys_=keys;
  return Core::Status::Success();
}

Core::Result<std::shared_ptr<const PlanarShockSheath>> PlanarShockSheath::Create(const SheathParameters& p) {
  using R=Core::Result<std::shared_ptr<const PlanarShockSheath>>;
  const double first=p.thicknessAtEpochM+p.thicknessRateMPerS*(p.beginTimeS-p.shock.epochS);
  const double last=p.thicknessAtEpochM+p.thicknessRateMPerS*(p.endTimeS-p.shock.epochS);
  if (!Valid(p.shock) || p.shock.identity.kind!=TransportSurfaceKind::Shock || p.shock.identity.generation==0 ||
      !(p.endTimeS>p.beginTimeS && first>0.0 && last>0.0 && p.auditAreaM2>0.0 &&
        p.conservationTolerance>0.0 && p.conservationTolerance<1.0 &&
        p.maximumPassiveWavePressureFraction>0.0 && p.maximumPassiveWavePressureFraction<1.0) ||
      !std::isfinite(p.beginTimeS+p.endTimeS+first+last+p.auditAreaM2+p.conservationTolerance+p.maximumPassiveWavePressureFraction) ||
      (p.magneticSector!=MagneticSector::Positive && p.magneticSector!=MagneticSector::Negative) ||
      !(p.outwardWaveEnergyJPerM3>=0.0 && p.inwardWaveEnergyJPerM3>=0.0) ||
      !std::isfinite(p.outwardWaveEnergyJPerM3+p.inwardWaveEnergyJPerM3))
    return R::Failure(Core::StatusCode::InvalidConfiguration,"invalid finite moving planar sheath geometry/coverage");
  const auto solved=SolveObliqueFastShock(p.upstream,p.shock.normal,p.shock.normalSpeedMPerS,
      p.gammaAdiabatic,p.conservationTolerance);
  if (!solved.ok()) return R::Failure(solved.status.code,solved.status.message);
  auto result=std::shared_ptr<PlanarShockSheath>(new PlanarShockSheath);
  result->parameters_=p; result->jump_=solved.value;
  const auto& d=solved.value.downstream;
  const double w1=Dot(p.upstream.velocityMPerS,p.shock.normal)-p.shock.normalSpeedMPerS;
  const double w2=Dot(d.velocityMPerS,p.shock.normal)-p.shock.normalSpeedMPerS;
  const double a1=static_cast<int>(p.magneticSector)*Dot(p.upstream.magneticFieldT,p.shock.normal)/std::sqrt(Constants::kVacuumPermeabilityHPerM*p.upstream.massDensityKgM3);
  const double a2=static_cast<int>(p.magneticSector)*Dot(d.magneticFieldT,p.shock.normal)/std::sqrt(Constants::kVacuumPermeabilityHPerM*d.massDensityKgM3);
  if (std::abs(w2+a2)<1.0e-12 || std::abs(w2-a2)<1.0e-12)
    return R::Failure(Core::StatusCode::UnsupportedCapability,"passive wave transmission has a zero-speed characteristic");
  auto waves=BuildDirectionalWaveState(p.outwardWaveEnergyJPerM3*std::abs(w1+a1)/std::abs(w2+a2),
      p.inwardWaveEnergyJPerM3*std::abs(w1-a1)/std::abs(w2-a2),
      p.magneticSector);
  if (!waves.ok()) return R::Failure(waves.status.code,waves.status.message);
  result->waves_=waves.value;
  const double upstreamWaveFraction=(p.outwardWaveEnergyJPerM3+p.inwardWaveEnergyJPerM3)/
      (p.upstream.pressurePa+p.upstream.massDensityKgM3*w1*w1);
  const double downstreamWaveFraction=waves.value.totalJPerM3/(d.pressurePa+d.massDensityKgM3*w2*w2);
  if (std::max(upstreamWaveFraction,downstreamWaveFraction)>p.maximumPassiveWavePressureFraction)
    return R::Failure(Core::StatusCode::UnsupportedCapability,"waves exceed declared passive-backreaction bound");
  const auto audit=result->Audit(32,p.beginTimeS);
  if (!audit.ok() || !audit.value.passed)
    return R::Failure(Core::StatusCode::NumericalFailure,"downstream sheath failed independent global conservation gate");
  return R::Success(result);
}
Core::Result<MhdPrimitiveState> PlanarShockSheath::Evaluate(Vec3 x,double time,bool onFront,InterfaceSide side) const {
  const auto& p=parameters_;
  if (!Finite(x)||!std::isfinite(time)||time<p.beginTimeS||time>p.endTimeS)
    return Core::Result<MhdPrimitiveState>::Failure(Core::StatusCode::OutOfDomain,"sheath query outside immutable time coverage");
  const double d=Distance(p.shock,x,time),thickness=p.thicknessAtEpochM+p.thicknessRateMPerS*(time-p.shock.epochS);
  if (d < -thickness) return Core::Result<MhdPrimitiveState>::Failure(Core::StatusCode::OutOfDomain,"no downstream extrapolation behind finite sheath");
  if (onFront && std::abs(d)>1.0e-8*std::max(1.0,p.thicknessAtEpochM))
    return Core::Result<MhdPrimitiveState>::Failure(Core::StatusCode::InvalidState,"one-sided front flag supplied away from shock");
  if (d==0.0 && !onFront) return Core::Result<MhdPrimitiveState>::Failure(Core::StatusCode::InvalidState,"shock query requires a one-sided state");
  return Core::Result<MhdPrimitiveState>::Success((d>0.0 || (onFront && side==InterfaceSide::Plus))?p.upstream:jump_.downstream);
}
Core::Result<SheathAudit> PlanarShockSheath::Audit(int cells,double time) const {
  using R=Core::Result<SheathAudit>;
  const auto& p=parameters_; const auto& a=p.upstream; const auto& b=jump_.downstream;
  if (cells<2 || !std::isfinite(time) || time<p.beginTimeS || time>p.endTimeS)
    return R::Failure(Core::StatusCode::OutOfDomain,"invalid sheath audit resolution/time");
  const double n1=Dot(a.magneticFieldT,p.shock.normal),n2=Dot(b.magneticFieldT,p.shock.normal);
  SheathAudit r;
  // Independent finite-volume face quadrature. Tangential side faces cancel
  // exactly for this planar family; no numerical div(B) stencil crosses it.
  double normalFlux=0.0;
  for (int i=0;i<cells;++i) normalFlux+=(n1-n2)*p.auditAreaM2/cells;
  r.normalFluxRelative=std::abs(normalFlux)/(p.auditAreaM2*std::max({Norm(a.magneticFieldT),Norm(b.magneticFieldT),1.0e-300}));
  r.divergenceRelative=r.normalFluxRelative;
  const double va=Dot(a.velocityMPerS,p.shock.normal),vb=Dot(b.velocityMPerS,p.shock.normal),vs=p.shock.normalSpeedMPerS;
  r.massRelative=Relative(a.massDensityKgM3*(va-vs),b.massDensityKgM3*(vb-vs));
  const double ea=EnergyDensity(a,p.gammaAdiabatic),eb=EnergyDensity(b,p.gammaAdiabatic);
  r.energyRelative=Relative(EnergyFlux(a,p.shock.normal,p.gammaAdiabatic)-vs*ea,
                           EnergyFlux(b,p.shock.normal,p.gammaAdiabatic)-vs*eb);
  const double rearSpeed=vs-p.thicknessRateMPerS;
  const double massNet=p.auditAreaM2*(-b.massDensityKgM3*(vb-vs)+b.massDensityKgM3*(vb-rearSpeed));
  const double energyNet=p.auditAreaM2*(-(EnergyFlux(b,p.shock.normal,p.gammaAdiabatic)-vs*eb)+
                                     (EnergyFlux(b,p.shock.normal,p.gammaAdiabatic)-rearSpeed*eb));
  r.massInventoryRateKgPerS=p.auditAreaM2*b.massDensityKgM3*p.thicknessRateMPerS;
  r.energyInventoryRateW=p.auditAreaM2*eb*p.thicknessRateMPerS;
  // Normalize inventory cancellation by the boundary throughput so a fixed
  // width (zero inventory rate) is not subjected to a zero-denominator ratio.
  r.massRelative=std::max(r.massRelative,std::abs(massNet-r.massInventoryRateKgPerS)/std::max(std::abs(p.auditAreaM2*b.massDensityKgM3*(vb-vs)),1.0e-300));
  r.energyRelative=std::max(r.energyRelative,std::abs(energyNet-r.energyInventoryRateW)/std::max(std::abs(p.auditAreaM2*(EnergyFlux(b,p.shock.normal,p.gammaAdiabatic)-vs*eb)),1.0e-300));
  const double an1=static_cast<int>(p.magneticSector)*n1/std::sqrt(Constants::kVacuumPermeabilityHPerM*a.massDensityKgM3),an2=static_cast<int>(p.magneticSector)*n2/std::sqrt(Constants::kVacuumPermeabilityHPerM*b.massDensityKgM3);
  r.waveFluxRelative=std::max(Relative(p.outwardWaveEnergyJPerM3*std::abs(va-vs+an1),waves_.outwardJPerM3*std::abs(vb-vs+an2)),
      Relative(p.inwardWaveEnergyJPerM3*std::abs(va-vs-an1),waves_.inwardJPerM3*std::abs(vb-vs-an2)));
  r.passed=std::max({r.divergenceRelative,r.normalFluxRelative,r.massRelative,r.energyRelative,r.waveFluxRelative})<=p.conservationTolerance;
  return R::Success(r);
}
Core::Result<ShockCrossingState> PlanarShockSheath::Cross(const std::string& id,
    const CrossingEvent& event,const FourMomentum& inertial,CrossingLedger* ledger) const {
  using R=Core::Result<ShockCrossingState>;
  const auto& p=parameters_;
  if (ledger==nullptr || id.empty() || event.surface.kind!=TransportSurfaceKind::Shock ||
      event.surface.stableId!=p.shock.identity.stableId || event.surface.generation!=p.shock.identity.generation ||
      event.incoming==event.outgoing || event.timeS<p.beginTimeS || event.timeS>p.endTimeS ||
      !std::isfinite(event.timeS) || !Finite(event.positionM) ||
      std::abs(Distance(p.shock,event.positionM,event.timeS))>1.0e-8*std::max(1.0,p.thicknessAtEpochM) ||
      !Finite(inertial.momentumKgMPerS) || !std::isfinite(inertial.totalEnergyJ) ||
      !(FourMomentumInvariant(inertial)>0.0) || !std::isfinite(FourMomentumInvariant(inertial)))
    return R::Failure(Core::StatusCode::InvalidState,"invalid shock crossing identity, root, or timelike momentum");
  const auto& outgoing=event.outgoing==InterfaceSide::Plus?p.upstream:jump_.downstream;
  const auto local=BoostFourMomentum(inertial,-1.0*outgoing.velocityMPerS);
  const auto shock=BoostFourMomentum(inertial,-p.shock.normalSpeedMPerS*p.shock.normal);
  if (!local.ok() || !shock.ok()) return R::Failure(Core::StatusCode::InvalidState,"shock/plasma frame transformation failed");
  ShockCrossingState result;
  result.inertial=inertial; result.localPlasma=local.value;
  const double magnitude=Norm(local.value.momentumKgMPerS),field=Norm(outgoing.magneticFieldT);
  result.pitchAngleCosine=magnitude>0.0&&field>0.0?Dot(local.value.momentumKgMPerS,outgoing.magneticFieldT)/(magnitude*field):0.0;
  result.shockFrameEnergyJ=shock.value.totalEnergyJ;
  const CrossingLedger::Key key{id,event.surface.stableId,event.surface.generation};
  if (!ledger->Contains(key)) { ledger->keys_.insert(key); result.applied=true; }
  // No impulsive energy gain is invented at a zero-potential mathematical
  // front. Energy changes in the local frame are the explicit Lorentz boost;
  // subsequent scattering uses that local state and transmitted wave profile.
  return R::Success(result);
}
Core::Result<SheathOrbitResult> PlanarShockSheath::Advance(const std::string& id,const OrbitState& initial,
    double mass,double charge,double duration,CrossingLedger* ledger,const OrbitControls& controls) const {
  using R=Core::Result<SheathOrbitResult>;
  const auto& p=parameters_;
  if (id.empty() || ledger==nullptr || !Finite(initial.positionM) || !Finite(initial.momentumKgMPerS) ||
      !(mass>0.0) || !std::isfinite(mass+charge+duration+initial.timeS) ||
      initial.timeS<p.beginTimeS || initial.timeS>p.endTimeS || initial.timeS+duration<p.beginTimeS || initial.timeS+duration>p.endTimeS ||
      !(controls.maximumGyroAngleRad>0.0 && controls.maximumGyroAngleRad<1.0 &&
        controls.maximumThicknessFraction>0.0 && controls.maximumThicknessFraction<1.0 && controls.eventTimeToleranceS>0.0) ||
      !std::isfinite(controls.maximumGyroAngleRad+controls.maximumThicknessFraction+controls.eventTimeToleranceS) ||
      controls.fieldGridSpacingM!=0.0 || controls.maximumSubsteps==0)
    return R::Failure(Core::StatusCode::InvalidConfiguration,"invalid unscattered sheath orbit/coverage; shock fields cannot be grid-blended");
  SheathOrbitResult result; result.state=initial;
  if (duration==0.0) return R::Success(result);
  // All mutations are staged. A rear-boundary, coverage or numerical failure
  // must not leave a committed shock event beside an unpublished particle step.
  CrossingLedger staged=*ledger;
  const double field=std::max(Norm(p.upstream.magneticFieldT),Norm(jump_.downstream.magneticFieldT));
  const double omega=std::abs(charge)*field/mass;
  const double electric=std::max(Norm(SEP::CoronalCME::Cross(p.upstream.velocityMPerS,p.upstream.magneticFieldT)),
      Norm(SEP::CoronalCME::Cross(jump_.downstream.velocityMPerS,jump_.downstream.magneticFieldT)));
  const double speedBound=std::min(Constants::kSpeedOfLightMPerS,
      (Norm(initial.momentumKgMPerS)+std::abs(charge)*electric*std::abs(duration))/mass)+std::abs(p.shock.normalSpeedMPerS);
  const double width=std::min(p.thicknessAtEpochM+p.thicknessRateMPerS*(p.beginTimeS-p.shock.epochS),
      p.thicknessAtEpochM+p.thicknessRateMPerS*(p.endTimeS-p.shock.epochS));
  const double maximumDt=std::min(omega>0.0?controls.maximumGyroAngleRad/omega:std::numeric_limits<double>::infinity(),
      speedBound>0.0?controls.maximumThicknessFraction*width/speedBound:std::numeric_limits<double>::infinity());
  const double needed=std::max(1.0,std::ceil(std::abs(duration)/maximumDt));
  if (!std::isfinite(needed) || needed>static_cast<double>(controls.maximumSubsteps))
    return R::Failure(Core::StatusCode::NumericalFailure,"sheath orbit exceeds declared substep budget");
  const auto count=static_cast<std::uint64_t>(needed);
  const double dt=duration/static_cast<double>(count);
  for (std::uint64_t k=0;k<count;++k) {
    const auto start=result.state;
    const double a=Distance(p.shock,start.positionM,start.timeS);
    const double leaving=Dot(Velocity(start.momentumKgMPerS,mass),p.shock.normal)-p.shock.normalSpeedMPerS;
    const auto side=(a>0.0 || (a==0.0 && leaving*duration>0.0))?InterfaceSide::Plus:InterfaceSide::Minus;
    const auto incoming=Evaluate(start.positionM,start.timeS,a==0.0,side);
    if (!incoming.ok()) return R::Failure(incoming.status.code,incoming.status.message);
    auto end=UniformBoris(start,incoming.value,mass,charge,dt);
    const double b=Distance(p.shock,end.positionM,end.timeS);
    if (a!=0.0 && ((a<0.0 && b>=0.0)||(a>0.0 && b<=0.0))) {
      double lo=0.0,hi=1.0;
      for (int iteration=0;iteration<100 && (hi-lo)*std::abs(dt)>controls.eventTimeToleranceS;++iteration) {
        const double f=0.5*(lo+hi);
        const auto trial=UniformBoris(start,incoming.value,mass,charge,f*dt);
        const double d=Distance(p.shock,trial.positionM,trial.timeS);
        if ((a<0.0 && d<0.0)||(a>0.0 && d>0.0)) lo=f; else hi=f;
      }
      const double fraction=0.5*(lo+hi);
      auto at=UniformBoris(start,incoming.value,mass,charge,fraction*dt);
      at.positionM=at.positionM-Distance(p.shock,at.positionM,at.timeS)*p.shock.normal;
      CrossingEvent event{p.shock.identity,at.timeS,fraction,at.positionM,
          a<0.0?InterfaceSide::Minus:InterfaceSide::Plus,a<0.0?InterfaceSide::Plus:InterfaceSide::Minus};
      FourMomentum momentum{Gamma(at.momentumKgMPerS,mass)*mass*Constants::kSpeedOfLightMPerS*Constants::kSpeedOfLightMPerS,at.momentumKgMPerS};
      const auto crossed=Cross(id,event,momentum,&staged);
      if (!crossed.ok()) return R::Failure(crossed.status.code,crossed.status.message);
      result.events.push_back(event); result.crossingStates.push_back(crossed.value);
      const auto outgoing=Evaluate(at.positionM,at.timeS,true,event.outgoing);
      if (!outgoing.ok()) return R::Failure(outgoing.status.code,outgoing.status.message);
      end=UniformBoris(at,outgoing.value,mass,charge,(1.0-fraction)*dt);
    }
    // Anchor the clock to the requested interval instead of accumulating dt.
    // Reverse stepping to the first covered epoch must not land just outside
    // immutable coverage because of floating-point summation roundoff.
    end.timeS=(k+1==count)?initial.timeS+duration:
        initial.timeS+static_cast<double>(k+1)*dt;
    if (!Finite(end.positionM)||!Finite(end.momentumKgMPerS))
      return R::Failure(Core::StatusCode::NumericalFailure,"sheath orbit produced a nonfinite state");
    const auto covered=Evaluate(end.positionM,end.timeS,
        Distance(p.shock,end.positionM,end.timeS)==0.0,side);
    if (!covered.ok()) return R::Failure(covered.status.code,covered.status.message);
    result.state=end; ++result.substeps;
  }
  result.state.timeS=initial.timeS+duration;
  *ledger=std::move(staged);
  return R::Success(result);
}

Core::Status ValidateDiscontinuityCapabilities(bool hcsRequested,bool sheathRequested,
    const std::shared_ptr<const FiniteHcsSheet>& hcs,const std::shared_ptr<const PlanarShockSheath>& sheath,
    bool compositeRequested) {
  if (compositeRequested) return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,
      "exterior finite HCS does not qualify the PFSS/SCS transition sheet");
  if (hcsRequested && !hcs) return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"finite-HCS operator/provider is absent");
  if (sheathRequested && !sheath) return Core::Status::Failure(Core::StatusCode::UnsupportedCapability,"independent globally conservative sheath provider is absent");
  return Core::Status::Success();
}
} } // namespace SEP::CoronalCME
