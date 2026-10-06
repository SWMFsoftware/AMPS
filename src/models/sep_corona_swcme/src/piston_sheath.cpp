#include "sep_corona_swcme/piston_sheath.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP { namespace CoronaSwcme { namespace {

using CoronalCME::KinematicValue;
using CoronalCME::Vec3;

double ShellCentroid(double left,double right) {
  return 0.75*(std::pow(right,4)-std::pow(left,4))/
      (std::pow(right,3)-std::pow(left,3));
}

Core::Result<PistonInitialState> SampleAmbientColumn(
    const PistonAmbientProjection& projection,
    const std::vector<double>& nodes) {
  using Return=Core::Result<PistonInitialState>;
  if(nodes.size()<2)return Return::Failure(
      Core::StatusCode::InvalidConfiguration,
      "piston ambient column requires at least one cell");
  PistonInitialState out;
  const std::size_t cells=nodes.size()-1;
  out.nodeVelocityMPerS.resize(nodes.size());
  out.cellDensityKgM3.resize(cells);out.cellPressurePa.resize(cells);
  out.cellTransverseMagneticFieldT.resize(cells);
  out.cellTransverseMagneticField2T.resize(cells);
  out.cellRadialMagneticFieldT.resize(cells);
  for(std::size_t i=0;i<nodes.size();++i) {
    const auto sample=projection.EvaluatePrimitive(nodes[i]);
    if(!sample.ok())return Return::Failure(
        sample.status.code,sample.status.message);
    out.nodeVelocityMPerS[i]=sample.value.radialVelocityMPerS;
  }
  for(std::size_t i=0;i<cells;++i) {
    const auto sample=projection.EvaluatePrimitive(
        ShellCentroid(nodes[i],nodes[i+1]));
    if(!sample.ok())return Return::Failure(
        sample.status.code,sample.status.message);
    out.cellDensityKgM3[i]=sample.value.densityKgM3;
    out.cellPressurePa[i]=sample.value.pressurePa;
    out.cellTransverseMagneticFieldT[i]=
        sample.value.transverseMagneticFieldT[0];
    out.cellTransverseMagneticField2T[i]=
        sample.value.transverseMagneticFieldT[1];
    out.cellRadialMagneticFieldT[i]=sample.value.radialMagneticFieldT;
  }
  return Return::Success(std::move(out));
}

} // namespace

struct PistonSheathModel::Ray {
  PistonRayInput input;
  PistonRayDisposition disposition=PistonRayDisposition::NoIntersection;
  std::shared_ptr<const PistonAmbientProjection> projection;
  std::shared_ptr<const PistonAmbientSourceTable> sourceTable;
  std::unique_ptr<PlanarPistonSolver> solver;
};

PistonSheathModel::~PistonSheathModel() = default;

namespace {

Core::Result<PlanarShockState> ExtractProductionShock(
    const PistonAmbientProjection& projection,
    const PlanarPistonSolver& solver) {
  using Return=Core::Result<PlanarShockState>;
  const auto raw=solver.DetectShock();
  if(!raw.ok())return Return::Failure(raw.status.code,raw.status.message);
  if(!raw.value.present)return raw;
  const auto cells=solver.Cells();
  if(!cells.ok())return Return::Failure(cells.status.code,cells.status.message);
  PlanarShockState out=raw.value;
  out.statesAvailable=false;
  // The production ambient changes by orders of magnitude across a handful
  // of low-coronal cells, so a broad absolute plateau average confounds the
  // radial background gradient with shock compression.  Estimate the local
  // downstream *perturbation relative to its own ambient*, then transport
  // that dimensionless perturbation to the Q-centroid.  The upstream limit is
  // the independently maintained ambient at that same radius.  No RH value is
  // used in this extraction.
  const auto upstream=projection.EvaluatePrimitive(out.radiusM);
  if(!upstream.ok())return Return::Failure(
      upstream.status.code,upstream.status.message);
  // The first cell immediately behind (but outside) the detected Q/P zone is
  // the specified downstream one-sided state.  Its perturbation is normalized
  // to the ambient at that cell and transported to the shock radius below;
  // using an absolute plateau average would confuse the steep coronal radial
  // gradient with compression.  A wider gap/window was tested during
  // CSWC0628 and increased the nonconvergent compression, so it is retained as
  // failed evidence rather than as a tuned production sampler.  No RH value
  // or compression cap enters this extraction.
  const int end=out.firstShockCell;
  const int begin=std::max(0,end-1);
  if(end<=begin)return Return::Success(out);
  double rhoRatio=0,pRatio=0,velocityDelta=0,btRatio=0;
  for(int cell=begin;cell<end;++cell) {
    const auto ambient=projection.EvaluatePrimitive(cells.value[cell].centerM);
    if(!ambient.ok())return Return::Failure(
        ambient.status.code,ambient.status.message);
    rhoRatio+=cells.value[cell].densityKgM3/ambient.value.densityKgM3;
    pRatio+=cells.value[cell].pressurePa/ambient.value.pressurePa;
    velocityDelta+=cells.value[cell].velocityMPerS-
        ambient.value.radialVelocityMPerS;
    const double ambientBt=std::hypot(
        ambient.value.transverseMagneticFieldT[0],
        ambient.value.transverseMagneticFieldT[1]);
    const double numericalBt=std::hypot(
        cells.value[cell].transverseMagneticFieldT,
        cells.value[cell].transverseMagneticField2T);
    if(!(ambientBt>0))return Return::Failure(
        Core::StatusCode::NumericalFailure,
        "production shock limit encountered zero transverse ambient field");
    btRatio+=numericalBt/ambientBt;
  }
  const double count=end-begin;
  rhoRatio/=count;pRatio/=count;velocityDelta/=count;btRatio/=count;
  out.upstreamDensityKgM3=upstream.value.densityKgM3;
  out.upstreamPressurePa=upstream.value.pressurePa;
  out.upstreamVelocityMPerS=upstream.value.radialVelocityMPerS;
  out.upstreamTransverseMagneticFieldT=
      upstream.value.transverseMagneticFieldT[0];
  out.upstreamTransverseMagneticField2T=
      upstream.value.transverseMagneticFieldT[1];
  out.upstreamRadialMagneticFieldT=upstream.value.radialMagneticFieldT;
  out.downstreamDensityKgM3=rhoRatio*out.upstreamDensityKgM3;
  out.downstreamPressurePa=pRatio*out.upstreamPressurePa;
  out.downstreamVelocityMPerS=out.upstreamVelocityMPerS+velocityDelta;
  out.downstreamTransverseMagneticFieldT=
      btRatio*out.upstreamTransverseMagneticFieldT;
  out.downstreamTransverseMagneticField2T=
      btRatio*out.upstreamTransverseMagneticField2T;
  out.downstreamRadialMagneticFieldT=out.upstreamRadialMagneticFieldT;
  out.compressionRatio=rhoRatio;
  const double denominator=out.downstreamDensityKgM3-out.upstreamDensityKgM3;
  if(!(denominator>0))return Return::Success(out);
  out.speedMPerS=(out.downstreamDensityKgM3*out.downstreamVelocityMPerS-
      out.upstreamDensityKgM3*out.upstreamVelocityMPerS)/denominator;
  out.statesAvailable=std::isfinite(out.speedMPerS)&&out.speedMPerS>
      out.upstreamVelocityMPerS&&rhoRatio>1;
  return Return::Success(out);
}

// Locate the outermost material cell whose total pressure differs from the
// frozen ambient authority by the configured dimensionless amount
// |P-P_a|/P_a.  The ambient denominator is part of the physical definition of
// the leading disturbance radius: using the evolved pressure instead would
// move the boundary merely because a cell became compressed.  This helper is
// shared by append admission, receipts and queries so all three interfaces
// classify the same finite material inventory.  No pressure floor is used;
// a non-positive ambient total pressure is a typed authority failure.
Core::Result<int> LastDisturbedCell(
    const PistonAmbientProjection& projection,
    const std::vector<PlanarCellState>& cells,double magneticPermeability,
    double relativeThreshold) {
  using Return=Core::Result<int>;
  int disturbed=-1;
  for(int cell=0;cell<static_cast<int>(cells.size());++cell) {
    const auto ambient=projection.EvaluatePrimitive(cells[cell].centerM);
    if(!ambient.ok())return Return::Failure(
        ambient.status.code,ambient.status.message);
    const double ambientBt2=
        ambient.value.transverseMagneticFieldT[0]*
            ambient.value.transverseMagneticFieldT[0]+
        ambient.value.transverseMagneticFieldT[1]*
            ambient.value.transverseMagneticFieldT[1];
    const double ambientTotal=ambient.value.pressurePa+
        ambientBt2/(2*magneticPermeability);
    if(!(ambientTotal>0)||!std::isfinite(ambientTotal))return Return::Failure(
        Core::StatusCode::InvalidState,
        "piston disturbance classifier requires positive ambient pressure");
    if(std::abs(cells[cell].totalPressurePa-ambientTotal)/ambientTotal>
        relativeThreshold)disturbed=cell;
  }
  return Return::Success(disturbed);
}

} // namespace

const char* Name(PistonBackgroundRegion region) noexcept {
  switch(region) {
    case PistonBackgroundRegion::UnsupportedFlank:return "unsupported-flank";
    case PistonBackgroundRegion::EjectaNotOwned:return "ejecta-not-owned";
    case PistonBackgroundRegion::Compression:return "compression";
    case PistonBackgroundRegion::Sheath:return "sheath";
    case PistonBackgroundRegion::Ambient:return "ambient";
  }
  return "unknown";
}

Core::Result<std::unique_ptr<PistonSheathModel>> PistonSheathModel::Create(
    std::shared_ptr<const EventConfiguration> event) {
  using Return=Core::Result<std::unique_ptr<PistonSheathModel>>;
  const PistonNumericsInput controls=event?event->pistonNumerics:
      PistonNumericsInput{};
  if(!event||!event->pistonContact.enabled||!event->pistonRays.enabled||
      event->sheathModel!="per-ray-lagrangian-piston-v1"||
      !controls.enabled||controls.artificialViscosity!="vnr"||
      controls.wellBalancedSources!="on"||
      controls.initialCells<32||!(controls.initialBufferM>0)||
      controls.sourceTablePoints<65||!(controls.trajectoryMaximumStepS>0)||
      !(controls.bufferCheckIntervalS>0)||
      controls.minimumBufferCells<4||controls.appendCells<4||
      !(controls.disturbanceRelativeThreshold>0&&
        controls.disturbanceRelativeThreshold<1)||
      !(controls.cfl>0&&controls.cfl<=0.8)||
      !(controls.quadraticViscosity>=0)||!(controls.linearViscosity>=0)||
      !(controls.shockThreshold>0))return Return::Failure(
          Core::StatusCode::InvalidConfiguration,
          "piston sheath controls or event closure are incomplete");
  const auto ambient=AmbientModel::Create(event);
  if(!ambient.ok())return Return::Failure(
      ambient.status.code,ambient.status.message);
  const auto contact=PistonContactModel::Create(event);
  if(!contact.ok())return Return::Failure(
      contact.status.code,contact.status.message);
  std::unique_ptr<PistonSheathModel> out(new PistonSheathModel);
  out->event_=std::move(event);out->ambient_=ambient.value;
  out->contact_=contact.value;out->controls_=controls;
  out->epochS_=out->event_->support.startS;out->generation_=1;

  for(const auto& rayInput:out->event_->pistonRays.rays) {
    std::unique_ptr<Ray> ray(new Ray);ray->input=rayInput;
    const auto contactState=out->contact_->EvaluateRay(
        rayInput.direction,out->epochS_);
    if(!contactState.ok())return Return::Failure(
        contactState.status.code,contactState.status.message);
    ray->disposition=contactState.value.disposition;
    if(!contactState.value.Supported()) {
      out->rays_.push_back(std::move(ray));continue;
    }
    const auto projection=PistonAmbientProjection::Create(out->ambient_,
        rayInput.direction,out->epochS_);
    if(!projection.ok())return Return::Failure(
        projection.status.code,projection.status.message);
    ray->projection=projection.value;
    const double inner=contactState.value.radiusM;
    const double outer=inner+controls.initialBufferM;
    if(!(outer<out->event_->support.coverageRadiusM))return Return::Failure(
        Core::StatusCode::OutOfDomain,
        "initial piston buffer exceeds ambient coverage");
    const auto sourceTable=PistonAmbientSourceTable::Create(projection.value,
        out->event_->support.firstValidPlasmaRadiusM,
        out->event_->support.coverageRadiusM*(1-1e-12),
        controls.sourceTablePoints,true);
    if(!sourceTable.ok())return Return::Failure(
        sourceTable.status.code,sourceTable.status.message);
    ray->sourceTable=sourceTable.value;
    const auto outerPath=PistonAmbientTrajectory::Create(projection.value,
        outer,out->epochS_,out->event_->support.endS,
        controls.trajectoryMaximumStepS);
    if(!outerPath.ok())return Return::Failure(
        outerPath.status.code,outerPath.status.message);
    std::vector<double> nodes(controls.initialCells+1);
    for(int i=0;i<=controls.initialCells;++i)
      nodes[i]=inner+(outer-inner)*i/controls.initialCells;
    auto state=SampleAmbientColumn(*projection.value,nodes);
    if(!state.ok())return Return::Failure(state.status.code,state.status.message);
    // The selected startup admits a small ambient/contact mismatch bounded by
    // epsilon_start.  The boundary node itself belongs to the piston and must
    // use the analytical contact velocity exactly.
    state.value.nodeVelocityMPerS.front()=contactState.value.radialSpeedMPerS;
    const auto piston=[contact=out->contact_,direction=rayInput.direction](double t) {
      const auto value=contact->EvaluateRay(direction,t);
      if(!value.ok())return Core::Result<KinematicValue>::Failure(
          value.status.code,value.status.message);
      if(!value.value.Supported())return Core::Result<KinematicValue>::Failure(
          Core::StatusCode::OutOfDomain,
          "contact-driven tube became unsupported");
      return Core::Result<KinematicValue>::Success({value.value.radiusM,
          value.value.radialSpeedMPerS,value.value.radialAccelerationMPerS2});
    };
    const auto source=[table=sourceTable.value](double radius,double) {
      return table->Evaluate(radius);
    };
    const auto outerHistory=[path=outerPath.value](double t) {
      return path->Evaluate(t);
    };
    const auto first=projection.value->EvaluatePrimitive(inner),
        last=projection.value->EvaluatePrimitive(outer);
    if(!first.ok()||!last.ok())return Return::Failure(
        Core::StatusCode::OutOfDomain,"piston tube endpoints lack ambient state");
    PlanarPistonInput input;
    input.geometry=PistonTubeGeometry::RadialSpherical;
    input.startS=out->epochS_;input.endS=out->event_->support.endS;
    input.leftPositionM=inner;input.columnLengthM=outer-inner;
    input.solidAngleSr=rayInput.solidAngleSr;
    input.initialDensityKgM3=first.value.densityKgM3;
    input.initialPressurePa=first.value.pressurePa;
    input.initialVelocityMPerS=last.value.radialVelocityMPerS;
    input.initialTransverseMagneticFieldT=first.value.transverseMagneticFieldT[0];
    input.initialTransverseMagneticField2T=first.value.transverseMagneticFieldT[1];
    input.initialRadialMagneticFieldT=first.value.radialMagneticFieldT;
    input.gammaAdiabatic=out->event_->composition.gammaAdiabatic;
    input.cells=controls.initialCells;input.cfl=controls.cfl;
    input.quadraticViscosity=controls.quadraticViscosity;
    input.linearViscosity=controls.linearViscosity;
    input.shockThreshold=controls.shockThreshold;
    auto solver=PlanarPistonSolver::CreateInitialized(input,piston,
        std::move(state.value),source,outerHistory);
    if(!solver.ok())return Return::Failure(
        solver.status.code,solver.status.message);
    ray->solver=std::move(solver.value);
    out->rays_.push_back(std::move(ray));
  }
  return Return::Success(std::move(out));
}

Core::Status PistonSheathModel::AdvanceTo(double epoch) {
  if(!std::isfinite(epoch)||epoch<epochS_||epoch>event_->support.endS)
    return Core::Status::Failure(Core::StatusCode::OutOfDomain,
        "piston sheath epoch is outside forward event support");
  std::vector<std::unique_ptr<PlanarPistonSolver>> candidates;
  candidates.reserve(rays_.size());
  for(const auto& ray:rays_)
    candidates.push_back(ray->solver?ray->solver->Clone():nullptr);

  // Buffer maintenance is scheduled on an absolute event-time lattice, not
  // at caller publication epochs.  Otherwise AdvanceTo(1800) and
  // AdvanceTo(600),AdvanceTo(1200),AdvanceTo(1800) append at different times
  // and evolve physically different outer domains.  All lattice substeps stay
  // private; a later ray/append failure still rolls back the entire requested
  // epoch.  A non-lattice target is advanced after the final maintenance time
  // without introducing a caller-dependent append decision.
  const double start=event_->support.startS;
  const double interval=controls_.bufferCheckIntervalS;
  long long nextIndex=static_cast<long long>(
      std::floor((epochS_-start)/interval+1e-12))+1;
  std::vector<std::pair<double,bool>> stops;
  for(;;++nextIndex) {
    const double scheduled=start+nextIndex*interval;
    if(scheduled>epoch+1e-10*std::max(1.0,std::abs(epoch)))break;
    stops.push_back({std::min(scheduled,epoch),true});
    if(scheduled>=epoch)break;
  }
  if(stops.empty()||stops.back().first<epoch)stops.push_back({epoch,false});

  for(const auto& stop:stops) {
    for(std::size_t index=0;index<rays_.size();++index) {
      const auto& ray=rays_[index];
      auto& candidate=candidates[index];
      if(!candidate)continue;
      const auto advanced=candidate->AdvanceTo(stop.first);
      if(!advanced.ok())return advanced;
      if(!stop.second||stop.first>=event_->support.endS)continue;
      const auto cells=candidate->Cells();
      if(!cells.ok())return cells.status;
      const auto disturbedResult=LastDisturbedCell(*ray->projection,cells.value,
          candidate->Input().magneticPermeabilityNPerA2,
          controls_.disturbanceRelativeThreshold);
      if(!disturbedResult.ok())return disturbedResult.status;
      const int disturbed=disturbedResult.value;
      const std::size_t ahead=disturbed<0?cells.value.size():
          cells.value.size()-static_cast<std::size_t>(disturbed+1);
      if(ahead>=static_cast<std::size_t>(controls_.minimumBufferCells))continue;

      const auto oldNodes=candidate->NodePositionsM();
      const auto& oldVelocities=candidate->NodeVelocitiesMPerS();
      const double width=oldNodes.back()-oldNodes[oldNodes.size()-2];
      const double newOuter=oldNodes.back()+controls_.appendCells*width;
      if(!(newOuter<event_->support.coverageRadiusM))return Core::Status::Failure(
          Core::StatusCode::OutOfDomain,
          "piston sheath cannot append the required ambient buffer inside coverage");
      std::vector<double> nodes(controls_.appendCells+1);
      for(int node=0;node<=controls_.appendCells;++node)
        nodes[node]=oldNodes.back()+(newOuter-oldNodes.back())*node/
            controls_.appendCells;
      auto state=SampleAmbientColumn(*ray->projection,nodes);
      if(!state.ok())return state.status;
      state.value.nodeVelocityMPerS.front()=oldVelocities.back();
      const auto outerPath=PistonAmbientTrajectory::Create(ray->projection,
          newOuter,stop.first,event_->support.endS,
          controls_.trajectoryMaximumStepS);
      if(!outerPath.ok())return outerPath.status;
      PistonAppendState append;
      append.nodePositionM=std::move(nodes);append.state=std::move(state.value);
      append.outerBoundary=[path=outerPath.value](double t) {
        return path->Evaluate(t);
      };
      const auto appended=candidate->AppendAmbient(std::move(append));
      if(!appended.ok())return appended;
    }
  }
  for(std::size_t i=0;i<rays_.size();++i)
    if(candidates[i])rays_[i]->solver=std::move(candidates[i]);
  epochS_=epoch;++generation_;
  return Core::Status::Success();
}

Core::Result<std::vector<PistonRayReceipt>> PistonSheathModel::Receipts() const {
  using Return=Core::Result<std::vector<PistonRayReceipt>>;
  std::vector<PistonRayReceipt> result;
  result.reserve(rays_.size());
  for(const auto& ray:rays_) {
    PistonRayReceipt receipt;
    receipt.rayId=ray->input.id;receipt.disposition=ray->disposition;
    receipt.solidAngleSr=ray->input.solidAngleSr;
    const auto contact=contact_->EvaluateRay(ray->input.direction,epochS_);
    if(!contact.ok())return Return::Failure(
        contact.status.code,contact.status.message);
    if(contact.value.Supported())receipt.contactRadiusM=contact.value.radiusM;
    if(ray->solver) {
      const auto cells=ray->solver->Cells();
      const auto nodes=ray->solver->NodePositionsM();
      const auto shock=ExtractProductionShock(*ray->projection,*ray->solver);
      if(!cells.ok()||!shock.ok())return Return::Failure(
          Core::StatusCode::InvalidState,"piston ray receipt is unavailable");
      receipt.cells=cells.value.size();receipt.outerRadiusM=nodes.back();
      receipt.shock=shock.value;
      const auto disturbedResult=LastDisturbedCell(*ray->projection,cells.value,
          ray->solver->Input().magneticPermeabilityNPerA2,
          controls_.disturbanceRelativeThreshold);
      if(!disturbedResult.ok())return Return::Failure(
          disturbedResult.status.code,disturbedResult.status.message);
      const int disturbed=disturbedResult.value;
      if(disturbed>=0)receipt.disturbanceRadiusM=cells.value[disturbed].centerM;
      receipt.ambientBufferCells=disturbed<0?cells.value.size():
          cells.value.size()-static_cast<std::size_t>(disturbed+1);
    }
    result.push_back(receipt);
  }
  return Return::Success(std::move(result));
}

Core::Result<PistonBackgroundSample> PistonSheathModel::QueryRay(
    std::uint64_t rayId,double radius) const {
  using Return=Core::Result<PistonBackgroundSample>;
  const auto found=std::find_if(rays_.begin(),rays_.end(),[&](const auto& ray) {
    return ray->input.id==rayId;
  });
  if(found==rays_.end()||!std::isfinite(radius))return Return::Failure(
      Core::StatusCode::OutOfDomain,"piston ray query is outside identity/support");
  const Ray& ray=**found;
  PistonBackgroundSample out;
  const auto contact=contact_->EvaluateRay(ray.input.direction,epochS_);
  if(!contact.ok())return Return::Failure(contact.status.code,contact.status.message);
  out.contact=contact.value;
  if(!ray.solver) {
    out.region=PistonBackgroundRegion::UnsupportedFlank;
    return Return::Success(out);
  }
  if(radius<contact.value.radiusM) {
    out.region=PistonBackgroundRegion::EjectaNotOwned;
    return Return::Success(out);
  }
  const auto nodes=ray.solver->NodePositionsM();
  const auto shock=ExtractProductionShock(*ray.projection,*ray.solver);
  if(!shock.ok())return Return::Failure(shock.status.code,shock.status.message);
  out.shock=shock.value;
  if(radius<=nodes.back()) {
    auto upper=std::upper_bound(nodes.begin(),nodes.end(),radius);
    const std::size_t cell=upper==nodes.begin()?0:
        std::min(static_cast<std::size_t>(upper-nodes.begin()-1),nodes.size()-2);
    const auto cells=ray.solver->Cells();
    if(!cells.ok())return Return::Failure(cells.status.code,cells.status.message);
    const auto disturbed=LastDisturbedCell(*ray.projection,cells.value,
        ray.solver->Input().magneticPermeabilityNPerA2,
        controls_.disturbanceRelativeThreshold);
    if(!disturbed.ok())return Return::Failure(
        disturbed.status.code,disturbed.status.message);
    out.plasma=cells.value[cell];out.plasmaAvailable=true;
    const double disturbanceRadius=disturbed.value>=0?
        cells.value[disturbed.value].centerM:contact.value.radiusM;
    // R_sh is a Q-weighted shock centroid, not the outer edge of a zero-width
    // discontinuity.  Cells between that centroid and r_d are the resolved
    // leading compression/shock transition.  Cells beyond r_d remain evolved
    // material (so conservation is retained) but are typed ambient because
    // they satisfy the same authority-relative disturbance criterion used by
    // outer-buffer admission.
    if(shock.value.present&&radius<=shock.value.radiusM)
      out.region=PistonBackgroundRegion::Sheath;
    else if(radius<=disturbanceRadius)
      out.region=PistonBackgroundRegion::Compression;
    else out.region=PistonBackgroundRegion::Ambient;
    return Return::Success(out);
  }
  const auto ambient=ambient_->Evaluate(radius*ray.input.direction,epochS_);
  if(!ambient.ok())return Return::Failure(ambient.status.code,ambient.status.message);
  out.region=PistonBackgroundRegion::Ambient;out.ambient=ambient.value;
  out.plasmaAvailable=true;
  return Return::Success(out);
}

} } // namespace SEP::CoronaSwcme
