#include "sep_corona_swcme/cross_ray_fluid.h"

#include "sep_coronal_cme/configuration_parser.h"

#include <cmath>
#include <locale>
#include <sstream>

namespace SEP { namespace CoronaSwcme { namespace CrossRay { namespace {

using CoronalCME::Vec3;

bool Finite(Vec3 v) {
  return std::isfinite(v.x)&&std::isfinite(v.y)&&std::isfinite(v.z);
}

// Wider accumulators reduce summation error without altering any physical
// value or closing a budget by subtraction. All public ledgers retain raw SI
// residuals; native/global qualification must apply its own frozen bounds.
struct Accumulator {
  long double mass=0,momentum[3]{},energy=0;
  void Add(double rho,Vec3 p,double e,long double multiplier) {
    mass+=multiplier*rho;momentum[0]+=multiplier*p.x;
    momentum[1]+=multiplier*p.y;momentum[2]+=multiplier*p.z;
    energy+=multiplier*e;
  }
  ConservedInventory Value() const {
    return {static_cast<double>(mass),{static_cast<double>(momentum[0]),
        static_cast<double>(momentum[1]),static_cast<double>(momentum[2])},
        static_cast<double>(energy)};
  }
};

ConservedInventory Inventory(const StaticSphericalMesh& mesh,
    const std::vector<FluidDensity>& cells) {
  Accumulator sum;
  for(std::size_t c=0;c<cells.size();++c)sum.Add(cells[c].massKgM3,
      cells[c].momentumKgM2S,cells[c].totalEnergyJPerM3,mesh.Cells()[c].measure.volumeM3);
  return sum.Value();
}

bool Finite(const ConservedInventory& value) {
  return std::isfinite(value.massKg)&&Finite(value.momentumKgMPerS)&&
      std::isfinite(value.totalEnergyJ);
}

} // namespace

Core::Result<std::shared_ptr<const StaticFluidState>> StaticFluidState::Create(
    std::shared_ptr<const StaticMagneticFlux> magnetic,
    std::vector<FluidDensity> cells,Thermodynamics gas) {
  using Return=Core::Result<std::shared_ptr<const StaticFluidState>>;
  if(!magnetic||cells.size()!=magnetic->Mesh()->Cells().size())return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"fluid cells do not match their authoritative CT mesh");
  const auto reconstructed=ReconstructCellField(*magnetic);
  if(!reconstructed.ok())return Return::Failure(
      reconstructed.status.code,reconstructed.status.message);
  std::shared_ptr<StaticFluidState> out(new StaticFluidState);
  out->magnetic_=std::move(magnetic);out->cells_=std::move(cells);out->gas_=gas;
  for(std::size_t c=0;c<out->cells_.size();++c) {
    const auto& value=out->cells_[c];
    const auto primitive=ToPrimitive({value.massKgM3,value.momentumKgM2S,
        value.totalEnergyJPerM3,reconstructed.value[c].fieldHciT},gas);
    if(!primitive.ok())return Return::Failure(primitive.status.code,
        "fluid/CT candidate rejected: "+primitive.status.message);
    out->primitive_.push_back(primitive.value);
  }
  const auto inventory=Inventory(*out->magnetic_->Mesh(),out->cells_);
  if(!Finite(inventory))return Return::Failure(Core::StatusCode::NumericalFailure,
      "fluid extensive inventory overflows SI ledger");
  out->ledger_.initial=inventory;out->ledger_.final=inventory;
  return Return::Success(std::move(out));
}

Core::Result<std::shared_ptr<const StaticFluidState>> StaticFluidState::Advance(
    const FluidFluxInterval& interval) const {
  using Return=Core::Result<std::shared_ptr<const StaticFluidState>>;
  const auto& mesh=*magnetic_->Mesh();
  if(interval.faces.size()!=mesh.Faces().size()||
      interval.sources.size()!=mesh.Cells().size())return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"fluid face/source inventory is incomplete");
  // CT validates event/mesh ownership, exact epoch/generation, all edge values
  // and the positive finite interval before any private fluid candidate exists.
  const auto nextField=magnetic_->Advance(interval.magnetic);
  if(!nextField.ok())return Return::Failure(nextField.status.code,nextField.status.message);
  const long double dt=interval.magnetic.endS-interval.magnetic.startS;
  std::vector<Accumulator> change(cells_.size());Accumulator boundary,source;
  for(std::size_t f=0;f<mesh.Faces().size();++f) {
    const auto& flux=interval.faces[f];const auto& face=mesh.Faces()[f];
    if(!std::isfinite(flux.massKgPerS)||!Finite(flux.momentumN)||
        !std::isfinite(flux.totalEnergyW))return Return::Failure(
            Core::StatusCode::InvalidState,"nonfinite integrated fluid face flux");
    if(face.leftCell!=NoMeshIndex)change[face.leftCell].Add(
        flux.massKgPerS,flux.momentumN,flux.totalEnergyW,-dt);
    if(face.rightCell!=NoMeshIndex)change[face.rightCell].Add(
        flux.massKgPerS,flux.momentumN,flux.totalEnergyW,dt);
    // Domain outward is +canonical at a missing right neighbor and -canonical
    // at a missing left neighbor. Internal transfers cancel identically and
    // must not be counted as external mass/work supplied to the CME region.
    if(face.leftCell==NoMeshIndex||face.rightCell==NoMeshIndex)
      boundary.Add(flux.massKgPerS,flux.momentumN,flux.totalEnergyW,
          face.rightCell==NoMeshIndex?dt:-dt);
  }
  for(std::size_t c=0;c<cells_.size();++c) {
    const auto& value=interval.sources[c];
    if(!Finite(value.forceN)||!std::isfinite(value.totalPowerW))return Return::Failure(
        Core::StatusCode::InvalidState,"nonfinite prescribed volume force/power");
    change[c].Add(0,value.forceN,value.totalPowerW,dt);
    source.Add(0,value.forceN,value.totalPowerW,dt);
  }
  std::vector<FluidDensity> candidate;candidate.reserve(cells_.size());
  for(std::size_t c=0;c<cells_.size();++c) {
    const auto& old=cells_[c];const long double v=mesh.Cells()[c].measure.volumeM3;
    // Static volume only: U_new=U_old+Delta_extensive/V. A moving mesh cannot
    // enter through this formula; it requires its actual new volume and GCL.
    candidate.push_back({static_cast<double>(old.massKgM3+change[c].mass/v),
        {static_cast<double>(old.momentumKgM2S.x+change[c].momentum[0]/v),
         static_cast<double>(old.momentumKgM2S.y+change[c].momentum[1]/v),
         static_cast<double>(old.momentumKgM2S.z+change[c].momentum[2]/v)},
        static_cast<double>(old.totalEnergyJPerM3+change[c].energy/v)});
  }
  // Publication validates pressure with NEW B and the conservatively updated
  // E. There is intentionally no E += (|B_new|^2-|B_old|^2)/(2mu0) reset.
  const auto prepared=Create(nextField.value,std::move(candidate),gas_);
  if(!prepared.ok())return prepared;
  std::shared_ptr<StaticFluidState> out(new StaticFluidState(*prepared.value));
  out->ledger_.initial=Inventory(mesh,cells_);
  out->ledger_.outwardBoundaryTransfer=boundary.Value();out->ledger_.volumeAddition=source.Value();
  Accumulator residual;
  for(const auto& term:{std::make_pair(out->ledger_.final,1.0),
                       std::make_pair(out->ledger_.initial,-1.0),
                       std::make_pair(out->ledger_.outwardBoundaryTransfer,1.0),
                       std::make_pair(out->ledger_.volumeAddition,-1.0)})
    residual.Add(term.first.massKg,term.first.momentumKgMPerS,term.first.totalEnergyJ,term.second);
  out->ledger_.residual=residual.Value();
  if(!Finite(out->ledger_.outwardBoundaryTransfer)||!Finite(out->ledger_.volumeAddition)||
      !Finite(out->ledger_.residual))return Return::Failure(Core::StatusCode::NumericalFailure,
          "fluid boundary/source/residual budget overflows SI ledger");
  return Return::Success(std::move(out));
}

std::string StaticFluidState::Fingerprint() const {
  std::ostringstream bytes;bytes.imbue(std::locale::classic());
  bytes<<"static-conservative-fluid-CT-v1\n"<<magnetic_->Fingerprint()<<'\n'<<std::hexfloat<<
      gas_.gammaAdiabatic<<','<<gas_.magneticPermeabilityNPerA2<<'\n';
  for(const auto& cell:cells_)bytes<<cell.massKgM3<<','<<cell.momentumKgM2S.x<<','<<
      cell.momentumKgM2S.y<<','<<cell.momentumKgM2S.z<<','<<cell.totalEnergyJPerM3<<'\n';
  return CoronalCME::ComputeContentChecksum(bytes.str());
}

} } } // namespace SEP::CoronaSwcme::CrossRay
