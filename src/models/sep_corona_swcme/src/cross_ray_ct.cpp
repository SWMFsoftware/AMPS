#include "sep_corona_swcme/cross_ray_ct.h"

#include "sep_coronal_cme/configuration_parser.h"

#include <cmath>
#include <limits>
#include <locale>
#include <new>
#include <sstream>

namespace SEP { namespace CoronaSwcme { namespace CrossRay { namespace {

// Bound floating-point arithmetic, not physical divergence. The positive
// bound is rounded upward, including a possible subnormal final conversion;
// this does not replace or modify any magnetic state value. Actual raw div B
// must still meet the separate production physical/numerical acceptance gates.
double ErrorBound(long double magnitude,std::size_t operations) {
  const long double eps=std::numeric_limits<double>::epsilon();
  const long double n=operations+4;
  const long double bound=n*eps/(1-n*eps)*magnitude;
  if(bound==0)return 0;
  return std::nextafter(static_cast<double>(bound),
      std::numeric_limits<double>::infinity());
}

bool Finite(const std::vector<double>& values) {
  for(double value:values)if(!std::isfinite(value))return false;
  return true;
}

} // namespace

Core::Result<std::shared_ptr<const StaticMagneticFlux>> StaticMagneticFlux::Publish(
    std::shared_ptr<const StaticSphericalMesh> mesh,std::vector<double> flux,
    std::vector<double> uncertainty,double epoch,std::uint64_t generation) {
  using Return=Core::Result<std::shared_ptr<const StaticMagneticFlux>>;
  if(!mesh||flux.size()!=mesh->Faces().size()||uncertainty.size()!=flux.size()||
      !std::isfinite(epoch)||generation==0)return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"CT candidate has incomplete geometry/epoch ownership");
  if(!Finite(flux)||!Finite(uncertainty))return Return::Failure(
      Core::StatusCode::NumericalFailure,"CT flux or arithmetic uncertainty is nonfinite");
  for(double value:uncertainty)if(value<0)return Return::Failure(
      Core::StatusCode::InvalidState,"CT uncertainty cannot be negative");
  std::shared_ptr<StaticMagneticFlux> out(new StaticMagneticFlux);
  out->mesh_=std::move(mesh);out->flux_=std::move(flux);
  out->uncertainty_=std::move(uncertainty);out->epochS_=epoch;out->generation_=generation;
  for(const auto& cell:out->mesh_->Cells()) {
    long double closed=0,bound=0,magnitude=0;std::size_t terms=0;
    for(const auto& face:cell.faces)if(!face.collapsedGeometricFace) {
      closed+=face.outwardOrientation*static_cast<long double>(out->flux_[face.index]);
      bound+=out->uncertainty_[face.index];magnitude+=std::abs(out->flux_[face.index]);
      ++terms;
    }
    bound+=ErrorBound(magnitude,terms);
    if(std::abs(closed)>bound)return Return::Failure(
        Core::StatusCode::DataIntegrityFailure,
        "CT cell has a magnetic monopole beyond certified arithmetic uncertainty; no cleaning applied");
    const double divergence=static_cast<double>(closed/cell.measure.volumeM3);
    if(!std::isfinite(divergence))return Return::Failure(
        Core::StatusCode::NumericalFailure,"CT raw divergence is not representable");
    out->divergence_.push_back(divergence);
  }
  return Return::Success(std::move(out));
}

Core::Result<std::shared_ptr<const StaticMagneticFlux>> StaticMagneticFlux::FromFaceFlux(
    std::shared_ptr<const StaticSphericalMesh> mesh,std::vector<double> flux) {
  using Return=Core::Result<std::shared_ptr<const StaticMagneticFlux>>;
  if(!mesh||flux.size()!=mesh->Faces().size())return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"CT initial flux inventory does not match its mesh");
  if(!Finite(flux))return Return::Failure(Core::StatusCode::InvalidState,
      "CT initial face flux is nonfinite");
  std::vector<double> bound;bound.reserve(flux.size());
  for(double value:flux)bound.push_back(ErrorBound(std::abs(value),4));
  const double epoch=mesh->EpochS();const auto generation=mesh->Generation();
  return Publish(std::move(mesh),std::move(flux),std::move(bound),epoch,generation);
}

Core::Result<std::shared_ptr<const StaticMagneticFlux>> StaticMagneticFlux::FromEdgePotential(
    std::shared_ptr<const StaticSphericalMesh> mesh,const std::vector<double>& potential) {
  using Return=Core::Result<std::shared_ptr<const StaticMagneticFlux>>;
  if(!mesh||potential.size()!=mesh->Edges().size())return Return::Failure(
      Core::StatusCode::InvalidConfiguration,"CT initial edge-potential inventory is incomplete");
  if(!Finite(potential))return Return::Failure(Core::StatusCode::InvalidState,
      "CT initial edge potential is nonfinite");
  std::vector<double> flux,bound;
  for(const auto& face:mesh->Faces()) {
    long double curl=0,magnitude=0;
    for(const auto& edge:face.positiveEdgeLoop) {
      curl+=edge.orientation*static_cast<long double>(potential[edge.index]);
      magnitude+=std::abs(potential[edge.index]);
    }
    flux.push_back(static_cast<double>(curl));
    bound.push_back(ErrorBound(magnitude,face.positiveEdgeLoop.size()));
  }
  const double epoch=mesh->EpochS();const auto generation=mesh->Generation();
  return Publish(std::move(mesh),std::move(flux),std::move(bound),epoch,generation);
}

Core::Result<std::shared_ptr<const StaticMagneticFlux>> StaticMagneticFlux::Advance(
    const EdgeEmfInterval& e) const {
  using Return=Core::Result<std::shared_ptr<const StaticMagneticFlux>>;
  if(e.meshFingerprint!=mesh_->Fingerprint()||e.eventFingerprint!=mesh_->EventFingerprint()||
      e.fieldGeneration!=generation_||e.startS!=epochS_)return Return::Failure(
          Core::StatusCode::InvalidState,"stale or foreign CT edge-EMF authority/epoch/generation");
  if(!std::isfinite(e.endS)||!(e.endS>e.startS)||
      e.valuesVolt.size()!=mesh_->Edges().size())return Return::Failure(
          Core::StatusCode::InvalidConfiguration,"CT interval or edge-EMF inventory is incomplete");
  if(!Finite(e.valuesVolt))return Return::Failure(Core::StatusCode::InvalidState,
      "CT edge EMF is nonfinite");
  if(generation_==std::numeric_limits<std::uint64_t>::max())return Return::Failure(
      Core::StatusCode::OutOfDomain,"CT generation increment overflows");
  const double step=e.endS-e.startS;
  if(!std::isfinite(step))return Return::Failure(Core::StatusCode::NumericalFailure,
      "CT interval duration overflows");
  std::vector<double> candidate,bound;
  candidate.reserve(flux_.size());bound.reserve(flux_.size());
  for(std::size_t f=0;f<mesh_->Faces().size();++f) {
    const auto& face=mesh_->Faces()[f];long double circulation=0,magnitude=0;
    for(const auto& edge:face.positiveEdgeLoop) {
      circulation+=edge.orientation*static_cast<long double>(e.valuesVolt[edge.index]);
      magnitude+=std::abs(e.valuesVolt[edge.index]);
    }
    // Faraday: Phi^{n+1}=Phi^n-dt integral_boundary E.dl, SI T m^2.
    // The SAME canonical edge value enters every incident face. Real-arithmetic
    // cell divergence cancels by integer incidence; no component projection,
    // subtraction of a measured residual, or field/energy reset is permitted.
    candidate.push_back(static_cast<double>(static_cast<long double>(flux_[f])-
        static_cast<long double>(step)*circulation));
    const long double error=static_cast<long double>(uncertainty_[f])+
        ErrorBound(std::abs(flux_[f])+std::abs(static_cast<long double>(step))*magnitude,
            face.positiveEdgeLoop.size()+2);
    bound.push_back(error==0?0:std::nextafter(static_cast<double>(error),
        std::numeric_limits<double>::infinity()));
  }
  return Publish(mesh_,std::move(candidate),std::move(bound),e.endS,generation_+1);
}

std::string StaticMagneticFlux::Fingerprint() const {
  // No mutable hash cache: immutable states can be queried from independent
  // readers without introducing a hidden thread-dependent publication epoch.
  std::ostringstream bytes;bytes.imbue(std::locale::classic());
  bytes<<"static-CT-face-flux-v1\n"<<mesh_->Fingerprint()<<'\n'<<
      mesh_->EventFingerprint()<<'\n'<<generation_<<'\n'<<std::hexfloat<<epochS_<<'\n';
  for(std::size_t i=0;i<flux_.size();++i)bytes<<flux_[i]<<','<<uncertainty_[i]<<'\n';
  return CoronalCME::ComputeContentChecksum(bytes.str());
}

} } } // namespace SEP::CoronaSwcme::CrossRay
