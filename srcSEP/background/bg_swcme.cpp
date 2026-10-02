// Canonical SWCME sampling and numerical derivatives live here; no AMPS
// center-node writes or MPI operations are permitted in this translation unit.
#include "bg_swcme.h"
#include "bg_parker.h"
#include "background_snapshot.h"
#include "swcme3d_input.hpp"
#include "sep_background_snapshot.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>

namespace SEP3D { namespace Background {
namespace {
Core::Status Invalid(const std::string& text) {
  return Core::Status(Core::StatusCode::BackgroundInvalid,text);
}
std::uint64_t Digest(const std::string& text) {
  // Compact per-sample FNV tag for interpolation provenance. The complete
  // configuration fingerprint remains in metadata and is checked at publication;
  // this digest is not a substitute for a scientific configuration manifest.
  std::uint64_t value=1469598103934665603ULL;
  for (unsigned char c:text) { value^=c; value*=1099511628211ULL; }
  return value;
}
// Model-neutral scratch state: B [T], U [m/s], electron density n [m^-3],
// fixed-composition mass density rho [kg/m^3], total thermal pressure p [Pa].
struct Primitive { Core::Vec3 B,U; double n=0,rho=0,p=0; };
}
struct SwcmeBackgroundProvider::Implementation {
  swcme::input3d::ResolvedConfiguration configuration;
  ParkerConfiguration ambient;
  double dt,cadence;
  Core::Vec3 origin;
  swcme3d::Model model;
  // Only candidate Prepared objects are mutable. Once installed, evaluation
  // sees one immutable StepState/metadata/inner-provider tuple for the epoch.
  struct Prepared {
    swcme3d::StepState step;
    SnapshotMetadata metadata;
    std::shared_ptr<AnalyticParkerProvider> inner;
  };
  std::shared_ptr<const Prepared> prepared;
  Implementation(const swcme::input3d::ResolvedConfiguration& c,
      const ParkerConfiguration& a,double d,double interval,const Core::Vec3& o)
      :configuration(c),ambient(a),dt(d),cadence(interval),origin(o),model(c.model) {}
  Core::Status Primitives(const std::vector<Core::Vec3>& points,
                         std::vector<Primitive>* output) const {
    // Internal scratch operation, called only after successful Prepare.
    // Partition below-domain points to the explicit inner ambient provider;
    // preserve original order through indices for one canonical outer batch.
    // Its output may contain a prefix on failure, so callers must check status.
    output->resize(points.size());
    std::vector<double> x,y,z;
    std::vector<std::size_t> indices;
    for (std::size_t i=0;i<points.size();++i) {
      const Core::Vec3 local=points[i]-origin;
      if (local.Norm()<swcme::solarwind::MIN_RADIUS_M) {
        // The canonical CME evaluator starts at 1.05 Rs. The explicitly
        // declared inner ambient continuation covers the photospheric shell.
        // Prepare proves the entire CME lies above this handoff, so no heated
        // or compressed CME cell is ever replaced by an ambient fallback.
        const auto a=prepared->inner->Evaluate(local);
        if (!a.status.ok()) return a.status;
        // Recover the inner provider's composition-consistent rho from
        // v_A=|B|/sqrt(mu0*rho), rather than treating n_e as proton density.
        (*output)[i]={a.B,a.U,a.numberDensityM3,
            a.absB*a.absB/(4e-7*Core::Const::kPi*a.alfvenSpeedMpS*a.alfvenSpeedMpS),a.pressurePa};
      } else {
        x.push_back(local.x);y.push_back(local.y);z.push_back(local.z);indices.push_back(i);
      }
    }
    std::array<std::vector<double>,9> values;
    // Structure-of-arrays order follows the canonical API exactly:
    // n, Ux, Uy, Uz, Bx, By, Bz, rho, p. A primitive is assembled only after
    // the checked batch succeeds; no application buffer is exposed here.
    for (auto& v:values) v.resize(indices.size());
    const auto status=model.evaluate_cartesian_primitive_checked(prepared->step,
        x.data(),y.data(),z.data(),values[0].data(),values[1].data(),values[2].data(),
        values[3].data(),values[4].data(),values[5].data(),values[6].data(),
        values[7].data(),values[8].data(),indices.size());
    if (!status.ok()) {
      // Preserve the canonical reason (e.g. WRONG_BRANCH), failing batch
      // index and heliocentric SI position. A bare SHOCK_SOLVER_FAILURE hides
      // the direction needed to reproduce a rank-local initialization error.
      // This routine also handles derivative stencils, so report the actual
      // queried point rather than guessing its parent physical-cell center.
      std::ostringstream diagnostic; diagnostic.precision(17);
      diagnostic << "SWCME mesh primitives: " << status.summary();
      if (status.sample_index<indices.size()) {
        const std::size_t k=status.sample_index;
        diagnostic << "; epoch_s=" << prepared->metadata.epochS
                   << "; local_position_m=(" << x[k] << ',' << y[k] << ',' << z[k] << ')';
      }
      return Invalid(diagnostic.str());
    }
    for (std::size_t k=0;k<indices.size();++k)
      (*output)[indices[k]]={{values[4][k],values[5][k],values[6][k]},
          {values[1][k],values[2][k],values[3][k]},values[0][k],values[7][k],values[8][k]};
    return Core::Status::OK();
  }
};
SwcmeBackgroundProvider::SwcmeBackgroundProvider(
    const swcme::input3d::ResolvedConfiguration& model,
    const ParkerConfiguration& ambient,double dt,double cadence,const Core::Vec3& origin)
    :implementation_(new Implementation(model,ambient,dt,cadence,origin)) {}
SwcmeBackgroundProvider::~SwcmeBackgroundProvider()=default;
Core::Status SwcmeBackgroundProvider::Validate() const {
  // Application/parser and canonical resolver freeze model options before
  // construction. This adapter additionally requires a finite local clock,
  // an origin and spherical support for the current handoff/mesh contract.
  const auto& p=*implementation_;
  if (!std::isfinite(p.dt)||p.dt<=0||!std::isfinite(p.cadence)||p.cadence<p.dt||
      !std::isfinite(p.origin.x)||!std::isfinite(p.origin.y)||!std::isfinite(p.origin.z)||
      p.configuration.model.shape!=swcme3d::ShockShape::Sphere)
    return Invalid("SWCME background requires finite clock/origin and spherical geometry");
  return AnalyticParkerProvider(p.ambient).Validate();
}
Core::Status SwcmeBackgroundProvider::Prepare(double timeS) {
  auto& p=*implementation_;
  const auto valid=Validate();if(!valid.ok()) return valid;
  if (!std::isfinite(timeS)||timeS<p.configuration.valid_from_s||
      timeS>p.configuration.valid_until_s)
    return Invalid("SWCME background epoch is outside declared event coverage");
  const long double tick=(static_cast<long double>(timeS)-p.configuration.valid_from_s)/p.dt;
  // Derive generation from the declared clock, not the number of calls.
  // This keeps ranks and repeated same-epoch preparation consistent. Round
  // to the nearest tick in long double and reserve generation zero as invalid.
  if (tick<0||tick>static_cast<long double>(UINT64_MAX-1)) return Invalid("SWCME background tick overflow");
  try {
    std::shared_ptr<Implementation::Prepared> candidate(new Implementation::Prepared);
    candidate->step=p.model.prepare_step(timeS-p.configuration.launch_epoch_s);
    if (candidate->step.region_config.mode==swcme::regions::Mode::FullICME) {
      // Include half the trailing smoothing width when proving support does
      // not overlap the ambient handoff; checking only the ejecta edge would
      // incorrectly replace part of its transition by unheated ambient fields.
      const auto b=swcme::regions::make_boundaries(candidate->step.r_sh_m,candidate->step.region_config);
      if (b.R_te_m-0.5*b.smooth_te_width_m<=swcme::solarwind::MIN_RADIUS_M)
        return Invalid("CME overlaps inner ambient handoff; supply a model covering that shell");
    }
    candidate->inner.reset(new AnalyticParkerProvider(p.ambient));
    const auto inner=candidate->inner->Prepare(timeS);if(!inner.ok())return inner;
    auto& m=candidate->metadata;
    // The application keeps fields frozen until the next background cadence.
    // Advertise that interval, clipped to event coverage; do not extrapolate
    // a prepared CME merely because the enclosing mesh still exists.
    m.provider=ProviderKind::Swcme;m.ownership=StorageOwnership::ModelOwned;
    m.epochS=m.validFromS=timeS;
    m.validUntilS=std::min(p.configuration.valid_until_s,timeS+p.cadence);
    m.generation=1+static_cast<std::uint64_t>(std::floor(tick+0.5L));
    if (p.prepared && m.generation<=p.prepared->metadata.generation &&
        timeS!=p.prepared->metadata.epochS)
      return Invalid("SWCME background epochs must advance monotonically");
    m.coordinateFrame=p.ambient.coordinateFrame;m.providerIdentity=CanonicalName();
    m.configurationFingerprint=SEP::Background::FingerprintConfiguration(ResolvedManifest());
    // The single swap follows every validation and model call above. Failure
    // leaves the previous epoch usable; prepared state never owns live mesh bytes.
    p.prepared=candidate;return Core::Status::OK();
  } catch(const std::exception& e) { return Invalid(std::string("SWCME mesh preparation: ")+e.what()); }
}
const SnapshotMetadata* SwcmeBackgroundProvider::PreparedMetadata() const {
  return implementation_->prepared?&implementation_->prepared->metadata:nullptr;
}
ProviderCapabilities SwcmeBackgroundProvider::Capabilities() const {
  ProviderCapabilities c;c.hasGradB=c.hasDivBhat=c.hasCurvature=c.hasGradU=true;
  c.hasFieldAlignedStrain=c.hasPlasmaState=c.supportsBatchEval=true;return c;
}
std::string SwcmeBackgroundProvider::ResolvedManifest() const {
  // Numerical differentiation, origin, inner handoff and species-heating
  // closure change transport just as model inputs do. Include them in identity
  // so fields with different policies cannot share an interpolation generation.
  const auto& p=*implementation_;
  std::ostringstream out;out.precision(17);
  out<<CanonicalName()<<';'<<p.configuration.normalized_manifest
     <<";inner_ambient="<<AnalyticParkerProvider(p.ambient).ResolvedManifest()
     <<";dt_s="<<p.dt<<";cadence_s="<<p.cadence
     <<";origin_m="<<p.origin.x<<','<<p.origin.y<<','<<p.origin.z
     <<";shock_solver=fast-branch-refinement-switch-on-v1"
     <<";gradient=cartesian-central2-boundary-one-sided2;relative_step=1e-4;layer_step_cap=width/16"
     <<";heating_partition=fixed-configured-temperature-ratios;inner_handoff_m="
     <<swcme::solarwind::MIN_RADIUS_M;return out.str();
}
BackgroundSample SwcmeBackgroundProvider::Evaluate(const Core::Vec3& x) const {
  BackgroundSample out;Core::Status status;
  EvaluateBatchDetailed(&x.x,&x.y,&x.z,1,&out,&status);
  if(!status.ok())out.status=status;
  return out;
}
Core::Status SwcmeBackgroundProvider::EvaluateBatchDetailed(
    const double* x,const double* y,const double* z,std::size_t n,
    BackgroundSample* out,Core::Status* statuses) const {
  const auto& p=*implementation_;
  if(!p.prepared) return Invalid("SWCME background must be prepared before evaluation");
  if(n && (!x||!y||!z||!out||!statuses))return Invalid("SWCME mesh batch pointer is null");
  Core::Status aggregate;
  // Bound scratch memory independently of mesh size. One authenticated
  // canonical batch evaluates the seven-point vector stencil for each chunk.
  // The epoch is frozen once, never prepared inside a cell/particle loop.
  constexpr std::size_t chunkSize=256;
  for(std::size_t begin=0;begin<n;begin+=chunkSize) {
    const std::size_t count=std::min(chunkSize,n-begin);
    std::vector<Core::Vec3> points;std::vector<double> steps(count);
    std::vector<Core::Status> inputStatus(count);
    std::vector<std::array<int,3>> oneSided(count);
    const auto boundaries=swcme::regions::make_boundaries(p.prepared->step.r_sh_m,p.prepared->step.region_config);
    for(std::size_t i=0;i<count;++i) {
      Core::Vec3 at(x[begin+i],y[begin+i],z[begin+i]);const double r=(at-p.origin).Norm();
      if(!std::isfinite(r)||r<p.ambient.sourceRadiusM) {
        inputStatus[i]=Invalid("SWCME mesh position outside physical ambient domain");
        // Do not let one malformed point poison the independent valid points.
        at=p.origin+Core::Vec3(2*Core::Const::AU,0,0);
      }
      double h=std::max(1.0,1e-4*(at-p.origin).Norm());
      // The 1-m floor avoids vanishing Cartesian steps near the inner shell;
      // the shock-width cap subsequently resolves the declared finite layer.
      // Seven points per cell are ordered: center, x pair, y pair, z pair.
      // oneSided=0 means (+h,-h); +/-1 means outward (sign*h,2*sign*h).
      if (boundaries.smooth_shock_width_m>0)
        h=std::min(h,boundaries.smooth_shock_width_m/16);
      steps[i]=h;points.push_back(at);
      for(int j=0;j<3;++j) {
        Core::Vec3 d;j==0?d.x=h:(j==1?d.y=h:d.z=h);
        // Do not sample the excluded Sun or shrink h to a roundoff-sized
        // value. At a photospheric boundary cell use a second-order outward
        // stencil; all other cells retain the symmetric Cartesian stencil.
        const auto local=at-p.origin;
        if ((local+d).Norm()<p.ambient.sourceRadiusM ||
            (local-d).Norm()<p.ambient.sourceRadiusM) {
          const double component=j==0?local.x:(j==1?local.y:local.z);
          const int sign=component<0?-1:1;oneSided[i][j]=sign;
          points.push_back(at+d*sign);points.push_back(at+d*(2*sign));
        } else {
          points.push_back(at+d);points.push_back(at-d);
        }
      }
    }
    std::vector<Primitive> values;
    const auto batch=p.Primitives(points,&values);
    for(std::size_t i=0;i<count;++i) {
      auto status=inputStatus[i].ok()?batch:inputStatus[i];
      std::vector<Primitive> recovered;
      const Primitive* cellValues=values.data()+7*i;
      if (inputStatus[i].ok() && !batch.ok()) {
        // A canonical batch can fail after a valid prefix. On this exceptional
        // path recover point/stencil status independently, so another cell's
        // model error cannot overwrite a good point's detailed output.
        const std::vector<Core::Vec3> local(points.begin()+7*i,points.begin()+7*(i+1));
        status=p.Primitives(local,&recovered);cellValues=recovered.data();
      }
      BackgroundSample sample;
      if(status.ok()) {
        const auto& v=cellValues[0];sample.B=v.B;sample.U=v.U;
        sample.numberDensityM3=v.n;sample.pressurePa=v.p;
        const auto ambient=swcme::solarwind::thermodynamic_state(p.prepared->step.common.solar_wind,v.n);
        // RH gives total pressure, not species heating fractions. Preserve the
        // configured Te/Tp and Ta/Tp ratios and scale all temperatures together.
        sample.temperatureK=p.prepared->step.common.solar_wind.T_K*v.p/ambient.pressure_Pa;
        sample.alfvenSpeedMpS=v.B.Norm()/std::sqrt(4e-7*Core::Const::kPi*v.rho);
        for(int j=0;j<3;++j) {
          const auto& left=cellValues[1+2*j];const auto& right=cellValues[2+2*j];
          const int sign=oneSided[i][j];
          // Central derivative: (f(+h)-f(-h))/(2h).
          // Boundary derivative: (-3f(0)+4f(sign*h)-f(2sign*h))/(2sign*h).
          // Differentiate components in Cartesian axes, not scalar |B| or
          // radial speed, so tangential fields and compression enter transport.
          const auto dB=sign?(-3*v.B+4*left.B-right.B)/(2*steps[i]*sign):
              (left.B-right.B)/(2*steps[i]);
          const auto dU=sign?(-3*v.U+4*left.U-right.U)/(2*steps[i]*sign):
              (left.U-right.U)/(2*steps[i]);
          sample.gradB(0,j)=dB.x;sample.gradB(1,j)=dB.y;sample.gradB(2,j)=dB.z;
          sample.gradU(0,j)=dU.x;sample.gradU(1,j)=dU.y;sample.gradU(2,j)=dU.z;
        }
        CompleteVectorDerivatives(&sample);sample.valid=true;
        sample.generation=p.prepared->metadata.generation;
        sample.configurationDigest=Digest(p.prepared->metadata.configurationFingerprint);
        status=ValidateCompleteSample(sample,Capabilities());
      }
      // Publish point results only after complete primitive/derivative checks.
      // Per-point statuses allow valid neighbours to survive a rejected point;
      // SnapshotBuilder still rejects the complete mesh candidate on any error.
      statuses[begin+i]=status;if(status.ok())out[begin+i]=sample;
      else if(aggregate.ok())aggregate=status;
    }
  }
  return aggregate;
}
} }
