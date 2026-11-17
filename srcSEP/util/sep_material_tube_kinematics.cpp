#include "sep_material_tube_kinematics.h"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <exception>
#include <limits>
#include <set>
#include <utility>
#ifdef __FAST_MATH__
#error "material signed-flux certificates require strict floating arithmetic, not fast-math reassociation"
#endif

namespace SEP { namespace FieldLine { namespace MaterialTube { namespace {
static_assert(std::numeric_limits<double>::is_iec559&&std::numeric_limits<double>::radix==2&&
    std::numeric_limits<double>::digits==53&&std::numeric_limits<double>::max_exponent==1024&&
    std::numeric_limits<double>::min_exponent==-1021&&std::numeric_limits<long double>::is_iec559&&
    std::numeric_limits<long double>::radix==2,"material flux certificate requires binary64 storage and IEEE binary wide arithmetic");
using Transport::Status;using Transport::StatusCode;
using WideVector=std::array<long double,3>;
using WideMatrix=std::array<long double,9>;
Status Error(const std::string& text){return Status::Error(StatusCode::InvalidArgument,text);}
Status Unresolved(const std::string& text){return Status::Error(StatusCode::UnresolvedCoefficient,text);}
Status Unsupported(const std::string& text){return Status::Error(StatusCode::UnsupportedConfiguration,text);}
bool Positive(double x){return std::isfinite(x)&&x>0;}
template<class Values>bool Finite(const Values& values){for(auto x:values)if(!std::isfinite(x))return false;return true;}
long double RoundoffBound(long double scale){return 512*std::numeric_limits<double>::epsilon()*scale;}
bool Roundoff(long double residual,long double scale){return std::isfinite(residual)&&std::fabs(residual)<=RoundoffBound(scale);}
template<class A,class B>long double Dot(const A& a,const B& b){long double v=0;for(int d=0;d<3;++d)v+=static_cast<long double>(a[d])*b[d];return v;}
struct Expansion {std::array<long double,16> terms{{0}};std::size_t size=0;};
bool NormalOrZero(long double x){return x==0||std::fpclassify(x)==FP_NORMAL;}
Status Grow(Expansion* expansion,long double term) {
  // General TwoSum: Q+e=next+error EXACTLY when its arithmetic does not
  // overflow/underflow. The product exponent preflight below guarantees a
  // common dyadic grid safely inside the normal range. Zero elimination
  // leaves an ordered nonoverlapping expansion, so an empty sum is exactzero.
  Expansion candidate;long double q=term;
  for(std::size_t i=0;i<expansion->size;++i) {
    const auto e=expansion->terms[i],next=q+e,virtualE=next-q;
    const auto error=(q-(next-virtualE))+(e-virtualE);
    if(!std::isfinite(next)||!NormalOrZero(next)||!NormalOrZero(error))
      return Unresolved("material flux expansion addition leaves its normal finite arithmetic domain");
    if(error!=0) {
      if(candidate.size==candidate.terms.size())return Unresolved("material flux expansion capacity exceeded");
      candidate.terms[candidate.size++]=error;
    }
    q=next;
  }
  if(q!=0) {
    if(candidate.size==candidate.terms.size())return Unresolved("material flux expansion capacity exceeded");
    candidate.terms[candidate.size++]=q;
  }
  *expansion=candidate;return Status::Ok();
}
template<class A,class B>Status ExactDot(const A& a,const B& b,Expansion* expansion) {
  Expansion result;
  for(int i=0;i<3;++i) {
    const long double x=a[i],y=b[i];if(x==0||y==0)continue;
    // FMA TwoProduct is NOT generally error-free at underflow. Explicitly
    // bound every possible exact product bit and reserve addition headroom.
    // Rejection is an arithmetic capability, never a physical field/flux floor.
    const int exponent=std::ilogb(x)+std::ilogb(y);
    const int lowestBit=exponent-2*(std::numeric_limits<long double>::digits-1);
    if(lowestBit<std::numeric_limits<long double>::min_exponent+4||
        exponent+1>std::numeric_limits<long double>::max_exponent-5)
      return Unresolved("material signed-flux products cannot be certified without underflow/overflow");
    const long double product=x*y,error=std::fma(x,y,-product);
    if(!std::isfinite(product)||!NormalOrZero(product)||!NormalOrZero(error))
      return Unresolved("material signed-flux FMA residual is outside its error-free arithmetic domain");
    auto status=Grow(&result,error);if(!status.ok())return status;
    status=Grow(&result,product);if(!status.ok())return status;
  }
  *expansion=result;return Status::Ok();
}
struct Enclosure {long double low=0,high=0;};
Enclosure Enclose(const Expansion& expansion) {
  Enclosure result;
  for(std::size_t i=0;i<expansion.size;++i) {
    // One-ULP outward bounds on each addition enclose the exact expansion
    // under the required round-to-nearest arithmetic. The empty
    // expansion remains exact[0,0]; no nonzero floor is invented for zero flux.
    result.low=std::nextafter(result.low+expansion.terms[i],-std::numeric_limits<long double>::infinity());
    result.high=std::nextafter(result.high+expansion.terms[i],std::numeric_limits<long double>::infinity());
  }
  return result;
}
long double Rounded(const Expansion& expansion) {
  long double result=0;for(std::size_t i=0;i<expansion.size;++i)result+=expansion.terms[i];return result;
}
Status FluxDifferenceBound(const Expansion& current,const Expansion& reference,long double* bound) {
  auto difference=current;
  for(std::size_t i=0;i<reference.size;++i){const auto status=Grow(&difference,-reference.terms[i]);if(!status.ok())return status;}
  const auto interval=Enclose(difference);*bound=std::max(std::fabs(interval.low),std::fabs(interval.high));
  return std::isfinite(*bound)?Status::Ok():Unresolved("material signed-flux difference enclosure is nonfinite");
}
Status CheckFluxRepresentability(const Vector& referenceB,const Vector& referenceArea,
    const WideVector& wideB,const WideVector& wideArea,const Vector& encodedB,const Vector& encodedArea,
    long double* referenceFlux,long double* encodedFlux) {
  Expansion ref,wide,encoded;auto status=ExactDot(referenceB,referenceArea,&ref);if(!status.ok())return status;
  status=ExactDot(wideB,wideArea,&wide);if(!status.ok())return status;
  status=ExactDot(encodedB,encodedArea,&encoded);if(!status.ok())return status;
  long double scale=0;
  if(ref.size) {
    const auto enclosure=Enclose(ref);
    if(enclosure.low<=0&&enclosure.high>=0)return Unresolved("nonzero reference magnetic flux is not numerically sign-resolved");
    // Use the LOWER magnitude bound, so the allowance cannot grow from a
    // cancellation-prone transformed product sum or enclosure roundoff.
    scale=std::min(std::fabs(enclosure.low),std::fabs(enclosure.high));
  } else {
    // Exactly zero flux has a dimensionally defined absolute contract based
    // on the REFERENCE field/area norms. It is a roundoff accuracy allowance,
    // not a fabricated nonzero flux or an undefined relative division byzero.
    Expansion bb,aa;status=ExactDot(referenceB,referenceB,&bb);if(!status.ok())return status;
    status=ExactDot(referenceArea,referenceArea,&aa);if(!status.ok())return status;
    const auto b=Enclose(bb),a=Enclose(aa);
    if(!(b.low>0)||!(a.low>0))return Unresolved("zero-flux reference norm scale is unresolved");
    const auto bn=std::nextafter(std::sqrt(b.low),0.0L),an=std::nextafter(std::sqrt(a.low),0.0L);
    scale=std::nextafter(bn*an,0.0L);
  }
  const auto allowed=std::nextafter(512*std::numeric_limits<double>::epsilon()*scale,0.0L);
  long double wideError=0,encodedError=0;
  status=FluxDifferenceBound(wide,ref,&wideError);if(!status.ok())return status;
  status=FluxDifferenceBound(encoded,ref,&encodedError);if(!status.ok())return status;
  if(!std::isfinite(allowed)||!(allowed>0)||wideError>allowed||encodedError>allowed)
    return Unresolved("material wide/encoded B-area magnetic flux loses its reference-scale numerical accuracy; no state repaired");
  *referenceFlux=Rounded(ref);*encodedFlux=Rounded(encoded);return Status::Ok();
}
Status Represent(long double x,double* value) {
  const double rounded=static_cast<double>(x);
  if(!std::isfinite(x)||!std::isfinite(rounded)||(x!=0&&rounded==0))
    return Unresolved("material geometry value underflowed/overflowed finite binary64 representation");
  *value=rounded;return Status::Ok();
}
template<class Out,class In>Status RepresentArray(const In& input,Out* output) {
  for(std::size_t i=0;i<input.size();++i){const auto s=Represent(input[i],&(*output)[i]);if(!s.ok())return s;}
  return Status::Ok();
}
struct Map {
  WideMatrix f,rate,cofactor,g;
  long double j=0,jdot=0;
};
Status Resolve(const Deformation& input,Map* map) {
  if(!Finite(input.gradient)||!Finite(input.ratePerS)||!Positive(input.declaredJacobian))
    return Error("material deformation needs finite F/Fdot and explicit positive declared J");
  Map m;for(int i=0;i<9;++i){m.f[i]=input.gradient[i];m.rate[i]=input.ratePerS[i];}
  const auto& f=m.f;auto& c=m.cofactor;
  // Cofactors are formed directly, avoiding an inverse followed by J*F^-T.
  // Invertibility is still required for G=Fdot*F^-1 and orientation-preserving
  // material maps. Cancellation uncertainty is a numerical failure, not a
  // determinant floor or a repaired deformation.
  c={{f[4]*f[8]-f[5]*f[7],f[5]*f[6]-f[3]*f[8],f[3]*f[7]-f[4]*f[6],
      f[2]*f[7]-f[1]*f[8],f[0]*f[8]-f[2]*f[6],f[1]*f[6]-f[0]*f[7],
      f[1]*f[5]-f[2]*f[4],f[2]*f[3]-f[0]*f[5],f[0]*f[4]-f[1]*f[3]}};
  m.j=f[0]*c[0]+f[1]*c[1]+f[2]*c[2];
  const long double determinantScale=std::fabs(f[0]*f[4]*f[8])+std::fabs(f[0]*f[5]*f[7])+
      std::fabs(f[1]*f[3]*f[8])+std::fabs(f[1]*f[5]*f[6])+std::fabs(f[2]*f[3]*f[7])+std::fabs(f[2]*f[4]*f[6]);
  if(!std::isfinite(m.j)||m.j<=0)return Error("material deformation is singular or reverses orientation");
  if(m.j<=64*std::numeric_limits<long double>::epsilon()*determinantScale)
    return Unresolved("material determinant sign/magnitude is unresolved by its arithmetic cancellation bound");
  if(!Roundoff(m.j-input.declaredJacobian,std::fabs(m.j)+input.declaredJacobian))
    return Error("declared material Jacobian differs from det(F)");
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    long double identity=0,scale=0;
    for(int k=0;k<3;++k){identity+=f[3*i+k]*c[3*j+k];scale+=std::fabs(f[3*i+k]*c[3*j+k]);}
    const auto expected=i==j?m.j:0;
    if(!Roundoff(identity-expected,scale+std::fabs(expected)))return Unresolved("material cofactor identity is numerically unresolved");
    m.g[3*i+j]=0;for(int k=0;k<3;++k)m.g[3*i+j]+=m.rate[3*i+k]*c[3*j+k]/m.j;
    m.jdot+=c[3*i+j]*m.rate[3*i+j];
  }
  if(!Finite(m.g)||!std::isfinite(m.jdot))return Unresolved("material velocity gradient/Jacobian rate is nonfinite");
  *map=m;return Status::Ok();
}
WideVector Multiply(const WideMatrix& m,const Vector& a) {
  WideVector result{{0,0,0}};for(int i=0;i<3;++i)for(int j=0;j<3;++j)result[i]+=m[3*i+j]*a[j];return result;
}
Status SampleValid(const AffineSample& s,Map* map) {
  if(!std::isfinite(s.timeS)||!Finite(s.translationM)||!Finite(s.translationVelocityMPerS))
    return Error("affine material sample lacks finite time/translation/velocity");
  return Resolve(s.deformation,map);
}
Status ReferenceValid(const ReferenceCell& ref) {
  if(!Positive(ref.volumeM3)||!Positive(ref.materialMassKg)||ref.identity.empty()||ref.provenance.empty()||ref.faces.size()<3)
    return Error("material reference requires true positive volume/mass, complete surfaces and provenance");
  std::set<std::uint64_t> ids;WideVector area{{0,0,0}},areaScale{{0,0,0}};
  WideMatrix moment{{0,0,0,0,0,0,0,0,0}},momentScale=moment;
  long double flux=0,fluxScale=0;bool lower=false,upper=false,lateral=false;
  for(const auto& face:ref.faces) {
    if(!ids.insert(face.id).second||!Finite(face.orientedAreaM2)||!Finite(face.firstMomentM3)||!std::isfinite(face.magneticFluxWb))
      return Error("material reference face namespace/moments/flux are invalid");
    if(face.role==FaceRole::LowerCap)lower=true;
    else if(face.role==FaceRole::UpperCap)upper=true;
    else if(face.role==FaceRole::Lateral){lateral=true;if(face.magneticFluxWb!=0)return Error("material magnetic side wall carries declared nonzero normal magnetic flux");}
    else return Unsupported("unknown material surface role");
    if(face.role!=FaceRole::Lateral&&face.magneticFluxWb==0)return Error("material tube cap must carry declared nonzero magnetic flux");
    for(int d=0;d<3;++d){area[d]+=face.orientedAreaM2[d];areaScale[d]+=std::fabs(face.orientedAreaM2[d]);}
    for(int d=0;d<9;++d){moment[d]+=face.firstMomentM3[d];momentScale[d]+=std::fabs(face.firstMomentM3[d]);}
    flux+=face.magneticFluxWb;fluxScale+=std::fabs(face.magneticFluxWb);
  }
  if(!lower||!upper||!lateral)return Error("material reference omits a cap or the lateral surface");
  for(int d=0;d<3;++d)if(!Roundoff(area[d],areaScale[d]))return Error("material reference oriented surfaces are not closed");
  // Divergence theorem, independently supplied area/moment/volume data:
  // sum_faces integral a_i*n_j dA=delta_ij*V. Do not define V from this sum.
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    const auto expected=i==j?ref.volumeM3:0;
    if(!Roundoff(moment[3*i+j]-expected,momentScale[3*i+j]+std::fabs(expected)))
      return Error("material reference face first moments disagree with independent volume");
  }
  return Roundoff(flux,fluxScale)?Status::Ok():Error("material reference signed cap magnetic flux is not closed");
}
Status MapFace(const ReferenceFace& ref,const AffineSample& sample,const Map& map,Vector* area,Matrix* moment) {
  const auto a=Multiply(map.cofactor,ref.orientedAreaM2);auto s=RepresentArray(a,area);if(!s.ok())return s;
  WideMatrix m{{0,0,0,0,0,0,0,0,0}};
  // Integral x_i*n_j dA = F_ik*M_kl*cof(F)_jl + c_i*A_current,j.
  // This is a surface-moment map; it does not replace a curved surface by a
  // plane, nor assume that the cap normal and local magnetic field align.
  for(int i=0;i<3;++i)for(int j=0;j<3;++j) {
    m[3*i+j]=sample.translationM[i]*a[j];
    for(int k=0;k<3;++k)for(int l=0;l<3;++l)m[3*i+j]+=map.f[3*i+k]*ref.firstMomentM3[3*k+l]*map.cofactor[3*j+l];
  }
  return RepresentArray(m,moment);
}
long double SurfaceVolumeRate(const ReferenceFace& face,const AffineSample& sample,const Map& map) {
  long double rate=0;
  for(int i=0;i<3;++i)for(int k=0;k<3;++k) {
    rate+=sample.translationVelocityMPerS[i]*map.cofactor[3*i+k]*face.orientedAreaM2[k];
    for(int j=0;j<3;++j)rate+=map.rate[3*i+j]*map.cofactor[3*i+k]*face.firstMomentM3[3*j+k];
  }
  return rate;
}
} // namespace

PointResult PushForwardPoint(const PointInput& input) {
  PointResult result;
  // TwoSum's proof assumes nearest rounding; this routine never changes the
  // process/thread rounding mode to manufacture a certificate. Fast-math is
  // also rejected at compilation above, before its algebra can erase errors.
  if(std::fegetround()!=FE_TONEAREST) {
    result.status=Unresolved("material signed-flux certificate requires round-to-nearest arithmetic");return result;
  }
  Map map;result.status=Resolve(input.deformation,&map);if(!result.status.ok())return result;
  if(!Positive(input.referenceMassDensityKgM3)||!Finite(input.referenceMagneticFieldT)||!Finite(input.referenceAreaVectorM2)||
      !(Dot(input.referenceMagneticFieldT,input.referenceMagneticFieldT)>0)||!(Dot(input.referenceAreaVectorM2,input.referenceAreaVectorM2)>0)) {
    result.status=Error("material point needs explicit positive rho/nonzero finite B/area");return result;
  }
  auto b=Multiply(map.f,input.referenceMagneticFieldT),rate=Multiply(map.rate,input.referenceMagneticFieldT);
  const auto a=Multiply(map.cofactor,input.referenceAreaVectorM2);
  for(int d=0;d<3;++d){b[d]/=map.j;rate[d]=rate[d]/map.j-map.jdot/map.j*b[d];}
  const long double rho=input.referenceMassDensityKgM3/map.j,rhoRate=-rho*map.jdot/map.j;
  // Narrow into private scratch BEFORE certifying the physical invariant.
  // Valid wide values do not guarantee that independently encoded B and area
  // retain a tiny flux after huge cancellation. Failed checks publish no state.
  PointState state;
  const std::pair<long double,double*> scalarFields[]={{map.j,&state.jacobian},{map.jdot,&state.jacobianRatePerS},
      {rho,&state.massDensityKgM3},{rhoRate,&state.massDensityRateKgM3PerS}};
  for(const auto& item:scalarFields) {
    result.status=Represent(item.first,item.second);if(!result.status.ok())return result;
  }
  result.status=RepresentArray(map.cofactor,&state.cofactor);if(!result.status.ok())return result;
  result.status=RepresentArray(map.g,&state.velocityGradientPerS);if(!result.status.ok())return result;
  result.status=RepresentArray(b,&state.magneticFieldT);if(!result.status.ok())return result;
  result.status=RepresentArray(rate,&state.magneticFieldRateTPerS);if(!result.status.ok())return result;
  result.status=RepresentArray(a,&state.areaVectorM2);if(!result.status.ok())return result;
  long double referenceFlux=0,encodedFlux=0;
  result.status=CheckFluxRepresentability(input.referenceMagneticFieldT,input.referenceAreaVectorM2,b,a,
      state.magneticFieldT,state.areaVectorM2,&referenceFlux,&encodedFlux);if(!result.status.ok())return result;
  result.status=Represent(referenceFlux,&state.referenceFluxWb);if(!result.status.ok())return result;
  result.status=Represent(encodedFlux,&state.currentFluxWb);if(!result.status.ok())return result;
  result.state=state;return result;
}

GeometryResult PrepareAffineInterval(const ReferenceCell& ref,const IntervalInput& interval) {
  GeometryResult result;
  if(interval.capability!=MapCapability::GlobalAffine) {
    result.status=Unsupported("material cell kernel requires explicit global-affine map; nonaffine native interpolation is unavailable");return result;
  }
  result.status=ReferenceValid(ref);if(!result.status.ok())return result;
  if(!interval.budget.declared||interval.budget.provenance.empty()||interval.operationProvenance.empty()||interval.quadrature.empty()) {
    result.status=Unsupported("material interval requires explicit time quadrature/error budget and operation provenance");return result;
  }
  for(double value:{interval.budget.volumeM3,interval.budget.deformation,interval.budget.translationM})
    if(!std::isfinite(value)||value<0){result.status=Error("material numerical quadrature budgets must be finite nonnegative absolute allowances");return result;}
  try {
    Map first,last;result.status=SampleValid(interval.initial,&first);if(!result.status.ok())return result;
    result.status=SampleValid(interval.final,&last);if(!result.status.ok())return result;
    const long double duration=static_cast<long double>(interval.final.timeS)-interval.initial.timeS;
    if(!(duration>0)){result.status=Error("material interval has nonpositive duration");return result;}
    auto candidate=std::make_shared<GeometryCandidate>();candidate->fieldLine=ref.fieldLine;candidate->segment=ref.segment;
    candidate->initialTimeS=interval.initial.timeS;candidate->finalTimeS=interval.final.timeS;
    candidate->materialMassKg=ref.materialMassKg;candidate->referenceIdentity=ref.identity;candidate->referenceProvenance=ref.provenance;
    candidate->operationProvenance=interval.operationProvenance;candidate->quadratureProvenance=interval.budget.provenance;
    const long double v0=first.j*ref.volumeM3,v1=last.j*ref.volumeM3;
    const std::pair<long double,double*> physicalFields[]={{v0,&candidate->initialVolumeM3},{v1,&candidate->finalVolumeM3},
        {ref.materialMassKg/v0,&candidate->initialMeanDensityKgM3},{ref.materialMassKg/v1,&candidate->finalMeanDensityKgM3}};
    for(const auto& item:physicalFields) {
      result.status=Represent(item.first,item.second);if(!result.status.ok())return result;
    }
    std::vector<long double> swept(ref.faces.size(),0);WideMatrix integratedRate{{0,0,0,0,0,0,0,0,0}},rateScale=integratedRate;
    WideVector integratedVelocity{{0,0,0}},velocityScale=integratedVelocity;long double weights=0,previous=interval.initial.timeS;
    bool firstPoint=true;
    for(const auto& point:interval.quadrature) {
      if(!Positive(point.weightS)||point.sample.timeS<interval.initial.timeS||point.sample.timeS>interval.final.timeS||
          (!firstPoint&&point.sample.timeS<=previous)) {result.status=Error("material time quadrature nodes/weights are invalid or unordered");return result;}
      firstPoint=false;previous=point.sample.timeS;Map map;result.status=SampleValid(point.sample,&map);if(!result.status.ok())return result;
      weights+=point.weightS;long double rateSum=0,rateAbs=0;
      for(std::size_t face=0;face<ref.faces.size();++face) {
        const auto rate=SurfaceVolumeRate(ref.faces[face],point.sample,map);rateSum+=rate;rateAbs+=std::fabs(rate);
        swept[face]+=point.weightS*rate;
      }
      const auto expectedRate=map.jdot*ref.volumeM3;
      if(!Roundoff(rateSum-expectedRate,rateAbs+std::fabs(expectedRate))) {
        result.status=Unresolved("material cap+side instantaneous GCL disagrees with independent Jacobian rate");return result;
      }
      for(int i=0;i<9;++i){const auto change=point.weightS*map.rate[i];integratedRate[i]+=change;rateScale[i]+=std::fabs(change);}
      for(int i=0;i<3;++i){const auto change=point.weightS*point.sample.translationVelocityMPerS[i];integratedVelocity[i]+=change;velocityScale[i]+=std::fabs(change);}
    }
    if(!Roundoff(weights-duration,weights+duration)){result.status=Error("material time weights do not span the declared interval");return result;}
    for(int i=0;i<9;++i) {
      const auto delta=last.f[i]-first.f[i],residual=integratedRate[i]-delta;
      candidate->maximumDeformationResidual=std::max(candidate->maximumDeformationResidual,static_cast<double>(std::fabs(residual)));
      if(std::fabs(residual)>interval.budget.deformation+RoundoffBound(rateScale[i]+std::fabs(last.f[i])+std::fabs(first.f[i]))) {
        result.status=Unresolved("material Fdot time quadrature exceeds its declared endpoint-deformation error budget");return result;
      }
    }
    for(int i=0;i<3;++i) {
      const auto delta=static_cast<long double>(interval.final.translationM[i])-interval.initial.translationM[i],residual=integratedVelocity[i]-delta;
      candidate->maximumTranslationResidualM=std::max(candidate->maximumTranslationResidualM,static_cast<double>(std::fabs(residual)));
      if(std::fabs(residual)>interval.budget.translationM+RoundoffBound(velocityScale[i]+std::fabs(interval.final.translationM[i])+std::fabs(interval.initial.translationM[i]))) {
        result.status=Unresolved("material translation time quadrature exceeds its declared endpoint error budget");return result;
      }
    }
    long double cap=0,side=0,sweptAbs=0;
    for(std::size_t i=0;i<ref.faces.size();++i) {
      const auto& face=ref.faces[i];FaceGeometry image;image.id=face.id;image.role=face.role;image.materialMagneticFluxWb=face.magneticFluxWb;
      result.status=MapFace(face,interval.initial,first,&image.initialAreaM2,&image.initialFirstMomentM3);if(!result.status.ok())return result;
      result.status=MapFace(face,interval.final,last,&image.finalAreaM2,&image.finalFirstMomentM3);if(!result.status.ok())return result;
      result.status=Represent(swept[i],&image.outwardGridSweptVolumeM3);if(!result.status.ok())return result;
      if(face.role==FaceRole::Lateral)side+=swept[i];else cap+=swept[i];sweptAbs+=std::fabs(swept[i]);candidate->faces.push_back(image);
    }
    const auto residual=v1-v0-cap-side,allowed=interval.budget.volumeM3+RoundoffBound(std::fabs(v0)+std::fabs(v1)+sweptAbs);
    if(!std::isfinite(residual)||!std::isfinite(allowed)||std::fabs(residual)>allowed) {
      result.status=Unresolved("material endpoint volume and independently integrated cap+side sweeps exceed declared GCL budget");return result;
    }
    const std::pair<long double,double*> ledgerFields[]={{cap,&candidate->capSweptVolumeM3},{side,&candidate->lateralSweptVolumeM3},
        {residual,&candidate->gclResidualM3},{allowed,&candidate->gclAllowedResidualM3}};
    for(const auto& item:ledgerFields) {
      result.status=Represent(item.first,item.second);if(!result.status.ok())return result;
    }
    candidate->withinDeclaredQuadratureBudget=true;
    result.candidate=std::move(candidate);result.status=Status::Ok();return result;
  } catch(const std::exception& error) {
    result.status=Error(std::string("material geometry preparation failed before publication: ")+error.what());return result;
  }
}
} } }
