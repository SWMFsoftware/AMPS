#include "sep_material_tube_kinematics_validation.h"
#include "sep_material_tube_kinematics.h"
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace M=FieldLine::MaterialTube;
const double B0=5e-9;
// The curved reference uses B_phi=C/R: C has units T*m, unlike the uniform
// point-test B0[T]. Naming both prevents an implicit dimensionless radius.
const double ToroidalFieldRadiusProductTM=5e-9;
M::Matrix Identity(){return {{1,0,0,0,1,0,0,0,1}};}
bool Near(double actual,double expected,double relative=2e-12) {
  return std::isfinite(actual)&&std::isfinite(expected)&&std::fabs(static_cast<long double>(actual)-expected)<=
      relative*std::max(std::fabs(static_cast<long double>(actual)),std::fabs(static_cast<long double>(expected)));
}
bool Small(double value,double scale){return std::isfinite(value)&&std::fabs(value)<=2e-12*scale;}
M::Vector Cross(const M::Vector& a,const M::Vector& b){return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};}
M::Vector Multiply(const M::Matrix& a,const M::Vector& v) {
  M::Vector result{{0,0,0}};for(int i=0;i<3;++i)for(int j=0;j<3;++j)result[i]+=a[3*i+j]*v[j];return result;
}
Result Finish(bool pass,const std::string& text) {
  Result result;result.status=pass?Status::Pass:Status::Fail;result.message=text;
  result.configuration.push_back("capability=neutral-material-kinematics-only;native_publication=false;production_qualification=false");
  result.configuration.push_back("point=Cauchy-B/rho-and-Nanson-area;cell=explicit-global-affine-map;all-cap+side-GCL");
  result.configuration.push_back("SI-inputs;true-material-mass;explicit-time-quadrature-and-numerical-budget;no-volume-repair");return result;
}
M::ReferenceCell Box() {
  // Exact rectangular material reference: lengths(2,3,4)m, V24m3.
  // Integrated mass is explicitly specified independently of any point rho.
  M::ReferenceCell r;r.volumeM3=24;r.materialMassKg=7e-19;r.identity="exact-box-2x3x4";
  r.provenance="analytic oriented face areas and centroid first moments of a centered rectangular box";
  const double area[3]={12,8,6};
  for(int d=0;d<3;++d)for(int side=0;side<2;++side) {
    M::ReferenceFace face;face.id=2*d+side;face.role=d==0?(side?M::FaceRole::UpperCap:M::FaceRole::LowerCap):M::FaceRole::Lateral;
    face.orientedAreaM2[d]=(side?1:-1)*area[d];face.firstMomentM3[3*d+d]=12;
    face.magneticFluxWb=d==0?(side?1:-1)*area[d]*B0:0;r.faces.push_back(face);
  }
  return r;
}
M::AffineSample LinearSample(const M::Matrix& rate,const M::Vector& velocity,double time) {
  // Test maps are upper triangular, so the independently declared determinant
  // is the diagonal product; production's general cofactor routine is unused.
  M::AffineSample sample;sample.timeS=time;sample.deformation.gradient=Identity();sample.deformation.ratePerS=rate;
  for(int i=0;i<9;++i)sample.deformation.gradient[i]+=time*rate[i];
  sample.deformation.declaredJacobian=(1+time*rate[0])*(1+time*rate[4])*(1+time*rate[8]);
  sample.translationVelocityMPerS=velocity;for(int d=0;d<3;++d)sample.translationM[d]=time*velocity[d];return sample;
}
M::IntervalInput LinearInterval(const M::Matrix& rate,double duration,const M::Vector& velocity={{0,0,0}}) {
  M::IntervalInput in;in.capability=M::MapCapability::GlobalAffine;
  in.initial=LinearSample(rate,velocity,0);in.final=LinearSample(rate,velocity,duration);
  // Two-node Gauss exactly integrates quadratic cofactor*constant Fdot for a
  // linear-in-time affine map. This is an analytic numerical-control choice,
  // not a production timestep or domain limit.
  for(int sign:{-1,1}) {
    M::TimeQuadraturePoint point;point.weightS=duration/2;
    point.sample=LinearSample(rate,velocity,duration/2*(1+sign/std::sqrt(3.0)));in.quadrature.push_back(point);
  }
  in.budget.declared=true;in.budget.provenance="analytic degree-two time integrands;2-point Gauss;zero allowance beyond arithmetic roundoff";
  in.operationProvenance="explicit manufactured global-affine material kinematics; no autonomous momentum/thermal solve";return in;
}
Result PointwiseCauchy() {
  M::PointInput in;in.deformation.gradient={{2,.5,-.2,0,3,.4,0,0,4}};in.deformation.declaredJacobian=24;
  in.referenceMagneticFieldT={{B0,2*B0,-B0}};in.referenceMassDensityKgM3=8e-21;
  const M::Vector tangent0{{1,2,0}},tangent1{{0,1,3}};in.referenceAreaVectorM2=Cross(tangent0,tangent1);
  const M::Matrix g{{.1,.2,0,0,-.05,.03,0,0,.04}};
  for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)
    in.deformation.ratePerS[3*i+j]+=g[3*i+k]*in.deformation.gradient[3*k+j];
  const auto result=M::PushForwardPoint(in);if(!result.status.ok())return Finish(false,result.status.message);
  const auto& s=result.state;const M::Vector expectedB{{3.2*B0/24,5.6*B0/24,-4*B0/24}};
  // Area oracle uses actual transformed spanning tangents and their cross
  // product, independently of the implementation's cofactor determinant.
  const auto expectedArea=Cross(Multiply(in.deformation.gradient,tangent0),Multiply(in.deformation.gradient,tangent1));
  bool pass=Near(s.jacobian,24)&&Near(s.massDensityKgM3,8e-21/24)&&Near(s.jacobianRatePerS,24*.09);
  for(int d=0;d<3;++d) {
    pass=pass&&Near(s.magneticFieldT[d],expectedB[d])&&Near(s.areaVectorM2[d],expectedArea[d]);
    double expectedRate=-.09*expectedB[d];for(int k=0;k<3;++k)expectedRate+=g[3*d+k]*expectedB[k];
    pass=pass&&Near(s.magneticFieldRateTPerS[d],expectedRate);
  }
  for(int i=0;i<9;++i)pass=pass&&(g[i]==0?Small(s.velocityGradientPerS[i],.2):Near(s.velocityGradientPerS[i],g[i]));
  pass=pass&&Near(s.currentFluxWb,s.referenceFluxWb)&&Near(s.massDensityRateKgM3PerS,-.09*8e-21/24);
  // Exact90-degree rigid rotation tests inverse-transpose orientation as well
  // as stretch/shear. B and area both rotate, with no density/J change.
  in.deformation.gradient={{0,-1,0,1,0,0,0,0,1}};in.deformation.ratePerS={{-.3,0,0,0,-.3,0,0,0,0}};
  in.deformation.declaredJacobian=1;in.referenceMagneticFieldT={{B0,0,0}};in.referenceAreaVectorM2={{2,0,0}};
  const auto rotation=M::PushForwardPoint(in);pass=pass&&rotation.status.ok();
  if(rotation.status.ok())pass=pass&&Near(rotation.state.magneticFieldT[1],B0)&&Near(rotation.state.areaVectorM2[1],2)&&
      Near(rotation.state.magneticFieldRateTPerS[0],-.3*B0)&&rotation.state.jacobianRatePerS==0&&Near(rotation.state.massDensityKgM3,8e-21);
  return Finish(pass,"general3D stretch/shear and rigid rotation satisfy independent tangent-area, Cauchy induction, density and full velocity-gradient/rate references");
}
Result IsotropicDilation() {
  const auto ref=Box();const double a=.1,t=2,lambda=1+a*t,j=lambda*lambda*lambda;
  const auto result=M::PrepareAffineInterval(ref,LinearInterval({{a,0,0,0,a,0,0,0,a}},t));
  if(!result.status.ok())return Finish(false,result.status.message);
  const auto& c=*result.candidate;bool pass=!c.qualifiedForProduction&&c.withinDeclaredQuadratureBudget&&
      Near(c.initialVolumeM3,24)&&Near(c.finalVolumeM3,24*j)&&c.materialMassKg==ref.materialMassKg&&
      Near(c.finalMeanDensityKgM3,ref.materialMassKg/(24*j))&&Near(c.capSweptVolumeM3,24*(j-1)/3)&&
      Near(c.lateralSweptVolumeM3,2*24*(j-1)/3);
  for(const auto& face:c.faces)pass=pass&&Near(face.outwardGridSweptVolumeM3,4*(j-1));
  auto r=Finish(pass,"exact isotropic material dilation retains true integrated mass, scales density byJ^-1 and closes endpoint volume against independently integrated cap AND side sweeps");
  r.metrics.push_back({"absolute_GCL_residual",std::fabs(c.gclResidualM3),c.gclAllowedResidualM3,"<=","m3"});return r;
}
Result ShearAndTranslation() {
  const auto ref=Box();const double t=1.3,shear=.2;const M::Vector velocity{{2,-1,.5}};
  const auto result=M::PrepareAffineInterval(ref,LinearInterval({{0,shear,0,0,0,0,0,0,0}},t,velocity));
  if(!result.status.ok())return Finish(false,result.status.message);
  const auto& c=*result.candidate;bool pass=Near(c.finalVolumeM3,24)&&c.materialMassKg==ref.materialMassKg;
  const double xSweep=12*(2*t+shear*t*t/2),ySweep=-8*t,zSweep=3*t;
  const double expected[6]={-xSweep,xSweep,-ySweep,ySweep,-zSweep,zSweep};
  for(std::size_t f=0;f<6;++f)pass=pass&&Near(c.faces[f].outwardGridSweptVolumeM3,expected[f]);
  pass=pass&&Near(c.faces[1].finalAreaM2[0],12)&&Near(c.faces[1].finalAreaM2[1],-12*shear*t)&&
      Small(c.capSweptVolumeM3,24)&&Small(c.lateralSweptVolumeM3,24);
  // Material cap normals tilt while a reference B along x remains along x.
  // A source bridge cannot infer an aligned-cap kinetic domain from this map.
  M::PointInput point;point.deformation=LinearInterval({{0,shear,0,0,0,0,0,0,0}},t).final.deformation;
  point.referenceMagneticFieldT={{B0,0,0}};point.referenceAreaVectorM2={{12,0,0}};point.referenceMassDensityKgM3=ref.materialMassKg/ref.volumeM3;
  const auto pushed=M::PushForwardPoint(point);pass=pass&&pushed.status.ok();
  if(pushed.status.ok())pass=pass&&pushed.state.magneticFieldT[1]==0&&pushed.state.areaVectorM2[1]!=0&&Near(pushed.state.currentFluxWb,12*B0);
  return Finish(pass,"isochoric shear plus full translation reproduces analytic signed face sweeps and tilted material cap normals while preserving magnetic flux and true mass");
}
M::ReferenceCell CurvedSector() {
  // Exact annular sector extrusion: R in[2,3]m, phi in[-.4,.4],
  // z in[-.75,.75]m. This has a curved centerline and two curved radial sides.
  // Areas/moments below are analytic surface integrals, not planar fitting.
  const double a=2,b=3,h=1.5,alpha=.4,base=(b*b-a*a)*alpha,span=(b-a)*h,mean=(a+b)/2;
  M::ReferenceCell r;r.volumeM3=base*h;
  // rho_ref(R)=(1+R/(1m))*1e-21kg/m3. Integrate with the actual R*dR*dphi*dz
  // volume measure; midpoint rho*V is not substituted for true material mass.
  r.materialMassKg=2*alpha*h*((b*b-a*a)/2+(b*b*b-a*a*a)/3)*1e-21;
  r.identity="exact-annular-sector-extrusion";
  r.provenance="analytic cylindrical surface integrals; curved reference toroidal magnetic tube Bphi=constant/R";
  for(int side=0;side<2;++side) {
    const double sign=side?1:-1,phi=sign*alpha;M::ReferenceFace f;f.id=side;
    f.role=side?M::FaceRole::UpperCap:M::FaceRole::LowerCap;
    f.orientedAreaM2={{-sign*std::sin(phi)*span,sign*std::cos(phi)*span,0}};
    const M::Vector centroid{{mean*std::cos(phi),mean*std::sin(phi),0}};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)f.firstMomentM3[3*i+j]=centroid[i]*f.orientedAreaM2[j];
    f.magneticFluxWb=sign*ToroidalFieldRadiusProductTM*h*std::log(b/a);r.faces.push_back(f);
  }
  for(int side=0;side<2;++side) {
    const double radius=side?b:a,sign=side?1:-1;M::ReferenceFace f;f.id=2+side;
    f.orientedAreaM2={{sign*2*radius*h*std::sin(alpha),0,0}};
    f.firstMomentM3[0]=sign*radius*radius*h*(alpha+std::sin(2*alpha)/2);
    f.firstMomentM3[4]=sign*radius*radius*h*(alpha-std::sin(2*alpha)/2);r.faces.push_back(f);
  }
  for(int side=0;side<2;++side) {
    const double sign=side?1:-1;M::ReferenceFace f;f.id=4+side;f.orientedAreaM2[2]=sign*base;
    f.firstMomentM3[2]=sign*2*std::sin(alpha)*(b*b*b-a*a*a)/3;f.firstMomentM3[8]=h*base/2;r.faces.push_back(f);
  }
  return r;
}
Result CurvedReference() {
  const auto ref=CurvedSector();const double ax=.2,ay=-.1,az=.05,s=.3,t=.8;
  const auto interval=LinearInterval({{ax,s,0,0,ay,0,0,0,az}},t);const auto result=M::PrepareAffineInterval(ref,interval);
  if(!result.status.ok())return Finish(false,result.status.message);
  const auto& c=*result.candidate;const double lx=1+ax*t,ly=1+ay*t,lz=1+az*t,j=lx*ly*lz;
  auto integral=[&](double a,double b){return t+(a+b)*t*t/2+a*b*t*t*t/3;};
  const double dx=ax*integral(ay,az),dy=ay*integral(ax,az),dz=az*integral(ax,ay),cross=s*(t+az*t*t/2);
  const double alpha=.4,h=1.5,mean=2.5,m=mean*h*std::sin(alpha)*std::cos(alpha),q=mean*h*std::sin(alpha)*std::sin(alpha);
  // Independently expanded cylindrical surface integrals and elementary time
  // antiderivatives. Production uses general tensors/moments/Gauss nodes.
  const double expected[6]={-dx*m+dy*m+cross*q,-dx*m+dy*m-cross*q,
      -4*h*(dx*(alpha+std::sin(2*alpha)/2)+dy*(alpha-std::sin(2*alpha)/2)),
       9*h*(dx*(alpha+std::sin(2*alpha)/2)+dy*(alpha-std::sin(2*alpha)/2)),
       dz*ref.volumeM3/2,dz*ref.volumeM3/2};
  bool pass=Near(c.finalVolumeM3,ref.volumeM3*j)&&Near(ref.materialMassKg,1.06e-20)&&
      Near(c.materialMassKg,1.06e-20)&&Near(c.finalMeanDensityKgM3,1.06e-20/(3*j));
  for(std::size_t f=0;f<6;++f)pass=pass&&Near(c.faces[f].outwardGridSweptVolumeM3,expected[f]);
  // Nonuniform toroidal Bphi=C/R is genuinely tangent to both curved radial
  // sides and z sides. Exact expanded mapped components check local flux,
  // without replacing it by a midpoint field times integrated cap area.
  for(double phi:{-.4,.4})for(double radius:{2.0,2.5,3.0}) {
    M::PointInput in;in.deformation=interval.final.deformation;in.referenceMassDensityKgM3=(1+radius)*1e-21;
    in.referenceMagneticFieldT={{-ToroidalFieldRadiusProductTM*std::sin(phi)/radius,ToroidalFieldRadiusProductTM*std::cos(phi)/radius,0}};
    in.referenceAreaVectorM2={{-.7*std::sin(phi),.7*std::cos(phi),0}};
    const auto point=M::PushForwardPoint(in);pass=pass&&point.status.ok();
    if(point.status.ok())pass=pass&&Near(point.state.magneticFieldT[0],ToroidalFieldRadiusProductTM*(-lx*std::sin(phi)+s*t*std::cos(phi))/(radius*j))&&
        Near(point.state.magneticFieldT[1],ToroidalFieldRadiusProductTM*ly*std::cos(phi)/(radius*j))&&
        Near(point.state.currentFluxWb,.7*ToroidalFieldRadiusProductTM/radius);
  }
  auto r=Finish(pass,"curved toroidal reference surfaces use exact cylindrical area/moment integrals; general affine stretch/shear matches independent cap/side antiderivatives and nonuniform-B Cauchy flux");
  r.metrics.push_back({"absolute_curved_GCL_residual",std::fabs(c.gclResidualM3),c.gclAllowedResidualM3,"<=","m3"});return r;
}
Result NonAffinePointLaw() {
  bool pass=true;const double q=.15,qdot=.04;
  for(double x:{-2.0,-.3,0.0,.7,2.0}) {
    // X=(a_x,a_y+q*a_x^2,a_z) has exact J1 and curved material field
    // lines y=q*x^2. This tests LOCAL deformation data only; its spatially
    // varying F cannot enter the affine-cell moment integration API.
    M::PointInput in;in.deformation.gradient=Identity();in.deformation.gradient[3]=2*q*x;
    in.deformation.ratePerS[3]=2*qdot*x;in.deformation.declaredJacobian=1;
    in.referenceMagneticFieldT={{B0,0,0}};in.referenceAreaVectorM2={{1,0,0}};in.referenceMassDensityKgM3=9e-21;
    const auto p=M::PushForwardPoint(in);pass=pass&&p.status.ok();
    if(p.status.ok())pass=pass&&Near(p.state.magneticFieldT[0],B0)&&Small(p.state.magneticFieldT[1]-2*q*x*B0,B0)&&
        Small(p.state.magneticFieldRateTPerS[1]-2*qdot*x*B0,B0)&&p.state.massDensityKgM3==9e-21&&Near(p.state.currentFluxWb,B0);
  }
  auto interval=LinearInterval({{0,0,0,0,0,0,0,0,0}},1);interval.capability=M::MapCapability::NonAffineUnavailable;
  const auto rejected=M::PrepareAffineInterval(Box(),interval);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  return Finish(pass,"nonaffine parabolic material lines satisfy independently expanded local induction/rates, while whole-cell nonaffine interpolation explicitly rejects");
}
Result AtomicRejection() {
  const auto reference=Box();const auto interval=LinearInterval({{.1,0,0,0,.1,0,0,0,.1}},2);
  bool pass=true;auto ref=reference;ref.faces[0].firstMomentM3[0]*=1.01;
  auto rejected=M::PrepareAffineInterval(ref,interval);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  ref=reference;ref.faces.back().orientedAreaM2[0]=1;
  rejected=M::PrepareAffineInterval(ref,interval);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  ref=reference;ref.faces.back().magneticFluxWb=B0;
  rejected=M::PrepareAffineInterval(ref,interval);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  auto in=interval;in.final.deformation.declaredJacobian=1;
  rejected=M::PrepareAffineInterval(reference,in);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  // A midpoint rule integrates Fdot exactly but not the quadratic volume rate.
  // The zero declared GCL budget must reject rather than repair a side sweep.
  in=interval;in.quadrature.resize(1);in.quadrature[0].weightS=2;
  in.quadrature[0].sample=LinearSample({{.1,0,0,0,.1,0,0,0,.1}},{{0,0,0}},1);
  rejected=M::PrepareAffineInterval(reference,in);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  // An explicit allowance may accept this known quadrature error, but its
  // measured0.048m3 residual and original sweeps/endpoint MUST be retained.
  // This is a numerical-contract test, not adoption of a production tolerance.
  in.budget.volumeM3=.05;in.budget.provenance="controlled midpoint error0.048m3; explicitly declared0.05m3 numerical-only allowance";
  const auto allowed=M::PrepareAffineInterval(reference,in);pass=pass&&allowed.status.ok();
  if(allowed.status.ok())pass=pass&&Near(allowed.candidate->gclResidualM3,.048)&&Near(allowed.candidate->finalVolumeM3,41.472)&&
      Near(allowed.candidate->capSweptVolumeM3+allowed.candidate->lateralSweptVolumeM3,17.424)&&
      !allowed.candidate->qualifiedForProduction;
  in=interval;in.budget.declared=false;
  rejected=M::PrepareAffineInterval(reference,in);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  in=interval;in.quadrature[1].sample.deformation.ratePerS[0]*=1.1;
  rejected=M::PrepareAffineInterval(reference,in);pass=pass&&!rejected.status.ok()&&!rejected.candidate;
  M::PointInput point;point.deformation.gradient={{-1,0,0,0,1,0,0,0,1}};point.deformation.declaredJacobian=1;
  point.referenceMagneticFieldT={{B0,0,0}};point.referenceAreaVectorM2={{1,0,0}};point.referenceMassDensityKgM3=1e-20;
  pass=pass&&!M::PushForwardPoint(point).status.ok();point.deformation.gradient[0]=0;
  pass=pass&&!M::PushForwardPoint(point).status.ok();point.deformation.gradient=Identity();point.referenceMassDensityKgM3=0;
  pass=pass&&!M::PushForwardPoint(point).status.ok();
  // Independent exact binary64 counterexample: rounded Bx=1e11, By=1e-9,
  // A=(1,-1e20,0) have true dot≈-6.22815914578e-6Wb, whereas referenceflux
  // is1e-9Wb. A relative bound on huge cancelling products can hide the entire
  // invariant. The numerical capability must reject, never repair B or area.
  point.deformation.gradient={{1,1e20,0,0,1,0,0,0,1}};point.deformation.ratePerS={{0,0,0,0,0,0,0,0,0}};
  point.deformation.declaredJacobian=1;point.referenceMagneticFieldT={{1e-9,1e-9,0}};
  point.referenceAreaVectorM2={{1,0,0}};point.referenceMassDensityKgM3=1;
  const auto lostNonzero=M::PushForwardPoint(point);pass=pass&&!lostNonzero.status.ok();
  point.referenceAreaVectorM2={{1,-1,0}};
  const auto lostZero=M::PushForwardPoint(point);pass=pass&&!lostZero.status.ok();
  // The SAME huge shear is allowed when this B/area pair retains its flux:
  // rejection depends on representability, never an imposed physical shearcap.
  point.referenceMagneticFieldT={{1e-9,0,0}};point.referenceAreaVectorM2={{1,0,0}};
  const auto resolvedShear=M::PushForwardPoint(point);pass=pass&&resolvedShear.status.ok()&&
      resolvedShear.state.currentFluxWb==resolvedShear.state.referenceFluxWb;
  // Zero referenceflux has an absolute reference-scale numerical criterion,
  // not an undefined relative division or a nonzero physical flux floor.
  point.deformation.gradient=Identity();point.referenceMagneticFieldT={{1e-9,0,0}};point.referenceAreaVectorM2={{0,1,0}};
  const auto exactZero=M::PushForwardPoint(point);pass=pass&&exactZero.status.ok()&&exactZero.state.referenceFluxWb==0&&exactZero.state.currentFluxWb==0;
  point.referenceMagneticFieldT={{1e-9,std::nextafter(1e-9,std::numeric_limits<double>::infinity()),0}};
  point.referenceAreaVectorM2={{1,-1,0}};
  const auto tinyNonzero=M::PushForwardPoint(point);pass=pass&&tinyNonzero.status.ok()&&tinyNonzero.state.referenceFluxWb!=0&&
      tinyNonzero.state.referenceFluxWb==tinyNonzero.state.currentFluxWb;
  // IEEE subnormal physical flux is retained if representable. A product
  // smaller than binary64's range and an overflowing physical flux explicitly
  // reject; the error-free FMA domain must not silently substitute zero/inf.
  point.referenceMagneticFieldT={{std::numeric_limits<double>::denorm_min(),0,0}};point.referenceAreaVectorM2={{1,0,0}};
  const auto subnormal=M::PushForwardPoint(point);pass=pass&&subnormal.status.ok()&&
      subnormal.state.currentFluxWb==std::numeric_limits<double>::denorm_min();
  point.referenceAreaVectorM2[0]=std::numeric_limits<double>::denorm_min();
  const auto belowRange=M::PushForwardPoint(point);pass=pass&&!belowRange.status.ok()&&belowRange.state.jacobian==0;
  point.referenceMagneticFieldT[0]=point.referenceAreaVectorM2[0]=std::numeric_limits<double>::max()/4;
  const auto aboveRange=M::PushForwardPoint(point);pass=pass&&!aboveRange.status.ok()&&aboveRange.state.jacobian==0;
  // Changing the global rounding mode is outside this exact-expansion API.
  // Restore it before evaluating any further fixture or returning a verdict.
  point.referenceMagneticFieldT={{1e-9,0,0}};point.referenceAreaVectorM2={{1,0,0}};
  const int savedRounding=std::fegetround();const bool changed=std::fesetround(FE_UPWARD)==0;
  const auto unsupportedRounding=M::PushForwardPoint(point);
  const bool restored=std::fesetround(savedRounding)==0;
  pass=pass&&changed&&restored&&!unsupportedRounding.status.ok()&&unsupportedRounding.state.jacobian==0;
  pass=pass&&reference.volumeM3==24&&reference.materialMassKg==7e-19&&reference.faces[0].firstMomentM3[0]==12&&
      interval.quadrature.size()==2&&interval.budget.volumeM3==0;
  return Finish(pass,"malformed closed geometry/flux, incorrect Jacobian, unresolved time GCL, missing budget, inconsistent rate and singular/reversing/invalid physical points reject without candidate or input mutation");
}
Descriptor Make(const char* id,const char* name,TestCallback callback) {
  Descriptor d;d.id=id;d.name=name;d.description=name;d.group="material-tube";d.initialization=InitializationLevel::None;
  d.runtime=RuntimeClass::Routine;d.seedCapability=SeedCapability::Deterministic;d.repeatSupported=true;
  d.seedPolicy="deterministic exact material-map/geometry references";d.supportedBuildModes="shared native and neutral strict sanitizer callbacks";
  d.stateIsolation="private explicit reference geometry, immutable prepared candidates; no native/global publication";d.callback=callback;return d;
}
} // namespace
std::vector<Descriptor> MaterialTubeKinematicsDescriptors() {
  return {Make("MTG01","Pointwise Cauchy induction and Nanson area",PointwiseCauchy),
      Make("MTG02","Exact material dilation mass and cap-side GCL",IsotropicDilation),
      Make("MTG03","Shear translation and tilted material caps",ShearAndTranslation),
      Make("MTG04","Independent curved-reference material GCL",CurvedReference),
      Make("MTG05","Nonaffine local law and explicit cell limitation",NonAffinePointLaw),
      Make("MTG06","Atomic rejection of invalid material geometry",AtomicRejection)};
}
} }
