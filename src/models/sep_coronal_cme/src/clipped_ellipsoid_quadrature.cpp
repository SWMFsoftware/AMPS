#include "sep_coronal_cme/clipped_ellipsoid_quadrature.h"

#include "sep_coronal_cme/constants.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SEP { namespace CoronalCME { namespace {
using Real=long double;
using Interval=std::pair<Real,Real>;
constexpr Real Pi=3.141592653589793238462643383279502884L;
constexpr unsigned MaximumDepth=20;
constexpr std::uint64_t MaximumEvaluations=50000000;
constexpr std::uint64_t MaximumCells=1048576;

struct Work {
  std::uint64_t evaluations=0;
  void Count() {
    if(++evaluations>MaximumEvaluations)
      throw std::runtime_error("solar cut-cell quadrature exhausted its integrand work bound");
  }
};

// All integration is dimensionless. Besides avoiding SI squares/products
// overflowing, long-double construction preserves the small difference between
// a solar radius and a just-emerging front before any canonical double sample
// is chosen. The eventual sample must independently be representably exposed.
struct Geometry {
  Real d,a,b,c,radius,scale;
  explicit Geometry(const FixedOrientationEllipsoid& shape,double solarRadius) {
    const auto& k=shape.Kinematics();
    scale=std::max({static_cast<Real>(k.centerDistanceM.value),
        static_cast<Real>(k.radialSemiAxisM.value),static_cast<Real>(k.firstLateralSemiAxisM.value),
        static_cast<Real>(k.secondLateralSemiAxisM.value),static_cast<Real>(solarRadius)});
    d=k.centerDistanceM.value/scale;a=k.radialSemiAxisM.value/scale;
    b=k.firstLateralSemiAxisM.value/scale;c=k.secondLateralSemiAxisM.value/scale;
    radius=solarRadius/scale;
  }
  Real LateralSquared(Real phi) const {
    const Real co=std::cos(phi),si=std::sin(phi);
    return b*b*co*co+c*c*si*si;
  }
  Real SignedCut(Real u,Real lateral) const {
    const Real radial=d+a*u;
    return (radial-radius)*(radial+radius)+lateral*(1-u)*(1+u);
  }
  Real AreaJacobian(Real u,Real phi) const {
    const Real co=std::cos(phi),si=std::sin(phi);
    const Real value=u*u/(a*a)+(1-u)*(1+u)*(co*co/(b*b)+si*si/(c*c));
    if(!(std::isfinite(value)&&value>0))
      throw std::runtime_error("solar cut-cell Jacobian is not finite and positive");
    return a*b*c*std::sqrt(value);
  }
  std::vector<Interval> Exposed(Real lo,Real hi,Real phi) const {
    const Real lateral=LateralSquared(phi),aa=a*a-lateral,bb=2*d*a;
    const Real cc=(d-radius)*(d+radius)+lateral;
    std::vector<Real> cuts{lo,hi};
    const auto append=[&](Real root){if(std::isfinite(root)&&root>lo&&root<hi)cuts.push_back(root);};
    if(aa==0)append(-cc/bb);
    else {
      const Real discriminant=std::fma(bb,bb,-4*aa*cc);
      if(discriminant>=0) {
        // The q form avoids subtracting two almost equal numbers for an
        // emerging cap's near-pole root. No root is moved to create area.
        const Real q=-.5L*(bb+std::copysign(std::sqrt(discriminant),bb));
        if(q==0)append(-bb/(2*aa));
        else {append(q/aa);append(cc/q);}
      }
    }
    std::sort(cuts.begin(),cuts.end());
    cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    std::vector<Interval> result;
    for(std::size_t i=1;i<cuts.size();++i) {
      const Real middle=cuts[i-1]+(cuts[i]-cuts[i-1])/2;
      if(SignedCut(middle,lateral)>0)result.emplace_back(cuts[i-1],cuts[i]);
    }
    return result;
  }
};

// Positive-weight nested estimates are evaluated separately, rather than
// reusing a coarse midpoint's inside/outside classification. Explicit topology
// splits below prevent an unseen exposed azimuthal island between nodes.
template<std::size_t N,class Function>
Real Gauss(const Function& function,Real lo,Real hi,
    const std::array<Real,N>& nodes,const std::array<Real,N>& weights,Work& work) {
  const Real middle=lo+(hi-lo)/2,half=(hi-lo)/2;
  Real sum=0;
  for(std::size_t i=0;i<N;++i) {
    work.Count();const Real left=function(middle-half*nodes[i]);
    work.Count();const Real right=function(middle+half*nodes[i]);
    if(!(std::isfinite(left)&&std::isfinite(right)&&left>=0&&right>=0))
      throw std::runtime_error("solar cut-cell area integrand is nonfinite or negative");
    sum+=weights[i]*(left+right);
  }
  return half*sum;
}
template<class Function>
Real Integrate(const Function& function,Real lo,Real hi,Real tolerance,Work& work,unsigned depth=0) {
  static constexpr std::array<Real,4> x8{{.183434642495649804939476142360184L,.525532409916328985817739049189246L,.796666477413626739591553936475831L,.960289856497536231683560868569473L}};
  static constexpr std::array<Real,4> w8{{.362683783378361982965150449277195L,.313706645877887287337962201986601L,.222381034453374470544355994426240L,.101228536290376259152531354309962L}};
  static constexpr std::array<Real,8> x16{{.095012509837637440185319335424958L,.281603550779258913230460501460496L,.458016777657227386342419442983578L,.617876244402643748446671764048791L,.755404408355003033895101194847442L,.865631202387831743880467897712393L,.944575023073232576077988415534608L,.989400934991649932596154173450333L}};
  static constexpr std::array<Real,8> w16{{.189450610455068496285396723208283L,.182603415044923588866763667969219L,.169156519395002538189312079030359L,.149595988816576732081501730547479L,.124628971255533872052476282192016L,.095158511682492784809925107602246L,.062253523938647892862843836994378L,.027152459411754094851780572456018L}};
  if(hi<=lo)return 0;
  const Real coarse=Gauss(function,lo,hi,x8,w8,work),fine=Gauss(function,lo,hi,x16,w16,work);
  if(std::abs(fine-coarse)<=tolerance*std::abs(fine))return fine;
  const Real middle=lo+(hi-lo)/2;
  if(depth>=MaximumDepth||middle==lo||middle==hi)
    throw std::runtime_error("solar cut-cell area quadrature did not converge within its subdivision bound");
  return Integrate(function,lo,middle,tolerance,work,depth+1)+
      Integrate(function,middle,hi,tolerance,work,depth+1);
}

void AddLateralEvents(const Geometry& g,Real lateral,Real lo,Real hi,std::vector<Real>* cuts) {
  if(g.b==g.c)return;
  const Real ratio=(lateral-g.b*g.b)/(g.c*g.c-g.b*g.b);
  if(!(ratio>=0&&ratio<=1))return;
  const Real angle=std::asin(std::sqrt(ratio));
  for(Real phi:{angle,Pi-angle,Pi+angle,2*Pi-angle})
    if(phi>lo&&phi<hi)cuts->push_back(phi);
}
std::vector<Real> AzimuthCuts(const Geometry& g,Real uLo,Real uHi,Real lo,Real hi) {
  std::vector<Real> cuts{lo,hi};
  for(Real phi:{Pi/2,Pi,3*Pi/2})if(phi>lo&&phi<hi)cuts.push_back(phi);
  // A cut enters/leaves a polar cell only at its u boundary or a double
  // quadratic root. Include the harmless quadratic-to-linear event too.
  for(Real u:{uLo,uHi})if(u>-1&&u<1) {
    const Real radial=g.d+g.a*u;
    AddLateralEvents(g,(g.radius-radial)*(g.radius+radial)/((1-u)*(1+u)),lo,hi,&cuts);
  }
  AddLateralEvents(g,g.a*g.a,lo,hi,&cuts);
  const Real coefficient=g.d*g.d-g.a*g.a-g.radius*g.radius;
  const Real discriminant=std::fma(coefficient,coefficient,-4*g.a*g.a*g.radius*g.radius);
  if(discriminant>=0) {
    const Real q=-.5L*(coefficient+std::copysign(std::sqrt(discriminant),coefficient));
    if(q==0)AddLateralEvents(g,-coefficient/2,lo,hi,&cuts);
    else {AddLateralEvents(g,q,lo,hi,&cuts);AddLateralEvents(g,g.a*g.a*g.radius*g.radius/q,lo,hi,&cuts);}
  }
  std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
  return cuts;
}
} // namespace

Core::Result<std::vector<SurfacePatch>> TessellateSolarExposed(
    const FixedOrientationEllipsoid& shape,double solarRadiusM,int polarCells,int azimuthCells) {
  using Return=Core::Result<std::vector<SurfacePatch>>;
  if(!std::isfinite(solarRadiusM)||solarRadiusM<0||polarCells<2||azimuthCells<3||
      static_cast<std::uint64_t>(polarCells)*static_cast<std::uint64_t>(azimuthCells)>MaximumCells)
    return Return::Failure(Core::StatusCode::InvalidConfiguration,
        "solar quadrature requires finite nonnegative radius, >=2x3 cells and <=1048576 cells");
  const auto validated=FixedOrientationEllipsoid::FromCenter(shape.Basis(),shape.Kinematics(),solarRadiusM);
  if(!validated.ok())return Return::Failure(validated.status.code,validated.status.message);
  try {
    const Geometry geometry(shape,solarRadiusM);Work work;std::vector<SurfacePatch> patches;
    for(int polar=0;polar<polarCells;++polar) {
      const Real uHi=polar==0?1:std::cos(Pi*polar/polarCells);
      const Real uLo=polar+1==polarCells?-1:std::cos(Pi*(polar+1)/polarCells);
      for(int azimuth=0;azimuth<azimuthCells;++azimuth) {
        const Real phiLo=2*Pi*azimuth/azimuthCells,phiHi=2*Pi*(azimuth+1)/azimuthCells;
        const auto cuts=AzimuthCuts(geometry,uLo,uHi,phiLo,phiHi);
        Real area=0,bestClearance=-1,bestU=0,bestPhi=0;
        for(std::size_t i=1;i<cuts.size();++i) {
          const Real lower=cuts[i-1],width=cuts[i]-lower,phiMiddle=lower+width/2;
          // Topology is fixed on this OPEN azimuth interval. Its midpoint
          // therefore cannot miss a positive-width island; select an exposed
          // representative privately, independently of the area quadrature.
          for(const auto& interval:geometry.Exposed(uLo,uHi,phiMiddle)) {
            const Real u=interval.first+(interval.second-interval.first)/2;
            const Real clearance=geometry.SignedCut(u,geometry.LateralSquared(phiMiddle));
            if(clearance>bestClearance){bestClearance=clearance;bestU=u;bestPhi=phiMiddle;}
          }
          const auto outer=[&](Real t) {
            // The sin^2 change of variable removes the endpoint square-root
            // opening of a tangential cut. No step function remains for the
            // adaptive area rule to discover by luck at a quadrature node.
            const Real sine=std::sin(Pi*t/2),phi=lower+width*sine*sine;
            Real slice=0;
            for(const auto& interval:geometry.Exposed(uLo,uHi,phi)) {
              const auto inner=[&](Real u){return geometry.AreaJacobian(u,phi);};
              if(interval.first<0&&interval.second>0)
                slice+=Integrate(inner,interval.first,0,2e-13L,work)+Integrate(inner,0,interval.second,2e-13L,work);
              else slice+=Integrate(inner,interval.first,interval.second,2e-13L,work);
            }
            return slice*width*Pi/2*std::sin(Pi*t);
          };
          area+=Integrate(outer,0,1,2e-11L,work);
        }
        if(area==0)continue;
        const Real physicalArea=area*geometry.scale*geometry.scale;
        if(!(std::isfinite(physicalArea)&&physicalArea>0&&physicalArea<=std::numeric_limits<double>::max())||bestClearance<=0)
          return Return::Failure(Core::StatusCode::NumericalFailure,"positive solar cut-cell area/sample is not representable");
        const double theta=static_cast<double>(std::acos(bestU)),phi=static_cast<double>(bestPhi);
        const Vec3 point=shape.Point(theta,phi);
        if(!(std::isfinite(Norm(point))&&Norm(point)>=solarRadiusM))
          return Return::Failure(Core::StatusCode::NumericalFailure,"positive exposed cell lacks a representably exposed canonical sample");
        const auto normal=shape.Evaluate(point);
        if(!normal.ok())return Return::Failure(normal.status.code,normal.status.message);
        // The public payload is double precision even though the cut/integral
        // staging is wider. Reject an underflowed area or a nonfinite canonical
        // normal instead of publishing a zero-area/undefined physical patch.
        const double storedArea=static_cast<double>(physicalArea);
        const double normalMagnitude=Norm(normal.value.outwardNormal);
        if(!(std::isfinite(storedArea)&&storedArea>0&&std::isfinite(normalMagnitude)&&normalMagnitude>0))
          return Return::Failure(Core::StatusCode::NumericalFailure,"solar cut-cell double payload is not finite and positive");
        patches.push_back({1+static_cast<std::uint64_t>(polar)*static_cast<std::uint64_t>(azimuthCells)+
            static_cast<std::uint64_t>(azimuth),point,normal.value.outwardNormal,
            storedArea,theta,phi});
      }
    }
    return Return::Success(std::move(patches));
  }catch(const std::exception& error) {
    return Return::Failure(Core::StatusCode::NumericalFailure,error.what());
  }
}

} }
