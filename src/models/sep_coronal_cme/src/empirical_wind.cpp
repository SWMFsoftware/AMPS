#include "sep_coronal_cme/empirical_wind.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace SEP { namespace CoronalCME { namespace {
// Polynomial helpers work in the normalized interval coordinate t in [0,1].
// Root isolation is recursive: derivative roots partition a polynomial into
// monotone intervals, within which a sign change brackets exactly one root.
double Poly(const std::vector<double>& c,double x) {
  double value=0.0;for(auto i=c.rbegin();i!=c.rend();++i)value=value*x+*i;return value;
}
std::vector<double> Derivative(const std::vector<double>& c) {
  std::vector<double> d;if(c.size()<2)return d;d.reserve(c.size()-1);
  for(std::size_t i=1;i<c.size();++i) d.push_back(i*c[i]);
  return d;
}
std::vector<double> Roots(const std::vector<double>& polynomial,double a,double b) {
  std::vector<double> roots;if(polynomial.size()<=1)return roots;
  if(polynomial.size()==2) { const double root=-polynomial[0]/polynomial[1];if(root>a&&root<b)roots.push_back(root);return roots; }
  std::vector<double> cuts{a};const auto stationary=Roots(Derivative(polynomial),a,b);
  cuts.insert(cuts.end(),stationary.begin(),stationary.end());cuts.push_back(b);
  for(std::size_t i=0;i+1<cuts.size();++i) {
    double left=cuts[i],right=cuts[i+1],fl=Poly(polynomial,left),fr=Poly(polynomial,right);
    if(std::abs(fl)<1e-13&&left>a)roots.push_back(left);
    if(fl*fr<0.0) { for(int n=0;n<100;++n){const double mid=0.5*(left+right),fm=Poly(polynomial,mid);if(fl*fm<=0){right=mid;fr=fm;}else{left=mid;fl=fm;}}roots.push_back(0.5*(left+right)); }
  }
  std::sort(roots.begin(),roots.end());roots.erase(std::unique(roots.begin(),roots.end(),
      [](double x,double y){return std::abs(x-y)<1e-10;}),roots.end());return roots;
}
bool Finite(double x){return std::isfinite(x);}
}

Core::Result<CertifiedPositiveProfile> CertifiedPositiveProfile::Create(
    std::vector<ProfileNode> nodes,double reference,bool noOvershoot) {
  if(nodes.size()<2||!(reference>0.0))return Core::Result<CertifiedPositiveProfile>::Failure(
      Core::StatusCode::InvalidConfiguration,"certified profile requires nodes and positive reference");
  CertifiedPositiveProfile result;result.reference_=reference;result.nodes_=nodes;
  for(std::size_t n=0;n+1<nodes.size();++n) {
    const auto& a=nodes[n];const auto& b=nodes[n+1];const double h=b.x-a.x;
    if(!(h>0&&a.value>0&&b.value>0&&Finite(a.first)&&Finite(a.second)&&Finite(b.first)&&Finite(b.second)))
      return Core::Result<CertifiedPositiveProfile>::Failure(Core::StatusCode::InvalidConfiguration,
          "profile nodes must increase and carry finite positive states/derivatives");
    const double za=std::log(a.value/reference),zb=std::log(b.value/reference);
    const double da=a.first/a.value,db=b.first/b.value;
    const double dda=a.second/a.value-da*da,ddb=b.second/b.value-db*db;
    // c0--c2 satisfy value/first/second derivative at the left node. The
    // remaining 3x3 system imposes those same physical constraints at t=1.
    std::vector<double> c(6);c[0]=za;c[1]=h*da;c[2]=0.5*h*h*dda;
    double matrix[3][4]={{1,1,1,zb-c[0]-c[1]-c[2]},
      {3,4,5,h*db-c[1]-2*c[2]},
      {6,12,20,h*h*ddb-2*c[2]}};
    for(int p=0;p<3;++p){int pivot=p;for(int r=p+1;r<3;++r)if(std::abs(matrix[r][p])>std::abs(matrix[pivot][p]))pivot=r;
      for(int q=p;q<4;++q)std::swap(matrix[p][q],matrix[pivot][q]);
      if(std::abs(matrix[p][p])<1e-14)return Core::Result<CertifiedPositiveProfile>::Failure(Core::StatusCode::NumericalFailure,"singular quintic system");
      const double scale=matrix[p][p];for(int q=p;q<4;++q)matrix[p][q]/=scale;
      for(int r=0;r<3;++r)if(r!=p){const double factor=matrix[r][p];for(int q=p;q<4;++q)matrix[r][q]-=factor*matrix[p][q];}}
    c[3]=matrix[0][3];c[4]=matrix[1][3];c[5]=matrix[2][3];
    if(noOvershoot) {
      const auto roots=Roots(Derivative(c),0.0,1.0);const double lo=std::min(za,zb)-1e-12,hi=std::max(za,zb)+1e-12;
      for(double root:roots)if(Poly(c,root)<lo||Poly(c,root)>hi)
        return Core::Result<CertifiedPositiveProfile>::Failure(Core::StatusCode::InvalidConfiguration,
            "quintic profile has a certified adjacent-node overshoot");
    }
    result.coefficients_.push_back(std::move(c));
  }
  return Core::Result<CertifiedPositiveProfile>::Success(std::move(result));
}

Core::Result<ProfileNode> CertifiedPositiveProfile::Evaluate(double x) const {
  if(nodes_.empty()||x<nodes_.front().x||x>nodes_.back().x)return Core::Result<ProfileNode>::Failure(
      Core::StatusCode::OutOfDomain,"profile extrapolation is forbidden");
  std::size_t segment=nodes_.size()-2;for(std::size_t i=0;i+1<nodes_.size();++i)if(x<=nodes_[i+1].x){segment=i;break;}
  const double h=nodes_[segment+1].x-nodes_[segment].x,t=(x-nodes_[segment].x)/h;
  const auto& c=coefficients_[segment];const auto d=Derivative(c),dd=Derivative(d);
  const double z=Poly(c,t),zt=Poly(d,t),ztt=Poly(dd,t),value=reference_*std::exp(z);
  return Core::Result<ProfileNode>::Success({x,value,value*zt/h,
      value*(zt*zt+ztt)/(h*h)});
}

Core::Result<TwoZoneState> BlendTwoZoneWind(double r,double ra,double rb,double inner,
    double outerSpeed,double field,double eta,double reference) {
  if(!(ra<rb&&inner>0&&outerSpeed>0&&field>0&&eta>0&&reference>0))
    return Core::Result<TwoZoneState>::Failure(Core::StatusCode::InvalidConfiguration,
        "two-zone wind requires ordered join radii and positive physical inputs");
  // One mass-per-flux value owns both zones. The outer density is derived,
  // never adjusted to make the independently observed inner product agree.
  const double outer=eta*field/outerSpeed;const double mismatch=std::log(inner/outer);
  double rho;
  if(r<=ra)rho=inner;else if(r>=rb)rho=outer;else{const double x=(r-ra)/(rb-ra);const double w=10*x*x*x-15*x*x*x*x+6*x*x*x*x*x;
    rho=reference*std::exp((1-w)*std::log(inner/reference)+w*std::log(outer/reference));}
  return Core::Result<TwoZoneState>::Success({rho,eta*field/rho,mismatch});
}

Core::Result<double> ConvertToCorotatingFieldAlignedSpeed(double stored,
    VelocityComponent component,VelocityFrame frame,double projection,double minimum,
    Vec3 omega,Vec3 position,Vec3 tangent) {
  if(!(stored>0&&Finite(stored)))return Core::Result<double>::Failure(
      Core::StatusCode::InvalidConfiguration,"velocity profile must be finite and positive");
  tangent=Unit(tangent);if(Norm(tangent)==0)return Core::Result<double>::Failure(
      Core::StatusCode::InvalidConfiguration,"profile tangent is undefined");
  if(component==VelocityComponent::Radial) {
    if(frame!=VelocityFrame::Inertial)return Core::Result<double>::Failure(
        Core::StatusCode::UnsupportedCapability,"radial/corotating is not schema-5 canonical");
    if(!(minimum>0&&minimum<=1&&projection>=minimum&&projection<=1))return Core::Result<double>::Failure(
        Core::StatusCode::InvalidConfiguration,"radial velocity projection guard failed");
    return Core::Result<double>::Success(stored/projection);
  }
  if(minimum!=0.0)return Core::Result<double>::Failure(Core::StatusCode::InvalidConfiguration,
      "field-aligned channels require an inactive zero radial guard");
  const double rotation=Dot(Cross(omega,position),tangent);
  return Core::Result<double>::Success(frame==VelocityFrame::Inertial?stored-rotation:stored);
}

Core::Result<CoverageCensus> BuildCoverageCensus(const std::vector<ConsumerMeasure>& consumers,
    bool eventNominal) {
  if(consumers.empty())return Core::Result<CoverageCensus>::Failure(
      Core::StatusCode::InvalidConfiguration,"coverage census requires consumers");
  double totals[5]={},covered[5]={};CoverageCensus result;
  for(const auto& c:consumers){const double values[5]={c.openFluxWb,c.openAreaM2,c.sourceNumberRate,c.observerExposureM2S,c.exportLengthM};
    for(int i=0;i<5;++i){if(values[i]<0||!Finite(values[i]))return Core::Result<CoverageCensus>::Failure(Core::StatusCode::InvalidConfiguration,"consumer measure is invalid");totals[i]+=values[i];if(c.covered)covered[i]+=values[i];}
    if(!c.covered)result.rejectedIds.push_back(c.stableId+":"+c.rejectionReason);
  }
  double* fractions[5]={&result.coveredFluxFraction,&result.coveredAreaFraction,&result.coveredSourceFraction,&result.coveredObserverFraction,&result.coveredExportFraction};
  for(int i=0;i<5;++i)*fractions[i]=totals[i]>0?covered[i]/totals[i]:1.0;
  if(eventNominal&&!result.rejectedIds.empty())return Core::Result<CoverageCensus>::Failure(
      Core::StatusCode::DataIntegrityFailure,"event-nominal consumer support is incomplete");
  return Core::Result<CoverageCensus>::Success(std::move(result));
}
} }
