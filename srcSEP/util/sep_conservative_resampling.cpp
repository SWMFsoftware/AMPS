#include "sep_conservative_resampling.h"
#include <cmath>
#include <limits>
#include <set>

namespace SEP { namespace PopulationControl { namespace Conservative {
namespace {
UnpreservedMoments Moments(const std::vector<Particle>& particles,
                           const std::vector<double>& weights,double c) {
  long double current=0,parallel=0,perpendicular=0;
  for(std::size_t i=0;i<particles.size();++i) {
    const auto& p=particles[i];
    const long double vp=p.parallelMPerS,vn=p.normalMPerS;
    const long double gamma=1.0L/std::sqrt(1.0L-(vp*vp+vn*vn)/(static_cast<long double>(c)*c));
    current+=weights[i]*static_cast<long double>(p.chargeC)*vp;
    parallel+=weights[i]*gamma*p.massKg*vp*vp;
    perpendicular+=0.5L*weights[i]*gamma*p.massKg*vn*vn;
  }
  UnpreservedMoments result;
  result.parallelCurrentCMPerS=static_cast<double>(current);
  result.parallelStressJ=static_cast<double>(parallel);
  result.perpendicularStressJ=static_cast<double>(perpendicular);
  return result;
}
} // namespace
Transport::Status EvaluateConstraints(const std::vector<Particle>& particles,double c,
    std::vector<std::vector<double>>* outputRows,std::vector<double>* outputWeights) {
  const auto invalid=[](){return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
      "resampling requires finite same-species intrinsic particles and subluminal velocities");};
  if(!outputRows||!outputWeights||particles.empty()||!std::isfinite(c)||c<=0)return invalid();
  std::set<std::uint64_t> identities;
  std::vector<double> weights;
  std::vector<std::vector<double>> rows(6);
  for(const auto& p:particles) {
    if(p.identity==0||!identities.insert(p.identity).second||p.species<0||
       p.species!=particles[0].species||p.massKg!=particles[0].massKg||
       p.chargeC!=particles[0].chargeC||!std::isfinite(p.massKg)||p.massKg<=0||
       !std::isfinite(p.chargeC)||!std::isfinite(p.localFraction)||
       p.localFraction<0||p.localFraction>1||!std::isfinite(p.parallelMPerS)||
       !std::isfinite(p.normalMPerS)||p.normalMPerS<0||
       !std::isfinite(p.weight)||p.weight<=0)return invalid();
    const long double vp=p.parallelMPerS,vn=p.normalMPerS,v2=vp*vp+vn*vn;
    const long double c2=static_cast<long double>(c)*c;
    if(v2>=c2)return invalid();
    const long double gamma=1.0L/std::sqrt(1.0L-v2/c2);
    // gamma-1=gamma^2*beta^2/(gamma+1), avoiding low-speed cancellation.
    const long double kinetic=gamma*gamma*p.massKg*v2/(gamma+1.0L);
    const double mu=v2?static_cast<double>(vp/std::sqrt(v2)):0.0;
    rows[0].push_back(1.0);rows[1].push_back(static_cast<double>(kinetic));
    rows[2].push_back(static_cast<double>(gamma*p.massKg*vp));
    rows[3].push_back(mu);rows[4].push_back(mu*mu);rows[5].push_back(p.localFraction);
    weights.push_back(p.weight);
  }
  outputRows->swap(rows);outputWeights->swap(weights);
  return Transport::Status::Ok();
}
double MaximumRelativeResidual(const std::vector<std::vector<double>>& rows,
    const std::vector<double>& before,const std::vector<double>& after) {
  if(before.size()!=after.size()||rows.empty())return std::numeric_limits<double>::infinity();
  double maximum=0;
  for(const auto& row:rows) {
    if(row.size()!=before.size())return std::numeric_limits<double>::infinity();
    long double a=0,b=0,scale=0;
    for(std::size_t i=0;i<before.size();++i) {
      if(!std::isfinite(row[i])||!std::isfinite(before[i])||!std::isfinite(after[i])||
         before[i]<0||after[i]<0)return std::numeric_limits<double>::infinity();
      const long double contribution=static_cast<long double>(before[i])*row[i];
      a+=contribution;b+=static_cast<long double>(after[i])*row[i];scale+=std::fabs(contribution);
    }
    const long double error=std::fabs(a-b);
    maximum=std::max(maximum,scale?static_cast<double>(error/scale):
        (error==0?0.0:std::numeric_limits<double>::infinity()));
  }
  return maximum;
}
UnpreservedMoments EvaluateUnpreserved(const std::vector<Particle>& particles,
    const std::vector<double>& weights,double c) {
  return Moments(particles,weights,c);
}
Endpoints Prepare(const std::vector<Particle>& particles,double tolerance,double c) {
  Endpoints result;
  std::vector<double> weights;
  std::vector<std::vector<double>> rows;
  const auto status=EvaluateConstraints(particles,c,&rows,&weights);
  if(!status.ok()) {result.weights.message=status.message;return result;}
  result.weights=PIC::ConservativeResampling::Prepare(weights,rows,tolerance);
  if(result.weights.status==PIC::ConservativeResampling::Status::Success) {
    result.before=Moments(particles,weights,c);
    result.positive=Moments(particles,result.weights.positive,c);
    result.negative=Moments(particles,result.weights.negative,c);
  }
  return result;
}
Transport::Status Select(const Endpoints& endpoints,double uniform,
                         std::vector<double>* weights,UnpreservedMoments* diagnostics) {
  if(!weights||endpoints.weights.status!=PIC::ConservativeResampling::Status::Success||
     !std::isfinite(uniform)||uniform<=0.0||uniform>=1.0)
    return Transport::Status::Error(Transport::StatusCode::InvalidArgument,
        "conservative resampling requires validated endpoints and an open uniform selector");
  const bool positive=uniform<endpoints.weights.positiveProbability;
  *weights=positive?endpoints.weights.positive:endpoints.weights.negative;
  if(diagnostics)*diagnostics=positive?endpoints.positive:endpoints.negative;
  return Transport::Status::Ok();
}
}}} // namespace SEP::PopulationControl::Conservative
