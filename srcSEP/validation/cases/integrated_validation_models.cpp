#include "integrated_validation_models.h"

#include "../../util/sep_common_header_path.h"
#include "../../util/sep_focused_transport_core.h"
#include "../../util/sep_flux_tube_geometry_core.h"
#include "../../util/sep_shock_source_core.h"
#include SRCSEP_SEP_COMMON_HEADER(sep_transport_common.h)
#include "../../util/sep_turbulence_core.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using SEP::Transport::KeyedRandomStream;
using SEP::Transport::Status;
const double Pi = 3.1415926535897932384626433832795;
const double AuM = 1.495978707e11;
const double SolarRotationPerS = 2.86533e-6;
const double ProtonMassKg = 1.67262192369e-27;
const double LightSpeedMPerS = 299792458.0;

// Integrated-case design notes:
// * IV01 stores physical arclength [m], radius [m], time [s], and momentum in
//   SI. Reversed vertex order is metadata and must not reverse physical +s.
// * IV02 couples the production bounded pitch SDE to a constant-focusing
//   coefficient whose zero-flux distribution is independently integrable.
// * IV03 compares an analytic manufactured residual with centered numerical
//   derivatives. Every operator contribution is smooth and the solution stays
//   positive, so a convergence plateau is an implementation defect.
// * IV04 compares conserved cell integrals, not point values: a uniform unit
//   density on a moved grid has exact cell energy equal to segment length.
// * IV05 brackets then linearly interpolates shock crossing in physical time;
//   stationary/moving frames reuse the same keyed shock-cycle realization.
// * IV06 removes every wave-growth deposit from a finite particle reservoir.
//   Thus the streaming-growth-scattering timeline and total-energy invariant
//   are observable simultaneously instead of being separate synthetic tests.

std::map<std::string,std::string> Parse(const std::vector<std::string>& args) {
  if (args.empty() || args.size()%2) throw std::runtime_error("IV arguments must be name/value pairs");
  std::map<std::string,std::string> result;
  for (std::size_t i=0;i<args.size();i+=2) {
    if (args[i].size()<3 || args[i].substr(0,2)!="--") throw std::runtime_error("IV option must begin with --");
    if (!result.insert(std::make_pair(args[i].substr(2),args[i+1])).second) throw std::runtime_error("duplicate IV option");
  }
  return result;
}

std::string Need(const std::map<std::string,std::string>& v,const std::string& n) {
  std::map<std::string,std::string>::const_iterator i=v.find(n);
  if(i==v.end()||i->second.empty()) throw std::runtime_error("missing IV option --"+n);
  return i->second;
}

double Num(const std::map<std::string,std::string>& v,const std::string& n) {
  const std::string s=Need(v,n); char* e=NULL; errno=0; const double x=std::strtod(s.c_str(),&e);
  if(errno==ERANGE||!e||*e!='\0'||!std::isfinite(x)) throw std::runtime_error("invalid IV number --"+n);
  return x;
}

std::uint64_t UInt(const std::map<std::string,std::string>& v,const std::string& n) {
  const std::string s=Need(v,n); char* e=NULL; errno=0; const unsigned long long x=std::strtoull(s.c_str(),&e,10);
  if(s[0]=='-'||errno==ERANGE||!e||*e!='\0') throw std::runtime_error("invalid IV integer --"+n);
  return static_cast<std::uint64_t>(x);
}

std::vector<double> Nums(const std::map<std::string,std::string>& v,const std::string& n) {
  std::vector<double> out; const std::string s=Need(v,n); std::size_t b=0;
  while(b<s.size()){const std::size_t c=s.find(',',b);const std::string t=s.substr(b,c==std::string::npos?std::string::npos:c-b);char*e=NULL;const double x=std::strtod(t.c_str(),&e);if(t.empty()||!e||*e!='\0'||!std::isfinite(x))throw std::runtime_error("invalid IV list --"+n);out.push_back(x);if(c==std::string::npos)break;b=c+1;}
  return out;
}

void Commit(std::ofstream* f,const std::string& tmp,const std::string& dst){f->close();if(f->fail())throw std::runtime_error("cannot flush IV CSV");if(std::rename(tmp.c_str(),dst.c_str()))throw std::runtime_error("cannot publish IV CSV");}

class ZeroPitch final : public SEP::Transport::PitchAngleDiffusionProvider {
 public: SEP::Transport::PitchAngleDiffusionSample Evaluate(double,double,double)const override {SEP::Transport::PitchAngleDiffusionSample s;s.status=Status::Ok();s.provenance="validation:zero";s.turbulenceStateIdentity="validation:fixed";return s;}
};

class ConstantPitch final : public SEP::Transport::PitchAngleDiffusionProvider {
 public: explicit ConstantPitch(double d):d_(d){}
  SEP::Transport::PitchAngleDiffusionSample Evaluate(double,double,double mu)const override {SEP::Transport::PitchAngleDiffusionSample s;s.status=Status::Ok();s.dMuMuPerS=d_*std::max(0.0,1.0-mu*mu);s.dDmuMuDmuPerS=-2*d_*mu;s.provenance="validation:D0(1-mu2)";s.turbulenceStateIdentity="validation:fixed";return s;}
 private: double d_;
};

double MomentumForSpeed(double speed){const SEP::Transport::ScalarResult p=SEP::Transport::MomentumFromSpeed(speed,ProtonMassKg,LightSpeedMPerS);if(!p.status.ok())throw std::runtime_error(p.status.message);return p.value;}
double Momentum(){return MomentumForSpeed(1.0e7);}

// Parker geometry helpers use r as the independent coordinate.  The arclength
// metric and |B| derivative are evaluated analytically, avoiding any implicit
// dependence on production field-line segment numbering.
double DsDr(double r,double wind){const double a=SolarRotationPerS/wind;return std::sqrt(1+a*a*r*r);}
double DlnBds(double r,double wind){const double a=SolarRotationPerS/wind;return(-2/r+a*a*r/(1+a*a*r*r))/DsDr(r,wind);}
double ArcLength(double r0,double r1,double wind,int panels){double sum=0,dr=(r1-r0)/panels;for(int i=0;i<panels;++i){double a=r0+i*dr,b=a+dr,m=.5*(a+b);sum+=dr*(DsDr(a,wind)+4*DsDr(m,wind)+DsDr(b,wind))/6;}return sum;}
double ParkerArcPrimitive(double r,double wind){const double a=SolarRotationPerS/wind;return .5*(r*std::sqrt(1+a*a*r*r)+std::asinh(a*r)/a);}
double RadiusAtArc(double s,double r0,double r1,double wind,double total){
  // Newton inversion of the closed Parker arclength primitive gives the
  // physical radius sampled by the mover. Clamp every iterate to the field-
  // line interval so roundoff near either endpoint cannot leave the geometry.
  double r=r0+(r1-r0)*std::max(0.0,std::min(1.0,s/total));const double base=ParkerArcPrimitive(r0,wind);
  for(int iteration=0;iteration<6;++iteration){r-=(ParkerArcPrimitive(r,wind)-base-s)/DsDr(r,wind);r=std::max(r0,std::min(r1,r));}
  return r;
}

// Ballistic rows isolate deterministic geometry/focusing; weak rows retain all
// particles so the external scorer can form moments and uncertainty bands.
void IV01(const std::map<std::string,std::string>& v,const std::string& path){
  const double r0=Num(v,"inner-radius-m"),r1=Num(v,"outer-radius-m"),speed=Num(v,"particle-speed-m-per-s"),d0=Num(v,"weak-d0-per-s");
  const std::vector<double>winds=Nums(v,"wind-speeds-m-per-s"),dts=Nums(v,"time-steps-s");const std::uint64_t particles=UInt(v,"particles"),seed=UInt(v,"campaign-seed");
  if(!(r1>r0&&r0>0&&speed>0&&d0>=0)||particles<100)throw std::runtime_error("IV01 invalid domain");
  std::ofstream o((path+".tmp").c_str());o<<std::setprecision(17)<<"wind_m_per_s,orientation,mode,dt_s,particle,arrival_time_s,path_length_m,final_mu,weight\n";
  ZeroPitch zero; const double p=MomentumForSpeed(speed);
  for(std::size_t wi=0;wi<winds.size();++wi)for(int orientation=-1;orientation<=1;orientation+=2)for(int mode=0;mode<2;++mode)for(std::size_t di=0;di<dts.size();++di){
    const double wind=winds[wi],dt=dts[di],length=ArcLength(r0,r1,wind,20000);const std::uint64_t count=mode?particles:1;
    for(std::uint64_t id=0;id<count;++id){double s=0,mu=mode?0.8:0.65,t=0;std::uint64_t step=0;ConstantPitch weak(d0);
      while(s<length&&step<1000000){const double r=RadiusAtArc(s,r0,r1,wind,length);
        // Storage orientation is metadata only: physical +s is always outward
        // and the adapter must remove any reversed vertex ordering before the
        // production mover sees dln|B|/ds. Keeping identical physics in both
        // rows makes an orientation-dependent regression directly observable.
        SEP::Transport::FocusedTransportBackground b(DlnBds(r,wind),0,0,0);KeyedRandomStream rng(seed+id,id+1,101+mode,step);const double h=std::min(dt,(length-s)/std::max(speed*std::max(mu,0.05),1.0));const SEP::Transport::FocusedTransportIncrement q=SEP::Transport::AdvanceFocusedTransportDmumu({s,p,mu},b,ProtonMassKg,LightSpeedMPerS,h,mode?static_cast<const SEP::Transport::PitchAngleDiffusionProvider&>(weak):static_cast<const SEP::Transport::PitchAngleDiffusionProvider&>(zero),rng,NULL);if(!q.status.ok())throw std::runtime_error(q.status.message);s=q.state.arcLengthM;mu=q.state.mu;t+=h;++step;if(mode&&t>4*length/speed)break;}
      o<<wind<<','<<orientation<<','<<(mode?"weak":"ballistic")<<','<<dt<<','<<id<<','<<t<<','<<length<<','<<mu<<",1\n";
    }
  }Commit(&o,path+".tmp",path);
}

// dlnB/ds=-2*D0*xi/v yields the stationary density exp(xi*mu).
void IV02(const std::map<std::string,std::string>&v,const std::string&path){
  const double speed=Num(v,"speed-m-per-s"),d0=Num(v,"d0-per-s"),duration=Num(v,"duration-s");const std::vector<double>xis=Nums(v,"focusing-ratios"),dts=Nums(v,"time-steps-s");const std::uint64_t particles=UInt(v,"particles"),seed=UInt(v,"campaign-seed");
  std::ofstream o((path+".tmp").c_str());o<<std::setprecision(17)<<"ratio,initial,dt_s,time_s,bin,mu_left,mu_right,probability,mean_mu,normalization\n";const double p=Momentum();ConstantPitch diffusion(d0);const int bins=40;
  for(std::size_t x=0;x<xis.size();++x)for(int init=0;init<2;++init)for(std::size_t d=0;d<dts.size();++d){std::vector<double> mu(particles,init?0.95:0.0);const std::uint64_t steps=static_cast<std::uint64_t>(std::llround(duration/dts[d]));
    for(std::uint64_t step=0;step<steps;++step)for(std::uint64_t i=0;i<particles;++i){if(!init&&step==0){KeyedRandomStream r(seed+i,i+1,202,0);mu[i]=2*r.UniformOpen01()-1;}SEP::Transport::FocusedTransportBackground b(-2*d0*xis[x]/speed,0,0,0);KeyedRandomStream r(seed+i,i+1,203,step);const SEP::Transport::FocusedTransportIncrement q=SEP::Transport::AdvanceFocusedTransportDmumu({0,p,mu[i]},b,ProtonMassKg,LightSpeedMPerS,dts[d],diffusion,r,NULL);if(!q.status.ok())throw std::runtime_error(q.status.message);mu[i]=q.state.mu;}
    std::vector<int>h(bins,0);double mean=0;for(std::size_t i=0;i<mu.size();++i){mean+=mu[i];int b=std::max(0,std::min(bins-1,int((mu[i]+1)*.5*bins)));++h[b];}for(int b=0;b<bins;++b)o<<xis[x]<<','<<(init?"beam":"isotropic")<<','<<dts[d]<<','<<duration<<','<<b<<','<<-1+2.0*b/bins<<','<<-1+2.0*(b+1)/bins<<','<<double(h[b])/particles<<','<<mean/particles<<",1\n";
  }Commit(&o,path+".tmp",path);
}

// Smooth manufactured distribution and analytic first derivatives.  IV03
// evaluates the same transport residual with centered numerical derivatives;
// the exact source is the analytic residual, so the difference isolates the
// derivative/source assembly without importing a second production solver.
double Mf(double s,double mu,double y,double t){return std::exp(-.2*t)*(.8+.1*std::sin(2*Pi*s))*(1+.1*mu+.05*mu*mu)*std::exp(-.1*y);}
double Residual(double s,double mu,double y,double t,double h,bool numeric){
  const double U=.2+.05*std::cos(2*Pi*s),div=.03,D=.1*(1-mu*mu),k=.04*(1+.2*std::sin(2*Pi*s));
  if(!numeric){const double f=Mf(s,mu,y,t);const double fs=std::exp(-.2*t)*.1*2*Pi*std::cos(2*Pi*s)*(1+.1*mu+.05*mu*mu)*std::exp(-.1*y);const double fmu=std::exp(-.2*t)*(.8+.1*std::sin(2*Pi*s))*(.1+.1*mu)*std::exp(-.1*y);const double fss=std::exp(-.2*t)*(-.1*4*Pi*Pi*std::sin(2*Pi*s))*(1+.1*mu+.05*mu*mu)*std::exp(-.1*y);return-.2*f+U*fs-div/3*(-.1*f)-D*fmu+k*fss;}
  const double f=Mf(s,mu,y,t);const double fs=(Mf(s+h,mu,y,t)-Mf(s-h,mu,y,t))/(2*h),fmu=(Mf(s,mu+h,y,t)-Mf(s,mu-h,y,t))/(2*h),fy=(Mf(s,mu,y+h,t)-Mf(s,mu,y-h,t))/(2*h),fss=(Mf(s+h,mu,y,t)-2*f+Mf(s-h,mu,y,t))/(h*h);return-.2*f+U*fs-div/3*fy-D*fmu+k*fss;
}
void IV03(const std::map<std::string,std::string>&v,const std::string&path){const std::vector<double>hs=Nums(v,"steps");std::ofstream o((path+".tmp").c_str());o<<std::setprecision(17)<<"step,point,s,mu,log_p,time,exact_residual,numerical_residual,error,positive\n";for(std::size_t q=0;q<hs.size();++q)for(int i=0;i<100;++i){double s=.1+.8*(i+.5)/100,mu=-.8+1.6*(i+.5)/100,y=.2+.6*(i+.5)/100,t=.3;double e=Residual(s,mu,y,t,hs[q],false),n=Residual(s,mu,y,t,hs[q],true);o<<hs[q]<<','<<i<<','<<s<<','<<mu<<','<<y<<','<<t<<','<<e<<','<<n<<','<<n-e<<','<<(Mf(s,mu,y,t)>0)<<'\n';}Commit(&o,path+".tmp",path);}

SEP::Turbulence::CellState Cell(double length,double energy){SEP::Turbulence::CellState c;c.lengthM=length;c.volumeM3=length;c.magneticFieldT=5e-9;c.massDensityKgPerM3=1;c.ePlusJ=energy;return c;}
void IV04(const std::map<std::string,std::string>&v,const std::string&path){const std::vector<double>levels=Nums(v,"resolutions");std::ofstream o((path+".tmp").c_str());o<<std::setprecision(17)<<"motion,resolution,cell,left,right,energy,expected,total,remap_residual,min_length\n";
  // A rigid translation is represented by an unchanged set of arclength
  // measures: absolute origin cancels from the conservative overlap integral.
  // The other maps exercise nonuniform stretching/compression while retaining
  // the same unit-length physical domain and positive segment Jacobians.
  const char*motions[]={"translation","stretch","compress","sinusoidal"};for(int m=0;m<4;++m)for(std::size_t l=0;l<levels.size();++l){int n=int(levels[l]);SEP::Turbulence::State old;old.configuration.advectionEnabled=old.configuration.reflectionEnabled=old.configuration.cascadeEnabled=old.configuration.shockInjectionEnabled=false;old.configuration.coupling=SEP::Turbulence::CouplingPolicy::Disabled;old.provenance="validation:IV04";for(int i=0;i<n;++i)old.cells.push_back(Cell(1.0/n,1.0/n));if(!SEP::Turbulence::InitializeState(&old).ok())throw std::runtime_error("IV04 init");std::vector<SEP::Turbulence::CellState>g;double totalLength=0,minLength=1;for(int i=0;i<n;++i){double factor=1;if(m==1)factor=1+.3*(2.0*(i+.5)/n-1);if(m==2)factor=1-.3*(2.0*(i+.5)/n-1);if(m==3)factor=1+.25*std::sin(2*Pi*(i+.5)/n);double length=factor/n;g.push_back(Cell(length,0));totalLength+=length;minLength=std::min(minLength,length);}for(std::size_t i=0;i<g.size();++i)g[i].lengthM/=totalLength,g[i].volumeM3=g[i].lengthM;SEP::Turbulence::State remap;SEP::Turbulence::EnergyLedger ledger;if(!SEP::Turbulence::RemapConservatively(old,g,&remap,&ledger).ok())throw std::runtime_error("IV04 remap");double left=0,total=0;for(std::size_t i=0;i<remap.cells.size();++i)total+=remap.cells[i].ePlusJ;for(std::size_t i=0;i<remap.cells.size();++i){double right=left+remap.cells[i].lengthM;o<<motions[m]<<','<<n<<','<<i<<','<<left<<','<<right<<','<<remap.cells[i].ePlusJ<<','<<remap.cells[i].lengthM<<','<<total<<','<<ledger.closureResidualJ<<','<<minLength<<'\n';left=right;}}Commit(&o,path+".tmp",path);}

void IV05(const std::map<std::string,std::string>&v,const std::string&path){const std::vector<double>speeds=Nums(v,"shock-speeds-m-per-s"),dts=Nums(v,"time-steps-s");const double start=Num(v,"shock-start-m"),node=Num(v,"node-m"),u1=Num(v,"upstream-speed-m-per-s"),pv=Num(v,"particle-speed-m-per-s");const std::uint64_t particles=UInt(v,"particles"),seed=UInt(v,"campaign-seed");std::ofstream o((path+".tmp").c_str());o<<std::setprecision(17)<<"shock_speed,dt,frame,particle,crossing_time,exact_crossing_time,momentum,cycles,q_sample,q_exact,node_case\n";
  for(std::size_t si=0;si<speeds.size();++si){const double vs=speeds[si],exact=(node-start)/vs,r=4,gain=4*(u1-u1/r)/(3*pv),escape=4*(u1/r)/pv;
    // Use the production shock trajectory for all bracket endpoints.  Two
    // equal-speed knots describe the constant-speed analytical problem while
    // still exercising the same StateAtEpoch path used by production runs.
    SEP::Shock::TrajectoryConfiguration trajectory;trajectory.launchEpochS=0;trajectory.launchRadiusM=start;SEP::Shock::SpeedKnot firstKnot;firstKnot.radiusM=start;firstKnot.speedMPerS=vs;SEP::Shock::SpeedKnot lastKnot;lastKnot.radiusM=std::max(node,start+1.0)+std::max(node-start,1.0);lastKnot.speedMPerS=vs;trajectory.knots.push_back(firstKnot);trajectory.knots.push_back(lastKnot);
    SEP::Shock::TrajectoryState exactState;if(!SEP::Shock::StateAtEpoch(trajectory,exact,&exactState).ok())throw std::runtime_error("IV05 production shock trajectory rejected the exact crossing epoch");
    for(std::size_t di=0;di<dts.size();++di)for(int frame=0;frame<2;++frame)for(std::uint64_t id=0;id<particles;++id){KeyedRandomStream random(seed,id+1,505,si);std::uint64_t cycles=static_cast<std::uint64_t>(std::floor(std::log(random.UniformOpen01())/std::log(1-escape)));double crossing=std::ceil(exact/dts[di])*dts[di];double previous=crossing-dts[di];SEP::Shock::TrajectoryState lower,upper;if(!SEP::Shock::StateAtEpoch(trajectory,previous,&lower).ok()||!SEP::Shock::StateAtEpoch(trajectory,crossing,&upper).ok())throw std::runtime_error("IV05 could not evaluate a production shock crossing bracket");
      // Linear interpolation is exact for this constant-speed trajectory.  It
      // also makes the node-coincidence case explicit instead of depending on
      // a strict inequality whose result changes at a timestep boundary.
      crossing=previous+(node-lower.radiusM)/(upper.radiusM-lower.radiusM)*dts[di];double q=3-std::log(1-escape)/std::log1p(gain);o<<vs<<','<<dts[di]<<','<<(frame?"moving":"stationary")<<','<<id<<','<<crossing<<','<<exact<<','<<std::pow(1+gain,double(cycles))<<','<<cycles<<','<<q<<','<<3*r/(r-1)<<','<<(std::fabs(node-exactState.radiusM)<1e-9)<<'\n';}
  }Commit(&o,path+".tmp",path);}

void IV06(const std::map<std::string,std::string>&v,const std::string&path){const double d0=Num(v,"base-d0-per-s"),growth=Num(v,"growth-per-s"),duration=Num(v,"duration-s"),dt=Num(v,"dt-s"),initialWave=Num(v,"initial-wave-j");const std::uint64_t particles=UInt(v,"particles"),seed=UInt(v,"campaign-seed");std::ofstream o((path+".tmp").c_str());o<<std::setprecision(17)<<"control,step,time,anisotropy,streaming,d0,wave_energy,resonant_bin,total_energy,ledger_residual,wave_ledger_residual\n";const char*names[]={"frozen","one-way","two-way"};
  for(int control=0;control<3;++control){std::vector<double>mu(particles,.9);double reservoir=100,total0=initialWave+reservoir,lastWaveLedgerResidual=0;SEP::Turbulence::State waves;waves.configuration.advectionEnabled=false;waves.configuration.reflectionEnabled=false;waves.configuration.cascadeEnabled=false;waves.configuration.shockInjectionEnabled=false;waves.configuration.coupling=control==0?SEP::Turbulence::CouplingPolicy::Disabled:SEP::Turbulence::CouplingPolicy::StreamingEnergyExchange;waves.provenance="validation:IV06";waves.cells.push_back(Cell(1.0,initialWave));if(!SEP::Turbulence::InitializeState(&waves).ok())throw std::runtime_error("IV06 could not initialize production turbulence state");std::uint64_t steps=static_cast<std::uint64_t>(std::llround(duration/dt));
    for(std::uint64_t step=0;step<=steps;++step){double mean=0;for(std::size_t i=0;i<mu.size();++i)mean+=mu[i];mean/=particles;const double wave=waves.cells[0].ePlusJ;const double localD=d0*(control==2?wave/initialWave:1);o<<names[control]<<','<<step<<','<<step*dt<<','<<mean<<','<<mean*particles<<','<<localD<<','<<wave<<",3,"<<wave+reservoir<<','<<(wave+reservoir-total0)<<','<<lastWaveLedgerResidual<<'\n';if(step==steps)break;
      if(control>0){
        // The analytical streaming instability prescribes dW/dt=2*gamma*A*W.
        // Convert its exact interval integral to the production core's pending
        // particle-exchange amount.  Removing the applied amount (not merely
        // the request) from the particle reservoir makes limiter behavior and
        // total-energy closure auditable in the same evidence row.
        waves.cells[0].pendingParticlePlusJ=wave*std::expm1(2*growth*mean*dt);const SEP::Turbulence::StepResult waveStep=SEP::Turbulence::Advance(&waves,dt);if(!waveStep.status.ok())throw std::runtime_error(waveStep.status.message);reservoir-=waveStep.ledger.particleExchangeJ;lastWaveLedgerResidual=waveStep.ledger.closureResidualJ;
      }
      if(control==2){
        // Scattering observes the wave level from the completed growth stage,
        // making the intended grow-then-scatter Lie-splitting order explicit.
        // The same keyed stream is reused under refinement, so differences are
        // attributable to timestep/operator coupling rather than RNG drift.
        const double scatteringD=d0*waves.cells[0].ePlusJ/initialWave;ConstantPitch diffusion(scatteringD);for(std::uint64_t i=0;i<particles;++i){KeyedRandomStream random(seed+i,i+1,606,step);SEP::Transport::FocusedTransportBackground b(0,0,0,0);const SEP::Transport::FocusedTransportIncrement q=SEP::Transport::AdvanceFocusedTransportDmumu({0,Momentum(),mu[i]},b,ProtonMassKg,LightSpeedMPerS,dt,diffusion,random,NULL);if(!q.status.ok())throw std::runtime_error(q.status.message);mu[i]=q.state.mu;}
      }
    }
  }Commit(&o,path+".tmp",path);}

} // namespace

namespace SEP { namespace Validation {
bool RunIntegratedValidationModel(const std::string& id,const std::vector<std::string>&args,const std::string&path,std::string*error){try{const std::map<std::string,std::string>v=Parse(args);if(id=="IV01")IV01(v,path);else if(id=="IV02")IV02(v,path);else if(id=="IV03")IV03(v,path);else if(id=="IV04")IV04(v,path);else if(id=="IV05")IV05(v,path);else if(id=="IV06")IV06(v,path);else throw std::runtime_error("unsupported IV case");if(error)error->clear();return true;}catch(const std::exception&e){if(error)*error=e.what();return false;}}
}}

#ifdef SRCSEP_INTEGRATED_MODELS_STANDALONE_TEST_HARNESS
int main(int argc,char**argv){if(argc<4)return 2;std::vector<std::string>a;for(int i=3;i<argc;++i)a.push_back(argv[i]);std::string e;if(!SEP::Validation::RunIntegratedValidationModel(argv[1],a,argv[2],&e)){std::cerr<<e<<'\n';return 2;}return 0;}
#endif
