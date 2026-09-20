// V02 1-D producer. This executable links the real srcSEP Parker and focused
// cores; it does not call a shared parity helper or any srcSEP3D routine.
#include "sep_parker_core.h"
#include "sep_focused_transport_core.h"

#include <cmath>
#include <fstream>
#include <iomanip>

namespace {
class Kappa final : public SEP::Transport::SpatialDiffusionProvider {
 public:
  SEP::Transport::SpatialDiffusionSample Evaluate(double, double) const override {
    SEP::Transport::SpatialDiffusionSample s; s.status=SEP::Transport::Status::Ok();
    s.valueState=SEP::Transport::CoefficientPhysics::ValueState::Finite;
    s.kappaParallelM2PerS=4.0; s.provenance="V02 constant kappa"; return s;
  }
};
class Dmumu final : public SEP::Transport::PitchAngleDiffusionProvider {
 public:
  SEP::Transport::PitchAngleDiffusionSample Evaluate(double,double,double) const override {
    SEP::Transport::PitchAngleDiffusionSample s; s.status=SEP::Transport::Status::Ok();
    s.valueState=SEP::Transport::CoefficientPhysics::ValueState::Finite;
    s.dMuMuPerS=0.25; s.dDmuMuDmuPerS=0.0; s.provenance="V02 constant Dmumu";
    s.turbulenceStateIdentity="v02"; return s;
  }
};
}

int main(int argc,char** argv) {
  if (argc!=2) return 2;
  constexpr std::uint64_t n=20000; constexpr double dt=0.5;
  long double sum=0,sum2=0,muSum=0,mu2=0; Kappa k; Dmumu d;
  for(std::uint64_t i=0;i<n;++i){
    SEP::Transport::KeyedRandomStream r(41,i,1,0);
    const auto p=SEP::Transport::AdvanceParker({0,2.0e-20},{2.0,0.0},1.0,dt,k,r);
    if(!p.status.ok()) return 3; sum+=p.state.arcLengthM; sum2+=p.state.arcLengthM*p.state.arcLengthM;
    SEP::Transport::KeyedRandomStream q(41,i,2,0);
    SEP::Transport::FocusedTransportBackground bg; bg.equationMode=SEP::Transport::FocusedEquationMode::FullGyrotropic;
    const double initialMu=-1.0+2.0*(i+0.5)/n;
    const auto f=SEP::Transport::AdvanceFocusedTransportDmumu(
        {0,2.0e-20,initialMu},bg,1.67262192369e-27,2.99792458e8,dt,d,q,nullptr);
    if(!f.status.ok()) return 4; muSum+=f.state.mu; mu2+=f.state.mu*f.state.mu;
  }
  const long double mean=sum/n, meanMu=muSum/n;
  std::ofstream out(argv[1]); if(!out) return 5;
  out<<std::setprecision(17)<<"{\n\"schema\":\"sep-cross-model-v1\",\n"
     <<"\"application\":\"srcSEP\",\n\"parker_mean_m\":"<<(double)mean<<",\n"
     <<"\"parker_variance_m2\":"<<(double)(sum2/n-mean*mean)<<",\n"
     <<"\"focused_mu_mean\":"<<(double)meanMu<<",\n"
     <<"\"focused_mu_variance\":"<<(double)(mu2/n-meanMu*meanMu)<<"\n}\n";
  return out ? 0 : 6;
}
