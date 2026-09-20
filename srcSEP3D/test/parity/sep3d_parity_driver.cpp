// V02 3-D producer. This is a distinct executable linked to srcSEP3D's real
// Parker/focused cores. Only the field-aligned projection is compared with the
// 1-D application, as required by the dimensional-reduction contract.
#include "parker_transport.h"
#include "focused_transport.h"

#include <fstream>
#include <iomanip>

int main(int argc,char** argv) {
  if(argc!=2) return 2;
  namespace T=SEP3D::Transport;
  constexpr std::uint64_t n=20000; constexpr double dt=0.5;
  long double sum=0,sum2=0,muSum=0,mu2=0;
  for(std::uint64_t i=0;i<n;++i){
    T::RandomKey key{41,i,0,0,T::RandomPurpose::ParkerParallel};
    T::KeyedRandomStream r(key); T::ParkerLocalState pl;
    pl.bHat={1,0,0}; pl.bulkVelocityMPerS={2,0,0}; pl.kappaParallelM2PerS=4.0;
    const auto p=T::AdvanceParker({{0,0,0},2.0e-20},pl,dt,&r);
    if(!p.status.ok()) return 3; sum+=p.state.positionM.x; sum2+=p.state.positionM.x*p.state.positionM.x;
    key.purpose=T::RandomPurpose::FocusedPitch; T::KeyedRandomStream q(key);
    T::FocusedLocalState fl; fl.bHat={1,0,0}; fl.dMuMuPerS=0.25;
    T::FocusedParticleState fs; fs.momentumKgMPerS=2.0e-20; fs.mu=-1.0+2.0*(i+0.5)/n;
    const auto f=T::AdvanceFocused(fs,fl,SEP3D::Core::Const::m_p,dt,&q);
    if(!f.status.ok()) return 4; muSum+=f.state.mu; mu2+=f.state.mu*f.state.mu;
  }
  const long double mean=sum/n,meanMu=muSum/n;
  std::ofstream out(argv[1]); if(!out) return 5;
  out<<std::setprecision(17)<<"{\n\"schema\":\"sep-cross-model-v1\",\n"
     <<"\"application\":\"srcSEP3D\",\n\"parker_mean_m\":"<<(double)mean<<",\n"
     <<"\"parker_variance_m2\":"<<(double)(sum2/n-mean*mean)<<",\n"
     <<"\"focused_mu_mean\":"<<(double)meanMu<<",\n"
     <<"\"focused_mu_variance\":"<<(double)(mu2/n-meanMu*meanMu)<<"\n}\n";
  return out ? 0 : 6;
}
