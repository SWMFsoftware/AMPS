#include "sep_standard_batsrus_wave_validation.h"
#include "sep_standard_batsrus_waves.h"
#include <cmath>
#include <cstring>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace W=Waves::StandardBatsrus;
namespace T=Transport;
constexpr double ReferenceMu0=1.256637062e-6;
struct Checks {
  unsigned failures=0;
  void Require(bool value){if(!value)++failures;}
  void Near(double actual,double expected){Require(std::isfinite(actual)&&std::abs(actual-expected)<=8.0*std::numeric_limits<double>::epsilon()*std::abs(expected));}
  Result Finish(const char* message)const {
    Result result;result.status=failures?Status::Fail:Status::Pass;
    result.message=failures?"standard BATSRUS wave-view assertion failed":message;
    result.configuration={"evidence_level=dependency-light-standard-producer-contract",
        "raw=total-Alfven-J/m3;local-view=magnetic-J/m3;factor=1/2",
        "Wave1=along-B;Wave2=opposite-B;no-polarity-or-line-order-swap;L/R=unresolved",
        "full-SWMF-import-gradient-composition-provider-qualification=false"};
    result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0,"<=","count"});return result;
  }
};
bool Same(const W::MagneticBranchView& a,const W::MagneticBranchView& b) {
  return a.plusMagneticJPerM3==b.plusMagneticJPerM3&&a.minusMagneticJPerM3==b.minusMagneticJPerM3&&
      a.plusVarianceT2==b.plusVarianceT2&&a.minusVarianceT2==b.minusVarianceT2&&a.totalVarianceT2==b.totalVarianceT2&&
      a.generation==b.generation&&a.checksum==b.checksum&&a.branchDirectionsKnown==b.branchDirectionsKnown&&
      a.circularPolarizationsKnown==b.circularPolarizationsKnown;
}
W::MagneticBranchView Sentinel(){W::MagneticBranchView view;view.plusMagneticJPerM3=17;view.checksum=29;return view;}
Result Conversion() {
  Checks check;W::RawTotalAlfvenEndpoint first{{{8e-9,2e-9}}},last{{{4e-9,10e-9}}};
  const auto beforeFirst=first,beforeLast=last;W::MagneticBranchView view;
  check.Require(W::InterpolateMagneticView(first,last,0.25,3,&view).ok());
  check.Near(view.plusMagneticJPerM3,3.5e-9);check.Near(view.minusMagneticJPerM3,2.0e-9);
  // Independent producer reference: total Alfven energy contains equal
  // kinetic/magnetic parts, so deltaB^2=mu0*(raw Wave1+raw Wave2).
  check.Near(view.plusVarianceT2,ReferenceMu0*7e-9);check.Near(view.minusVarianceT2,ReferenceMu0*4e-9);
  check.Near(view.totalVarianceT2,ReferenceMu0*11e-9);
  check.Require(view.totalVarianceT2!=2.0*ReferenceMu0*11e-9);
  const auto local=T::CoefficientPhysics::MagneticVarianceFromEnergyDensity(7e-9,4e-9);
  check.Require(local.status.ok());check.Near(local.totalT2,2.0*ReferenceMu0*11e-9);
  check.Require(std::memcmp(&first,&beforeFirst,sizeof(first))==0&&std::memcmp(&last,&beforeLast,sizeof(last))==0);
  return check.Finish("standard raw totals convert to half-energy magnetic views while local magnetic normalization remains unchanged");
}
Result Ordering() {
  Checks check;W::RawTotalAlfvenEndpoint first{{{9e-9,1e-9}}},last{{{3e-9,5e-9}}};
  W::MagneticBranchView direct,reversed;
  check.Require(W::InterpolateMagneticView(first,last,0.25,7,&direct).ok()&&
      W::InterpolateMagneticView(last,first,0.75,7,&reversed).ok());
  check.Require(direct.plusVarianceT2==reversed.plusVarianceT2&&direct.minusVarianceT2==reversed.minusVarianceT2&&
      direct.plusVarianceT2>direct.minusVarianceT2&&direct.branchDirectionsKnown&&reversed.branchDirectionsKnown);
  // No B-sign argument exists at this conversion boundary: its branch labels
  // refer to actual B, never permanently radial or geometric directions.
  check.Require(std::string(W::ContractIdentity()).find("Wave1=along-B;Wave2=opposite-B")!=std::string::npos);
  return check.Finish("reversing endpoint order and fraction preserves physical B-directed branch amplitudes without a branch swap");
}
Result EndpointRejection() {
  Checks check;const W::RawTotalAlfvenEndpoint valid{{{1e-8,2e-8}}};
  for(double value:{-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN(),
      std::numeric_limits<double>::denorm_min(),1e-320}) {
    for(unsigned branch=0;branch<2;++branch)for(unsigned endpoint=0;endpoint<2;++endpoint) {
      auto first=valid,last=valid;(endpoint?last:first).waveJPerM3[branch]=value;
      const auto beforeFirst=first,beforeLast=last;const auto sentinel=Sentinel();auto output=sentinel;
      // The bad endpoint is assigned exactly zero interpolation weight. Its
      // invalid physical state must not be hidden by a valid weighted mean.
      const double fraction=endpoint?0.0:1.0;
      check.Require(!W::InterpolateMagneticView(first,last,fraction,1,&output).ok()&&Same(output,sentinel));
      check.Require(std::memcmp(&first,&beforeFirst,sizeof(first))==0&&std::memcmp(&last,&beforeLast,sizeof(last))==0);
    }
  }
  return check.Finish("invalid raw endpoints and positive half-energy/variance underflow reject atomically even at zero interpolation weight");
}
Result NumericLimits() {
  Checks check;const W::RawTotalAlfvenEndpoint first{{{1e-8,2e-8}}};const auto sentinel=Sentinel();
  for(double fraction:{-1.0,2.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    auto output=sentinel;check.Require(!W::InterpolateMagneticView(first,first,fraction,1,&output).ok()&&Same(output,sentinel));
  }
  auto output=sentinel;check.Require(!W::InterpolateMagneticView(first,first,0.5,0,&output).ok()&&Same(output,sentinel));
  check.Require(!W::InterpolateMagneticView(first,first,0.5,1,nullptr).ok());
  W::RawTotalAlfvenEndpoint empty,positive{{{1e-8,0}}};
  output=sentinel;check.Require(W::InterpolateMagneticView(empty,positive,std::numeric_limits<double>::denorm_min(),1,&output).code==T::StatusCode::StepUnderflow&&Same(output,sentinel));
  const double maximum=std::numeric_limits<double>::max();W::RawTotalAlfvenEndpoint largest{{{maximum,maximum}}};
  check.Require(W::InterpolateMagneticView(largest,largest,0.5,1,&output).ok());
  check.Near(output.totalVarianceT2,static_cast<double>(2.0L*ReferenceMu0*maximum));
  return check.Finish("closed-unit sampling and generation contracts reject extrapolation; wider convex conversion avoids spurious raw-sum overflow");
}
Result Identity() {
  Checks check;W::RawTotalAlfvenEndpoint first{{{4e-9,2e-9}}},last=first;W::MagneticBranchView a,b,c;
  check.Require(W::InterpolateMagneticView(first,last,0.0,3,&a).ok()&&W::InterpolateMagneticView(first,last,0.0,3,&b).ok()&&Same(a,b));
  check.Require(W::InterpolateMagneticView(first,last,0.0,4,&c).ok()&&c.checksum!=a.checksum);
  last.waveJPerM3[0]*=2.0;check.Require(W::InterpolateMagneticView(first,last,0.0,3,&c).ok()&&c.checksum!=a.checksum&&c.plusVarianceT2==a.plusVarianceT2);
  check.Require(a.generation==3&&a.checksum!=0&&a.branchDirectionsKnown&&!a.circularPolarizationsKnown);
  check.Require(std::string(W::Representation()).find("total-to-magnetic-half")!=std::string::npos&&
      std::string(W::Representation()).find("L/R=unspecified")!=std::string::npos);
  W::RawTotalAlfvenEndpoint zero;check.Require(W::InterpolateMagneticView(zero,zero,0.5,3,&c).ok()&&c.totalVarianceT2==0.0&&!c.circularPolarizationsKnown);
  return check.Finish("contract/generation/complete endpoint identity are pinned; no circular polarization is fabricated from two integrated waves");
}
} // namespace
std::vector<Descriptor> StandardBatsrusWaveDescriptors() {
  const char* names[]={"Standard BATSRUS total-to-magnetic energy conversion","B-directed branch labels under endpoint reversal",
      "Atomic raw endpoint and underflow rejection","Imported wave interpolation numerical domain","Imported wave contract identity and absent polarizations"};
  const TestCallback callbacks[]={Conversion,Ordering,EndpointRejection,NumericLimits,Identity};std::vector<Descriptor> result;
  for(unsigned i=0;i<5;++i) {
    Descriptor descriptor;descriptor.id="SBW0"+std::to_string(i+1);descriptor.name=names[i];descriptor.description=names[i];
    descriptor.group="standard-batsrus-waves";descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Routine;
    descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;descriptor.seedPolicy="deterministic producer-contract values; no draws";
    descriptor.stateIsolation="private typed endpoints and output values; no native imported array mutation";
    descriptor.supportedBuildModes="dependency-light/native identical callbacks; full imported provider remains unqualified";
    descriptor.callback=callbacks[i];result.push_back(descriptor);
  }
  return result;
}
} }
