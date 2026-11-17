#include "sep_flux_tube_state_validation.h"
#include "sep_flux_tube_state.h"
#include <cmath>
#include <limits>

namespace SEP { namespace Testing { namespace {
namespace F=FieldLine::FluxTubeState;
struct Checks {
  unsigned failures=0;
  void Require(bool value){if(!value)++failures;}
  Result Finish(const char* message)const {
    Result result;result.status=failures?Status::Fail:Status::Pass;
    result.message=failures?"flux-owner value-state assertion failed":message;
    result.configuration={"evidence_level=dependency-light-owner-values-and-fieldwise-codec",
        "units=reference-area-m2;magnetic-flux-Wb;opaque-callbacks-not-serialized",
        "native-checkpoint-and-fresh-process-gates=separate"};
    result.metrics.push_back({"assertion_failures",static_cast<double>(failures),0,"<=","count"});return result;
  }
};
F::State Example() {
  F::State state;state.referenceAreaM2=3.25;state.records.resize(3);
  state.records[0].magneticFluxWb=0.125;state.records[0].generation=17;state.records[0].provenance="explicit source flux";
  state.records[2].magneticFluxWb=0.75;state.records[2].generation=31;state.records[2].provenance="material-line lineage";
  return state;
}
Result ExactRoundTrip() {
  Checks check;const auto original=Example();std::string encoded,again;F::State restored;
  check.Require(F::Serialize(original,&encoded).ok()&&F::Deserialize(encoded,&restored).ok()&&F::Serialize(restored,&again).ok());
  check.Require(encoded==again&&restored.referenceAreaM2==3.25&&restored.records.size()==3&&
      restored.records[0].magneticFluxWb==0.125&&restored.records[0].generation==17&&
      restored.records[0].provenance=="explicit source flux"&&F::Uninitialized(restored.records[1])&&
      restored.records[2].magneticFluxWb==0.75&&restored.records[2].generation==31);
  return check.Finish("fieldwise flux owner preserves actual values, generations, provenance and explicit uninitialized slots exactly");
}
Result PhysicalStateRejection() {
  Checks check;const auto original=Example();
  check.Require(F::Validate(original,3,3.25).ok());
  check.Require(!F::Validate(original,2,3.25).ok()&&!F::Validate(original,3,3.5).ok());
  for(double bad:{0.0,-1.0,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()}) {
    auto state=original;state.records[0].magneticFluxWb=bad;
    check.Require(!F::Validate(state,3,3.25).ok());
    std::string sentinel="unchanged";check.Require(!F::Serialize(state,&sentinel).ok()&&sentinel=="unchanged");
  }
  auto state=original;state.schema=2;check.Require(!F::Validate(state,3,3.25).ok());
  state=original;state.records[1].generation=1;
  check.Require(!F::Validate(state,3,3.25).ok());
  state=original;state.records[1].magneticFluxWb=-0.0;
  check.Require(!F::Validate(state,3,3.25).ok());
  return check.Finish("unknown schemas, missing provenance, malformed flux and partial lazy slots reject without inventing normalization");
}
Result CodecRejection() {
  Checks check;const auto original=Example();std::string good,baseline;check.Require(F::Serialize(original,&good).ok());
  F::State output=original;check.Require(F::Serialize(output,&baseline).ok());
  for(const auto malformed:{good.substr(0,good.size()-1),good+"trailing",std::string("legacy")}) {
    check.Require(!F::Deserialize(malformed,&output).ok());std::string after;
    check.Require(F::Serialize(output,&after).ok()&&after==baseline);
  }
  auto oversized=good;for(unsigned i=24;i<32;++i)oversized[i]=static_cast<char>(255);
  check.Require(!F::Deserialize(oversized,&output).ok());
  return check.Finish("truncation, trailing bytes and oversized counts reject before output replacement or unbounded allocation");
}
Result AreaNormalization() {
  Checks check;auto state=Example();std::string encoded;F::State restored;
  check.Require(F::Serialize(state,&encoded).ok()&&F::Deserialize(encoded,&restored).ok());
  // Independent SI reference: A=Phi/B, using an explicit nondefault owner
  // rather than deriving flux again from the current first-vertex field.
  const double fieldT=0.025;
  const auto area=FieldLine::FluxTubeGeometryCore::AreaFromMagneticFlux(restored.records[0],Units::MagneticFieldT(fieldT));
  check.Require(area.Value()==5.0&&restored.records[0].magneticFluxWb!=state.referenceAreaM2*fieldT);
  F::State lazy;lazy.referenceAreaM2=3.25;
  check.Require(F::Validate(lazy,3,3.25).ok()&&F::Serialize(lazy,&encoded).ok()&&F::Deserialize(encoded,&restored).ok()&&restored.records.empty());
  return check.Finish("explicit conserved flux controls SI area and an observational capture/codec does not materialize absent lazy entries");
}
} // namespace
std::vector<Descriptor> FluxTubeStateDescriptors() {
  const char* names[]={"Exact flux-owner fieldwise continuation","Flux-owner physical/state validation",
      "Atomic flux-owner codec rejection","Restored explicit flux controls area normalization"};
  const TestCallback callbacks[]={ExactRoundTrip,PhysicalStateRejection,CodecRejection,AreaNormalization};
  std::vector<Descriptor> result;
  for(unsigned i=0;i<4;++i) {
    Descriptor descriptor;descriptor.id="FLCP0"+std::to_string(i+1);descriptor.name=names[i];descriptor.description=names[i];
    descriptor.group="flux-checkpoint";descriptor.initialization=InitializationLevel::None;descriptor.runtime=RuntimeClass::Routine;
    descriptor.seedCapability=SeedCapability::Deterministic;descriptor.repeatSupported=true;
    descriptor.seedPolicy="deterministic owner values and malformed binary records; no draws";
    descriptor.supportedBuildModes="dependency-light and native callback factory; owner/codec values only";
    descriptor.stateIsolation="private value records and strings; no native owner mutation";
    descriptor.callback=callbacks[i];result.push_back(descriptor);
  }
  return result;
}
} }
