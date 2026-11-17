#include "sep_flux_tube_state.h"
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace SEP { namespace FieldLine { namespace FluxTubeState { namespace {
Transport::Status Error(const std::string& text) {
  return Transport::Status::Error(Transport::StatusCode::InvalidParticleState,text);
}
const std::uint64_t Magic=UINT64_C(0x534550464c555831);
void Word(std::uint64_t value,std::string* out) {
  for(unsigned i=0;i<8;++i)out->push_back(static_cast<char>(value>>(8*i)));
}
void Real(double value,std::string* out) {
  static_assert(sizeof(double)==8,"flux state requires binary64 storage");
  std::uint64_t bits;std::memcpy(&bits,&value,8);Word(bits,out);
}
struct Reader {
  const std::string& bytes;std::size_t at=0;
  void Need(std::uint64_t count)const {
    if(count>bytes.size()-at)throw std::runtime_error("truncated flux-owner state");
  }
  std::uint64_t Word() {
    Need(8);std::uint64_t value=0;
    for(unsigned i=0;i<8;++i)value|=static_cast<std::uint64_t>(static_cast<unsigned char>(bytes[at++]))<<(8*i);
    return value;
  }
  double Real(){const auto bits=Word();double value;std::memcpy(&value,&bits,8);return value;}
  std::string Text(){const auto count=Word();Need(count);auto value=bytes.substr(at,count);at+=count;return value;}
};
} // namespace
bool Uninitialized(const FluxTubeGeometryCore::MagneticFluxRecord& record) {
  return record.magneticFluxWb==0.0&&!std::signbit(record.magneticFluxWb)&&
      record.generation==0&&record.provenance.empty();
}
Transport::Status Validate(const State& state,std::size_t lines,double configuredArea) {
  if(state.schema!=1||!std::isfinite(state.referenceAreaM2)||state.referenceAreaM2<=0.0||
      !std::isfinite(configuredArea)||configuredArea<=0.0||state.referenceAreaM2!=configuredArea||
      state.records.size()>lines)
    return Error("flux-owner schema/seed area or fixed physical line namespace is incompatible");
  for(const auto& record:state.records) {
    // Only the exact zero/empty triple means an uninitialized lazy slot. A
    // malformed partially filled record cannot be hidden by deriving new flux.
    if(Uninitialized(record))continue;
    const auto valid=FluxTubeGeometryCore::ValidateMagneticFluxRecord(record);
    if(!valid.ok())return Error("invalid flux-owner record; missing history cannot be reconstructed implicitly");
  }
  return Transport::Status::Ok();
}
Transport::Status Serialize(const State& state,std::string* output) {
  if(!output)return Error("flux-owner serialized output is absent");
  const auto valid=Validate(state,state.records.size(),state.referenceAreaM2);if(!valid.ok())return valid;
  try {
    std::string bytes;Word(Magic,&bytes);Word(state.schema,&bytes);Real(state.referenceAreaM2,&bytes);
    Word(state.records.size(),&bytes);
    for(const auto& record:state.records) {
      Real(record.magneticFluxWb,&bytes);Word(record.generation,&bytes);
      Word(record.provenance.size(),&bytes);bytes+=record.provenance;
    }
    output->swap(bytes);return Transport::Status::Ok();
  }catch(const std::exception& error){return Error(error.what());}
}
Transport::Status Deserialize(const std::string& encoded,State* output) {
  if(!output)return Error("flux-owner decoded output is absent");
  try {
    Reader in{encoded};if(in.Word()!=Magic)return Error("incompatible flux-owner state magic");
    State candidate;candidate.schema=in.Word();candidate.referenceAreaM2=in.Real();
    const auto count=in.Word();
    // Each record has at least three words. Bound allocation by the received
    // payload before constructing a vector, including malicious count fields.
    if(count>(encoded.size()-in.at)/24)return Error("flux-owner record count exceeds encoded payload");
    candidate.records.resize(static_cast<std::size_t>(count));
    for(auto& record:candidate.records) {
      record.magneticFluxWb=in.Real();record.generation=in.Word();record.provenance=in.Text();
    }
    if(in.at!=encoded.size())return Error("trailing flux-owner state bytes");
    const auto valid=Validate(candidate,candidate.records.size(),candidate.referenceAreaM2);
    if(!valid.ok())return valid;
    *output=std::move(candidate);return Transport::Status::Ok();
  }catch(const std::exception& error){return Error(error.what());}
}
} } }
