#ifndef SRCSEP_UTIL_SEP_FLUX_TUBE_STATE_H
#define SRCSEP_UTIL_SEP_FLUX_TUBE_STATE_H
#include "sep_flux_tube_geometry_core.h"
#include <string>
#include <vector>

namespace SEP { namespace FieldLine { namespace FluxTubeState {

// This is the actual owner state, not a reconstruction from the current last
// cell/field value. The reference seed area [m^2] and per-line conserved flux
// [Wb], physical generation and provenance all affect future volume/source
// normalization. Missing table slots are explicit not-yet-initialized state.
struct State {
  std::uint64_t schema=1;
  double referenceAreaM2=0.0;
  std::vector<FluxTubeGeometryCore::MagneticFluxRecord> records;
};
bool Uninitialized(const FluxTubeGeometryCore::MagneticFluxRecord& record);
Transport::Status Validate(const State& state,std::size_t physicalLineCount,
                            double configuredSeedAreaM2);
// Fieldwise binary representation, explicitly little endian. Decode stages
// every string/value and rejects malformed/trailing data before output changes;
// no pointer, object padding or opaque area-profile callback is serialized.
Transport::Status Serialize(const State& state,std::string* output);
Transport::Status Deserialize(const std::string& encoded,State* output);

} } }
#endif
