#ifndef SEP3D_OUTPUT_REDUCED_FRONT_OUTPUT_H
#define SEP3D_OUTPUT_REDUCED_FRONT_OUTPUT_H

// Dependency-light presentation of the shared reduced front.  Physics and
// classification stay in sep_corona_swcme; this layer only maps immutable SI
// records to a Tecplot surface that can be overlaid on AMPS volume output.
#include "provider.h"

#include <string>

namespace SEP3D { namespace Output {

// Serialize the provider-owned surface as FETRIANGLE with nodal HCI geometry
// and cell-centred physical state. Each element is the same exact-curved-area
// patch classified by the shared provider; the chord triangulation is not a
// second shock model and its smaller planar area is diagnostic only.
// Tecplot cannot reliably ingest NaN on every supported installation, so a
// non-shock record uses a finite *display encoding*: Mach=0, both compression
// ratios=1, and the downstream display columns repeat the upstream ambient
// primitive.  shock_accepted=0 and downstream_valid=0 remain authoritative;
// the repeated values are not an RH state or downstream CME plasma.  Separate
// theta_Bn_valid and magnetic_compression_valid columns protect the rarer
// magnetic-null cases from acquiring an apparently valid zero/one diagnostic.
// The exact support rim remains the sole open boundary, while one shared apex
// vertex and periodic azimuth connectivity eliminate the former pole hole.
std::string SerializeReducedFrontTecplot(
    const SEP::CoronaSwcme::ShockFront::Epoch& epoch,
    const SEP::CoronaSwcme::ShockFront::Configuration& configuration);

} } // namespace SEP3D::Output

#endif
