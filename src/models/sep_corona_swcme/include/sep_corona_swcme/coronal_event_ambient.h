#ifndef SEP_CORONA_SWCME_CORONAL_EVENT_AMBIENT_H
#define SEP_CORONA_SWCME_CORONAL_EVENT_AMBIENT_H

#include "sep_corona_swcme/ambient_state.h"
#include "sep_coronal_cme/background_event.h"

namespace SEP {
namespace CoronaSwcme {

// Bind the existing immutable canonical ambient to one coronal background
// event epoch. The dependency points from this composite to its coronal core;
// no coronal source/header links back to a composite, application or MPI API.
// Acquisition requires the expected actual ambient identity, HCI frame, finite
// covered time, explicit epoch origin and nonzero generation before returning
// an owning const handle. No source plan or regional closure is constructed.
Core::Result<std::shared_ptr<const CoronalCME::AmbientEpoch>> BindAmbientEpoch(
    std::shared_ptr<const AmbientModel> ambient, double timeS, std::uint64_t generation,
    const std::string& epochUtc, const std::string& expectedAmbientIdentity);

}  // namespace CoronaSwcme
}  // namespace SEP
#endif
