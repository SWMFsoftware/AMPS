#ifndef _GENERIC_COMET_RADIATIVE_COOLING_H_
#define _GENERIC_COMET_RADIATIVE_COOLING_H_

#include <string>

#if __has_include("models/cometary_physics/RadiativeCooling.h")
#include "models/cometary_physics/RadiativeCooling.h"
#elif __has_include("../models/cometary_physics/RadiativeCooling.h")
#include "../models/cometary_physics/RadiativeCooling.h"
#else
#error "Cannot locate models/cometary_physics/RadiativeCooling.h"
#endif

namespace Comet {
namespace RadiativeCooling {

// One input-selected model code is stored for every compiled species.  Zero
// disables cooling for that species.  The generic input parser resolves names
// to these indices and rejects use of a species-specific law on another
// species.
extern int SpeciesModel[];

bool Validate(std::string *error=NULL);

// Apply the selected cooling laws after one AMPS particle iteration and
// before the next.  Macroscopic density and rotational temperature come from
// the most recently completed sampling window; the active collection buffer
// is never read.
void Apply();

}
}

#endif
