// Compatibility interface for the unchanged Moon source-level API.
// The actual implementation and data semantics now belong to the shared,
// PIC/MPI-independent chemical database module. No rate, source selection,
// domain, unit, stochastic state or physical activation is changed here.
// The generated tree places this header in build/main/species; standalone
// source adapters compile the original srcMoon/species copy. Both resolve
// the same maintained shared header without an absolute installation path.
#ifndef SRCMOON_SHARED_PHOTOCHEMISTRY_COMPAT_H
#define SRCMOON_SHARED_PHOTOCHEMISTRY_COMPAT_H

#if __has_include("../../models/chemical_database/Photochemistry.h")
#include "../../models/chemical_database/Photochemistry.h"
#else
#include "../../src/models/chemical_database/Photochemistry.h"
#endif

namespace Moon { namespace Species {
namespace ReactionData = ::AMPS::ChemicalDatabase::ReactionData;
namespace Photochemistry = ::AMPS::ChemicalDatabase::Photochemistry;
} }
#endif
