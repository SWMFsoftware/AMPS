/*
 * Process-wide controls shared by the three canonical field-line movers.
 *
 * Step 14 split these retained configuration values out of the former
 * monolithic mover.cpp.  Implementations now live in parker_mover.cpp,
 * focused_transport_dmumu.cpp, and focused_transport_mfp.cpp; keeping only
 * genuine shared state here makes the production object list mirror that
 * architecture and prevents dead historical movers from being linked.
 */
#include "sep.h"

// Cap the modeled turbulence fluctuation when the explicitly configured guard
// is enabled.  The value is consumed by the coefficient helpers in sep.h.
double SEP::MaxTurbulenceLevel = 0.1;
bool SEP::MaxTurbulenceEnforceLimit = false;

// Preserve the supported lower bound that prevents a modeled parallel mean
// free path from becoming smaller than the local Larmor radius.
bool SEP::LimitMeanFreePath = false;

// Canonical movers consult this switch before applying the common exact
// plasma-frame adiabatic momentum update.
bool SEP::AccountAdiabaticCoolingFlag = true;
