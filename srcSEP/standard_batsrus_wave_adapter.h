#ifndef SRCSEP_STANDARD_BATSRUS_WAVE_ADAPTER_H
#define SRCSEP_STANDARD_BATSRUS_WAVE_ADAPTER_H
#include "util/sep_standard_batsrus_waves.h"
namespace SEP { namespace Transport { namespace PICAdapter {
struct ParticleContext;
// The actual production coefficient accessor. Requires an active immutable
// read-only SWMF snapshot matching this context's generation, valid native
// coordinate/segment and allocated two-component imported vertex wave datums.
// It copies raw totals into a typed private view; PIC buffers remain unchanged.
Status ReadStandardBatsrusWaveView(const ParticleContext& context,
    Waves::StandardBatsrus::MagneticBranchView* output);
} } }
#endif
