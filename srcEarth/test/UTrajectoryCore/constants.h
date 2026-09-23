#ifndef _EARTH_UTRAJECTORYCORE_CONSTANTS_H_
#define _EARTH_UTRAJECTORYCORE_CONSTANTS_H_

// GridlessParticleMovers.h normally receives this generated AMPS constant from the
// configured build.  The standalone test supplies the same exact SI value so it can
// compile the production mover implementation without the full AMPS dependency graph.
constexpr double SpeedOfLight=299792458.0;

#endif
