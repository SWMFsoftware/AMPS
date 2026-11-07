#ifndef AMPS_CORE_TEST_MESH_H
#define AMPS_CORE_TEST_MESH_H

#include "pic.h"

namespace AMPS {
namespace CoreTests {

/* Lazily construct and return the owner leaf of the shared one-rank fixture.
 * The domain is [-1,1]^3 m with a radius-0.2 m sphere at the origin. The
 * fixture is process-local: the native adapter starts every case in a fresh
 * MPI process, while individual cases still restore any callback they replace.
 */
cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* MoverFixtureNode(
    const double* position);
cInternalSphericalData* MoverFixtureSphere();

}  // namespace CoreTests
}  // namespace AMPS

#endif
