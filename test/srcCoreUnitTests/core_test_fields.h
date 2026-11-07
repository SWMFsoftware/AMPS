#ifndef AMPS_CORE_TEST_FIELDS_H
#define AMPS_CORE_TEST_FIELDS_H

/*
 * Uniform-field provider used only by the configured core-test application.
 *
 * AMPS' non-relativistic Boris mover expresses the Lorentz term as
 *
 *   dv/dt = a + Omega x v,
 *   a = (q/m) E,  Omega = -(q/m) B,
 *
 * with E in V/m, B in tesla, q in coulombs, and m in kilograms.  The test
 * fixture supplies those two vectors directly so analytic constant-field
 * solutions can validate the production pusher without a plasma coupler or a
 * field-solver interpolation error.  The relativistic hook returns the same
 * SI fields.  Both arrays default to zero, preserving the neutral/ballistic
 * cases and making every test responsible for installing and restoring its
 * own state.
 */

#include "pic.h"

namespace AMPS {
namespace CoreTests {

extern double UniformElectricField[3];
extern double UniformMagneticField[3];

void BorisUniformFieldAcceleration(
    double* acceleration, double* rotation, int species, long int particle,
    double* position, double* velocity,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

void RelativisticBorisUniformFields(double* electric_field,
                                    double* magnetic_field);

}  // namespace CoreTests
}  // namespace AMPS

#endif
