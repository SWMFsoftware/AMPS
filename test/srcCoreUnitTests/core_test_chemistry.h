#ifndef AMPS_CORE_TEST_CHEMISTRY_H
#define AMPS_CORE_TEST_CHEMISTRY_H

/*
 * Application-level photochemistry hook for the dedicated T6 daughter build.
 *
 * The generic AMPS whole-mesh chemistry driver owns traversal of each native
 * cell list and delegates one parent record at a time through
 * _PIC_PHOTOLYTIC_REACTIONS__REACTION_PROCESSOR_.  The normal core-test deck
 * deliberately retains AMPS' no-daughter default.  A separate profile defines
 * _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE_=ON and redirects that production
 * hook to T6PhotochemistryDaughterProcessor below.  Keeping the declaration in
 * a physical-model header makes it visible while pic_photolytic_reactions.cpp
 * is compiled, rather than replacing the whole-mesh driver with test code.
 *
 * The callback represents an already-admitted H2O + photon -> OH + H event.
 * Event-time/rate sampling is tested separately; this hook validates daughter
 * insertion, individual statistical weights, kinematics, and intrusive-list
 * ownership.  It is test-application code and must not be enabled in physical
 * applications.
 */

#include "pic.h"

#define _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE__OFF_ 0
#define _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE__ON_ 1

#ifndef _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE_
#define _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE_ \
  _AMPS_CORE_TEST_T6_DAUGHTER_PROCESSOR_MODE__OFF_
#endif

namespace AMPS {
namespace CoreTests {

enum T6DaughterProcessorStatus {
  T6_DAUGHTER_PROCESSOR_OK = 0,
  T6_DAUGHTER_PROCESSOR_INVALID_PARENT_WEIGHT = 1,
  T6_DAUGHTER_PROCESSOR_INVALID_PRODUCT_WEIGHT = 2,
  T6_DAUGHTER_PROCESSOR_INVALID_CHANNEL = 3
};

extern int T6DaughterProcessorStatusCode;
extern int T6DaughterProcessorEventCount;

void ResetT6DaughterProcessorState();

void T6PhotochemistryDaughterProcessor(
    long int parent, long int& first_particle_cell,
    cTreeNodeAMR<PIC::Mesh::cDataBlockAMR>* node);

}  // namespace CoreTests
}  // namespace AMPS

#endif
