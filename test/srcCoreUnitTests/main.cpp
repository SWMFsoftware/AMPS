#include "pic.h"
#include "src/test/native/test_harness.h"

#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>

namespace AMPS {
namespace CoreTests {

// Application objects are collected in a static archive by the standard AMPS
// build.  Referencing one externally visible symbol ensures the linker loads
// particle_buffer.o and therefore runs its static test registrations.
void LinkParticleBufferTests();
void LinkPhotochemistryTests();
void LinkCollisionTests();
void LinkDistributionTests();
void LinkInternalDegreesOfFreedomTests();
void LinkMoverTests();
void LinkT4Tests();
void LinkT5Tests();
void LinkT6PhotochemistryTests();
void LinkT6ReactionTests();
void LinkSpeciesTests();
void LinkUnitTests();
void LinkWeightTimeT2Tests();
void LinkParticleLayoutT2Tests();
void LinkSplitMergeT2Tests();

}  // namespace CoreTests
}  // namespace AMPS

int main(int argc, char** argv) {
  AMPS::CoreTests::LinkParticleBufferTests();
  AMPS::CoreTests::LinkPhotochemistryTests();
  AMPS::CoreTests::LinkCollisionTests();
  AMPS::CoreTests::LinkDistributionTests();
  AMPS::CoreTests::LinkInternalDegreesOfFreedomTests();
  AMPS::CoreTests::LinkMoverTests();
  AMPS::CoreTests::LinkT4Tests();
  AMPS::CoreTests::LinkT5Tests();
  // T6 source files are always archived, but charge-exchange/electron-impact
  // registrations remain compile-profile guarded.  Referencing both anchors
  // makes the linked registry explicit and allows the dedicated reaction
  // profile to activate its optional-module cases without a second main.
  AMPS::CoreTests::LinkT6PhotochemistryTests();
  AMPS::CoreTests::LinkT6ReactionTests();
  AMPS::CoreTests::LinkSpeciesTests();
  AMPS::CoreTests::LinkUnitTests();
  AMPS::CoreTests::LinkWeightTimeT2Tests();
  // Runtime extension offsets must be registered before any production
  // initializer can allocate the particle buffer.  Split/merge only needs an
  // archive anchor, but keeping both T2 anchors here makes every generated
  // profile's linked case inventory explicit and deterministic.
  AMPS::CoreTests::LinkParticleLayoutT2Tests();
  AMPS::CoreTests::LinkSplitMergeT2Tests();

  // Use the production pre-parser initialization rather than constructing a
  // test-only particle layout.  With this application's input deck it installs
  // the LB internal-energy fields, collision sampling state, mover callbacks,
  // and Huebner photochemical tables.  Tests therefore exercise the same
  // compile-time configuration and byte offsets as an AMPS model run.
  PIC::Init_BeforeParser();

#if _PIC_PHOTOLYTIC_REACTIONS_MODE_ == _PIC_PHOTOLYTIC_REACTIONS_MODE_ON_
  // Huebner table initialization and the generic exponential-decay scheduler
  // are separate production layers.  PIC::Init_BeforeParser initializes the
  // former; initialize the latter explicitly so its constant-lifetime API can
  // be tested without building a full model mesh.
  PIC::ChemicalReactions::PhotolyticReactions::Init();
#endif

  // Block allocation in the mover fixture requires all configured sampling
  // offsets to be frozen first.  No mesh is allocated here; the mover suite
  // creates a single uniform block lazily after process-local tests run.
  PIC::Mesh::initCellSamplingDataBuffer();

  std::string filter;
  std::vector<std::string> exact_cases;
  bool list_only = false;
  AMPS::Test::CaseOrder order = AMPS::Test::CASE_ORDER_NORMAL;
  unsigned int shuffle_seed = 0;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--list") == 0) {
      list_only = true;
    } else if (std::strcmp(argv[i], "--filter") == 0 && i + 1 < argc) {
      filter = argv[++i];
    } else if (std::strcmp(argv[i], "--case") == 0 && i + 1 < argc) {
      exact_cases.push_back(argv[++i]);
    } else if (std::strcmp(argv[i], "--order") == 0 && i + 1 < argc) {
      if (!AMPS::Test::ParseCaseOrder(argv[++i], order)) {
        if (PIC::ThisThread == 0)
          std::fprintf(stderr,
                       "invalid --order; use normal, reverse, or shuffle\n");
        MPI_Finalize();
        return 2;
      }
    } else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
      char* end = NULL;
      const unsigned long parsed = std::strtoul(argv[++i], &end, 10);
      if (end == argv[i] || *end != '\0') {
        if (PIC::ThisThread == 0)
          std::fprintf(stderr, "--seed requires an unsigned integer\n");
        MPI_Finalize();
        return 2;
      }
      shuffle_seed = static_cast<unsigned int>(parsed);
    } else {
      if (PIC::ThisThread == 0) {
        std::fprintf(stderr,
            "usage: %s [--list] [--filter SUBSTRING | --case ID ...] "
            "[--order normal|reverse|shuffle] [--seed UINT]\n", argv[0]);
      }
      MPI_Finalize();
      return 2;
    }
  }

  // These first module tests exercise process-local particle state.  Execute
  // them once on rank zero while every rank participates in initialization and
  // the final reduction.  MPI/ownership tests will use dedicated collective
  // cases rather than multiplying process-local results in the log.
  int local_status = 0;
  if (PIC::ThisThread == 0) {
    local_status = AMPS::Test::RunAll(filter, list_only, exact_cases, order,
                                      shuffle_seed);
  }
  // MAX is valid because the harness exit vocabulary is ordered 0=pass,
  // 1=assertion failure, 2=usage/no-selection error.  Any rank-local error
  // therefore dominates success before every process finalizes consistently.
  int global_status = 0;
  MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Finalize();
  return global_status;
}
