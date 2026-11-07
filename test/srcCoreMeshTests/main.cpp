#include "mesh_test_fixture.h"

#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  using namespace AMPS::CoreMeshTests;
  LinkMeshTopologyTests();
  LinkMeshAdaptationTests();
  LinkMeshCoarsenDataTests();
  LinkPICInterpolationTests();
  LinkMeshCutCellTests();

  // Reserve the T3 fixture's application-owned block datum before AMPS fixes
  // cDataBlockAMR's associated-data stride and allocates any mesh block.
  RegisterMeshCoarsenData();

  // Use AMPS' configured MPI, mesh, associated-buffer, and diagnostics setup.
  // The fixture delays actual tree construction until a selected case needs
  // it, which keeps --list side-effect free and allows the outer runner to
  // isolate each mesh mutation in a fresh process.
  PIC::Init_BeforeParser();
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
          std::fprintf(stderr, "invalid --order; use normal, reverse, or shuffle\n");
        MPI_Finalize();
        return 2;
      }
    } else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
      char* end = NULL;
      const unsigned long value = std::strtoul(argv[++i], &end, 10);
      if (end == argv[i] || *end != '\0') {
        if (PIC::ThisThread == 0)
          std::fprintf(stderr, "--seed requires an unsigned integer\n");
        MPI_Finalize();
        return 2;
      }
      shuffle_seed = static_cast<unsigned int>(value);
    } else {
      if (PIC::ThisThread == 0)
        std::fprintf(stderr,
            "usage: %s [--list] [--filter TEXT | --case ID ...] "
            "[--order normal|reverse|shuffle] [--seed UINT]\n", argv[0]);
      MPI_Finalize();
      return 2;
    }
  }

  int local_status = 0;
  if (PIC::ThisThread == 0)
    local_status = AMPS::Test::RunAll(filter, list_only, exact_cases, order,
                                      shuffle_seed);
  int global_status = 0;
  MPI_Allreduce(&local_status, &global_status, 1, MPI_INT, MPI_MAX,
                MPI_GLOBAL_COMMUNICATOR);
  MPI_Finalize();
  return global_status;
}
