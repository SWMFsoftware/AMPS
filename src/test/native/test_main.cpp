#include "test_harness.h"

#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  // Keep this dependency-light executable's CLI identical to the configured
  // AMPS application.  A test can therefore move across that boundary without
  // inventing another selection or listing protocol.
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
        std::fprintf(stderr, "invalid --order; use normal, reverse, or shuffle\n");
        return 2;
      }
    } else if (std::strcmp(argv[i], "--seed") == 0 && i + 1 < argc) {
      char* end = NULL;
      const unsigned long parsed = std::strtoul(argv[++i], &end, 10);
      if (end == argv[i] || *end != '\0') {
        std::fprintf(stderr, "--seed requires an unsigned integer\n");
        return 2;
      }
      shuffle_seed = static_cast<unsigned int>(parsed);
    } else {
      std::fprintf(stderr,
          "usage: %s [--list] [--filter SUBSTRING | --case ID ...] "
          "[--order normal|reverse|shuffle] [--seed UINT]\n", argv[0]);
      return 2;
    }
  }
  return AMPS::Test::RunAll(filter, list_only, exact_cases, order, shuffle_seed);
}
