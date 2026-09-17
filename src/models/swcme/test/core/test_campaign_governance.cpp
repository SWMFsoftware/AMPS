#include "test_framework.hpp"

#include <cstdlib>
#include <iostream>

namespace {

// REP01/EVT01 exercise the real Python campaign manager instead of duplicating
// its JSON rules in C++.  Each command selects one deterministic unittest and
// propagates its process status into the common validation registry.
void run_python_contract(swcme_test::Context& context,const char* test_name) {
  const std::string command=std::string("python3 -m unittest ")+
      "python.test_run_tests.CampaignRunnerTests."+test_name;
  context.expect_true(std::system(command.c_str())==0,
                      std::string(test_name)+" campaign contract");
}

}  // namespace

void test_rep01(swcme_test::Context& context) {
  std::cout << "REP01 fixture and campaign reproducibility\n";
  run_python_contract(context,"test_rep01_manifest_fingerprint_is_canonical");
}

void test_evt01(swcme_test::Context& context) {
  std::cout << "EVT01 campaign schema and completeness\n";
  run_python_contract(context,"test_evt01_schema_rejects_false_completeness");
}

