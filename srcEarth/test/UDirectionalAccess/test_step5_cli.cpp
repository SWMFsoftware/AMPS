#include "../../util/cutoff_cli.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failures=0;

void Check(bool condition,const std::string& message) {
  if (!condition) {
    ++failures;
    std::cerr << "FAIL: " << message << "\n";
  }
}

EarthUtil::CliOptions Parse(std::vector<std::string> values) {
  std::vector<char*> argv;
  argv.reserve(values.size());
  for (std::string& value:values) argv.push_back(&value[0]);
  return EarthUtil::ParseCli(static_cast<int>(argv.size()),argv.data());
}

} // namespace

int main() {
  const EarthUtil::CliOptions options=Parse({
      "amps",
      "--cutoff-direct-access-adaptive", "T",
      "--cutoff-direct-access-adaptive-max-depth", "12",
      "--cutoff-direct-access-adaptive-guard-depth", "2",
      "--cutoff-direct-access-adaptive-tolerance-gv", "0.0025",
      "--cutoff-direct-access-adaptive-relative-tolerance", "0.0002",
      "--cutoff-direct-access-adaptive-max-samples", "96"});

  Check(options.cutoffDirectAccessAdaptive==1,
        "adaptive boolean CLI value is preserved");
  Check(options.cutoffDirectAccessAdaptiveMaxDepth==12,
        "adaptive max depth CLI value is preserved");
  Check(options.cutoffDirectAccessAdaptiveGuardDepth==2,
        "adaptive guard depth CLI value is preserved");
  Check(std::fabs(options.cutoffDirectAccessAdaptiveTolerance_GV-0.0025)<1.0e-15,
        "adaptive absolute tolerance CLI value is preserved");
  Check(std::fabs(options.cutoffDirectAccessAdaptiveRelativeTolerance-0.0002)<1.0e-15,
        "adaptive relative tolerance CLI value is preserved");
  Check(options.cutoffDirectAccessAdaptiveMaxSamples==96,
        "adaptive sample cap CLI value is preserved");

  const EarthUtil::CliOptions defaults=Parse({"amps"});
  Check(defaults.cutoffDirectAccessAdaptiveTolerance_GV<0.0 &&
        defaults.cutoffDirectAccessAdaptiveRelativeTolerance<0.0 &&
        defaults.cutoffDirectAccessAdaptiveMaxSamples<0,
        "absent Step-5 options retain no-override sentinels");

  const std::string help=EarthUtil::HelpMessage("amps");
  Check(help.find("--cutoff-direct-access-adaptive-tolerance-gv")!=std::string::npos &&
        help.find("--cutoff-direct-access-adaptive-relative-tolerance")!=std::string::npos &&
        help.find("--cutoff-direct-access-adaptive-max-samples")!=std::string::npos,
        "help text exposes all Step-5 convergence controls");

  bool negativeSamplesRejected=false;
  try {
    (void)Parse({"amps","--cutoff-direct-access-adaptive-max-samples","-1"});
  }
  catch (const std::runtime_error&) {
    negativeSamplesRejected=true;
  }
  Check(negativeSamplesRejected,"negative adaptive sample cap is rejected");

  bool nonfiniteToleranceRejected=false;
  try {
    (void)Parse({"amps","--cutoff-direct-access-adaptive-tolerance-gv","nan"});
  }
  catch (const std::runtime_error&) {
    nonfiniteToleranceRejected=true;
  }
  Check(nonfiniteToleranceRejected,"non-finite adaptive tolerance is rejected");

  if (failures!=0) {
    std::cerr << "UDirectionalAccess CLI: " << failures << " failure(s)\n";
    return EXIT_FAILURE;
  }
  std::cout << "UDirectionalAccess CLI: PASS (U-F17)\n";
  return EXIT_SUCCESS;
}
