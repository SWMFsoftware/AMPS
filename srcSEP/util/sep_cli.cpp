#include "sep_cli.h"
#include "sep_configuration_matrix.h"
#ifndef SEP_CLI_PARSE_ONLY
#include "../sep.h"
#endif

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace SEP {
namespace Util {
namespace CLI {

namespace {

// Convert a string to lower case for case-insensitive command-line values.  The
// unsigned-char cast is intentional; std::tolower has undefined behavior for
// negative signed char values.
std::string ToLower(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return value;
}

// Parse a boolean switch value.  Accepting several common spellings makes the
// CLI convenient for shell scripts and interactive use while still rejecting
// ambiguous inputs with a clear error message.
bool ParseBoolValue(const std::string& raw_value, bool& value) {
  const std::string v = ToLower(raw_value);

  if (v == "on" || v == "true" || v == "yes" || v == "1" || v == "enable" || v == "enabled") {
    value = true;
    return true;
  }

  if (v == "off" || v == "false" || v == "no" || v == "0" || v == "disable" || v == "disabled") {
    value = false;
    return true;
  }

  return false;
}


// Parse the turbulence model selector.  This deliberately accepts several
// spellings used in notes/scripts so users do not need to remember one exact
// token.  The model name controls only the representation of wave energy; the
// independent switches --coupling, --cascade, and --reflection still turn the
// corresponding physics on/off.
bool ParseTurbulenceModelValue(const std::string& raw_value, Options::TurbulenceModel& value) {
  const std::string v = ToLower(raw_value);

  if (v == "integrated" || v == "legacy" || v == "old" || v == "branch-integrated") {
    value = Options::TurbulenceModel::Integrated;
    return true;
  }

  if (v == "wave-number-resolved" || v == "wavenumber-resolved" ||
      v == "k-resolved" || v == "k" || v == "spectral" || v == "new") {
    value = Options::TurbulenceModel::WaveNumberResolved;
    return true;
  }

  return false;
}

// Split options of the form "--option=value".  If there is no '=' character,
// option_name receives the complete argument and option_value is left empty.
void SplitOption(const std::string& arg, std::string& option_name, std::string& option_value) {
  const std::string::size_type pos = arg.find('=');
  if (pos == std::string::npos) {
    option_name = arg;
    option_value.clear();
  }
  else {
    option_name = arg.substr(0, pos);
    option_value = arg.substr(pos + 1);
  }
}

// Read an option value either from "--option=value" or from the following
// command-line token in "--option value" form.  The index i is advanced only
// when the following token is consumed.
bool GetOptionValue(int argc, char** argv, int& i,
                    const std::string& option_name,
                    const std::string& value_from_equals,
                    std::string& value,
                    std::ostream& err) {
  if (!value_from_equals.empty()) {
    value = value_from_equals;
    return true;
  }

  if (i + 1 >= argc) {
    err << "ERROR: option '" << option_name << "' requires a value.\n";
    return false;
  }

  value = argv[++i];
  return true;
}

// Shared implementation for switches that set one boolean model flag.  The
// function handles both --option=value and --option value syntaxes and reports a
// model-specific error if the value cannot be interpreted as boolean.
bool ParseBooleanOption(int argc, char** argv, int& i,
                        const std::string& option_name,
                        const std::string& value_from_equals,
                        const char* description,
                        bool& destination,
                        std::ostream& err) {
  std::string raw_value;
  if (!GetOptionValue(argc, argv, i, option_name, value_from_equals, raw_value, err)) return false;

  bool parsed_value = false;
  if (!ParseBoolValue(raw_value, parsed_value)) {
    err << "ERROR: invalid value '" << raw_value << "' for " << description
        << ". Use on/off, true/false, yes/no, or 1/0.\n";
    return false;
  }

  destination = parsed_value;
  return true;
}

// Parse a non-negative integer option.  This is used for diagnostic-output
// cadence controls where zero has a useful meaning: disable that diagnostic.
// strtol is used instead of atoi so malformed values such as "10abc" are
// rejected rather than silently truncated.
bool ParseNonNegativeIntegerOption(int argc, char** argv, int& i,
                                   const std::string& option_name,
                                   const std::string& value_from_equals,
                                   const char* description,
                                   int& destination,
                                   std::ostream& err) {
  std::string raw_value;
  if (!GetOptionValue(argc, argv, i, option_name, value_from_equals, raw_value, err)) return false;

  errno = 0;
  char* end_ptr = nullptr;
  const long parsed_value = std::strtol(raw_value.c_str(), &end_ptr, 10);

  if (errno != 0 || end_ptr == raw_value.c_str() || (end_ptr && *end_ptr != '\0') ||
      parsed_value < 0 || parsed_value > std::numeric_limits<int>::max()) {
    err << "ERROR: invalid value '" << raw_value << "' for " << description
        << ". Use a non-negative integer; 0 disables the diagnostic.\n";
    return false;
  }

  destination = static_cast<int>(parsed_value);
  return true;
}

// Parse a strictly positive integer option.  This is separate from the
// non-negative diagnostic-cadence parser above because injection-particle counts
// are used in denominators in field_line.cpp when computing particle-weight
// correction factors.  Accepting zero would therefore create either a singular
// weight or a run with a silently disabled injection source.
bool ParsePositiveIntegerOption(int argc, char** argv, int& i,
                                const std::string& option_name,
                                const std::string& value_from_equals,
                                const char* description,
                                int& destination,
                                std::ostream& err) {
  std::string raw_value;
  if (!GetOptionValue(argc, argv, i, option_name, value_from_equals, raw_value, err)) return false;

  errno = 0;
  char* end_ptr = nullptr;
  const long parsed_value = std::strtol(raw_value.c_str(), &end_ptr, 10);

  if (errno != 0 || end_ptr == raw_value.c_str() || (end_ptr && *end_ptr != '\0') ||
      parsed_value <= 0 || parsed_value > std::numeric_limits<int>::max()) {
    err << "ERROR: invalid value '" << raw_value << "' for " << description
        << ". Use a positive integer greater than zero.\n";
    return false;
  }

  destination = static_cast<int>(parsed_value);
  return true;
}

// Strict floating-point parsing is shared by the Step 11 SI-valued controls.
// The non_negative flag is used for coefficients and boundary data; strictly
// positive lengths/wave numbers reject zero as a physically singular input.
bool ParseDoubleOption(int argc, char** argv, int& i,
                       const std::string& option_name,
                       const std::string& value_from_equals,
                       const char* description, bool non_negative,
                       double& destination, std::ostream& err) {
  std::string raw_value;
  if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                      raw_value, err)) return false;
  errno = 0;
  char* end_ptr = nullptr;
  const double parsed = std::strtod(raw_value.c_str(), &end_ptr);
  const bool range_ok = non_negative ? parsed >= 0.0 : parsed > 0.0;
  if (errno != 0 || end_ptr == raw_value.c_str() ||
      (end_ptr && *end_ptr != '\0') || !std::isfinite(parsed) || !range_ok) {
    err << "ERROR: invalid value '" << raw_value << "' for " << description
        << ". Use a finite " << (non_negative ? "non-negative" : "positive")
        << " number.\n";
    return false;
  }
  destination = parsed;
  return true;
}

// Divergence is signed, so it cannot use the positive/non-negative helper.
// This parser still rejects suffixes, overflow, NaN, and infinity.
bool ParseFiniteDoubleOption(int argc, char** argv, int& i,
                             const std::string& option_name,
                             const std::string& value_from_equals,
                             const char* description, double& destination,
                             std::ostream& err) {
  std::string raw_value;
  if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                      raw_value, err)) return false;
  errno = 0;
  char* end_ptr = nullptr;
  const double parsed = std::strtod(raw_value.c_str(), &end_ptr);
  if (errno != 0 || end_ptr == raw_value.c_str() ||
      (end_ptr && *end_ptr != '\0') || !std::isfinite(parsed)) {
    err << "ERROR: invalid value '" << raw_value << "' for " << description
        << ". Use a finite number.\n";
    return false;
  }
  destination = parsed;
  return true;
}

} // anonymous namespace

void PrintHelp(const char* program_name, std::ostream& out) {
  const char* exe = (program_name && program_name[0] != '\0') ? program_name : "amps";

  out
      << "SEP standalone driver command-line options\n"
      << "\n"
      << "Usage:\n"
      << "  " << exe << " [options]\n"
      << "\n"
      << "Runtime initialization:\n"
      << "  --input <path>              Parse the SI Parker/domain/mesh input before\n"
      << "                               AMPS initialization. If omitted, the legacy\n"
      << "                               hard-coded srcSEP mesh is preserved.\n"
      << "  --input=<path>              Equivalent equals-sign form.\n"
      << "\n"
      << "Turbulence physics switches:\n"
      << "  --coupling <on|off>          Enable/disable SEP particle coupling to the\n"
      << "                               Alfven turbulence wave energy.\n"
      << "  --coupling-mode <on|off>     Alias for --coupling.\n"
      << "  --particle-coupling <on|off> Alias for --coupling.\n"
      << "  --no-coupling                Shortcut for --coupling off.\n"
      << "\n"
      << "  --cascade <on|off>           Enable/disable nonlinear turbulence cascade.\n"
      << "  --cascase <on|off>           Accepted alias for --cascade.\n"
      << "  --no-cascade                 Shortcut for --cascade off.\n"
      << "  --no-cascase                 Accepted alias for --no-cascade.\n"
      << "\n"
      << "  --reflection <on|off>        Enable/disable large-scale-gradient reflection\n"
      << "                               between W+ and W- waves.\n"
      << "  --no-reflection              Shortcut for --reflection off.\n"
      << "\n"
      << "Turbulence representation:\n"
      << "  --turbulence-model <integrated|wave-number-resolved>\n"
      << "                               Select the turbulence-energy representation.\n"
      << "                               integrated: legacy/default E+,E- only.\n"
      << "                               wave-number-resolved: store and advect E±(k_j);\n"
      << "                               particle coupling modifies the resonant k-bin.\n"
      << "  --turbulence-model=...      Same option using --option=value syntax.\n"
      << "  --wave-number-resolved      Shortcut for --turbulence-model wave-number-resolved.\n"
      << "  --integrated-turbulence     Shortcut for --turbulence-model integrated.\n"
      << "  --turbulence-source <name>  prescribed, self-consistent-integrated,\n"
      << "                               self-consistent-spectral, swmf-read-only, or\n"
      << "                               swmf-initial-then-local. Source ownership is\n"
      << "                               independent of the particle mover.\n"
      << "  --turbulence-coupling-policy <disabled|streaming-energy-exchange>\n"
      << "  --turbulence-inner-boundary <specified-incoming-energy|specified-incoming-flux|\n"
      << "                               transparent-outflow|fixed-reservoir>\n"
      << "  --turbulence-outer-boundary <policy>  Set the outer boundary object.\n"
      << "  --turbulence-inner-value <SI> Boundary energy [J] or inward flux [W].\n"
      << "  --turbulence-outer-value <SI> Boundary energy [J] or inward flux [W].\n"
      << "  --reflection-coefficient <C> Non-negative dimensionless reflection factor.\n"
      << "  --cascade-coefficient <C>  Non-negative dimensionless cascade factor.\n"
      << "  --turbulence-advection <on|off> Enable/disable conservative wave advection.\n"
      << "  --shock-injection <on|off> Enable/disable the shock turbulence source.\n"
      << "  --turbulence-cfl <value>   Set the finite-volume CFL safety in (0,1].\n"
      << "  --turbulence-operator-safety <value> Shared local-rate accuracy safety in (0,1].\n"
      << "  --turbulence-max-source-fraction <value> Maximum source/energy ratio per step.\n"
      << "  --turbulence-max-cascade-fraction <value> Maximum cascade fraction per substep.\n"
      << "  --turbulence-min-substep <s> Reject smaller combined operator steps.\n"
      << "  --turbulence-max-substeps <N> Reject plans requiring more than N stages.\n"
      << "  --turbulence-correlation-length <m> Positive perpendicular scale [m].\n"
      << "  --spectral-k-min <1/m> --spectral-k-max <1/m> --spectral-bins <N>\n"
      << "                               Configure the logarithmic spectral authority.\n"
      << "  --turbulence-conservation-tolerance <relative>\n"
      << "                               Set the non-negative signed-ledger tolerance.\n"
      << "\n"
      << "Particle mover selection:\n"
      << "  --particle-mover <name>     Select one production field-line mover.\n"
      << "                               Choices: parker, fte-dmumu, fte-mfp.\n"
      << "                               Default: fte-dmumu.\n"
      << "  --particle-mover=...        Same option using --option=value syntax.\n"
      << "  --mover <name>              Alias for --particle-mover.\n"
      << "  --sep-mover <name>          Alias for --particle-mover.\n"
      << "  --list-movers               List the same three canonical movers and\n"
      << "                               capability metadata, then exit.\n"
      << "\n"
      << "Transport coefficient registry:\n"
      << "  --coefficient-source <name> Source authority: prescribed, self-consistent,\n"
      << "                               or swmf. Default: prescribed.\n"
      << "  --spatial-diffusion-provider <name>\n"
      << "                               kappa provider: from-dmumu or from-mfp.\n"
      << "  --pitch-angle-diffusion-provider <name>\n"
      << "                               Dmumu provider: configured, constant,\n"
      << "                               jokipii-1966, or florinskiy.\n"
      << "  --mean-free-path-provider <name>\n"
      << "                               lambda provider: qlt, qlt1, tenishev-2005,\n"
      << "                               chen-2024, or from-spatial.\n"
      << "  --invalid-coefficient-policy <fail|ballistic>\n"
      << "                               Reject invalid coefficients, or use the exact\n"
      << "                               lambda=+infinity zero-event-rate limit.\n"
      << "                               Incompatible conversion cycles and source/\n"
      << "                               provider combinations are rejected at parse.\n"
      << "  --resonance-gap-policy <reject|ballistic>\n"
      << "                               Declare the physical 90-degree-gap behavior.\n"
      << "  --turbulence-amplitude-policy <reject|limit-to-mean-field>\n"
      << "                               Reject deltaB/B>1 or ledger an explicit limit.\n"
      << "  --constant-dmumu <1/s>      Non-negative constant provider value.\n"
      << "  --prescribed-delta-b-over-b <ratio>\n"
      << "  --coefficient-correlation-length <m>\n"
      << "  --coefficient-k-min <1/m> --coefficient-k-max <1/m>\n"
      << "  --coefficient-reference-radius <m>\n"
      << "  --coefficient-k-min-radial-exponent <q>\n"
      << "  --coefficient-k-max-radial-exponent <q>\n"
      << "  --coefficient-quadrature-absolute <m2/s>\n"
      << "  --coefficient-quadrature-relative <fraction>\n"
      << "                               Configure named prescribed scales and adaptive\n"
      << "                               Dmumu-to-kappa integration tolerances.\n"
      << "\n"
      << "Transport numerical error controls:\n"
      << "  --transport-geometry-fraction <fraction>\n"
      << "  --transport-deterministic-tolerance <fraction>\n"
      << "  --transport-stochastic-mu-rms <fraction>\n"
      << "  --transport-cooling-log-change <fraction>\n"
      << "  --transport-focusing-mu-change <fraction>\n"
      << "  --transport-shock-fraction <fraction>\n"
      << "  --transport-min-step <s>    Set the shared mover error budget.\n"
      << "\n"
      << "Field-line SEP injection controls:\n"
      << "  --particles-per-iteration <N>\n"
      << "                               Set SEP::FieldLine::InjectionParameters::\n"
      << "                               nParticlesPerIteration.  Default: N=300.\n"
      << "                               This controls the target number of injected\n"
      << "                               SEP macroparticles per active injection event;\n"
      << "                               the code adjusts statistical weights when the\n"
      << "                               physical injection estimate differs from N.\n"
      << "                               N must be a positive integer.\n"
      << "  --particles-per-iteration=... Same option using --option=value syntax.\n"
      << "  --n-particles-per-iteration <N>\n"
      << "                               Alias for --particles-per-iteration.\n"
      << "\n"
      << "Authoritative run controls (WP30):\n"
      << "  --total-iterations <N>       Positive standalone driver iteration count.\n"
      << "  --shock-model <analytical|swcme1d>\n"
      << "  --cme-scenario <fast|slow>   Select the frozen SWCME parameter set.\n"
      << "  --swcme-failure-policy <strict|clamp-radius|diagnostic-fallback>\n"
      << "                               strict rejects every invalid query (default);\n"
      << "                               recovery modes are counted and fingerprinted.\n"
      << "  --swcme-fallback-density <m^-3> --swcme-fallback-speed <m/s>\n"
      << "  --swcme-fallback-divergence <s^-1>\n"
      << "                               Explicit sample used only by diagnostic-fallback.\n"
      << "  --swcme-override <key=value> Repeatable canonical SWCME/shock/source override.\n"
      << "                               Values carry explicit units, e.g.\n"
      << "                               --swcme-override 'cme.launch_speed=1400 km/s'.\n"
      << "  --field-line-seed-area <m2> Positive seed surface area used once to form Phi_i.\n"
      << "  --shock-turbulence-efficiency <0..1>\n"
      << "  --shock-turbulence-plus-fraction <0..1>\n"
      << "  --merge-minimum <N> --merge-maximum <N>\n"
      << "                               Configure particle merge/split population bounds.\n"
      << "  --injection-particles-per-iteration <N>\n"
      << "                               Alias for --particles-per-iteration.\n"
      << "\n"
      << "Wave-number-resolved spectrum output:\n"
      << "  --spectrum-output-interval <N>\n"
      << "                               Write the Tecplot 2-D W+/- (s,k), sigma_c(s,k)\n"
      << "                               diagnostic every N main-loop iterations when\n"
      << "                               --turbulence-model wave-number-resolved is active.\n"
      << "                               Default: N=100.  Use N=0 to disable.\n"
      << "  --spectrum-output-interval=... Same option using --option=value syntax.\n"
      << "  --wave-number-output-interval <N>\n"
      << "                               Alias for --spectrum-output-interval.\n"
      << "  --turbulence-spectrum-output-interval <N>\n"
      << "                               Alias for --spectrum-output-interval.\n"
      << "\n"
      << "Diagnostics and development tests:\n"
      << "  --test-manager <on|off>      Enable/disable the standalone SEP TestManager()\n"
      << "                               diagnostics after AMPS/field-line initialization.\n"
      << "  --testmanager <on|off>       Alias for --test-manager.\n"
      << "  --run-test-manager           Shortcut for --test-manager on.\n"
      << "  --no-test-manager            Shortcut for --test-manager off.\n"
      << "\n"
      << "Selectable standalone component tests (test-only execution):\n"
      << "  --list-tests                 List stable test IDs and metadata, then exit\n"
      << "                               before model initialization.\n"
      << "  --test <ID>                  Run one test; repeat to select several tests.\n"
      << "  --test=<ID>                  Equivalent equals-sign form.\n"
      << "  --test-group <GROUP>         Run every test in a group; repeatable.\n"
      << "  --test-group=<GROUP>         Equivalent equals-sign form.\n"
      << "  --all-tests                  Run the bounded routine component-test set.\n"
      << "                               Extended tests remain individually selectable.\n"
      << "                               IDs/groups are case-insensitive and overlaps\n"
      << "                               execute once in stable ID order.\n"
      << "  --test-json <path>           Write the complete result/metric/evidence\n"
      << "                               summary as srcsep-component-tests-v1 JSON.\n"
      << "  --test-junit <path>          Write the same outcomes as JUnit XML.\n"
      << "  --test-input <path>          Reviewed native input for one end-to-end\n"
      << "                               validation test selected by explicit ID.\n"
      << "  --test-output-dir <path>     Isolated directory for that test's native\n"
      << "                               model artifacts.\n"
      << "\n"
      << "General:\n"
      << "  -h, --help                   Print this help message and exit before AMPS\n"
      << "                               initialization.\n"
      << "\n"
      << "Accepted boolean values are: on/off, true/false, yes/no, 1/0,\n"
      << "enable/disable, and enabled/disabled.\n"
      << "\n"
      << "Defaults:\n"
      << "  coupling=on, cascade=on, reflection=on, turbulence-model=integrated, particle-mover=fte-dmumu,\n"
      << "  coefficient-source=prescribed, spatial=from-dmumu, pitch-angle=configured,\n"
      << "  mean-free-path=tenishev-2005, invalid-coefficient-policy=fail,\n"
      << "  particles-per-iteration=300, test-manager=off, spectrum-output-interval=100.\n"
      << "\n"
      << "Examples:\n"
      << "  " << exe << " --coupling off --cascade off --reflection off\n"
      << "  " << exe << " --coupling-mode=on --no-cascade --reflection=on\n"
      << "  " << exe << " --run-test-manager\n"
      << "  " << exe << " --list-tests\n"
      << "  " << exe << " --test DXX01 --test=TURB01\n"
      << "  " << exe << " --test-group parker\n"
      << "  " << exe << " --all-tests --test-json results.json --test-junit results.xml\n"
      << "  " << exe << " --test CV01 --test-input CV01-native.args --test-output-dir evidence/CV01\n"
      << "  " << exe << " --all-tests\n"
      << "  " << exe << " --turbulence-model wave-number-resolved --coupling on\n"
      << "  " << exe << " --wave-number-resolved --particle-mover fte-mfp --coupling on\n"
      << "  " << exe << " --particle-mover fte-mfp --mean-free-path-provider qlt1\n"
      << "  " << exe << " --coefficient-source swmf --mean-free-path-provider from-spatial\n"
      << "  " << exe << " --particles-per-iteration 1000\n"
      << "  " << exe << " --wave-number-resolved --spectrum-output-interval 25\n";
}

bool ParseCommandLine(int argc, char** argv, Options& options,
                      std::ostream& out, std::ostream& err) {

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i] ? argv[i] : "";

    if (arg == "-h" || arg == "--help") {
      options.printHelp = true;
      return true;
    }

    std::string option_name, value_from_equals;
    SplitOption(arg, option_name, value_from_equals);

    if (option_name == "--input") {
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                          options.inputPath, err)) return false;
      if (options.inputPath.empty()) {
        err << "ERROR: option '--input' requires a non-empty path.\n";
        return false;
      }
      continue;
    }

    if (option_name == "--list-tests") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '--list-tests' does not take a value.\n";
        return false;
      }
      options.listTests = true;
      continue;
    }

    if (option_name == "--list-movers") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '--list-movers' does not take a value.\n";
        return false;
      }
      options.listMovers = true;
      continue;
    }

    if (option_name == "--all-tests") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '--all-tests' does not take a value.\n";
        return false;
      }
      options.runAllTests = true;
      continue;
    }

    if (option_name == "--test" || option_name == "--test-group") {
      std::string selector;
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                          selector, err)) {
        return false;
      }
      if (selector.empty()) {
        err << "ERROR: option '" << option_name
            << "' requires a non-empty selector.\n";
        return false;
      }

      if (option_name == "--test") options.testIds.push_back(selector);
      else options.testGroups.push_back(selector);
      continue;
    }

    if (option_name == "--test-json" || option_name == "--test-junit" ||
        option_name == "--test-input" || option_name == "--test-output-dir") {
      std::string path;
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                          path, err)) return false;
      if (path.empty()) {
        err << "ERROR: option '" << option_name
            << "' requires a non-empty output path.\n";
        return false;
      }
      if (option_name == "--test-json") options.testJsonPath = path;
      else if (option_name == "--test-junit") options.testJunitPath = path;
      else if (option_name == "--test-input") options.testInputPath = path;
      else options.testArtifactDirectory = path;
      continue;
    }

    if (option_name == "--coupling" || option_name == "--coupling-mode" ||
        option_name == "--particle-coupling") {
      if (!ParseBooleanOption(argc, argv, i, option_name, value_from_equals,
                              "particle/turbulence coupling", options.particleCouplingMode, err)) {
        return false;
      }
      continue;
    }

    if (option_name == "--cascade" || option_name == "--cascase") {
      if (!ParseBooleanOption(argc, argv, i, option_name, value_from_equals,
                              "nonlinear turbulence cascade", options.cascadeActive, err)) {
        return false;
      }
      continue;
    }

    if (option_name == "--reflection") {
      if (!ParseBooleanOption(argc, argv, i, option_name, value_from_equals,
                              "turbulence reflection", options.reflectionActive, err)) {
        return false;
      }
      continue;
    }

    if (option_name == "--turbulence-advection" ||
        option_name == "--shock-injection") {
      bool* destination = option_name == "--turbulence-advection"
          ? &options.turbulence.advectionEnabled
          : &options.turbulence.shockInjectionEnabled;
      if (!ParseBooleanOption(argc, argv, i, option_name, value_from_equals,
                              option_name.c_str(), *destination, err))
        return false;
      continue;
    }

    if (option_name == "--turbulence-model") {
      std::string raw_value;
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals, raw_value, err)) return false;

      Options::TurbulenceModel parsed_model = Options::TurbulenceModel::Integrated;
      if (!ParseTurbulenceModelValue(raw_value, parsed_model)) {
        err << "ERROR: invalid value '" << raw_value << "' for --turbulence-model. "
            << "Use integrated or wave-number-resolved.\n";
        return false;
      }

      options.turbulenceModel = parsed_model;
      options.turbulence.representation =
          parsed_model == Options::TurbulenceModel::WaveNumberResolved
              ? Turbulence::Representation::Spectral
              : Turbulence::Representation::Integrated;
      if (options.turbulence.source == Turbulence::Source::SelfConsistentIntegrated ||
          options.turbulence.source == Turbulence::Source::SelfConsistentSpectral)
        options.turbulence.source =
            parsed_model == Options::TurbulenceModel::WaveNumberResolved
                ? Turbulence::Source::SelfConsistentSpectral
                : Turbulence::Source::SelfConsistentIntegrated;
      continue;
    }

    if (option_name == "--turbulence-source" ||
        option_name == "--turbulence-coupling-policy" ||
        option_name == "--turbulence-inner-boundary" ||
        option_name == "--turbulence-outer-boundary") {
      std::string raw_value;
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                          raw_value, err)) return false;
      const std::string canonical = ToLower(raw_value);
      bool parsed = false;
      if (option_name == "--turbulence-source") {
        parsed = Turbulence::ParseSource(canonical, &options.turbulence.source);
        if (parsed && options.turbulence.source ==
                          Turbulence::Source::SelfConsistentSpectral) {
          options.turbulence.representation = Turbulence::Representation::Spectral;
          options.turbulenceModel = Options::TurbulenceModel::WaveNumberResolved;
        }
        else if (parsed && options.turbulence.source ==
                               Turbulence::Source::SelfConsistentIntegrated) {
          options.turbulence.representation = Turbulence::Representation::Integrated;
          options.turbulenceModel = Options::TurbulenceModel::Integrated;
        }
      }
      else if (option_name == "--turbulence-coupling-policy")
        parsed = Turbulence::ParseCouplingPolicy(canonical,
                                                &options.turbulence.coupling);
      else if (option_name == "--turbulence-inner-boundary")
        parsed = Turbulence::ParseBoundaryPolicy(
            canonical, &options.turbulence.innerBoundary.policy);
      else
        parsed = Turbulence::ParseBoundaryPolicy(
            canonical, &options.turbulence.outerBoundary.policy);
      if (!parsed) {
        err << "ERROR: invalid value '" << raw_value << "' for "
            << option_name << ". Run with --help for canonical names.\n";
        return false;
      }
      if (option_name == "--turbulence-coupling-policy")
        options.particleCouplingMode =
            options.turbulence.coupling != Turbulence::CouplingPolicy::Disabled;
      continue;
    }

    if (option_name == "--turbulence-inner-value" ||
        option_name == "--turbulence-outer-value" ||
        option_name == "--reflection-coefficient" ||
        option_name == "--cascade-coefficient" ||
        option_name == "--turbulence-cfl" ||
        option_name == "--turbulence-operator-safety" ||
        option_name == "--turbulence-max-source-fraction" ||
        option_name == "--turbulence-max-cascade-fraction" ||
        option_name == "--turbulence-min-substep" ||
        option_name == "--turbulence-correlation-length" ||
        option_name == "--spectral-k-min" || option_name == "--spectral-k-max" ||
        option_name == "--turbulence-conservation-tolerance") {
      double* destination = nullptr;
      bool non_negative = true;
      if (option_name == "--turbulence-inner-value")
        destination = &options.turbulence.innerBoundary.value;
      else if (option_name == "--turbulence-outer-value")
        destination = &options.turbulence.outerBoundary.value;
      else if (option_name == "--reflection-coefficient")
        destination = &options.turbulence.reflectionCoefficient;
      else if (option_name == "--cascade-coefficient")
        destination = &options.turbulence.cascadeCoefficient;
      else if (option_name == "--turbulence-cfl") {
        destination = &options.turbulence.cflSafety;
        non_negative = false;
      }
      else if (option_name == "--turbulence-operator-safety") {
        destination = &options.turbulence.operatorAccuracySafety;
        non_negative = false;
      }
      else if (option_name == "--turbulence-max-source-fraction") {
        destination = &options.turbulence.maximumSourceFraction;
        non_negative = false;
      }
      else if (option_name == "--turbulence-max-cascade-fraction") {
        destination = &options.turbulence.maximumCascadeFraction;
        non_negative = false;
      }
      else if (option_name == "--turbulence-min-substep") {
        destination = &options.turbulence.minimumSubstepS;
        non_negative = false;
      }
      else if (option_name == "--turbulence-correlation-length") {
        destination = &options.turbulence.perpendicularCorrelationLengthM;
        non_negative = false;
      }
      else if (option_name == "--spectral-k-min") {
        destination = &options.turbulence.spectralKMinPerM;
        non_negative = false;
      }
      else if (option_name == "--spectral-k-max") {
        destination = &options.turbulence.spectralKMaxPerM;
        non_negative = false;
      }
      else destination = &options.turbulence.conservationRelativeTolerance;
      if (!ParseDoubleOption(argc, argv, i, option_name, value_from_equals,
                             option_name.c_str(), non_negative,
                             *destination, err)) return false;
      continue;
    }

    if (option_name == "--spectral-bins") {
      int bins = 0;
      if (!ParsePositiveIntegerOption(argc, argv, i, option_name,
                                      value_from_equals, "spectral bin count",
                                      bins, err)) return false;
      options.turbulence.spectralBins = static_cast<std::size_t>(bins);
      continue;
    }

    if (option_name == "--turbulence-max-substeps") {
      int count = 0;
      if (!ParsePositiveIntegerOption(argc, argv, i, option_name,
                                      value_from_equals,
                                      "turbulence maximum substeps",
                                      count, err)) return false;
      options.turbulence.maximumSubsteps = static_cast<std::uint64_t>(count);
      continue;
    }

    if (option_name == "--particle-mover" || option_name == "--mover" ||
        option_name == "--sep-mover") {
      std::string raw_value;
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals, raw_value, err)) return false;

      Mover::ProductionMover parsed_mover =
          Mover::ProductionMover::FocusedTransportDiffusion;
      std::string alias_warning;
      if (!Mover::ParseProductionMover(raw_value, parsed_mover, alias_warning)) {
        err << "ERROR: invalid value '" << raw_value << "' for " << option_name << ".\n"
            << "       Accepted production movers: parker, fte-dmumu, fte-mfp.\n";
        return false;
      }

      options.particleMover = parsed_mover;
      options.particleMoverProvided = true;
      if (!alias_warning.empty()) out << "WARNING: " << alias_warning << ".\n";
      continue;
    }

    if (option_name == "--coefficient-source" ||
        option_name == "--spatial-diffusion-provider" ||
        option_name == "--pitch-angle-diffusion-provider" ||
        option_name == "--mean-free-path-provider" ||
        option_name == "--invalid-coefficient-policy" ||
        option_name == "--resonance-gap-policy" ||
        option_name == "--turbulence-amplitude-policy") {
      std::string raw_value;
      if (!GetOptionValue(argc, argv, i, option_name, value_from_equals,
                          raw_value, err)) return false;
      bool parsed = false;
      if (option_name == "--coefficient-source") {
        parsed = Transport::Coefficient::ParseSource(
            raw_value, &options.coefficients.source);
      }
      else if (option_name == "--spatial-diffusion-provider") {
        parsed = Transport::Coefficient::ParseSpatial(
            raw_value, &options.coefficients.spatial);
      }
      else if (option_name == "--pitch-angle-diffusion-provider") {
        parsed = Transport::Coefficient::ParsePitchAngle(
            raw_value, &options.coefficients.pitchAngle);
      }
      else if (option_name == "--mean-free-path-provider") {
        parsed = Transport::Coefficient::ParseMeanFreePath(
            raw_value, &options.coefficients.meanFreePath);
      }
      else if (option_name == "--invalid-coefficient-policy") {
        parsed = Transport::Coefficient::ParseInvalidPolicy(
            raw_value, &options.coefficients.invalidPolicy);
      }
      else if (option_name == "--resonance-gap-policy") {
        parsed = Transport::Coefficient::ParseResonanceGapPolicy(
            raw_value, &options.coefficients.resonanceGapPolicy);
      }
      else {
        parsed = Transport::Coefficient::ParseTurbulenceAmplitudePolicy(
            raw_value, &options.coefficients.amplitudePolicy);
      }
      if (!parsed) {
        err << "ERROR: invalid value '" << raw_value << "' for "
            << option_name << ". Run with --help for canonical names.\n";
        return false;
      }
      continue;
    }

    if (option_name == "--constant-dmumu" ||
        option_name == "--prescribed-delta-b-over-b" ||
        option_name == "--coefficient-correlation-length" ||
        option_name == "--coefficient-k-min" ||
        option_name == "--coefficient-k-max" ||
        option_name == "--coefficient-reference-radius" ||
        option_name == "--coefficient-k-min-radial-exponent" ||
        option_name == "--coefficient-k-max-radial-exponent" ||
        option_name == "--coefficient-quadrature-absolute" ||
        option_name == "--coefficient-quadrature-relative" ||
        option_name == "--transport-geometry-fraction" ||
        option_name == "--transport-deterministic-tolerance" ||
        option_name == "--transport-stochastic-mu-rms" ||
        option_name == "--transport-cooling-log-change" ||
        option_name == "--transport-focusing-mu-change" ||
        option_name == "--transport-shock-fraction" ||
        option_name == "--transport-min-step") {
      double* destination = NULL;
      bool nonNegative = false;
      if (option_name == "--constant-dmumu") {
        destination = &options.constantDmumuPerS;
        nonNegative = true;
        options.constantDmumuProvided = true;
      }
      else if (option_name == "--prescribed-delta-b-over-b") {
        destination = &options.coefficients.prescribedDeltaBOverB;
        nonNegative = true;
      }
      else if (option_name == "--coefficient-correlation-length")
        destination = &options.coefficients.correlationLengthAt1AuM;
      else if (option_name == "--coefficient-k-min")
        destination = &options.coefficients.spectrum.kMinAtReferencePerM;
      else if (option_name == "--coefficient-k-max")
        destination = &options.coefficients.spectrum.kMaxAtReferencePerM;
      else if (option_name == "--coefficient-reference-radius")
        destination = &options.coefficients.spectrum.referenceRadiusM;
      else if (option_name == "--coefficient-k-min-radial-exponent") {
        destination = &options.coefficients.spectrum.kMinRadialExponent;
        nonNegative = true;
      }
      else if (option_name == "--coefficient-k-max-radial-exponent") {
        destination = &options.coefficients.spectrum.kMaxRadialExponent;
        nonNegative = true;
      }
      else if (option_name == "--coefficient-quadrature-absolute") {
        destination = &options.coefficients.spatialQuadrature.absoluteToleranceM2PerS;
        nonNegative = true;
      }
      else if (option_name == "--coefficient-quadrature-relative")
        destination = &options.coefficients.spatialQuadrature.relativeTolerance;
      else if (option_name == "--transport-geometry-fraction")
        destination = &options.numericalTolerances.geometryFraction;
      else if (option_name == "--transport-deterministic-tolerance")
        destination = &options.numericalTolerances.deterministicRelativeTolerance;
      else if (option_name == "--transport-stochastic-mu-rms")
        destination = &options.numericalTolerances.stochasticPitchRms;
      else if (option_name == "--transport-cooling-log-change")
        destination = &options.numericalTolerances.coolingLogChange;
      else if (option_name == "--transport-focusing-mu-change")
        destination = &options.numericalTolerances.focusingPitchChange;
      else if (option_name == "--transport-shock-fraction")
        destination = &options.numericalTolerances.shockFraction;
      else {
        destination = &options.numericalTolerances.minimumStepS;
        nonNegative = true;
      }
      if (!ParseDoubleOption(argc, argv, i, option_name, value_from_equals,
                             option_name.c_str(), nonNegative,
                             *destination, err)) return false;
      if (option_name == "--constant-dmumu") {
        // Keep the pure registry value synchronized with the compatibility
        // callback value.  Provider provenance fingerprints the registry
        // record, while legacy diagnostics still read the public scalar.
        options.coefficients.constantDmumuPerS =
            options.constantDmumuPerS;
      }
      continue;
    }

    if (option_name == "--particles-per-iteration" ||
        option_name == "--n-particles-per-iteration" ||
        option_name == "--injection-particles-per-iteration" ||
        option_name == "--field-line-injection-particles") {
      if (!ParsePositiveIntegerOption(argc, argv, i, option_name, value_from_equals,
                                      "field-line injection particles per iteration",
                                      options.injectionParticlesPerIteration, err)) {
        return false;
      }
      continue;
    }

    if (option_name=="--total-iterations" || option_name=="--merge-minimum" ||
        option_name=="--merge-maximum") {
      int* destination=option_name=="--total-iterations" ? &options.totalIterations
          : (option_name=="--merge-minimum" ? &options.mergeMinimum
                                             : &options.mergeMaximum);
      if (!ParsePositiveIntegerOption(argc,argv,i,option_name,value_from_equals,
                                      "WP30 positive run control",*destination,err))
        return false;
      if (option_name=="--total-iterations") options.totalIterationsProvided=true;
      else if (option_name=="--merge-minimum") options.mergeMinimumProvided=true;
      else options.mergeMaximumProvided=true;
      continue;
    }

    if (option_name=="--field-line-seed-area" ||
        option_name=="--shock-turbulence-efficiency" ||
        option_name=="--shock-turbulence-plus-fraction") {
      double* destination=option_name=="--field-line-seed-area"
          ? &options.fieldLineSeedAreaM2
          : (option_name=="--shock-turbulence-efficiency"
              ? &options.shockTurbulenceEfficiency
              : &options.shockTurbulencePlusFraction);
      if (!ParseDoubleOption(argc,argv,i,option_name,value_from_equals,
                             "WP30 SI/fraction run control",
                             option_name!="--field-line-seed-area",
                             *destination,err))
        return false;
      if (option_name=="--field-line-seed-area") options.fieldLineSeedAreaProvided=true;
      else if (option_name=="--shock-turbulence-efficiency")
        options.shockTurbulenceEfficiencyProvided=true;
      else options.shockTurbulencePlusFractionProvided=true;
      continue;
    }

    if (option_name=="--shock-model" || option_name=="--cme-scenario") {
      std::string value;
      if (!GetOptionValue(argc,argv,i,option_name,value_from_equals,value,err)) return false;
      value=ToLower(value);
      if (option_name=="--shock-model") {
        if (value=="analytical" || value=="analytic") options.analyticalShock=true;
        else if (value=="swcme1d" || value=="swcme") options.analyticalShock=false;
        else { err<<"ERROR: --shock-model requires analytical or swcme1d.\n";return false; }
        options.shockModelProvided=true;
      }
      else {
        if (value=="slow") options.slowCmeScenario=true;
        else if (value=="fast") options.slowCmeScenario=false;
        else { err<<"ERROR: --cme-scenario requires fast or slow.\n";return false; }
        options.cmeScenarioProvided=true;
      }
      continue;
    }

    if (option_name=="--swcme-failure-policy") {
      std::string value;
      if (!GetOptionValue(argc,argv,i,option_name,value_from_equals,value,err))
        return false;
      value=ToLower(value);
      if (value=="strict")
        options.swcmeFailurePolicy=Options::SwcmeFailurePolicy::Strict;
      else if (value=="clamp-radius" || value=="clamp")
        options.swcmeFailurePolicy=Options::SwcmeFailurePolicy::ClampRadius;
      else if (value=="diagnostic-fallback" || value=="fallback")
        options.swcmeFailurePolicy=
            Options::SwcmeFailurePolicy::DiagnosticFallback;
      else {
        err<<"ERROR: --swcme-failure-policy requires strict, clamp-radius, "
              "or diagnostic-fallback.\n";
        return false;
      }
      options.swcmeFailurePolicyProvided=true;
      continue;
    }

    if (option_name=="--swcme-override") {
      std::string assignment;
      if (!GetOptionValue(argc,argv,i,option_name,value_from_equals,
                          assignment,err)) return false;
      const std::string::size_type separator=assignment.find('=');
      if (separator==std::string::npos || separator==0 ||
          separator+1>=assignment.size()) {
        err<<"ERROR: --swcme-override requires non-empty key=value.\n";
        return false;
      }
      options.swcmeOverrides.push_back(assignment);
      continue;
    }

    if (option_name=="--swcme-fallback-density" ||
        option_name=="--swcme-fallback-speed" ||
        option_name=="--swcme-fallback-divergence") {
      double* destination=option_name=="--swcme-fallback-density"
          ? &options.swcmeFallbackDensityM3
          : (option_name=="--swcme-fallback-speed"
              ? &options.swcmeFallbackSpeedMPerS
              : &options.swcmeFallbackDivergencePerS);
      const bool parsed=option_name=="--swcme-fallback-divergence"
          ? ParseFiniteDoubleOption(argc,argv,i,option_name,value_from_equals,
                                    "SWCME fallback divergence [s^-1]",
                                    *destination,err)
          : ParseDoubleOption(argc,argv,i,option_name,value_from_equals,
                              "SWCME fallback density/speed in SI",false,
                              *destination,err);
      if (!parsed) return false;
      if (option_name=="--swcme-fallback-density")
        options.swcmeFallbackDensityProvided=true;
      else if (option_name=="--swcme-fallback-speed")
        options.swcmeFallbackSpeedProvided=true;
      else options.swcmeFallbackDivergenceProvided=true;
      continue;
    }

    if (option_name == "--spectrum-output-interval" ||
        option_name == "--wave-number-output-interval" ||
        option_name == "--turbulence-spectrum-output-interval" ||
        option_name == "--spectral-output-interval") {
      if (!ParseNonNegativeIntegerOption(argc, argv, i, option_name, value_from_equals,
                                         "wave-number-resolved spectrum output interval",
                                         options.spectralOutputInterval, err)) {
        return false;
      }
      continue;
    }

    if (option_name == "--test-manager" || option_name == "--testmanager") {
      if (!ParseBooleanOption(argc, argv, i, option_name, value_from_equals,
                              "SEP TestManager diagnostics", options.runTestManager, err)) {
        return false;
      }
      continue;
    }

    if (option_name == "--wave-number-resolved" || option_name == "--k-resolved" || option_name == "--spectral-turbulence") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.turbulenceModel = Options::TurbulenceModel::WaveNumberResolved;
      options.turbulence.representation = Turbulence::Representation::Spectral;
      if (options.turbulence.source == Turbulence::Source::SelfConsistentIntegrated)
        options.turbulence.source = Turbulence::Source::SelfConsistentSpectral;
      continue;
    }

    if (option_name == "--integrated-turbulence" || option_name == "--legacy-turbulence") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.turbulenceModel = Options::TurbulenceModel::Integrated;
      options.turbulence.representation = Turbulence::Representation::Integrated;
      if (options.turbulence.source == Turbulence::Source::SelfConsistentSpectral)
        options.turbulence.source = Turbulence::Source::SelfConsistentIntegrated;
      continue;
    }

    if (option_name == "--no-coupling" || option_name == "--no-particle-coupling") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.particleCouplingMode = false;
      continue;
    }

    if (option_name == "--no-cascade" || option_name == "--no-cascase") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.cascadeActive = false;
      continue;
    }

    if (option_name == "--no-reflection") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.reflectionActive = false;
      continue;
    }

    if (option_name == "--run-test-manager") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.runTestManager = true;
      continue;
    }

    if (option_name == "--no-test-manager") {
      if (!value_from_equals.empty()) {
        err << "ERROR: option '" << option_name << "' does not take a value.\n";
        return false;
      }
      options.runTestManager = false;
      continue;
    }

    err << "ERROR: unknown SEP command-line option '" << arg << "'.\n"
        << "Run with -h or --help to see available options.\n";
    return false;
  }

  const bool executionRequested = IsComponentTestExecutionRequested(options);
  if (options.listMovers &&
      (options.listTests || executionRequested || options.printHelp)) {
    err << "ERROR: --list-movers cannot be combined with help or component-test selectors.\n";
    return false;
  }
  if (options.listTests && executionRequested) {
    err << "ERROR: --list-tests cannot be combined with --test, --test-group, "
        << "or --all-tests.\n";
    return false;
  }
  if (!options.inputPath.empty() &&
      (options.listTests || options.listMovers || executionRequested)) {
    err << "ERROR: --input configures a production initialization and cannot "
        << "be combined with test/list modes.\n";
    return false;
  }

  // --all-tests already denotes the complete bounded selection.  Rejecting a
  // simultaneous explicit selector prevents scripts from assuming that an
  // extended test was included in --all-tests when routine-only policy omits it.
  if (options.runAllTests &&
      (!options.testIds.empty() || !options.testGroups.empty())) {
    err << "ERROR: --all-tests cannot be combined with --test or --test-group.\n";
    return false;
  }

  if ((!options.testJsonPath.empty() || !options.testJunitPath.empty()) &&
      !executionRequested) {
    err << "ERROR: --test-json and --test-junit require --test, --test-group, "
        << "or --all-tests.\n";
    return false;
  }

  // Input and artifact paths form one indivisible end-to-end test contract.
  // Requiring one explicit ID prevents a group/all invocation from sharing a
  // case-specific file accidentally, and requiring both paths prevents a case
  // from falling back to an unreviewed default or the process working directory.
  const bool hasTestInput = !options.testInputPath.empty();
  const bool hasTestArtifacts = !options.testArtifactDirectory.empty();
  if (hasTestInput != hasTestArtifacts) {
    err << "ERROR: --test-input and --test-output-dir must be provided together.\n";
    return false;
  }
  if (hasTestInput &&
      (options.testIds.size() != 1 || !options.testGroups.empty() ||
       options.runAllTests)) {
    err << "ERROR: end-to-end test input requires exactly one explicit --test ID.\n";
    return false;
  }

  const Transport::Status coefficientStatus =
      Transport::Coefficient::ValidateConfiguration(options.coefficients);
  if (!coefficientStatus.ok()) {
    err << "ERROR: invalid coefficient configuration: "
        << coefficientStatus.message << ".\n";
    return false;
  }
  const Transport::Status compatibilityStatus =
      Transport::Coefficient::ValidateMoverCompatibility(
          options.coefficients,
          Mover::Describe(options.particleMover).canonicalName);
  if (!compatibilityStatus.ok()) {
    err << "ERROR: unsupported mover/coefficient combination: "
        << compatibilityStatus.message << ".\n";
    return false;
  }
  const Transport::Status toleranceStatus =
      Transport::ValidateNumericalTolerances(options.numericalTolerances);
  if (!toleranceStatus.ok()) {
    err << "ERROR: invalid transport error controls: "
        << toleranceStatus.message << ".\n";
    return false;
  }
  if (!std::isfinite(options.constantDmumuPerS) ||
      options.constantDmumuPerS < 0.0) {
    err << "ERROR: constant Dmumu must be finite and non-negative [1/s].\n";
    return false;
  }
  if (options.shockTurbulenceEfficiency>1.0 ||
      options.shockTurbulencePlusFraction>1.0 ||
      options.mergeMaximum<options.mergeMinimum) {
    err << "ERROR: shock fractions must be in [0,1] and merge maximum must not be below minimum.\n";
    return false;
  }

  // Compatibility switches feed the one authoritative Step 11 configuration;
  // they are not a second ownership path.
  options.turbulence.coupling = options.particleCouplingMode
      ? Turbulence::CouplingPolicy::StreamingEnergyExchange
      : Turbulence::CouplingPolicy::Disabled;
  options.turbulence.cascadeEnabled = options.cascadeActive;
  options.turbulence.reflectionEnabled = options.reflectionActive;
  options.turbulence.diagnosticCadence =
      static_cast<std::size_t>(options.spectralOutputInterval);
  const Transport::Status turbulenceStatus =
      Turbulence::ValidateConfiguration(options.turbulence);
  if (!turbulenceStatus.ok()) {
    err << "ERROR: invalid turbulence configuration: "
        << turbulenceStatus.message << ".\n";
    return false;
  }
  ConfigurationMatrix::Combination combination;
  combination.mover=options.particleMover;
  combination.coefficientSource=options.coefficients.source;
  combination.turbulenceSource=options.turbulence.source;
  combination.coupling=options.turbulence.coupling;
  const Transport::Status matrixStatus=
      ConfigurationMatrix::Preflight(combination,false);
  if (!matrixStatus.ok()) {
    err << "ERROR: unsupported production configuration: "
        << matrixStatus.message << ".\n";
    return false;
  }

  return true;
}

bool IsComponentTestExecutionRequested(const Options& options) {
  return options.runAllTests || !options.testIds.empty() ||
         !options.testGroups.empty();
}

#ifndef SEP_CLI_PARSE_ONLY
void ApplyTurbulenceOptions(const Options& options) {
  // These are the three runtime switches used by the turbulence operators.  The
  // assignment is centralized here so the standalone driver no longer hard-codes
  // the values in main.cpp; future input-file parsing can reuse this function or
  // the same Options structure.
  SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode = options.particleCouplingMode;
  SEP::AlfvenTurbulence_Kolmogorov::Cascade::active = options.cascadeActive;
  SEP::AlfvenTurbulence_Kolmogorov::Reflection::active = options.reflectionActive;

  // Configure the field-line SEP injection macroparticle count.
  //
  // This option writes the variable defined in field_line.cpp:
  //   SEP::FieldLine::InjectionParameters::nParticlesPerIteration
  //
  // The injection routine uses this number as the target Monte-Carlo particle
  // count and compensates differences between the physically estimated source
  // strength and the requested macroparticle count by adjusting statistical
  // weights.  Applying the CLI value here, after the optional post-compile
  // input file is parsed in main.cpp and before particle injection starts, lets
  // command-line runs override the hard-coded default without recompilation.
  SEP::FieldLine::InjectionParameters::nParticlesPerIteration =
      options.injectionParticlesPerIteration;

  // Runtime selection goes through the production registry; the CLI never
  // writes a raw mover function pointer.
  SEP::Mover::SelectProductionMover(options.particleMover);

  // ParseCommandLine has already rejected incompatible combinations.  Reapply
  // the same validation at the mutation boundary so non-CLI callers cannot
  // install a configuration that the production movers would interpret
  // recursively or against the wrong background authority.
  const Transport::Status compatibilityStatus =
      Transport::Coefficient::ValidateMoverCompatibility(
          options.coefficients,
          Mover::Describe(options.particleMover).canonicalName);
  if (!compatibilityStatus.ok())
    exit(__LINE__, __FILE__, compatibilityStatus.message.c_str());

  const Transport::Status coefficientStatus =
      Transport::Coefficient::SetActiveConfiguration(options.coefficients);
  if (!coefficientStatus.ok())
    exit(__LINE__, __FILE__, coefficientStatus.message.c_str());

  const Transport::Status toleranceStatus =
      Transport::SetActiveNumericalTolerances(options.numericalTolerances);
  if (!toleranceStatus.ok())
    exit(__LINE__, __FILE__, toleranceStatus.message.c_str());

  if (options.constantDmumuProvided)
    SEP::Diffusion::ConstPitchAngleDiffusionValue =
        options.coefficients.constantDmumuPerS;

  // A constant model selected in post-compile input may retain its input-file
  // value when the CLI does not override it.  Validate that effective value at
  // the final mutation boundary instead of silently resetting it to zero.
  if (options.coefficients.pitchAngle ==
          Transport::Coefficient::PitchAngleKind::Constant &&
      (!std::isfinite(SEP::Diffusion::ConstPitchAngleDiffusionValue) ||
       SEP::Diffusion::ConstPitchAngleDiffusionValue < 0.0)) {
    exit(__LINE__, __FILE__,
         "effective constant Dmumu must be finite and non-negative [1/s]");
  }
  // Explicit registry models also update the legacy compatibility callback so
  // diagnostics and any not-yet-migrated read-only sampling path observe the
  // same selection.  The production adapters themselves use source-bound pure
  // kernels and do not infer capability from these function-pointer values.
  switch (options.coefficients.pitchAngle) {
    case Transport::Coefficient::PitchAngleKind::Constant:
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient =
          SEP::Diffusion::Constant::GetPitchAngleDiffusionCoefficient;
      break;
    case Transport::Coefficient::PitchAngleKind::Jokipii1966:
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient =
          SEP::Diffusion::Jokopii1966AJ::GetPitchAngleDiffusionCoefficient;
      break;
    case Transport::Coefficient::PitchAngleKind::Florinskiy:
      SEP::Diffusion::GetPitchAngleDiffusionCoefficient =
          SEP::Diffusion::Florinskiy::GetPitchAngleDiffusionCoefficient;
      break;
    case Transport::Coefficient::PitchAngleKind::Configured:
      // Post-compile input owns the compatibility callback in this mode.
      break;
  }

  const Transport::Status turbulenceStatus =
      Turbulence::SetActiveConfiguration(options.turbulence);
  if (!turbulenceStatus.ok())
    exit(__LINE__, __FILE__, turbulenceStatus.message.c_str());

  // Select the wave-energy representation.  This affects only the turbulence
  // energy transport/coupling kernels.  All existing output and scattering code
  // still sees the integrated CellIntegratedWaveEnergy datum, which the new
  // spectral model keeps synchronized by summing over k-bins.
  SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::TurbulenceModelMode =
      (options.turbulenceModel == Options::TurbulenceModel::WaveNumberResolved)
          ? SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::ModelMode::WaveNumberResolved
          : SEP::AlfvenTurbulence_Kolmogorov::WaveNumberResolved::ModelMode::Integrated;
}
#endif

void PrintTurbulenceOptions(const Options& options, std::ostream& out) {
  out << "SEP turbulence CLI configuration:\n"
      << "  particle/turbulence coupling: " << (options.particleCouplingMode ? "on" : "off") << "\n"
      << "  nonlinear cascade:            " << (options.cascadeActive ? "on" : "off") << "\n"
      << "  wave reflection:               " << (options.reflectionActive ? "on" : "off") << "\n"
      << "  turbulence model:              "
      << (options.turbulenceModel == Options::TurbulenceModel::WaveNumberResolved
              ? "wave-number-resolved" : "integrated") << "\n"
      << "  turbulence source:             "
      << Turbulence::SourceName(options.turbulence.source) << "\n"
      << "  turbulence coupling policy:    "
      << Turbulence::CouplingPolicyName(options.turbulence.coupling) << "\n"
      << "  inner boundary:                "
      << Turbulence::BoundaryPolicyName(options.turbulence.innerBoundary.policy)
      << " (value=" << options.turbulence.innerBoundary.value << ")\n"
      << "  outer boundary:                "
      << Turbulence::BoundaryPolicyName(options.turbulence.outerBoundary.policy)
      << " (value=" << options.turbulence.outerBoundary.value << ")\n"
      << "  reflection/cascade factors:    "
      << options.turbulence.reflectionCoefficient << "/"
      << options.turbulence.cascadeCoefficient << "\n"
      << "  perpendicular scale [m]:       "
      << options.turbulence.perpendicularCorrelationLengthM << "\n"
      << "  spectral k range [1/m]:        "
      << options.turbulence.spectralKMinPerM << " .. "
      << options.turbulence.spectralKMaxPerM << " ("
      << options.turbulence.spectralBins << " bins)\n"
      << "  conservation tolerance:        "
      << options.turbulence.conservationRelativeTolerance << "\n"
      << "  advection/shock source:        "
      << (options.turbulence.advectionEnabled ? "on" : "off") << "/"
      << (options.turbulence.shockInjectionEnabled ? "on" : "off")
      << " (CFL=" << options.turbulence.cflSafety << ")\n"
      << "  operator safety/source/cascade: "
      << options.turbulence.operatorAccuracySafety << "/"
      << options.turbulence.maximumSourceFraction << "/"
      << options.turbulence.maximumCascadeFraction << "\n"
      << "  operator min step/max stages:  "
      << options.turbulence.minimumSubstepS << " s/"
      << options.turbulence.maximumSubsteps << "\n"
      << "  particle mover:                "
      << Mover::Describe(options.particleMover).canonicalName << "\n"
      << "  coefficient contract:          "
      << Mover::CoefficientContractName(
             Mover::Describe(options.particleMover).capabilities.coefficientContract)
      << "\n"
      << "  coefficient source:            "
      << Transport::Coefficient::SourceName(options.coefficients.source) << "\n"
      << "  spatial provider:              "
      << Transport::Coefficient::SpatialName(options.coefficients.spatial) << "\n"
      << "  pitch-angle provider:          "
      << Transport::Coefficient::PitchAngleName(
             options.coefficients.pitchAngle) << "\n"
      << "  mean-free-path provider:       "
      << Transport::Coefficient::MeanFreePathName(
             options.coefficients.meanFreePath) << "\n"
      << "  invalid-coefficient policy:    "
      << Transport::Coefficient::InvalidPolicyName(
             options.coefficients.invalidPolicy) << "\n"
      << "  resonance-gap policy:          "
      << Transport::Coefficient::ResonanceGapPolicyName(
             options.coefficients.resonanceGapPolicy) << "\n"
      << "  amplitude policy:              "
      << Transport::Coefficient::TurbulenceAmplitudePolicyName(
             options.coefficients.amplitudePolicy) << "\n"
      << "  constant Dmumu [1/s]:          "
      << options.constantDmumuPerS << "\n"
      << "  prescribed deltaB/B:           "
      << options.coefficients.prescribedDeltaBOverB << "\n"
      << "  correlation length at 1 AU [m]: "
      << options.coefficients.correlationLengthAt1AuM << "\n"
      << "  coefficient k range [1/m]:     "
      << options.coefficients.spectrum.kMinAtReferencePerM << " .. "
      << options.coefficients.spectrum.kMaxAtReferencePerM << "\n"
      << "  quadrature abs/rel tolerance:  "
      << options.coefficients.spatialQuadrature.absoluteToleranceM2PerS
      << " / " << options.coefficients.spatialQuadrature.relativeTolerance << "\n"
      << "  coefficient fingerprint:      "
      << Transport::Coefficient::ConfigurationFingerprint(
             options.coefficients) << "\n"
      << "  transport geometry fraction:   "
      << options.numericalTolerances.geometryFraction << "\n"
      << "  deterministic tolerance:       "
      << options.numericalTolerances.deterministicRelativeTolerance << "\n"
      << "  stochastic pitch RMS:          "
      << options.numericalTolerances.stochasticPitchRms << "\n"
      << "  cooling/focusing limits:       "
      << options.numericalTolerances.coolingLogChange << " / "
      << options.numericalTolerances.focusingPitchChange << "\n"
      << "  shock fraction/min step [s]:   "
      << options.numericalTolerances.shockFraction << " / "
      << options.numericalTolerances.minimumStepS << "\n"
      << "  injected particles/iteration:  " << options.injectionParticlesPerIteration << "\n"
      << "  spectrum output interval:      " << options.spectralOutputInterval
      << " iteration(s)" << (options.spectralOutputInterval == 0 ? " (disabled)" : "") << "\n"
      << "  TestManager diagnostics:       " << (options.runTestManager ? "on" : "off") << "\n";

  if (options.particleCouplingMode &&
      options.turbulenceModel == Options::TurbulenceModel::WaveNumberResolved &&
      !Mover::Describe(options.particleMover).capabilities.accumulatesWaveStreaming) {
    out << "  WARNING: wave-number-resolved particle coupling needs a mover that fills "
        << "G_+(k),G_-(k).\n";
  }
}

} // namespace CLI
} // namespace Util
} // namespace SEP
