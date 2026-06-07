#include "sep_cli.h"
#include "../sep.h"

#include <algorithm>
#include <cctype>
#include <iostream>
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
    err << "ERROR: option '" << option_name << "' requires a value: on/off, true/false, yes/no, or 1/0.\n";
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

} // anonymous namespace

void PrintHelp(const char* program_name, std::ostream& out) {
  const char* exe = (program_name && program_name[0] != '\0') ? program_name : "amps";

  out
      << "SEP standalone driver command-line options\n"
      << "\n"
      << "Usage:\n"
      << "  " << exe << " [options]\n"
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
      << "Diagnostics and development tests:\n"
      << "  --test-manager <on|off>      Enable/disable the standalone SEP TestManager()\n"
      << "                               diagnostics after AMPS/field-line initialization.\n"
      << "  --testmanager <on|off>       Alias for --test-manager.\n"
      << "  --run-test-manager           Shortcut for --test-manager on.\n"
      << "  --no-test-manager            Shortcut for --test-manager off.\n"
      << "\n"
      << "General:\n"
      << "  -h, --help                   Print this help message and exit before AMPS\n"
      << "                               initialization.\n"
      << "\n"
      << "Accepted boolean values are: on/off, true/false, yes/no, 1/0,\n"
      << "enable/disable, and enabled/disabled.\n"
      << "\n"
      << "Defaults:\n"
      << "  coupling=on, cascade=on, reflection=on, test-manager=off.\n"
      << "\n"
      << "Examples:\n"
      << "  " << exe << " --coupling off --cascade off --reflection off\n"
      << "  " << exe << " --coupling-mode=on --no-cascade --reflection=on\n"
      << "  " << exe << " --run-test-manager\n";
}

bool ParseCommandLine(int argc, char** argv, Options& options,
                      std::ostream& out, std::ostream& err) {
  (void)out; // Reserved for future informational parse-time messages.

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i] ? argv[i] : "";

    if (arg == "-h" || arg == "--help") {
      options.printHelp = true;
      return true;
    }

    std::string option_name, value_from_equals;
    SplitOption(arg, option_name, value_from_equals);

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

    if (option_name == "--test-manager" || option_name == "--testmanager") {
      if (!ParseBooleanOption(argc, argv, i, option_name, value_from_equals,
                              "SEP TestManager diagnostics", options.runTestManager, err)) {
        return false;
      }
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

  return true;
}

void ApplyTurbulenceOptions(const Options& options) {
  // These are the three runtime switches used by the turbulence operators.  The
  // assignment is centralized here so the standalone driver no longer hard-codes
  // the values in main.cpp; future input-file parsing can reuse this function or
  // the same Options structure.
  SEP::AlfvenTurbulence_Kolmogorov::ParticleCouplingMode = options.particleCouplingMode;
  SEP::AlfvenTurbulence_Kolmogorov::Cascade::active = options.cascadeActive;
  SEP::AlfvenTurbulence_Kolmogorov::Reflection::active = options.reflectionActive;
}

void PrintTurbulenceOptions(const Options& options, std::ostream& out) {
  out << "SEP turbulence CLI configuration:\n"
      << "  particle/turbulence coupling: " << (options.particleCouplingMode ? "on" : "off") << "\n"
      << "  nonlinear cascade:            " << (options.cascadeActive ? "on" : "off") << "\n"
      << "  wave reflection:               " << (options.reflectionActive ? "on" : "off") << "\n"
      << "  TestManager diagnostics:       " << (options.runTestManager ? "on" : "off") << "\n";
}

} // namespace CLI
} // namespace Util
} // namespace SEP
