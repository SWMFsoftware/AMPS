/*
 * shieldSimSourceTool
 * -------------------
 * A deliberately separate Stage-2 executable for validating and canonicalizing
 * published source products.  Keeping this operation outside shieldSim's
 * transport CLI enforces the stage boundary: successful ingestion does not
 * imply that a directional geosynchronous measurement or an incomplete
 * near-Earth GCR measurement is an approved transport boundary condition.
 */

#include "PublishedSourceAdapters.hh"

#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace {

void Usage(std::ostream& output) {
  output <<
      "Usage:\n"
      "  shieldSimSourceTool ams02 --info=<HAPI-info.json> --data=<HAPI-data.csv> \\\n"
      "      --interpretation=<interpretation.json> --output=<canonical.json|->\n"
      "  shieldSimSourceTool goes16 --csv=<event-timeseries.csv> \\\n"
      "      --definition=<event-definition.json> --output=<fluence.json|->\n\n"
      "The AMS command emits shieldSim.spectrum_dataset.v2. The GOES command\n"
      "emits shieldSim.goes_event_fluence.v1 with both directions retained.\n"
      "Neither command starts Geant4 transport. Unknown or duplicate options\n"
      "are fatal so a misspelled scientific input cannot be ignored.\n";
}

std::map<std::string,std::string> ParseOptions(int argc, char** argv, int first) {
  std::map<std::string,std::string> result;
  for (int i = first; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument.rfind("--", 0) != 0 || argument.find('=') == std::string::npos)
      throw std::runtime_error("option must use --name=value form: " + argument);
    const std::size_t equals = argument.find('=');
    const std::string name = argument.substr(2, equals - 2);
    const std::string value = argument.substr(equals + 1);
    if (name.empty() || value.empty()) throw std::runtime_error("empty option name/value: " + argument);
    if (!result.emplace(name, value).second) throw std::runtime_error("duplicate option --" + name);
  }
  return result;
}

std::string Take(std::map<std::string,std::string>& options, const std::string& name) {
  const auto found = options.find(name);
  if (found == options.end()) throw std::runtime_error("missing required option --" + name);
  const std::string result = found->second;
  options.erase(found);
  return result;
}

void Write(const std::string& path, const std::string& canonicalJson) {
  if (path == "-") {
    std::cout << canonicalJson << '\n';
    return;
  }
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot open output file '" + path + "'");
  output << canonicalJson << '\n';
  if (!output) throw std::runtime_error("failed while writing output file '" + path + "'");
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 1 || std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h") {
      Usage(std::cout);
      return argc == 1 ? 2 : 0;
    }
    const std::string command = argv[1];
    auto options = ParseOptions(argc, argv, 2);
    if (command == "ams02") {
      const std::string info = Take(options, "info");
      const std::string data = Take(options, "data");
      const std::string interpretation = Take(options, "interpretation");
      const std::string output = Take(options, "output");
      if (!options.empty()) throw std::runtime_error("unknown option --" + options.begin()->first);
      const auto dataset = shield_sim::AMS02HapiAdapter::Load(info, data, interpretation);
      Write(output, dataset.canonicalMetadata.DumpCanonical());
    } else if (command == "goes16") {
      const std::string csv = Take(options, "csv");
      const std::string definition = Take(options, "definition");
      const std::string output = Take(options, "output");
      if (!options.empty()) throw std::runtime_error("unknown option --" + options.begin()->first);
      const auto adapter = shield_sim::GOESRSGPSAdapter::Load(csv, definition);
      Write(output, adapter.Integrate().canonicalJson.DumpCanonical());
    } else {
      throw std::runtime_error("unknown command '" + command + "'");
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "shieldSimSourceTool: " << error.what() << '\n';
    return 2;
  }
}
