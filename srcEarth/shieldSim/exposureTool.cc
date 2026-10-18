/*
 * shieldSimExposureTool
 * ---------------------
 * Auditable Stage-7 operations that do not require a Geant4 process.  Inputs
 * and outputs are strict versioned JSON.  This tool cannot turn the incomplete
 * AMS component into a full GCR environment or enable GOES transport without
 * an explicit boundary mapping.
 */

#include "Stage7Qualification.hh"
#include "StrictJson.hh"
#include "BuildProvenance.hh"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void Usage(std::ostream& output) {
  output <<
      "Usage:\n"
      "  shieldSimExposureTool integrate INPUT.json OUTPUT.json\n"
      "  shieldSimExposureTool fold-spectrum INPUT.json OUTPUT.json\n"
      "  shieldSimExposureTool allocate INPUT.json OUTPUT.json\n\n"
      "integrate preserves explicit UTC intervals, gaps, configurations, anomaly\n"
      "timestamps, covariance, and sampled source ensembles. fold-spectrum uses\n"
      "adaptive quadrature and log-energy response interpolation without\n"
      "extrapolation. allocate performs deterministic pilot-based Neyman\n"
      "allocation and never enables Geant4 physics biasing.\n";
}

void Write(const std::string& path, const shield_sim::Json& document) {
  std::ofstream output(path, std::ios::binary | std::ios::trunc);
  if (!output) throw std::runtime_error("cannot write output file '" + path + "'");
  output << document.DumpCanonical() << '\n';
  if (!output) throw std::runtime_error("failed while writing output file '" + path + "'");
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
      Usage(std::cout);
      return 0;
    }
    if (argc != 4) {
      Usage(std::cerr);
      return 2;
    }
    const std::string command = argv[1];
    const shield_sim::Json input = shield_sim::Json::ParseFile(argv[2]);
    shield_sim::Json result;
    if (command == "integrate")
      result = shield_sim::IntegrateTimeResolvedExposure(input);
    else if (command == "fold-spectrum")
      result = shield_sim::FoldSpectrumDocument(input);
    else if (command == "allocate")
      result = shield_sim::AllocateHistoriesDocument(input);
    else
      throw std::runtime_error("unknown command '" + command + "'");
    // Stage-7 products can outlive the build directory in which they were
    // created.  Embed the same source/build identity used by transport output
    // so a time integral or response fold is never attributed by filename.
    auto& root = result.AsObject("Stage-7 result");
    root["build_provenance"] = shield_sim::Json(shield_sim::Json::Object{
        {"build_type", shield_sim_build::BuildType},
        {"code_revision", shield_sim_build::GitRevision},
        {"compiler", shield_sim_build::Compiler},
        {"compiler_flags", shield_sim_build::CompilerFlags},
        {"platform", shield_sim_build::Platform},
        {"source_tree_state", shield_sim_build::GitState}});
    Write(argv[3], result);
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "shieldSimExposureTool: " << error.what() << '\n';
    return 2;
  }
}
