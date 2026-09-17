#include "publication.h"

#include <filesystem>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

namespace SEP3D {
namespace Output {
namespace {

namespace fs = std::filesystem;

Core::Status Error(const std::string& message) {
  return Core::Status(Core::StatusCode::Error, message);
}

bool SafeToken(const std::string& value) {
  return !value.empty() && value.find_first_of("\r\n,=") == std::string::npos;
}

bool WriteCells(const fs::path& path, const SamplingSnapshot& snapshot) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out << "cell_id,species,represented_particles,number_density_m-3,"
         "flux_x_m-2_s-1,flux_y_m-2_s-1,flux_z_m-2_s-1,"
         "kinetic_energy_density_J_m-3,first_pitch_moment\n";
  out << std::scientific << std::setprecision(17);
  for (const CellMoment& m : snapshot.cellMoments)
    out << m.cellId << ',' << m.species << ',' << m.representedParticles << ','
        << m.numberDensityM3 << ',' << m.weightedFluxM2PerS.x << ','
        << m.weightedFluxM2PerS.y << ',' << m.weightedFluxM2PerS.z << ','
        << m.kineticEnergyDensityJPerM3 << ',' << m.firstPitchMoment << '\n';
  out.close(); return out.good();
}

bool WriteSpacecraft(const fs::path& path, const SamplingSnapshot& snapshot) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out << "spacecraft,species,energy_min_J,energy_max_J,"
         "represented_particles_J-1,dipole_anisotropy\n";
  out << std::scientific << std::setprecision(17);
  for (const VirtualSpacecraftProduct& p : snapshot.spacecraft)
    for (std::size_t i = 0; i < p.representedParticlesPerJ.size(); ++i)
      out << p.name << ',' << p.species << ',' << p.kineticEnergyEdgesJ[i]
          << ',' << p.kineticEnergyEdgesJ[i + 1] << ','
          << p.representedParticlesPerJ[i] << ',' << p.dipoleAnisotropy << '\n';
  out.close(); return out.good();
}

bool WriteFieldLines(const fs::path& path, const SamplingSnapshot& snapshot) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out << "projection,species,distance_min_m,distance_max_m,"
         "represented_particles_m-1\n";
  out << std::scientific << std::setprecision(17);
  for (const FieldLineProjection& p : snapshot.fieldLines)
    for (std::size_t i = 0; i < p.representedParticlesPerM.size(); ++i)
      out << p.name << ',' << p.species << ',' << p.distanceEdgesM[i]
          << ',' << p.distanceEdgesM[i + 1] << ','
          << p.representedParticlesPerM[i] << '\n';
  out.close(); return out.good();
}

bool WriteShocks(const fs::path& path, const SamplingSnapshot& snapshot) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  out << "step,species,injected_count,escaped_count,absorbed_count,"
         "failed_count,shock_crossings_count\n";
  for (const ShockDiagnostic& s : snapshot.shocks)
    out << s.step << ',' << s.species << ',' << s.injected << ',' << s.escaped
        << ',' << s.absorbed << ',' << s.failed << ',' << s.crossings << '\n';
  out.close(); return out.good();
}

Core::Status CheckHeader(const fs::path& path, const std::string& expected) {
  std::ifstream input(path, std::ios::binary);
  std::string header;
  if (!std::getline(input, header) || header != expected)
    return Error("publication CSV schema mismatch: " + path.string());
  return Core::Status::OK();
}

}  // namespace

std::string HashFileFNV1a64(const std::string& path, Core::Status* status) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    if (status) *status = Error("cannot open artifact for hashing: " + path);
    return {};
  }
  std::uint64_t hash = UINT64_C(14695981039346656037);
  char buffer[8192];
  while (input) {
    input.read(buffer, sizeof(buffer));
    const std::streamsize count = input.gcount();
    for (std::streamsize i = 0; i < count; ++i) {
      hash ^= static_cast<unsigned char>(buffer[i]);
      hash *= UINT64_C(1099511628211);
    }
  }
  if (!input.eof()) {
    if (status) *status = Error("artifact read failed while hashing: " + path);
    return {};
  }
  std::ostringstream text;
  text << std::hex << std::setw(16) << std::setfill('0') << hash;
  if (status) *status = Core::Status::OK();
  return text.str();
}

PublicationResult Publish(const std::string& outputDirectory,
                          const std::string& prefix,
                          const PublicationMetadata& metadata,
                          const SamplingSnapshot& snapshot) {
  PublicationResult result;
  if (!snapshot.status.ok() || !SafeToken(prefix) ||
      !SafeToken(metadata.configurationFingerprint) ||
      !SafeToken(metadata.codeIdentity) ||
      !SafeToken(metadata.snapshotFingerprint) ||
      !std::isfinite(metadata.simulationTimeS) ||
      metadata.snapshotGeneration == 0) {
    result.status = Error("publication metadata or sampling snapshot is invalid");
    return result;
  }
  std::ostringstream leaf;
  leaf << prefix << '.' << std::setw(8) << std::setfill('0') << metadata.sequence;
  const fs::path root(outputDirectory);
  const fs::path final = root / leaf.str();
  const fs::path staging = root / (leaf.str() + ".staging");
  std::error_code ec;
  fs::create_directories(root, ec);
  if (ec || fs::exists(staging) || fs::exists(final)) {
    result.status = Error("publication destination already exists or is unavailable");
    return result;
  }
  if (!fs::create_directory(staging, ec) || ec) {
    result.status = Error("cannot create publication staging directory");
    return result;
  }

  const std::vector<std::pair<std::string, bool>> written = {
      {"cells.csv", WriteCells(staging / "cells.csv", snapshot)},
      {"spacecraft.csv", WriteSpacecraft(staging / "spacecraft.csv", snapshot)},
      {"field_lines.csv", WriteFieldLines(staging / "field_lines.csv", snapshot)},
      {"shocks.csv", WriteShocks(staging / "shocks.csv", snapshot)}};
  for (const auto& file : written) {
    if (!file.second) {
      fs::remove_all(staging, ec);
      result.status = Error("failed writing staged output artifact: " + file.first);
      return result;
    }
    Core::Status hashStatus;
    const std::string hash = HashFileFNV1a64(
        (staging / file.first).string(), &hashStatus);
    if (!hashStatus.ok()) {
      fs::remove_all(staging, ec); result.status = hashStatus; return result;
    }
    result.artifactHashes[file.first] = hash;
  }

  std::ofstream manifest(staging / "manifest.txt",
                         std::ios::binary | std::ios::trunc);
  manifest << "schema=srcSEP3D-output-v1\n"
           << "sequence=" << metadata.sequence << '\n'
           << std::scientific << std::setprecision(17)
           << "simulation_time_s=" << metadata.simulationTimeS << '\n'
           << "snapshot_generation=" << metadata.snapshotGeneration << '\n'
           << "configuration_fingerprint="
           << metadata.configurationFingerprint << '\n'
           << "code_identity=" << metadata.codeIdentity << '\n'
           << "snapshot_fingerprint=" << metadata.snapshotFingerprint << '\n'
           << "sampling_completed=" << snapshot.nextState.completedSamplings << '\n'
           << "observations_processed="
           << snapshot.nextState.observationsProcessed << '\n';
  for (const auto& artifact : result.artifactHashes)
    manifest << "artifact=" << artifact.first << ',' << artifact.second << '\n';
  manifest.close();
  if (!manifest.good()) {
    fs::remove_all(staging, ec);
    result.status = Error("failed writing staged output manifest");
    return result;
  }
  fs::rename(staging, final, ec);
  if (ec) {
    fs::remove_all(staging, ec);
    result.status = Error("atomic publication rename failed");
    return result;
  }
  result.directory = final.string();
  result.status = Core::Status::OK();
  return result;
}

Core::Status ParseAndVerifyPublication(const std::string& directory,
                                       ParsedPublication* output) {
  if (output == nullptr) return Error("publication parser output is null");
  const fs::path root(directory);
  std::ifstream manifest(root / "manifest.txt", std::ios::binary);
  if (!manifest) return Error("publication manifest is absent");
  ParsedPublication candidate;
  std::map<std::string, std::string> values;
  std::string line;
  while (std::getline(manifest, line)) {
    const std::size_t equals = line.find('=');
    if (equals == std::string::npos) return Error("malformed manifest record");
    const std::string key = line.substr(0, equals);
    const std::string value = line.substr(equals + 1);
    if (key == "artifact") {
      const std::size_t comma = value.find(',');
      if (comma == std::string::npos ||
          !candidate.artifactHashes.emplace(
              value.substr(0, comma), value.substr(comma + 1)).second)
        return Error("malformed or duplicate artifact manifest record");
    } else if (!values.emplace(key, value).second) {
      return Error("duplicate manifest key");
    }
  }
  const char* required[] = {"schema", "sequence", "simulation_time_s",
      "snapshot_generation", "configuration_fingerprint", "code_identity",
      "snapshot_fingerprint", "sampling_completed", "observations_processed"};
  for (const char* key : required)
    if (values.find(key) == values.end()) return Error("manifest key is absent");
  if (values["schema"] != "srcSEP3D-output-v1" ||
      candidate.artifactHashes.size() != 4)
    return Error("publication schema or artifact count is unsupported");
  try {
    candidate.metadata.sequence = std::stoull(values["sequence"]);
    candidate.metadata.simulationTimeS = std::stod(values["simulation_time_s"]);
    candidate.metadata.snapshotGeneration =
        std::stoull(values["snapshot_generation"]);
  } catch (const std::exception&) {
    return Error("publication manifest number is invalid");
  }
  candidate.metadata.configurationFingerprint =
      values["configuration_fingerprint"];
  candidate.metadata.codeIdentity = values["code_identity"];
  candidate.metadata.snapshotFingerprint = values["snapshot_fingerprint"];
  for (const auto& artifact : candidate.artifactHashes) {
    Core::Status hashStatus;
    const std::string observed = HashFileFNV1a64(
        (root / artifact.first).string(), &hashStatus);
    if (!hashStatus.ok()) return hashStatus;
    if (observed != artifact.second) return Error("publication artifact hash mismatch");
  }
  const std::vector<std::pair<std::string, std::string>> headers = {
      {"cells.csv", "cell_id,species,represented_particles,number_density_m-3,flux_x_m-2_s-1,flux_y_m-2_s-1,flux_z_m-2_s-1,kinetic_energy_density_J_m-3,first_pitch_moment"},
      {"spacecraft.csv", "spacecraft,species,energy_min_J,energy_max_J,represented_particles_J-1,dipole_anisotropy"},
      {"field_lines.csv", "projection,species,distance_min_m,distance_max_m,represented_particles_m-1"},
      {"shocks.csv", "step,species,injected_count,escaped_count,absorbed_count,failed_count,shock_crossings_count"}};
  for (const auto& header : headers) {
    const Core::Status checked = CheckHeader(root / header.first, header.second);
    if (!checked.ok()) return checked;
  }
  *output = candidate;
  return Core::Status::OK();
}

}  // namespace Output
}  // namespace SEP3D
