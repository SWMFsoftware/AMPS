#include "sep_field_line_bundle_io.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <regex>
#include <sstream>
#include <system_error>

namespace SEP { namespace FieldLine { namespace {

namespace fs = std::filesystem;

std::uint32_t RotateRight(std::uint32_t value, unsigned count) {
  return (value >> count) | (value << (32U - count));
}

std::string Sha256(const std::string& bytes) {
  static constexpr std::array<std::uint32_t, 64> rounds{{
      0x428a2f98U,0x71374491U,0xb5c0fbcfU,0xe9b5dba5U,0x3956c25bU,
      0x59f111f1U,0x923f82a4U,0xab1c5ed5U,0xd807aa98U,0x12835b01U,
      0x243185beU,0x550c7dc3U,0x72be5d74U,0x80deb1feU,0x9bdc06a7U,
      0xc19bf174U,0xe49b69c1U,0xefbe4786U,0x0fc19dc6U,0x240ca1ccU,
      0x2de92c6fU,0x4a7484aaU,0x5cb0a9dcU,0x76f988daU,0x983e5152U,
      0xa831c66dU,0xb00327c8U,0xbf597fc7U,0xc6e00bf3U,0xd5a79147U,
      0x06ca6351U,0x14292967U,0x27b70a85U,0x2e1b2138U,0x4d2c6dfcU,
      0x53380d13U,0x650a7354U,0x766a0abbU,0x81c2c92eU,0x92722c85U,
      0xa2bfe8a1U,0xa81a664bU,0xc24b8b70U,0xc76c51a3U,0xd192e819U,
      0xd6990624U,0xf40e3585U,0x106aa070U,0x19a4c116U,0x1e376c08U,
      0x2748774cU,0x34b0bcb5U,0x391c0cb3U,0x4ed8aa4aU,0x5b9cca4fU,
      0x682e6ff3U,0x748f82eeU,0x78a5636fU,0x84c87814U,0x8cc70208U,
      0x90befffaU,0xa4506cebU,0xbef9a3f7U,0xc67178f2U}};
  std::vector<std::uint8_t> data(bytes.begin(), bytes.end());
  const std::uint64_t bitLength = static_cast<std::uint64_t>(data.size()) * 8U;
  data.push_back(0x80U);
  while (data.size() % 64U != 56U) data.push_back(0U);
  for (int shift = 56; shift >= 0; shift -= 8)
    data.push_back(static_cast<std::uint8_t>(bitLength >> shift));
  std::array<std::uint32_t, 8> state{{0x6a09e667U,0xbb67ae85U,
      0x3c6ef372U,0xa54ff53aU,0x510e527fU,0x9b05688cU,0x1f83d9abU,
      0x5be0cd19U}};
  for (std::size_t block = 0; block < data.size(); block += 64U) {
    std::array<std::uint32_t, 64> words{};
    for (std::size_t i = 0; i < 16; ++i) {
      const std::size_t p = block + 4U * i;
      words[i] = (static_cast<std::uint32_t>(data[p]) << 24U) |
          (static_cast<std::uint32_t>(data[p+1]) << 16U) |
          (static_cast<std::uint32_t>(data[p+2]) << 8U) | data[p+3];
    }
    for (std::size_t i = 16; i < 64; ++i) {
      const std::uint32_t s0 = RotateRight(words[i-15],7) ^
          RotateRight(words[i-15],18) ^ (words[i-15] >> 3U);
      const std::uint32_t s1 = RotateRight(words[i-2],17) ^
          RotateRight(words[i-2],19) ^ (words[i-2] >> 10U);
      words[i] = words[i-16] + s0 + words[i-7] + s1;
    }
    std::uint32_t a=state[0],b=state[1],c=state[2],d=state[3];
    std::uint32_t e=state[4],f=state[5],g=state[6],h=state[7];
    for (std::size_t i = 0; i < 64; ++i) {
      const std::uint32_t s1=RotateRight(e,6)^RotateRight(e,11)^RotateRight(e,25);
      const std::uint32_t choose=(e&f)^((~e)&g);
      const std::uint32_t t1=h+s1+choose+rounds[i]+words[i];
      const std::uint32_t s0=RotateRight(a,2)^RotateRight(a,13)^RotateRight(a,22);
      const std::uint32_t majority=(a&b)^(a&c)^(b&c);
      const std::uint32_t t2=s0+majority;
      h=g; g=f; f=e; e=d+t1; d=c; c=b; b=a; a=t1+t2;
    }
    state[0]+=a; state[1]+=b; state[2]+=c; state[3]+=d;
    state[4]+=e; state[5]+=f; state[6]+=g; state[7]+=h;
  }
  std::ostringstream output;
  output << std::hex << std::setfill('0');
  for (std::uint32_t word : state) output << std::setw(8) << word;
  return output.str();
}

bool SafeToken(const std::string& token) {
  return !token.empty() && std::all_of(token.begin(), token.end(), [](char c) {
    return std::isalnum(static_cast<unsigned char>(c)) || c == '_' ||
        c == '-' || c == '.';
  });
}

std::string SerializeLine(const LineRecord& line) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "L\t" << line.stableLineId << '\t' << line.open << '\t'
      << line.singleSmoothSector << '\t' << static_cast<int>(line.historyCoverage)
      << '\t' << static_cast<int>(line.intersectionPresence) << '\t'
      << static_cast<int>(line.sourceMeasureStatus) << '\t'
      << line.unsignedMagneticFluxWb << '\t' << line.quadratureWeight << '\t'
      << line.observerCoverage.beginTimeS << '\t'
      << line.observerCoverage.endTimeS << '\t'
      << line.observerCoverage.complete << '\t'
      << line.connection.geometricRootCount << '\t'
      << line.connection.sourceEligibleRootCount << '\t'
      << line.connection.firstGeometricTimeS << '\t'
      << line.connection.lastGeometricTimeS << '\t'
      << line.connection.firstSourceActiveTimeS << '\t'
      << line.connection.lastSourceActiveTimeS << '\t'
      << line.connection.evaluatedThroughTimeS << '\t'
      << line.connection.rejectionCauses << '\n';
  for (const NodeState& n : line.nodes) {
    out << "N\t" << n.arcLengthM << '\t'
        << n.positionM.x << '\t' << n.positionM.y << '\t' << n.positionM.z << '\t'
        << n.outwardTangent.x << '\t' << n.outwardTangent.y << '\t'
        << n.outwardTangent.z << '\t' << n.magneticFieldT.x << '\t'
        << n.magneticFieldT.y << '\t' << n.magneticFieldT.z << '\t'
        << n.plasmaVelocityMPerS.x << '\t' << n.plasmaVelocityMPerS.y << '\t'
        << n.plasmaVelocityMPerS.z << '\t' << n.massDensityKgM3 << '\t'
        << n.pressurePa << '\t' << n.temperatureK << '\t' << n.focusingLengthM
        << '\t' << n.outwardWaveEnergyJPerM3 << '\t'
        << n.inwardWaveEnergyJPerM3 << '\t' << n.tubeAreaM2 << '\t'
        << n.region << '\t' << n.primaryTopology << '\t'
        << n.secondaryTopology << '\t' << n.magneticSector << '\t'
        << n.sourceLabel << '\t' << n.forwardLongitudeJacobian << '\t'
        << n.inverseLongitudeJacobian << '\t' << n.mappingValid << '\t'
        << n.interfaceIdentity << '\n';
  }
  for (const FrontIntersection& i : line.intersections)
    out << "I\t" << i.stableId << '\t' << i.frontGeneration << '\t'
        << i.timeS << '\t' << i.arcLengthM << '\t'
        << i.geometricIntersection << '\t' << i.sourceEligible << '\t'
        << i.rejectionCauses << '\n';
  for (const LineObserverMapping& observer : line.observers) {
    out << "O\t" << observer.observerId << '\t' << observer.detectorFrame
        << '\t' << observer.valid << '\t' << observer.energyEdgesJ.size();
    for (double edge : observer.energyEdgesJ) out << '\t' << edge;
    out << '\n';
    for (const ObserverOverlapComponent& c : observer.components)
      out << "C\t" << observer.observerId << '\t' << c.stableComponentId
          << '\t' << c.beginTimeS << '\t' << c.endTimeS << '\t'
          << c.closestArcLengthM << '\t' << c.separationM << '\t'
          << c.representedVolumeM3 << '\t' << c.exposureM3S << '\n';
  }
  return out.str();
}

Core::Result<LineRecord> ParseLine(const std::string& bytes) {
  std::istringstream input(bytes);
  std::string row;
  LineRecord line;
  std::map<std::string, std::size_t> observerIndex;
  while (std::getline(input, row)) {
    std::istringstream values(row);
    std::string type;
    values >> type;
    if (type == "L") {
      int history=0, presence=0, measure=0;
      values >> line.stableLineId >> line.open >> line.singleSmoothSector >> history
          >> presence >> measure >> line.unsignedMagneticFluxWb
          >> line.quadratureWeight >> line.observerCoverage.beginTimeS
          >> line.observerCoverage.endTimeS >> line.observerCoverage.complete
          >> line.connection.geometricRootCount
          >> line.connection.sourceEligibleRootCount
          >> line.connection.firstGeometricTimeS
          >> line.connection.lastGeometricTimeS
          >> line.connection.firstSourceActiveTimeS
          >> line.connection.lastSourceActiveTimeS
          >> line.connection.evaluatedThroughTimeS
          >> line.connection.rejectionCauses;
      line.historyCoverage = static_cast<HistoryCoverage>(history);
      line.intersectionPresence = static_cast<IntersectionPresence>(presence);
      line.sourceMeasureStatus = static_cast<SourceMeasureStatus>(measure);
    } else if (type == "N") {
      NodeState n;
      values >> n.arcLengthM >> n.positionM.x >> n.positionM.y >> n.positionM.z
          >> n.outwardTangent.x >> n.outwardTangent.y >> n.outwardTangent.z
          >> n.magneticFieldT.x >> n.magneticFieldT.y >> n.magneticFieldT.z
          >> n.plasmaVelocityMPerS.x >> n.plasmaVelocityMPerS.y
          >> n.plasmaVelocityMPerS.z >> n.massDensityKgM3 >> n.pressurePa
          >> n.temperatureK >> n.focusingLengthM
          >> n.outwardWaveEnergyJPerM3 >> n.inwardWaveEnergyJPerM3
          >> n.tubeAreaM2 >> n.region >> n.primaryTopology
          >> n.secondaryTopology >> n.magneticSector >> n.sourceLabel
          >> n.forwardLongitudeJacobian >> n.inverseLongitudeJacobian
          >> n.mappingValid >> n.interfaceIdentity;
      line.nodes.push_back(n);
    } else if (type == "I") {
      FrontIntersection i;
      values >> i.stableId >> i.frontGeneration >> i.timeS >> i.arcLengthM
          >> i.geometricIntersection >> i.sourceEligible >> i.rejectionCauses;
      line.intersections.push_back(i);
    } else if (type == "O") {
      LineObserverMapping observer;
      std::size_t edgeCount = 0;
      values >> observer.observerId >> observer.detectorFrame >> observer.valid
          >> edgeCount;
      observer.energyEdgesJ.resize(edgeCount);
      for (double& edge : observer.energyEdgesJ) values >> edge;
      observerIndex[observer.observerId] = line.observers.size();
      line.observers.push_back(observer);
    } else if (type == "C") {
      std::string observerId;
      ObserverOverlapComponent c;
      values >> observerId >> c.stableComponentId >> c.beginTimeS >> c.endTimeS
          >> c.closestArcLengthM >> c.separationM >> c.representedVolumeM3
          >> c.exposureM3S;
      const auto found = observerIndex.find(observerId);
      if (found == observerIndex.end())
        return Core::Result<LineRecord>::Failure(
            Core::StatusCode::DataIntegrityFailure,
            "overlap component precedes its observer record");
      line.observers[found->second].components.push_back(c);
    } else {
      return Core::Result<LineRecord>::Failure(
          Core::StatusCode::DataIntegrityFailure, "unknown bundle member row");
    }
    if (!values)
      return Core::Result<LineRecord>::Failure(
          Core::StatusCode::DataIntegrityFailure, "malformed bundle member row");
  }
  auto status = ValidateLine(line);
  if (!status.ok())
    return Core::Result<LineRecord>::Failure(status.code, status.message);
  return Core::Result<LineRecord>::Success(line);
}

std::string IdentityBytes(const FieldLineSet& set) {
  std::vector<LineRecord> lines = set.lines;
  std::sort(lines.begin(), lines.end(), [](const LineRecord& a,
                                          const LineRecord& b) {
    return a.stableLineId < b.stableLineId;
  });
  std::ostringstream bytes;
  bytes << std::setprecision(17) << set.schemaVersion << '\n' << set.frame
        << '\n' << set.epochS << '\n' << set.historyBeginS << '\n'
        << set.historyEndS << '\n' << set.rotationProvenance << '\n';
  for (const LineRecord& line : lines) bytes << SerializeLine(line);
  return bytes.str();
}

std::string ExtractString(const std::string& json, const std::string& key) {
  std::smatch match;
  const std::regex expression("\\\"" + key + "\\\":\\\"([^\\\"]*)\\\"");
  return std::regex_search(json, match, expression) ? match[1].str() : "";
}
double ExtractNumber(const std::string& json, const std::string& key,
                     double missing = std::numeric_limits<double>::quiet_NaN()) {
  std::smatch match;
  const std::regex expression("\\\"" + key +
      "\\\":([-+0-9.eE]+)");
  return std::regex_search(json, match, expression)
      ? std::stod(match[1].str()) : missing;
}

}  // namespace

std::string BundleIdentity(const FieldLineSet& set) {
  return Sha256(IdentityBytes(set));
}

Core::Status WriteBundleTransactional(const FieldLineSet& set,
                                      const std::string& targetDirectory) {
  auto valid = ValidateFieldLineSet(set);
  if (!valid.ok()) return valid;
  if (set.bundleId != BundleIdentity(set))
    return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                "bundle ID does not match canonical content");
  for (const auto& line : set.lines)
    if (!SafeToken(line.stableLineId))
      return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                  "line ID is not a safe serialized token");
  fs::path target(targetDirectory);
  fs::path temporary = target;
  temporary += ".tmp-" + set.bundleId.substr(0, 12);
  std::error_code error;
  if (fs::exists(target) || fs::exists(temporary))
    return Core::Status::Failure(Core::StatusCode::InvalidState,
                                "bundle target/transaction path already exists");
  if (!fs::create_directories(temporary, error) || error)
    return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                "cannot create bundle transaction directory");
  std::vector<LineRecord> lines = set.lines;
  std::sort(lines.begin(), lines.end(), [](const LineRecord& a,
                                          const LineRecord& b) {
    return a.stableLineId < b.stableLineId;
  });
  std::ostringstream members;
  bool first = true;
  for (const auto& line : lines) {
    const std::string filename = "line_" + line.stableLineId + ".tsv";
    const std::string bytes = SerializeLine(line);
    std::ofstream output(temporary / filename, std::ios::binary);
    output << bytes;
    output.close();
    if (!output) {
      fs::remove_all(temporary, error);
      return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                  "failed to write a bundle member");
    }
    if (!first) members << ',';
    first = false;
    members << "{\"line_id\":\"" << line.stableLineId
            << "\",\"file\":\"" << filename << "\",\"sha256\":\""
            << Sha256(bytes) << "\"}";
  }
  std::ostringstream manifest;
  manifest << std::setprecision(17)
      << "{\"schema_version\":" << set.schemaVersion
      << ",\"bundle_id\":\"" << set.bundleId
      << "\",\"frame\":\"" << set.frame
      << "\",\"epoch_s\":" << set.epochS
      << ",\"history_begin_s\":" << set.historyBeginS
      << ",\"history_end_s\":" << set.historyEndS
      << ",\"rotation_provenance\":\"" << set.rotationProvenance
      << "\",\"members\":[" << members.str() << "]}\n";
  std::ofstream manifestOutput(temporary / "manifest.json", std::ios::binary);
  manifestOutput << manifest.str();
  manifestOutput.close();
  if (!manifestOutput) {
    fs::remove_all(temporary, error);
    return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                "failed to write canonical bundle manifest");
  }
  fs::rename(temporary, target, error);
  if (error) {
    fs::remove_all(temporary, error);
    return Core::Status::Failure(Core::StatusCode::DataIntegrityFailure,
                                "atomic bundle publication failed");
  }
  return Core::Status::Success();
}

Core::Result<FieldLineSet> ReadBundle(const std::string& directory) {
  const fs::path root(directory);
  std::ifstream manifestInput(root / "manifest.json", std::ios::binary);
  std::ostringstream manifestBytes;
  manifestBytes << manifestInput.rdbuf();
  if (!manifestInput)
    return Core::Result<FieldLineSet>::Failure(
        Core::StatusCode::DataIntegrityFailure, "bundle manifest is unreadable");
  const std::string manifest = manifestBytes.str();
  FieldLineSet set;
  set.schemaVersion = static_cast<int>(ExtractNumber(manifest, "schema_version"));
  set.bundleId = ExtractString(manifest, "bundle_id");
  set.frame = ExtractString(manifest, "frame");
  set.epochS = ExtractNumber(manifest, "epoch_s");
  set.historyBeginS = ExtractNumber(manifest, "history_begin_s");
  set.historyEndS = ExtractNumber(manifest, "history_end_s");
  set.rotationProvenance = ExtractString(manifest, "rotation_provenance");
  if (set.schemaVersion != 1)
    return Core::Result<FieldLineSet>::Failure(
        Core::StatusCode::UnsupportedCapability, "unsupported bundle schema");
  const std::regex memberExpression(
      "\\{\\\"line_id\\\":\\\"([^\\\"]+)\\\",\\\"file\\\":\\\"([^\\\"]+)\\\",\\\"sha256\\\":\\\"([0-9a-f]{64})\\\"\\}");
  for (std::sregex_iterator it(manifest.begin(), manifest.end(), memberExpression),
       end; it != end; ++it) {
    const std::string lineId = (*it)[1].str();
    const std::string filename = (*it)[2].str();
    const std::string checksum = (*it)[3].str();
    if (!SafeToken(lineId) || filename != "line_" + lineId + ".tsv")
      return Core::Result<FieldLineSet>::Failure(
          Core::StatusCode::DataIntegrityFailure, "unsafe/mismatched member identity");
    std::ifstream memberInput(root / filename, std::ios::binary);
    std::ostringstream memberBytes;
    memberBytes << memberInput.rdbuf();
    if (!memberInput || Sha256(memberBytes.str()) != checksum)
      return Core::Result<FieldLineSet>::Failure(
          Core::StatusCode::DataIntegrityFailure, "bundle member checksum failed");
    auto line = ParseLine(memberBytes.str());
    if (!line.ok() || line.value.stableLineId != lineId)
      return Core::Result<FieldLineSet>::Failure(
          Core::StatusCode::DataIntegrityFailure, "bundle member validation failed");
    set.lines.push_back(line.value);
  }
  if (set.lines.empty())
    return Core::Result<FieldLineSet>::Failure(
        Core::StatusCode::DataIntegrityFailure, "bundle manifest has no members");
  auto valid = ValidateFieldLineSet(set);
  if (!valid.ok() || BundleIdentity(set) != set.bundleId)
    return Core::Result<FieldLineSet>::Failure(
        Core::StatusCode::DataIntegrityFailure,
        "bundle validation or canonical identity failed");
  return Core::Result<FieldLineSet>::Success(set);
}

} }  // namespace SEP::FieldLine
