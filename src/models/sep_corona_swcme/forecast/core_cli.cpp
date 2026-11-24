#include "provider.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
namespace Front = SEP::CoronaSwcme::ShockFront;
using SEP::CoronaSwcme::AmbientModel;
using SEP::CoronalCME::Vec3;

std::string ReadFile(const std::filesystem::path& path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open " + path.string());
  std::ostringstream bytes;
  bytes << stream.rdbuf();
  if (stream.bad()) throw std::runtime_error("cannot read " + path.string());
  return bytes.str();
}

struct AmbientMember {
  std::shared_ptr<const Front::Configuration> configuration;
  std::shared_ptr<const AmbientModel> model;
};

AmbientMember LoadMember(const std::filesystem::path& path) {
  // Reuse the canonical event resolver, including the relative harmonic asset
  // checksum. The event's prescribed trajectory is independent of this replay.
  using Result = SEP::Core::Result<std::string>;
  auto config = Front::ResolveConfiguration(ReadFile(path), [&](const std::string& asset) {
    try {
      return Result::Success(ReadFile(path.parent_path() / asset));
    } catch (const std::runtime_error& error) {
      return Result::Failure(SEP::Core::StatusCode::InvalidConfiguration, error.what());
    }
  }, path.string());
  if (!config.ok()) throw std::runtime_error(config.status.message);
  auto ambient = AmbientModel::Create(config.value->ambient);
  if (!ambient.ok()) throw std::runtime_error(ambient.status.message);
  return {config.value, ambient.value};
}

void Sample(const std::vector<AmbientMember>& members, std::istream& input) {
  Front::SseKinematicState state;
  Vec3 ray;
  int member;
  double ambientEpoch;
  if (!(input >> state.apexRadiusM >> state.apexSpeedMPerS
        >> state.direction.x >> state.direction.y >> state.direction.z
        >> state.halfWidthRad >> ray.x >> ray.y >> ray.z >> member >> ambientEpoch) ||
      member < -1 || (member >= 0 && static_cast<std::size_t>(member) >= members.size()))
    throw std::runtime_error("invalid ray request");

  // The shared core finds the supported leading intersection in HCI metres
  // and takes the level-set time derivative in m/s. Never use Rdot projected
  // along a heliocentric ray as the speed of an oblique shock patch.
  auto geometry = Front::EvaluateSseRay(state, ray);
  if (!geometry.ok()) {
    if (geometry.status.code != SEP::Core::StatusCode::OutOfDomain)
      throw std::runtime_error(geometry.status.message);
    std::cout << Front::Name(Front::FrontStatus::OutsideFrontSupport) << '\n';
    return;
  }

  const auto& g = geometry.value;
  // SSE's mathematical flank may extend inside the photosphere at an early
  // epoch. Publish only the exterior map. Ambient members can declare larger
  // inner validity radii separately.
  constexpr double solarRadiusM = 695700000.0;
  if (SEP::CoronalCME::Norm(g.positionM) < solarRadiusM) {
    std::cout << Front::Name(Front::FrontStatus::BelowPhysicalInnerBoundary) << '\n';
    return;
  }
  const double missing = std::numeric_limits<double>::quiet_NaN();
  double mach = missing, compression = missing, obliquity = missing;
  std::string status = Front::Name(Front::FrontStatus::AmbientUnavailable);
  if (member >= 0) {
    const auto& ambient = members.at(member);
    // The time supplied here is elapsed from this *ambient member's* own UTC
    // reference epoch. The canonical provider owns its spatial/time domain.
    auto upstream = ambient.model->Evaluate(g.positionM, ambientEpoch);
    if (upstream.ok()) {
      const auto& p = upstream.value;
      const auto& c = *ambient.configuration;
      auto jump = Front::EvaluateLocalJump(
          {p.plasma.massDensityKgM3, p.plasma.pressurePa, p.velocityMPerS, p.magneticFieldT},
          g.outwardNormal, g.normalSpeedMPerS, c.ambient.composition.gammaAdiabatic,
          c.weakMachTolerance, c.rhResidualTolerance,
          c.ambient.ambient.minimumMagneticFieldT);
      if (!jump.ok()) throw std::runtime_error(jump.status.message);
      const auto& j = jump.value;
      status = Front::Name(j.status);
      mach = j.fastMach;
      if (j.magneticDirectionValid) obliquity = j.characteristics.obliquityRad;
      if (j.downstreamValid) compression = j.jump.compressionRatio;
    } else if (upstream.status.code != SEP::Core::StatusCode::OutOfDomain) {
      // Preserve the distinction between a gap in declared coverage and a
      // provider failure. Both remain unresolved when counting probabilities.
      status = "ambient-model-error";
      std::cerr << "ambient member " << member << ": " << upstream.status.message << '\n';
    }
  }
  // Missing plasma quantities stay missing through the wire format. The
  // caller can still issue a geometric map from a fully supported ray.
  std::cout << status << ' ' << g.positionM.x << ' ' << g.positionM.y << ' '
            << g.positionM.z << ' ' << g.outwardNormal.x << ' ' << g.outwardNormal.y
            << ' ' << g.outwardNormal.z << ' ' << g.normalSpeedMPerS << ' '
            << mach << ' ' << compression << ' ' << obliquity << '\n';
}

Front::TrajectoryState DragState(double elapsedS, double radiusM, double speedMPerS,
                                 double windMPerS, double gammaPerM) {
  // Arrival and time-series queries start from the same fitted handoff state
  // at elapsed zero. The shared core preserves continuity of R and dR/dt.
  auto state = Front::EvaluateQuadraticDrag(elapsedS, 0, radiusM, speedMPerS,
                                            windMPerS, gammaPerM);
  if (!state.ok()) throw std::runtime_error(state.status.message);
  return state.value;
}

void SampleDrag(std::istream& input) {
  double elapsed, r, v, wind, drag;
  if (!(input >> elapsed >> r >> v >> wind >> drag))
    throw std::runtime_error("invalid drag state request");
  const auto state = DragState(elapsed, r, v, wind, drag);
  std::cout << "DRAG_STATE " << state.apexRadiusM << ' ' << state.apexSpeedMPerS << '\n';
}

void Arrival(std::istream& input) {
  double r, v, wind, drag, target, horizon;
  if (!(input >> r >> v >> wind >> drag >> target >> horizon) ||
      !std::isfinite(target) || !std::isfinite(horizon) ||
      v <= 0 || wind <= 0 || target <= r || horizon <= 0)
    throw std::runtime_error("invalid arrival request");
  // The request is an apex-radius crossing derived from a supported ray with
  // fixed direction and half width. Positive apex/wind speeds make the shared
  // drag history monotonic. No search is performed beyond the declared horizon.
  auto at = [&](double t) { return DragState(t, r, v, wind, drag); };
  if (at(horizon).apexRadiusM < target) {
    std::cout << "BEYOND_FORECAST_HORIZON\n";
    return;
  }
  double lo = 0, hi = horizon;
  for (int i = 0; i < 64 && hi - lo > 1e-3; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (at(mid).apexRadiusM < target) lo = mid;
    else hi = mid;
  }
  std::cout << "ARRIVAL " << hi << ' ' << at(hi).apexSpeedMPerS << '\n';
}

void DescribeMember(const std::vector<AmbientMember>& members, std::istream& input) {
  std::size_t index;
  if (!(input >> index) || index >= members.size())
    throw std::runtime_error("invalid ambient member index");
  // The validated reader owns the UTC epoch string and verifies the hash of
  // supporting assets. The caller resolves the time offset and records this
  // same physical model identity in the output receipt.
  const auto& member = *members.at(index).configuration;
  std::cout << member.referenceEpoch << ' ' << member.physicsFingerprint << '\n';
}
} // namespace

int main(int argc, char** argv) {
  try {
    // All members are loaded once. Each request is pure and line-local, so
    // malformed or truncated input cannot shift the fields of a later query.
    std::vector<AmbientMember> members;
    for (int i = 1; i < argc; ++i) members.push_back(LoadMember(argv[i]));
    std::cout << std::setprecision(17);
    std::string line;
    while (std::getline(std::cin, line)) {
      std::istringstream input(line);
      char command;
      if (!(input >> command)) throw std::runtime_error("missing command");
      if (command == 'R') Sample(members, input);
      else if (command == 'A') Arrival(input);
      else if (command == 'T') SampleDrag(input);
      else if (command == 'M') DescribeMember(members, input);
      else throw std::runtime_error("unknown command");
      std::string extra;
      if (input >> extra) throw std::runtime_error("extra request fields");
    }
    if (!std::cin.eof() || !std::cout) throw std::runtime_error("forecast core I/O failure");
  } catch (const std::exception& error) {
    std::cerr << "shock forecast core: " << error.what() << '\n';
    return 1;
  }
}
