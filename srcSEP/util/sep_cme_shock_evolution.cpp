#include "sep_cme_shock_evolution.h"
#include "sep_transactional_output.h"
#include "provider.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace SEP::Testing {
namespace {
namespace Front = CoronaSwcme::ShockFront;
using CoronalCME::Vec3;
using std::filesystem::path;

constexpr int stepS = 1200;
constexpr int endS = 18000;
constexpr const char* legacySchema = "cme-shock-evolution-v1";
constexpr const char* eventSchema = "cme-shock-evolution-v2";

std::string ReadFile(const path& name) {
  std::ifstream stream(name, std::ios::binary);
  if (!stream) throw std::runtime_error("cannot open " + name.string());
  return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

struct Prepared {
  path ambientPath;
  std::string schema, referenceUtc, issueUtc, observationsSha256, observationsPath;
  double radiusM, speedMPerS, accelerationMPerS2, widthRad;
  Vec3 axis;
  double fitStartS, traceStartS;
  double handoffS, handoffMinM, windMPerS, dragPerM, horizonS;
  bool dragEnabled = false;
  std::vector<Vec3> rays;
};

Prepared LoadPrepared(const path& manifest, const std::string& bytes) {
  // This is a versioned Python-to-C++ test contract, in HCI and SI. It
  // carries coefficients from the fit and the upstream path, without any
  // sampled state or expected solution. The native process owns all queries.
  std::istringstream input(bytes);
  input.imbue(std::locale::classic());
  Prepared p{};
  std::string ambientName;
  std::size_t rays = 0;
  if (!(input >> p.schema >> p.referenceUtc) ||
      (p.schema != legacySchema && p.schema != eventSchema)) {
    throw std::runtime_error("invalid prepared CME event schema");
  }
  const bool legacy = p.schema == legacySchema;
  if (!legacy && (!(input >> p.issueUtc >> p.observationsSha256 >> std::quoted(p.observationsPath)) ||
      p.observationsPath.empty() || p.observationsSha256.size() != 64 ||
      p.observationsSha256.find_first_not_of("0123456789abcdef") != std::string::npos)) {
    throw std::runtime_error("invalid prepared observation provenance");
  }
  if (!(input >> std::quoted(ambientName)) || ambientName.empty() ||
      (!legacy && !(input >> p.fitStartS >> p.handoffS >> p.traceStartS)) ||
      !(input >> p.radiusM >> p.speedMPerS >> p.accelerationMPerS2
        >> p.axis.x >> p.axis.y >> p.axis.z >> p.widthRad)) {
    throw std::runtime_error("invalid prepared CME fit or upstream path");
  }
  if (legacy) {
    p.dragEnabled = true;
    if (!(input >> p.handoffS >> p.handoffMinM >> p.windMPerS >> p.dragPerM >> p.horizonS))
      throw std::runtime_error("invalid legacy CME drag record");
  } else {
    std::string propagation;
    input >> propagation;
    p.dragEnabled = propagation == "quadratic-drag";
    if ((!p.dragEnabled && propagation != "fitted-only") ||
        (p.dragEnabled && !(input >> p.windMPerS >> p.dragPerM >> p.horizonS >> p.handoffMinM))) {
      throw std::runtime_error("invalid prepared CME propagation record");
    }
  }
  if (!(input >> rays) || rays == 0 || (legacy && rays != 49)) {
    throw std::runtime_error("invalid prepared CME ray count");
  }
  p.ambientPath = manifest.parent_path() / ambientName;
  for (std::size_t i = 0; i < rays; ++i) {
    Vec3 ray;
    if (!(input >> ray.x >> ray.y >> ray.z) ||
        !(std::abs(CoronalCME::Norm(ray) - 1.0) <= (legacy ? 1e-12 : 1e-10))) {
      throw std::runtime_error("prepared CME event requires unit HCI sample directions");
    }
    p.rays.push_back(ray);
  }
  std::string trailing;
  // Negative acceleration is allowed if the measured branch advances at
  // both ends. Positive speed throughout follows from linearity of V(t).
  const double neededOuterS = p.traceStartS + endS - p.handoffS;
  const bool outerValid = p.dragEnabled && p.handoffMinM > 0 && p.windMPerS > 0 &&
      p.dragPerM >= 0 && p.horizonS > 0 && p.horizonS >= neededOuterS;
  if (input >> trailing || !(p.radiusM > 0 && p.speedMPerS > 0
      && std::isfinite(p.accelerationMPerS2) && p.fitStartS >= 0 && p.handoffS > p.fitStartS
      && p.traceStartS >= p.fitStartS && p.speedMPerS + p.accelerationMPerS2 * p.handoffS > 0
      && (neededOuterS <= 0 || outerValid)) || (legacy && !(p.handoffS < endS))) {
    throw std::runtime_error("prepared CME event does not cover the five-hour trace");
  }
  return p;
}

Front::SseKinematicState Trace(const Prepared& p, double timeS) {
  // The manufactured event tests the same constant-acceleration curve used
  // by the Python fitter. Evaluate it through the last measured epoch. Use
  // the analytically matched outer drag law afterward. R and V are
  // continuous at the handoff while acceleration may change.
  const double fitTime = std::min<double>(timeS, p.handoffS);
  const auto fitted = Front::EvaluateFittedApex(
      fitTime, p.radiusM, p.speedMPerS, p.accelerationMPerS2);
  Front::SseKinematicState state{};
  state.apexRadiusM = fitted.apexRadiusM;
  state.apexSpeedMPerS = fitted.apexSpeedMPerS;
  state.direction = p.axis;
  state.halfWidthRad = p.widthRad;
  if (timeS > p.handoffS) {
    if (state.apexRadiusM < p.handoffMinM) {
      throw std::runtime_error("fitted CME failed the declared minimum drag handoff radius");
    }
    const auto outer = Front::EvaluateQuadraticDrag(
        timeS, p.handoffS, state.apexRadiusM, state.apexSpeedMPerS, p.windMPerS, p.dragPerM);
    if (!outer.ok()) throw std::runtime_error(outer.status.message);
    state.apexRadiusM = outer.value.apexRadiusM;
    state.apexSpeedMPerS = outer.value.apexSpeedMPerS;
  }
  return state;
}

std::time_t UtcSeconds(const std::string& epoch) {
  // Whole-second civil UTC matches the ambient contract. Round trip to
  // reject calendar normalization, including invalid days and leap seconds.
  // As with the fitter, spanning a leap second requires an external table.
  std::istringstream input(epoch);
  std::tm value{};
  input.imbue(std::locale::classic());
  if (!(input >> std::get_time(&value, "%Y-%m-%dT%H:%M:%SZ"))) {
    throw std::runtime_error("invalid prepared UTC reference epoch");
  }
  const auto seconds = timegm(&value);
  std::tm normalized{};
  if (!gmtime_r(&seconds, &normalized)) throw std::runtime_error("UTC reference is out of range");
  std::ostringstream roundTrip;
  roundTrip.imbue(std::locale::classic());
  roundTrip << std::put_time(&normalized, "%Y-%m-%dT%H:%M:%SZ");
  if (roundTrip.str() != epoch) throw std::runtime_error("invalid prepared UTC reference epoch");
  return seconds;
}
}  // namespace

Result RunCmeShockEvolution() {
  Result result;
  const auto& context = GetExecutionContext();
  if (context.inputPath.empty() || context.artifactDirectory.empty()) {
    // Extended discovery can include this case in a mixed catalog. A real
    // trace needs its own prepared event and isolated output directory.
    result.status = Status::Skip;
    result.message = "select CME individually after preparing the event with regression.py";
    return result;
  }

  // Resolve the upstream archive's harmonic checksum and its UTC origin
  // using the canonical event reader.
  const auto preparedBytes = ReadFile(context.inputPath);
  const auto prepared = LoadPrepared(context.inputPath, preparedBytes);
  const auto cfg = Front::ResolveConfiguration(ReadFile(prepared.ambientPath),
      [&](const std::string& asset) {
        try {
          return Core::Result<std::string>::Success(
              ReadFile(prepared.ambientPath.parent_path() / asset));
        } catch (const std::exception& failure) {
          return Core::Result<std::string>::Failure(Core::StatusCode::DataIntegrityFailure,
                                                     failure.what());
        }
      }, prepared.ambientPath.string());
  if (!cfg.ok()) throw std::runtime_error(cfg.status.message);
  const auto& model = *cfg.value;
  const bool legacy = prepared.schema == legacySchema;
  if (legacy && model.referenceEpoch != prepared.referenceUtc) {
    throw std::runtime_error("fit and upstream must have the same UTC reference in this regression");
  }
  const auto ambient = CoronaSwcme::AmbientModel::Create(model.ambient);
  if (!ambient.ok()) throw std::runtime_error(ambient.status.message);
  const double ambientOffsetS = legacy ? 0.0 :
      std::difftime(UtcSeconds(prepared.referenceUtc), UtcSeconds(model.referenceEpoch));
  const double lastEpochS = prepared.traceStartS + endS;
  if (prepared.traceStartS + ambientOffsetS < model.ambient.support.startS ||
      lastEpochS + ambientOffsetS > model.ambient.support.endS) {
    throw std::runtime_error("upstream event does not cover the five-hour trace in UTC");
  }
  // Validate the handoff before publishing the first snapshot.
  Trace(prepared, lastEpochS);
  const auto fingerprint = Output::Checksum64(model.physicsFingerprint + '\n' + preparedBytes);
  result.configuration = {"schema=" + prepared.schema, "coordinate_frame=HCI",
      "reference_utc=" + prepared.referenceUtc, "ambient=" + ambient.value->Identity(),
      "prepared_event=" + fingerprint, "shock_step_s=" + std::to_string(stepS)};
  if (!legacy) {
    result.configuration.push_back("issue_utc=" + prepared.issueUtc);
    result.configuration.push_back("observations_sha256=" + prepared.observationsSha256);
    result.configuration.push_back("observation_data=" + prepared.observationsPath);
    result.configuration.push_back("ambient_epoch_offset_s=" + std::to_string(ambientOffsetS));
  }
  int totalFrames = 0;
  double maxJumpResidual = 0.0;
  // An initial map plus fifteen completed intervals. The integer clock gives
  // exact 20-minute epochs and never touches the PIC or particle-sampling clock.
  for (int offsetS = 0; offsetS <= endS; offsetS += stepS) {
    const double timeS = prepared.traceStartS + offsetS;
    const auto state = Trace(prepared, timeS);
    std::ostringstream csv;
    csv.imbue(std::locale::classic());
    csv << std::setprecision(17)
        << "elapsed_s,ray_id,x_m,y_m,z_m,n_x,n_y,n_z,normal_speed_m_s,status,"
           "fast_mach,density_compression,theta_bn_rad,max_rh_residual\n";
    for (std::size_t i = 0; i < prepared.rays.size(); ++i) {
      // Solve the leading root and normal speed on a fixed HCI ray. The jump
      // uses the local upstream flow relative to the moving shock. No plasma
      // property is carried over from an earlier snapshot.
      const auto geometry = Front::EvaluateSseRay(state, prepared.rays[i], i);
      if (!geometry.ok()) throw std::runtime_error("shock ray " + std::to_string(i) +
          " at elapsed s=" + std::to_string(timeS) + ": " + geometry.status.message);
      const auto& g = geometry.value;
      const auto upstream = ambient.value->Evaluate(g.positionM, timeS + ambientOffsetS);
      if (!upstream.ok()) throw std::runtime_error(upstream.status.message);
      const auto& p = upstream.value;
      const auto jump = Front::EvaluateLocalJump(
          {p.plasma.massDensityKgM3, p.plasma.pressurePa, p.velocityMPerS, p.magneticFieldT},
          g.outwardNormal, g.normalSpeedMPerS, model.ambient.composition.gammaAdiabatic,
          model.weakMachTolerance, model.rhResidualTolerance,
          model.ambient.ambient.minimumMagneticFieldT);
      if (!jump.ok()) throw std::runtime_error(jump.status.message);
      const auto& j = jump.value;
      if ((!j.downstreamValid && j.status != Front::FrontStatus::SubfastFront &&
           j.status != Front::FrontStatus::NonForwardInflow) ||
          (legacy && i == 0 && !j.downstreamValid)) {
        result.status = Status::Fail;
        result.message = "unexpected shock classification at t=" + std::to_string(timeS) +
            " ray=" + std::to_string(i) + " " + Front::Name(j.status) + " " + j.reason;
        return result;
      }

      const auto x = g.positionM;
      const auto n = g.outwardNormal;
      csv << timeS << ',' << i << ',' << x.x << ',' << x.y << ',' << x.z << ','
          << n.x << ',' << n.y << ',' << n.z << ',' << g.normalSpeedMPerS << ','
          << Front::Name(j.status) << ',' << j.fastMach << ',';
      // Blank downstream columns on physical no-shock branches preserve
      // the distinction between absent states and a valid zero.
      if (j.downstreamValid) csv << j.jump.compressionRatio;
      csv << ',';
      if (j.magneticDirectionValid) csv << j.characteristics.obliquityRad;
      csv << ',';
      if (j.downstreamValid) {
        csv << j.jump.residuals.maximum;
        maxJumpResidual = std::max(maxJumpResidual, j.jump.residuals.maximum);
      }
      csv << '\n';
    }
    std::ostringstream name;
    name << "shock-" << std::setw(5) << std::setfill('0') << offsetS << ".csv";
    const auto write = Output::WriteTransactional(context.artifactDirectory,
        name.str(), csv.str(), {prepared.schema, prepared.rays.size(), fingerprint});
    if (!write.status.ok()) throw std::runtime_error(write.status.message);
    result.artifacts.push_back(write.finalPath);
    ++totalFrames;
  }
  result.metrics.emplace_back("frames", totalFrames, endS / stepS + 1, "==", "count");
  result.metrics.emplace_back("last_epoch", lastEpochS, lastEpochS, "==", "s");
  result.metrics.emplace_back("max_jump_residual", maxJumpResidual, model.rhResidualTolerance, "<=", "1");
  result.status = Status::Pass;
  result.message = "wrote the initial shock and all fifteen 20-minute steps through five hours";
  return result;
}
}  // namespace SEP::Testing
