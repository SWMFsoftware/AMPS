#include "sep_coronal_cme/background_event.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

#include "sep_coronal_cme/clipped_ellipsoid_quadrature.h"
#include "sep_coronal_cme/constants.h"

namespace SEP {
namespace CoronalCME {
namespace {
using Core::Status;
using Core::StatusCode;
bool Finite(double x) {
  return std::isfinite(x);
}
Status Invalid(const std::string& text) {
  return Status::Failure(StatusCode::InvalidConfiguration, text);
}
const KinematicValue& Component(const EllipsoidKinematics& k, std::size_t i) {
  if (i == 0)
    return k.centerDistanceM;
  if (i == 1)
    return k.radialSemiAxisM;
  if (i == 2)
    return k.firstLateralSemiAxisM;
  return k.secondLateralSemiAxisM;
}
std::string Value(const std::map<std::string, std::string>& input,
                  const std::string& key) {
  const auto found = input.find(key);
  if (found == input.end())
    throw std::runtime_error("missing event field " + key);
  return found->second;
}
double Number(const std::map<std::string, std::string>& input, const std::string& key) {
  std::istringstream in(Value(input, key));
  in.imbue(std::locale::classic());
  double value;
  if (!(in >> value) || !Finite(value))
    throw std::runtime_error("event field must be finite SI: " + key);
  in >> std::ws;
  if (!in.eof())
    throw std::runtime_error("trailing event numeric field: " + key);
  return value;
}
// A schema-5 snapshot asset carries analytic value/rate pairs. Its explicitly
// selected monotone-C1 contract is checked over EACH polynomial, not inferred
// from ordered knot values alone. No overshoot is repaired by clipping slopes.
Status Monotone(const ScalarGeometryHistory& history) {
  for (std::size_t i = 1; i < history.knots.size(); ++i) {
    const auto& l = history.knots[i - 1];
    const auto& r = history.knots[i];
    const double h = r.timeS - l.timeS, delta = r.value - l.value;
    if (h <= 0)
      return Invalid("history times must be strictly increasing");
    if (delta == 0) {
      if (l.ratePerS != 0 || r.ratePerS != 0)
        return Invalid("flat monotone history has nonzero slope");
      continue;
    }
    const double sign = delta > 0 ? 1. : -1.;
    const double a =
        sign * (6 * l.value - 6 * r.value + 3 * h * (l.ratePerS + r.ratePerS)) / h;
    const double b =
        sign * (-6 * l.value + 6 * r.value - h * (4 * l.ratePerS + 2 * r.ratePerS)) / h;
    const double c = sign * l.ratePerS;
    double minimum = std::min(c, a + b + c);
    if (a > 0) {
      const double x = -b / (2 * a);
      if (x > 0 && x < 1)
        minimum = std::min(minimum, (a * x + b) * x + c);
    }
    if (!Finite(minimum) || minimum < 0)
      return Invalid("tabulated monotone-C1 history overshoots between knots");
  }
  return Status::Success();
}
Core::Result<GeometryHistory> ReadAsset(const GeometryAsset& asset,
                                        const std::string& frame) {
  using Return = Core::Result<GeometryHistory>;
  if (asset.expectedSha256.size() != 64 ||
      ComputeContentChecksum(asset.bytes) != asset.expectedSha256)
    return Return::Failure(StatusCode::DataIntegrityFailure,
                           "ellipsoid asset content checksum mismatch");
  // This versioned geometry-only record is shared by snapshot histories and
  // independent pistons. It has no ambient, regional plasma, source or native
  // configuration side effects. Unknown/trailing fields reject atomically.
  std::istringstream in(asset.bytes);
  in.imbue(std::locale::classic());
  std::string magic, key, declaredFrame, kind;
  unsigned version = 0;
  double reference = 0, lon = 0, lat = 0, tilt = 0;
  if (!(in >> magic >> version >> key >> std::quoted(declaredFrame)) ||
      magic != "SCCM_FIXED_ELLIPSOID" || version != 1 || key != "frame" ||
      declaredFrame != frame)
    return Return::Failure(StatusCode::DataIntegrityFailure,
                           "ellipsoid asset schema/frame mismatch");
  if (!(in >> key >> reference) || key != "reference_time_s" ||
      !(in >> key >> lon >> lat >> tilt) || key != "direction_rad" ||
      !(in >> key >> std::quoted(kind)) || key != "evolution")
    return Return::Failure(StatusCode::InvalidConfiguration,
                           "incomplete ellipsoid asset header");
  auto basis = BuildRadialPrincipalBasis(lat, lon, tilt);
  if (!basis.ok())
    return Return::Failure(basis.status.code, basis.status.message);
  GeometryHistory result;
  result.basis = basis.value;
  if (kind == "independent-component-laws") {
    const char* names[] = {"center_distance", "radial_semiaxis", "lateral_semiaxis_1",
                           "lateral_semiaxis_2"};
    for (std::size_t i = 0; i < 4; ++i) {
      auto& c = result.components[i];
      std::string name;
      if (!(in >> key >> std::quoted(name) >> c.referenceValueM >> c.initialRateMPerS >>
            c.finalRateMPerS >> c.transitionStartS >> c.transitionDurationS) ||
          key != "component" || name != names[i])
        return Return::Failure(
            StatusCode::InvalidConfiguration,
            "ellipsoid asset requires four ordered complete component laws");
      c.referenceTimeS = reference;
    }
  } else if (kind == "tabulated-snapshots") {
    std::size_t count = 0;
    if (!(in >> key >> count) || key != "samples" || count < 2 || count > 100000)
      return Return::Failure(
          StatusCode::InvalidConfiguration,
          "ellipsoid asset requires a bounded complete snapshot table");
    for (auto& c : result.components) {
      c.kind = HistoryKind::CubicHermiteC1;
      c.referenceTimeS = reference;
    }
    for (std::size_t i = 0; i < count; ++i) {
      double time;
      if (!(in >> time))
        return Return::Failure(StatusCode::InvalidConfiguration,
                               "truncated ellipsoid history time");
      for (auto& c : result.components) {
        HermiteKnot knot;
        knot.timeS = time;
        if (!(in >> knot.value >> knot.ratePerS))
          return Return::Failure(StatusCode::InvalidConfiguration,
                                 "truncated ellipsoid history value/rate");
        c.knots.push_back(knot);
      }
    }
    for (const auto& c : result.components) {
      auto valid = Monotone(c);
      if (!valid.ok())
        return Return::Failure(valid.code, valid.message);
    }
  } else
    return Return::Failure(StatusCode::UnsupportedCapability,
                           "unsupported ellipsoid asset evolution");
  in >> std::ws;
  if (!in.eof())
    return Return::Failure(
        StatusCode::InvalidConfiguration,
        "trailing ellipsoid asset fields (dynamic attitude is unsupported)");
  return Return::Success(std::move(result));
}
Status CertifyAttachment(const GeometryHistory& history, const EventSupport& support,
                         const EventNumerics& controls, double begin, double end,
                         unsigned depth) {
  const double middle = begin + (end - begin) / 2;
  const auto state = history.At(middle);
  if (!state.ok())
    return state.status;
  const auto extent = EvaluateEllipsoidRadialExtent(state.value, support.solarRadiusM);
  if (!extent.ok())
    return extent.status;
  if (!extent.value.intersectsSolarSurface)
    return Status::Failure(StatusCode::OutOfDomain,
                           "attached history detaches from the physical solar sphere");
  double centerBound = 0, axesBound = 0;
  for (std::size_t i = 0; i < 4; ++i) {
    const auto range = history.components[i].Range(begin, end);
    if (!range.ok())
      return range.status;
    const double value = Component(state.value, i).value;
    const double bound = std::max(std::abs(range.value.first - value),
                                  std::abs(range.value.second - value));
    if (i == 0)
      centerBound = bound;
    else
      axesBound = std::max(axesBound, bound);
  }
  // For any unit parameter vector, |delta X| <= |delta center| +
  // max_i|delta a_i|. Min/max radial support are therefore Lipschitz in this
  // SAME bound. This interval proof cannot miss an interior detachment even
  // when both endpoints/knots happen to be attached.
  const double bound = centerBound + axesBound;
  if (extent.value.minimumRadiusM + bound <= support.solarRadiusM &&
      extent.value.maximumRadiusM - bound >= support.solarRadiusM)
    return Status::Success();
  if (depth >= controls.maximumCertificateDepth ||
      end - begin <= controls.minimumCertificateIntervalS || middle == begin ||
      middle == end)
    return Status::Failure(StatusCode::NumericalFailure,
                           "solar attachment interval certificate exhausted its "
                           "explicit work/resolution bound");
  const auto left =
      CertifyAttachment(history, support, controls, begin, middle, depth + 1);
  return left.ok()
             ? CertifyAttachment(history, support, controls, middle, end, depth + 1)
             : left;
}
Status CertifyPiston(const GeometryHistory& front, const GeometryHistory& piston,
                     double separation, const EventNumerics& controls, double begin,
                     double end, unsigned depth) {
  std::array<std::pair<double, double>, 4> f, p;
  for (std::size_t i = 0; i < 4; ++i) {
    const auto fr = front.components[i].Range(begin, end),
               pr = piston.components[i].Range(begin, end);
    if (!fr.ok() || !pr.ok())
      return !fr.ok() ? fr.status : pr.status;
    f[i] = fr.value;
    p[i] = pr.value;
  }
  const Vec3 fb[] = {front.basis.radial, front.basis.firstLateral,
                     front.basis.secondLateral};
  const Vec3 pb[] = {piston.basis.radial, piston.basis.firstLateral,
                     piston.basis.secondLateral};
  // Transform the ENTIRE piston into the front's normalized quadratic space.
  // ||b+A u|| <= ||b||+sqrt(||A||_1 ||A||_inf) for all ||u||=1.
  // The interval coefficients bound the complete time span as well. This is
  // a sufficient containment/standoff certificate, never a finite patch scan
  // accepted as proof for all rays. Unresolved valid shapes reject explicitly.
  double centerSquared = 0, normOne = 0, normInfinity = 0;
  double matrix[3][3]{};
  for (std::size_t i = 0; i < 3; ++i) {
    const double dot = Dot(fb[i], piston.basis.radial);
    const double low =
        std::min(dot * p[0].first, dot * p[0].second) - (i == 0 ? f[0].second : 0);
    const double high =
        std::max(dot * p[0].first, dot * p[0].second) - (i == 0 ? f[0].first : 0);
    const double c = std::max(std::abs(low), std::abs(high)) / f[i + 1].first;
    centerSquared += c * c;
    for (std::size_t j = 0; j < 3; ++j)
      matrix[i][j] = std::abs(Dot(fb[i], pb[j])) * p[j + 1].second / f[i + 1].first;
  }
  for (std::size_t i = 0; i < 3; ++i) {
    double row = 0, column = 0;
    for (std::size_t j = 0; j < 3; ++j) {
      row += matrix[i][j];
      column += matrix[j][i];
    }
    normOne = std::max(normOne, column);
    normInfinity = std::max(normInfinity, row);
  }
  const double normalizedBound =
      std::sqrt(centerSquared) + std::sqrt(normOne * normInfinity);
  const double minimumAxis = std::min({f[1].first, f[2].first, f[3].first});
  if (Finite(normalizedBound) && normalizedBound < 1 &&
      minimumAxis * (1 - normalizedBound) >= separation)
    return Status::Success();
  const double middle = begin + (end - begin) / 2;
  if (depth >= controls.maximumCertificateDepth ||
      end - begin <= controls.minimumCertificateIntervalS || middle == begin ||
      middle == end)
    return Status::Failure(StatusCode::NumericalFailure,
                           "independent piston containment/standoff cannot be "
                           "certified over the complete history");
  const auto left =
      CertifyPiston(front, piston, separation, controls, begin, middle, depth + 1);
  return left.ok() ? CertifyPiston(front, piston, separation, controls, middle, end,
                                   depth + 1)
                   : left;
}
std::string SealedIdentity(const BackgroundEventDefinition& definition) {
  std::ostringstream text;
  text.imbue(std::locale::classic());
  text << "SCCM-background-event-1;" << std::hexfloat;
  for (const auto& value : {definition.inputIdentity, definition.ambientIdentity,
                            definition.support.frame, definition.support.epochUtc})
    text << value.size() << ':' << value;
  const auto& s = definition.support;
  text << s.beginS << ';' << s.endS << ';' << s.solarRadiusM << ';'
       << s.qualifiedSourceInnerRadiusM << ';' << s.plasmaMinimumRadiusM << ';'
       << s.plasmaMaximumRadiusM << ';' << static_cast<int>(s.attachment) << ';'
       << static_cast<int>(s.plasmaCoverage) << ';' << definition.gammaAdiabatic << ';'
       << definition.includeElectronMass << ';' << definition.pistonEnabled << ';'
       << definition.minimumPistonSeparationM << ';';
  for (const auto* geometry :
       {&definition.front,
        definition.pistonEnabled ? &definition.piston : &definition.front}) {
    for (const auto& axis : {geometry->basis.radial, geometry->basis.firstLateral,
                             geometry->basis.secondLateral})
      text << axis.x << ';' << axis.y << ';' << axis.z << ';';
    for (const auto& history : geometry->components) {
      text << static_cast<int>(history.kind) << ';' << history.referenceTimeS << ';'
           << history.referenceValueM << ';' << history.initialRateMPerS << ';'
           << history.finalRateMPerS << ';' << history.transitionStartS << ';'
           << history.transitionDurationS << ';' << history.knots.size() << ';';
      for (const auto& knot : history.knots)
        text << knot.timeS << ';' << knot.value << ';' << knot.ratePerS << ';';
    }
  }
  return ComputeContentChecksum(text.str());
}
bool SameEpoch(const AmbientEpochIdentity& first, const AmbientEpochIdentity& second) {
  return first.physicsIdentity == second.physicsIdentity &&
         first.frame == second.frame && first.epochUtc == second.epochUtc &&
         first.timeS == second.timeS && first.validFromS == second.validFromS &&
         first.validUntilS == second.validUntilS &&
         first.gammaAdiabatic == second.gammaAdiabatic &&
         first.generation == second.generation &&
         first.includeElectronMass == second.includeElectronMass;
}
Status ValidateAmbient(const EventAmbientSample& sample, const AmbientEpochIdentity& id,
                       QuerySide requested) {
  const bool closed = sample.topology == FieldLineTopology::ClosedBelowOuterBoundary;
  if (sample.region != AmbientRegion::PfssOpen &&
      sample.region != AmbientRegion::PfssClosed &&
      sample.region != AmbientRegion::FiniteShellScs &&
      sample.region != AmbientRegion::ParkerExterior)
    return Status::Failure(StatusCode::DataIntegrityFailure,
                           "unknown ambient region code");
  if (sample.topology == FieldLineTopology::Invalid ||
      sample.topology == FieldLineTopology::Separatrix ||
      (closed ? (sample.region != AmbientRegion::PfssClosed || sample.sector != 0)
              : (sample.region == AmbientRegion::PfssClosed ||
                 (sample.sector != 1 && sample.sector != -1))))
    return Status::Failure(
        StatusCode::DataIntegrityFailure,
        "ambient topology/region/signed sector is invalid or ambiguous");
  if ((!sample.oneSided && sample.side != QuerySide::Ordinary) ||
      (requested != QuerySide::Ordinary &&
       (!sample.oneSided || sample.side != requested)))
    return Status::Failure(
        StatusCode::DataIntegrityFailure,
        "ambient physical trace metadata is missing or inconsistent");
  if (sample.oneSided && (requested == QuerySide::Ordinary || sample.side != requested))
    return Status::Failure(
        StatusCode::DataIntegrityFailure,
        "ambient discontinuity requires the explicitly requested one-sided state");
  if (!sample.derivativesValid)
    return Status::Failure(StatusCode::UnsupportedCapability,
                           "ambient epoch lacks required one-sided spatial tensors");
  for (const auto* tensor : {&sample.gradientBPerM, &sample.gradientUPerS})
    for (double value : *tensor)
      if (!Finite(value))
        return Status::Failure(StatusCode::InvalidState,
                               "ambient spatial tensor is nonfinite");
  const auto eos = EvaluatePlasmaFromElectronDensity(
      sample.electronDensityM3, sample.electronTemperatureK, sample.ions,
      id.includeElectronMass, id.gammaAdiabatic, Norm(sample.primitive.magneticFieldT));
  if (!eos.ok())
    return eos.status;
  const auto consistent = [](double a, double b) {
    return Finite(a) && Finite(b) &&
           std::abs(a - b) <= 1e-10 * std::max(std::abs(a), std::abs(b));
  };
  if (!consistent(eos.value.massDensityKgM3, sample.primitive.massDensityKgM3) ||
      !consistent(eos.value.pressurePa, sample.primitive.pressurePa))
    return Status::Failure(StatusCode::DataIntegrityFailure,
                           "ambient primitive density/pressure does not match its "
                           "explicit composition/EOS");
  return Status::Success();
}
}  // namespace
// namespace
Core::Result<KinematicValue> ScalarGeometryHistory::At(double t) const {
  if (kind == HistoryKind::CubicHermiteC1)
    return EvaluateCubicHermiteHistory(knots, t);
  if (kind != HistoryKind::SmoothRateTransition)
    return Core::Result<KinematicValue>::Failure(StatusCode::UnsupportedCapability,
                                                 "unsupported geometry history kind");
  if (!Finite(transitionDurationS))
    return Core::Result<KinematicValue>::Failure(
        StatusCode::InvalidConfiguration, "geometry transition duration is not finite");
  return EvaluateSmoothKinematics(t, referenceTimeS, referenceValueM, transitionStartS,
                                  transitionDurationS, initialRateMPerS,
                                  finalRateMPerS);
}
Core::Result<std::pair<double, double>> ScalarGeometryHistory::Range(double begin,
                                                                     double end) const {
  using Return = Core::Result<std::pair<double, double>>;
  if (!Finite(begin) || !Finite(end) || end < begin)
    return Return::Failure(StatusCode::InvalidConfiguration,
                           "geometry range requires a finite ordered interval");
  auto first = At(begin), last = At(end);
  if (!first.ok() || !last.ok())
    return Return::Failure((!first.ok() ? first : last).status.code,
                           (!first.ok() ? first : last).status.message);
  double low = std::min(first.value.value, last.value.value),
         high = std::max(first.value.value, last.value.value);
  const auto include = [&](double time) {
    if (time < begin || time > end)
      return Status::Success();
    const auto point = At(time);
    if (!point.ok())
      return point.status;
    low = std::min(low, point.value.value);
    high = std::max(high, point.value.value);
    return Status::Success();
  };
  if (kind == HistoryKind::SmoothRateTransition) {
    for (double time : {transitionStartS, transitionStartS + transitionDurationS}) {
      const auto status = include(time);
      if (!status.ok())
        return Return::Failure(status.code, status.message);
    }
    if ((initialRateMPerS < 0 && finalRateMPerS > 0) ||
        (initialRateMPerS > 0 && finalRateMPerS < 0)) {
      double l = 0, r = 1;
      for (unsigned i = 0; i < 80; ++i) {
        const double m = (l + r) / 2;
        const double s = m * m * m * (10 + m * (-15 + 6 * m));
        const double v = initialRateMPerS + (finalRateMPerS - initialRateMPerS) * s;
        if ((v > 0) == (initialRateMPerS > 0))
          l = m;
        else
          r = m;
      }
      const auto status = include(transitionStartS + transitionDurationS * (l + r) / 2);
      if (!status.ok())
        return Return::Failure(status.code, status.message);
    }
  } else {
    for (std::size_t i = 1; i < knots.size(); ++i) {
      const auto& l = knots[i - 1];
      const auto& r = knots[i];
      if (r.timeS < begin || l.timeS > end)
        continue;
      const double h = r.timeS - l.timeS,
                   a = 6 * l.value - 6 * r.value + 3 * h * (l.ratePerS + r.ratePerS);
      const double b = -6 * l.value + 6 * r.value -
                       h * (4 * l.ratePerS + 2 * r.ratePerS),
                   c = h * l.ratePerS;
      // At(begin/end) above already validated the complete immutable knot
      // table. Include its exact stored endpoints directly. Re-validating
      // all N knots at each of N segments made supported100000-row assets
      // quadratic without adding any independent physical validation.
      if (l.timeS >= begin && l.timeS <= end) {
        low = std::min(low, l.value);
        high = std::max(high, l.value);
      }
      if (r.timeS >= begin && r.timeS <= end) {
        low = std::min(low, r.value);
        high = std::max(high, r.value);
      }
      std::vector<double> roots;
      if (a == 0) {
        if (b != 0)
          roots.push_back(-c / b);
      } else {
        const double disc = b * b - 4 * a * c;
        if (disc >= 0) {
          const double q = -.5 * (b + std::copysign(std::sqrt(disc), b));
          if (q == 0)
            roots.push_back(-b / (2 * a));
          else {
            roots.push_back(q / a);
            roots.push_back(c / q);
          }
        }
      }
      for (double x : roots)
        if (x > 0 && x < 1) {
          const double rootTime = l.timeS + h * x;
          if (rootTime < begin || rootTime > end)
            continue;
          // Reuse the canonical polynomial/derivative evaluator on the
          // already validated LOCAL pair. This preserves the same extremum
          // formula and arithmetic while keeping each root evaluation O(1).
          const auto extremum = EvaluateCubicHermiteHistory({l, r}, rootTime);
          if (!extremum.ok())
            return Return::Failure(extremum.status.code, extremum.status.message);
          low = std::min(low, extremum.value.value);
          high = std::max(high, extremum.value.value);
        }
    }
  }
  if (!Finite(low) || !Finite(high))
    return Return::Failure(StatusCode::NumericalFailure,
                           "geometry history range is nonfinite");
  return Return::Success({low, high});
}
Core::Result<EllipsoidKinematics> GeometryHistory::At(double time) const {
  KinematicValue values[4];
  for (std::size_t i = 0; i < 4; ++i) {
    const auto value = components[i].At(time);
    if (!value.ok())
      return Core::Result<EllipsoidKinematics>::Failure(value.status.code,
                                                        value.status.message);
    values[i] = value.value;
  }
  return Core::Result<EllipsoidKinematics>::Success(
      {values[0], values[1], values[2], values[3]});
}

Core::Result<BackgroundEventDefinition> MaterializeBackgroundEvent(
    const ModelConfiguration& config, const EventSupport& support,
    const EventNumerics& numerics, const std::string& ambientIdentity,
    const GeometryAssetReader& reader) {
  using Return = Core::Result<BackgroundEventDefinition>;
  try {
    const auto& a = config.assignments;
    const auto configured = ValidateConfiguration(config);
    if (!configured.ok())
      return Return::Failure(configured.code, configured.message);
    if (config.physicsFingerprint.empty() ||
        ComputePhysicsFingerprint(a) != config.physicsFingerprint ||
        Number(a, "plasma_eos.adiabatic_index") != config.gammaAdiabatic)
      return Return::Failure(
          StatusCode::DataIntegrityFailure,
          "materialized event configuration fingerprint/EOS is stale");
    if (a.count("source.enabled") && Value(a, "source.enabled") != "false")
      return Return::Failure(
          StatusCode::UnsupportedCapability,
          "background event assembly requires explicitly disabled particle sources");
    if (Value(a, "shock_geometry.model") != "triaxial-ellipsoid" ||
        Value(a, "shock_geometry.axis_convention") != "radial-principal-axis" ||
        Value(a, "shock_geometry.orientation_evolution") != "fixed" ||
        Value(a, "shock_geometry.snapshot_extrapolation") != "reject")
      return Return::Failure(StatusCode::UnsupportedCapability,
                             "background event requires maintained fixed-orientation "
                             "triaxial geometry and no extrapolation");
    if (Value(a, "shock_geometry.direction_frame") != support.frame ||
        Value(a, "domain.coordinate_frame") != support.frame ||
        Value(a, "domain.epoch_utc") != support.epochUtc)
      return Return::Failure(StatusCode::DataIntegrityFailure,
                             "event and support frame/epoch disagree");
    BackgroundEventDefinition out;
    out.support = support;
    out.numerics = numerics;
    out.ambientIdentity = ambientIdentity;
    out.gammaAdiabatic = config.gammaAdiabatic;
    out.includeElectronMass = config.includeElectronMassInDensity;
    if (config.startTimeS != support.beginS || config.endTimeS != support.endS ||
        config.solarRadiusM != support.solarRadiusM ||
        Number(a, "domain.qualified_source_inner_radius_m") !=
            support.qualifiedSourceInnerRadiusM)
      return Return::Failure(StatusCode::DataIntegrityFailure,
                             "event run/radial bounds disagree with explicit support");
    const auto basis =
        BuildRadialPrincipalBasis(Number(a, "shock_geometry.direction_latitude_rad"),
                                  Number(a, "shock_geometry.direction_longitude_rad"),
                                  Number(a, "shock_geometry.lateral_axis_tilt_rad"));
    if (!basis.ok())
      return Return::Failure(basis.status.code, basis.status.message);
    out.front.basis = basis.value;
    const auto evolution = Value(a, "shock_geometry.evolution");
    std::string assetIdentity;
    if (evolution == "independent-component-laws") {
      const char* names[] = {"center_distance", "radial_semiaxis", "lateral_semiaxis_1",
                             "lateral_semiaxis_2"};
      for (std::size_t i = 0; i < 4; ++i) {
        auto& c = out.front.components[i];
        const auto prefix = "shock_geometry." + std::string(names[i]) + ".";
        if (Value(a, prefix + "kinematics") != "smooth-rate-transition" ||
            Value(a, prefix + "after_transition") != "constant-final-rate")
          return Return::Failure(
              StatusCode::UnsupportedCapability,
              "geometry component requires complete maintained smooth-rate law");
        c.referenceTimeS = Number(a, "shock_geometry.reference_time_s");
        c.referenceValueM = Number(a, prefix + "initial_value_m");
        c.initialRateMPerS = Number(a, prefix + "initial_rate_m_per_s");
        c.finalRateMPerS = Number(a, prefix + "final_rate_m_per_s");
        c.transitionStartS = Number(a, prefix + "transition_start_time_s");
        c.transitionDurationS = Number(a, prefix + "transition_duration_s");
      }
      if (Value(a, "shock_geometry.snapshots_file") != "none" ||
          Value(a, "shock_geometry.snapshot_interpolation") != "none")
        return Return::Failure(
            StatusCode::InvalidConfiguration,
            "independent laws cannot also select a snapshot authority");
      const auto parameterization = Value(a, "shock_geometry.initial_parameterization");
      const double apex = Number(a, "shock_geometry.initial_apex_radius_m");
      if (parameterization == "apex-and-principal-axes") {
        if (out.front.components[0].referenceValueM != 0 || apex <= 0)
          return Return::Failure(
              StatusCode::InvalidConfiguration,
              "apex form conflicts with an independently supplied initial center");
        out.front.components[0].referenceValueM =
            apex - out.front.components[1].referenceValueM;
      } else if (parameterization != "center-and-principal-axes" || apex != 0)
        return Return::Failure(
            StatusCode::InvalidConfiguration,
            "exactly one initial radial parameterization is required");
    } else if (evolution == "tabulated-snapshots") {
      const char* names[] = {"center_distance", "radial_semiaxis", "lateral_semiaxis_1",
                             "lateral_semiaxis_2"};
      for (const auto* name : names) {
        const auto prefix = "shock_geometry." + std::string(name) + ".";
        if (Value(a, prefix + "kinematics") != "inactive" ||
            Value(a, prefix + "after_transition") != "none")
          return Return::Failure(
              StatusCode::InvalidConfiguration,
              "snapshot history has an active competing component law");
        for (const auto* number :
             {"initial_value_m", "initial_rate_m_per_s", "final_rate_m_per_s",
              "transition_start_time_s", "transition_duration_s"})
          if (Number(a, prefix + number) != 0)
            return Return::Failure(
                StatusCode::InvalidConfiguration,
                "snapshot component-law fields must be explicitly inactive/zero");
      }
      if (!reader ||
          Value(a, "shock_geometry.snapshot_interpolation") != "monotone-cubic-c1")
        return Return::Failure(StatusCode::UnsupportedCapability,
                               "snapshot history needs an explicit verified asset "
                               "reader and monotone-C1 contract");
      const auto asset = reader(Value(a, "shock_geometry.snapshots_file"));
      if (!asset.ok())
        return Return::Failure(asset.status.code, asset.status.message);
      const auto parsed = ReadAsset(asset.value, support.frame);
      if (!parsed.ok())
        return Return::Failure(parsed.status.code, parsed.status.message);
      for (std::size_t i = 0; i < 3; ++i) {
        const Vec3 p[] = {parsed.value.basis.radial, parsed.value.basis.firstLateral,
                          parsed.value.basis.secondLateral};
        const Vec3 f[] = {out.front.basis.radial, out.front.basis.firstLateral,
                          out.front.basis.secondLateral};
        if (Norm(p[i] - f[i]) > 1e-12)
          return Return::Failure(StatusCode::DataIntegrityFailure,
                                 "snapshot asset changes fixed parent orientation");
      }
      for (const auto& c : parsed.value.components)
        if (c.kind != HistoryKind::CubicHermiteC1)
          return Return::Failure(StatusCode::DataIntegrityFailure,
                                 "snapshot asset contains a competing evolution law");
      if (parsed.value.components[0].referenceTimeS !=
          Number(a, "shock_geometry.reference_time_s"))
        return Return::Failure(StatusCode::DataIntegrityFailure,
                               "snapshot reference time differs from parent geometry");
      const auto initial =
          parsed.value.At(Number(a, "shock_geometry.reference_time_s"));
      if (!initial.ok())
        return Return::Failure(initial.status.code, initial.status.message);
      const auto parameterization = Value(a, "shock_geometry.initial_parameterization");
      const double apex = Number(a, "shock_geometry.initial_apex_radius_m");
      if ((parameterization == "center-and-principal-axes" && apex != 0) ||
          (parameterization == "apex-and-principal-axes" &&
           apex != initial.value.centerDistanceM.value +
                       initial.value.radialSemiAxisM.value) ||
          (parameterization != "center-and-principal-axes" &&
           parameterization != "apex-and-principal-axes"))
        return Return::Failure(
            StatusCode::DataIntegrityFailure,
            "snapshot radial parameterization conflicts with complete history");
      out.front = parsed.value;
      assetIdentity += asset.value.expectedSha256;
    } else
      return Return::Failure(StatusCode::UnsupportedCapability,
                             "unsupported/research-only geometry evolution");
    const auto enabled = Value(a, "cme_piston.enabled");
    if (enabled != "true" && enabled != "false")
      return Return::Failure(StatusCode::InvalidConfiguration,
                             "piston enabled must be explicit true/false");
    out.pistonEnabled = enabled == "true";
    if (out.pistonEnabled) {
      // The enabled flag and an existing asset do not bind a physical
      // authority on their own. The named record must explicitly select the
      // independent piston record; accepting "none" here would silently
      // promote otherwise inactive geometry fields into a piston model.
      if (!reader || Value(a, "cme_piston.geometry_source") != "independent-record" ||
          Value(a, "cme_piston.geometry_record") != "piston_geometry" ||
          Value(a, "piston_geometry.model") != "triaxial-ellipsoid-versioned-file")
        return Return::Failure(
            StatusCode::UnsupportedCapability,
            "piston requires its own complete independently verified geometry asset");
      const auto asset = reader(Value(a, "piston_geometry.record_file"));
      if (!asset.ok())
        return Return::Failure(asset.status.code, asset.status.message);
      const auto parsed = ReadAsset(asset.value, support.frame);
      if (!parsed.ok())
        return Return::Failure(parsed.status.code, parsed.status.message);
      out.piston = parsed.value;
      out.minimumPistonSeparationM =
          Number(a, "piston_geometry.minimum_front_separation_m");
      assetIdentity += asset.value.expectedSha256;
    }
    if (!out.pistonEnabled &&
        (Value(a, "cme_piston.geometry_source") != "none" ||
         Value(a, "cme_piston.geometry_record") != "none" ||
         Value(a, "piston_geometry.model") != "none" ||
         Value(a, "piston_geometry.record_file") != "none" ||
         Number(a, "piston_geometry.minimum_front_separation_m") != 0))
      return Return::Failure(StatusCode::InvalidConfiguration,
                             "disabled piston has a competing active geometry record");
    // Length-prefix variable strings. Event identity is physical: surface
    // quadrature and certificate work limits have separate numerical meaning
    // and do not redefine the eruption when a convergence grid is refined.
    std::ostringstream identity;
    identity.imbue(std::locale::classic());
    for (const auto& value : {config.physicsFingerprint, assetIdentity, ambientIdentity,
                              support.frame, support.epochUtc})
      identity << value.size() << ':' << value;
    identity << std::hexfloat << support.plasmaMinimumRadiusM << ';'
             << support.plasmaMaximumRadiusM << ';'
             << static_cast<int>(support.attachment) << ';'
             << static_cast<int>(support.plasmaCoverage);
    out.inputIdentity = ComputeContentChecksum(identity.str());
    return Return::Success(std::move(out));
  } catch (const std::exception& e) {
    return Return::Failure(StatusCode::InvalidConfiguration, e.what());
  }
}

Core::Result<std::shared_ptr<BackgroundEvent>> BackgroundEvent::Create(
    BackgroundEventDefinition definition) {
  using Return = Core::Result<std::shared_ptr<BackgroundEvent>>;
  const auto& s = definition.support;
  const auto& n = definition.numerics;
  if (!(Finite(s.beginS) && Finite(s.endS) && s.endS > s.beginS &&
        Finite(s.solarRadiusM) && s.solarRadiusM > 0 &&
        Finite(s.qualifiedSourceInnerRadiusM) &&
        s.qualifiedSourceInnerRadiusM >= s.solarRadiusM &&
        Finite(s.plasmaMinimumRadiusM) && s.plasmaMinimumRadiusM >= s.solarRadiusM &&
        Finite(s.plasmaMaximumRadiusM) &&
        s.plasmaMaximumRadiusM > s.plasmaMinimumRadiusM && !s.frame.empty() &&
        !s.epochUtc.empty() && Finite(definition.gammaAdiabatic) &&
        definition.gammaAdiabatic > 1 && !definition.inputIdentity.empty() &&
        !definition.ambientIdentity.empty() && n.polarCells >= 2 &&
        n.azimuthCells >= 3 && n.maximumCertificateDepth > 0 &&
        n.maximumCertificateDepth <= 30 && Finite(n.minimumCertificateIntervalS) &&
        n.minimumCertificateIntervalS > 0))
    return Return::Failure(StatusCode::InvalidConfiguration,
                           "background event requires complete finite "
                           "physics/support/identity and explicit numerical controls");
  if (s.attachment != AttachmentPolicy::RequireSolarIntersection &&
      s.attachment != AttachmentPolicy::AllowDetachedFront)
    return Return::Failure(StatusCode::UnsupportedCapability,
                           "unknown attachment policy");
  if (s.plasmaCoverage != PlasmaCoveragePolicy::RequireEveryExposedPatch &&
      s.plasmaCoverage != PlasmaCoveragePolicy::GeometryOnlyOutsideCoverage)
    return Return::Failure(StatusCode::UnsupportedCapability,
                           "unknown plasma coverage policy");
  for (const auto* geometry :
       {&definition.front,
        definition.pistonEnabled ? &definition.piston : &definition.front}) {
    for (const auto& component : geometry->components) {
      const auto range = component.Range(s.beginS, s.endS);
      if (!range.ok())
        return Return::Failure(range.status.code, range.status.message);
      if (range.value.first <= 0)
        return Return::Failure(
            StatusCode::OutOfDomain,
            "nonpositive ellipsoid center/axis between history knots");
    }
    const auto state = geometry->At(s.beginS);
    if (!state.ok())
      return Return::Failure(state.status.code, state.status.message);
    const auto shape = FixedOrientationEllipsoid::FromCenter(
        geometry->basis, state.value, s.solarRadiusM);
    if (!shape.ok())
      return Return::Failure(shape.status.code, shape.status.message);
  }
  const auto initial = definition.front.At(s.beginS);
  if (initial.value.centerDistanceM.value + initial.value.radialSemiAxisM.value <
      s.qualifiedSourceInnerRadiusM)
    return Return::Failure(
        StatusCode::OutOfDomain,
        "initial front apex is below the declared source qualification radius");
  if (s.attachment == AttachmentPolicy::RequireSolarIntersection) {
    const auto valid = CertifyAttachment(definition.front, s, n, s.beginS, s.endS, 0);
    if (!valid.ok())
      return Return::Failure(valid.code, valid.message);
  }
  if (definition.pistonEnabled) {
    if (!Finite(definition.minimumPistonSeparationM) ||
        definition.minimumPistonSeparationM < 0)
      return Return::Failure(StatusCode::InvalidConfiguration,
                             "piston standoff must be finite and nonnegative");
    const auto valid =
        CertifyPiston(definition.front, definition.piston,
                      definition.minimumPistonSeparationM, n, s.beginS, s.endS, 0);
    if (!valid.ok())
      return Return::Failure(valid.code, valid.message);
  }
  definition.physicsIdentity = SealedIdentity(definition);
  std::shared_ptr<BackgroundEvent> result(new BackgroundEvent);
  result->definition_ = std::move(definition);
  return Return::Success(std::move(result));
}
Core::Result<std::shared_ptr<const BackgroundEventEpoch>> BackgroundEvent::Prepare(
    double time, std::shared_ptr<const AmbientEpoch> ambient, QuerySide side) {
  using Return = Core::Result<std::shared_ptr<const BackgroundEventEpoch>>;
  try {
    const auto& d = definition_;
    if (side != QuerySide::Ordinary && side != QuerySide::Minus &&
        side != QuerySide::Plus)
      return Return::Failure(StatusCode::InvalidConfiguration,
                             "unknown ambient query side");
    if (!Finite(time) || time < d.support.beginS || time > d.support.endS)
      return Return::Failure(StatusCode::OutOfDomain,
                             "event epoch is outside complete declared histories");
    if (!ambient)
      return Return::Failure(
          StatusCode::InvalidState,
          "background event requires an owning immutable ambient epoch");
    const auto id = ambient->Identity();
    if (id.generation == 0 || id.physicsIdentity != d.ambientIdentity ||
        id.frame != d.support.frame || id.epochUtc != d.support.epochUtc ||
        id.timeS != time || !Finite(id.validFromS) || !Finite(id.validUntilS) ||
        id.validFromS > time || id.validUntilS < time ||
        id.gammaAdiabatic != d.gammaAdiabatic ||
        id.includeElectronMass != d.includeElectronMass)
      return Return::Failure(
          StatusCode::DataIntegrityFailure,
          "ambient/event epoch, identity, frame, EOS or coverage mismatch");
    if (current_ && id.generation < current_->ambient->Identity().generation)
      return Return::Failure(StatusCode::DataIntegrityFailure,
                             "ambient epoch generation regressed");
    if (current_ && (time < current_->timeS ||
                     current_->generation == std::numeric_limits<std::uint64_t>::max()))
      return Return::Failure(StatusCode::InvalidState,
                             "background event cannot reverse its publication clock or "
                             "overflow generation");
    auto state = d.front.At(time);
    if (!state.ok())
      return Return::Failure(state.status.code, state.status.message);
    const auto shape = FixedOrientationEllipsoid::FromCenter(d.front.basis, state.value,
                                                             d.support.solarRadiusM);
    if (!shape.ok())
      return Return::Failure(shape.status.code, shape.status.message);
    const auto extent =
        EvaluateEllipsoidRadialExtent(state.value, d.support.solarRadiusM);
    if (!extent.ok())
      return Return::Failure(extent.status.code, extent.status.message);
    // Strict coverage is a CONTINUOUS surface admission, independent of the
    // tessellation. In particular an unsupported apex/solar strip between
    // quadrature centers cannot disappear merely by coarsening the grid.
    const double exposedMinimum =
        std::max(d.support.solarRadiusM, extent.value.minimumRadiusM);
    if (d.support.plasmaCoverage == PlasmaCoveragePolicy::RequireEveryExposedPatch &&
        (exposedMinimum < d.support.plasmaMinimumRadiusM ||
         extent.value.maximumRadiusM > d.support.plasmaMaximumRadiusM))
      return Return::Failure(
          StatusCode::OutOfDomain,
          "continuous exposed front radial extent exceeds declared plasma support");
    const auto patches =
        TessellateSolarExposed(shape.value, d.support.solarRadiusM,
                               d.numerics.polarCells, d.numerics.azimuthCells);
    if (!patches.ok())
      return Return::Failure(patches.status.code, patches.status.message);
    auto candidate = std::make_shared<BackgroundEventEpoch>();
    candidate->generation = current_ ? current_->generation + 1 : 1;
    candidate->timeS = time;
    candidate->eventIdentity = d.physicsIdentity;
    candidate->front = shape.value;
    candidate->pistonEnabled = d.pistonEnabled;
    candidate->ambient = ambient;
    candidate->frontRadialExtent = extent.value;
    if (d.pistonEnabled) {
      const auto pk = d.piston.At(time);
      if (!pk.ok())
        return Return::Failure(pk.status.code, pk.status.message);
      const auto p = FixedOrientationEllipsoid::FromCenter(d.piston.basis, pk.value,
                                                           d.support.solarRadiusM);
      if (!p.ok())
        return Return::Failure(p.status.code, p.status.message);
      candidate->piston = p.value;
      const auto pe = EvaluateEllipsoidRadialExtent(pk.value, d.support.solarRadiusM);
      if (!pe.ok())
        return Return::Failure(pe.status.code, pe.status.message);
      // An independently configured piston may initially be wholly buried.
      // Its complete geometry/history remains present, but no below-Sun patch
      // is generated or assigned ambient plasma by an invented continuation.
      if (pe.value.maximumRadiusM > d.support.solarRadiusM) {
        const auto pp =
            TessellateSolarExposed(p.value, d.support.solarRadiusM,
                                   d.numerics.polarCells, d.numerics.azimuthCells);
        if (!pp.ok())
          return Return::Failure(pp.status.code, pp.status.message);
        candidate->pistonPatches = pp.value;
      }
    }
    std::vector<ShockPatchInput> shocks;
    for (const auto& geometry : patches.value) {
      BackgroundEventPatch patch;
      patch.geometry = geometry;
      const double radius = Norm(geometry.centerM);
      // ||dX/dtheta|| and ||dX/dphi|| are bounded by the largest
      // semiaxis. Together with the exact global radial extent this encloses
      // the ENTIRE solar-exposed parameter cell around its exposed representative. The
      // representative can shift
      // anywhere inside a cut cell, so FULL angular widths bound its extent.
      // Boundary-crossing uncertainty remains explicit in diagnostic mode;
      // a center value is never advertised as a full-cell coverage proof.
      const double maximumAxis = std::max({state.value.radialSemiAxisM.value,
                                           state.value.firstLateralSemiAxisM.value,
                                           state.value.secondLateralSemiAxisM.value});
      const double cellBound =
          maximumAxis * (Constants::kPi / d.numerics.polarCells +
                         2 * Constants::kPi / d.numerics.azimuthCells);
      patch.parameterCellMinimumRadiusM = std::max(exposedMinimum, radius - cellBound);
      patch.parameterCellMaximumRadiusM =
          std::min(extent.value.maximumRadiusM, radius + cellBound);
      if (patch.parameterCellMinimumRadiusM >= d.support.plasmaMinimumRadiusM &&
          patch.parameterCellMaximumRadiusM <= d.support.plasmaMaximumRadiusM)
        patch.parameterCellCoverage = ParameterCellCoverage::InsidePlasmaSupport;
      else if (patch.parameterCellMaximumRadiusM < d.support.plasmaMinimumRadiusM ||
               patch.parameterCellMinimumRadiusM > d.support.plasmaMaximumRadiusM)
        patch.parameterCellCoverage = ParameterCellCoverage::OutsidePlasmaSupport;
      patch.sourceRadialSupport = radius >= d.support.qualifiedSourceInnerRadiusM
                                      ? SourceRadialSupport::QualifiedRadius
                                      : SourceRadialSupport::BelowQualifiedSourceRadius;
      const auto normal = shape.value.Evaluate(geometry.centerM);
      if (!normal.ok())
        return Return::Failure(normal.status.code, normal.status.message);
      patch.normalSpeedMPerS = normal.value.normalSpeedMPerS;
      if (radius < d.support.plasmaMinimumRadiusM)
        patch.plasmaSupport = PatchSupport::BelowPlasmaCoverage;
      else if (radius > d.support.plasmaMaximumRadiusM)
        patch.plasmaSupport = PatchSupport::AbovePlasmaCoverage;
      if (patch.plasmaSupport != PatchSupport::PlasmaAvailable) {
        if (d.support.plasmaCoverage == PlasmaCoveragePolicy::RequireEveryExposedPatch)
          return Return::Failure(StatusCode::OutOfDomain,
                                 "exposed front patch has no declared plasma support");
        candidate->frontPatches.push_back(std::move(patch));
        continue;
      }
      const auto sampled = ambient->At(geometry.centerM, side);
      if (!sampled.ok())
        return Return::Failure(sampled.status.code, sampled.status.message);
      const auto valid = ValidateAmbient(sampled.value, id, side);
      if (!valid.ok())
        return Return::Failure(valid.code, valid.message);
      patch.ambient = sampled.value;
      ShockPatchInput input;
      input.stableId = geometry.physicalId;
      input.parentId = geometry.physicalId;
      input.areaM2 = geometry.areaM2;
      input.outwardNormal = geometry.outwardNormal;
      input.upstream = sampled.value.primitive;
      input.shockNormalSpeedMPerS = patch.normalSpeedMPerS;
      input.sourceEnabled = false;
      input.interface = {InterfacePolicy::DiagnosticKinematic,
                         id.generation,
                         true,
                         true,
                         true,
                         false,
                         "background-only:" + id.physicsIdentity};
      shocks.push_back(input);
      candidate->frontPatches.push_back(std::move(patch));
    }
    // The child provider is COPIED into the private transaction. Its own
    // successful Prepare cannot advance a public generation if a later epoch
    // allocation/check fails. Both owning handles publish together below.
    auto shockCandidate = shockProvider_;
    if (!shocks.empty()) {
      ShockPreparationOptions options;
      options.gammaAdiabatic = d.gammaAdiabatic;
      const auto prepared =
          shockCandidate.Prepare(time, id.generation, shocks, options);
      if (!prepared.ok())
        return Return::Failure(prepared.status.code, prepared.status.message);
      candidate->shocks = prepared.value;
    }
    if (!SameEpoch(id, ambient->Identity()))
      return Return::Failure(
          StatusCode::DataIntegrityFailure,
          "ambient owner identity changed during complete epoch preparation");
    std::shared_ptr<const BackgroundEventEpoch> publication = std::move(candidate);
    shockProvider_ = std::move(shockCandidate);
    current_ = publication;
    return Return::Success(std::move(publication));
  } catch (const std::exception& e) {
    return Return::Failure(
        StatusCode::NumericalFailure,
        std::string("background event preparation failed without publication: ") +
            e.what());
  }
}

}  // namespace CoronalCME
}  // namespace SEP
