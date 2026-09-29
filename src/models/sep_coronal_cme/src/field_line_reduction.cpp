#include "sep_coronal_cme/field_line_reduction.h"

#include "sep_field_line_bundle_io.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <utility>

namespace SEP { namespace CoronalCME { namespace {

Vec3 FromNeutral(FieldLine::Vector3 value) {
  return {value.x, value.y, value.z};
}
FieldLine::Vector3 ToNeutral(Vec3 value) {
  return {value.x, value.y, value.z};
}

Core::Result<Vec3> Direction(Vec3 position, double time, int sign,
                             const FieldLineStateEvaluator& evaluator) {
  auto state = evaluator(position, time);
  if (!state.ok()) return Core::Result<Vec3>::Failure(
      state.status.code, state.status.message);
  const Vec3 magnetic = FromNeutral(state.value.magneticFieldT);
  if (Norm(magnetic) == 0.0)
    return Core::Result<Vec3>::Failure(Core::StatusCode::InvalidState,
                                      "field-line trace reached B=0");
  return Core::Result<Vec3>::Success(static_cast<double>(sign) * Unit(magnetic));
}

Core::Result<Vec3> Rk4Step(Vec3 position, double time, double step,
                           int sign,
                           const FieldLineStateEvaluator& evaluator) {
  auto k1 = Direction(position, time, sign, evaluator);
  if (!k1.ok()) return k1;
  auto k2 = Direction(position + 0.5 * step * k1.value, time, sign, evaluator);
  if (!k2.ok()) return k2;
  auto k3 = Direction(position + 0.5 * step * k2.value, time, sign, evaluator);
  if (!k3.ok()) return k3;
  auto k4 = Direction(position + step * k3.value, time, sign, evaluator);
  if (!k4.ok()) return k4;
  return Core::Result<Vec3>::Success(position + step / 6.0 *
      (k1.value + 2.0 * k2.value + 2.0 * k3.value + k4.value));
}

Vec3 SegmentSphereEvent(Vec3 begin, Vec3 end, double radius) {
  const Vec3 direction = end - begin;
  const double a = Dot(direction, direction);
  const double b = 2.0 * Dot(begin, direction);
  const double c = Dot(begin, begin) - radius * radius;
  const double discriminant = std::max(0.0, b * b - 4.0 * a * c);
  const double first = (-b - std::sqrt(discriminant)) / (2.0 * a);
  const double second = (-b + std::sqrt(discriminant)) / (2.0 * a);
  double fraction = first >= 0.0 && first <= 1.0 ? first : second;
  fraction = std::max(0.0, std::min(1.0, fraction));
  return begin + fraction * direction;
}

struct Branch {
  std::vector<Vec3> points;
  bool reachedSun = false;
  bool reachedOuter = false;
};

Core::Result<Branch> TraceBranch(const FieldLineTraceRequest& request,
                                 int sign,
                                 const FieldLineStateEvaluator& evaluator) {
  Branch branch;
  branch.points.push_back(request.seedM);
  for (int step = 0; step < request.maximumStepsPerBranch; ++step) {
    const Vec3 current = branch.points.back();
    auto next = Rk4Step(current, request.timeS, request.nominalStepM,
                        sign, evaluator);
    if (!next.ok())
      return Core::Result<Branch>::Failure(next.status.code,
                                          next.status.message);
    const double currentRadius = Norm(current);
    const double nextRadius = Norm(next.value);
    if (currentRadius > request.solarRadiusM &&
        nextRadius <= request.solarRadiusM) {
      branch.points.push_back(SegmentSphereEvent(
          current, next.value, request.solarRadiusM));
      branch.reachedSun = true;
      return Core::Result<Branch>::Success(branch);
    }
    if (currentRadius < request.outerRadiusM &&
        nextRadius >= request.outerRadiusM) {
      branch.points.push_back(SegmentSphereEvent(
          current, next.value, request.outerRadiusM));
      branch.reachedOuter = true;
      return Core::Result<Branch>::Success(branch);
    }
    branch.points.push_back(next.value);
  }
  return Core::Result<Branch>::Failure(Core::StatusCode::OutOfDomain,
                                      "field-line branch did not reach a boundary");
}

}  // namespace

Core::Result<FieldLine::LineRecord> TraceSunConnectedLine(
    const FieldLineTraceRequest& request,
    const FieldLineStateEvaluator& evaluator) {
  if (request.stableLineId.empty() ||
      !(Norm(request.seedM) > request.solarRadiusM &&
        Norm(request.seedM) < request.outerRadiusM) ||
      !(request.solarRadiusM > 0.0) ||
      !(request.outerRadiusM > request.solarRadiusM) ||
      !(request.nominalStepM > 0.0) || request.maximumStepsPerBranch <= 0 ||
      !(request.unsignedMagneticFluxWb > 0.0) || !evaluator)
    return Core::Result<FieldLine::LineRecord>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "field-line trace request is incomplete or seed is out of domain");
  auto plus = TraceBranch(request, 1, evaluator);
  auto minus = TraceBranch(request, -1, evaluator);
  if (!(plus.ok() && minus.ok()))
    return Core::Result<FieldLine::LineRecord>::Failure(
        Core::StatusCode::OutOfDomain,
        "both-sign trace did not produce Sun and outer branches");
  const Branch* sun = nullptr;
  const Branch* outer = nullptr;
  if (plus.value.reachedSun && minus.value.reachedOuter) {
    sun = &plus.value; outer = &minus.value;
  } else if (minus.value.reachedSun && plus.value.reachedOuter) {
    sun = &minus.value; outer = &plus.value;
  } else {
    return Core::Result<FieldLine::LineRecord>::Failure(
        Core::StatusCode::InvalidState,
        "trace does not contain one Sun-connected and one outer branch");
  }
  std::vector<Vec3> points(sun->points.rbegin(), sun->points.rend());
  points.insert(points.end(), outer->points.begin() + 1, outer->points.end());
  FieldLine::LineRecord line;
  line.stableLineId = request.stableLineId;
  line.open = true;
  line.singleSmoothSector = true;
  line.historyCoverage = FieldLine::HistoryCoverage::Complete;
  line.intersectionPresence = FieldLine::IntersectionPresence::None;
  line.sourceMeasureStatus = FieldLine::SourceMeasureStatus::FiniteTube;
  line.quadratureWeight = request.unsignedMagneticFluxWb;
  line.unsignedMagneticFluxWb = request.unsignedMagneticFluxWb;
  double arcLength = 0.0;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i) arcLength += Norm(points[i] - points[i - 1]);
    auto state = evaluator(points[i], request.timeS);
    if (!state.ok())
      return Core::Result<FieldLine::LineRecord>::Failure(
          state.status.code, state.status.message);
    state.value.positionM = ToNeutral(points[i]);
    state.value.arcLengthM = arcLength;
    Vec3 tangent;
    if (i == 0) tangent = points[1] - points[0];
    else if (i + 1 == points.size()) tangent = points[i] - points[i - 1];
    else tangent = points[i + 1] - points[i - 1];
    state.value.outwardTangent = ToNeutral(Unit(tangent));
    line.nodes.push_back(state.value);
  }
  auto measured = AssignFiniteFluxTubeMeasure(&line,
                                               request.unsignedMagneticFluxWb);
  if (!measured.ok())
    return Core::Result<FieldLine::LineRecord>::Failure(
        measured.code, measured.message);
  return Core::Result<FieldLine::LineRecord>::Success(line);
}

Core::Result<std::vector<FieldLine::FrontIntersection>>
FindFrontIntersections(const FieldLine::LineRecord& line,
                       const ImplicitSurface& surface, double time,
                       std::uint64_t generation, bool eligible,
                       double tolerance) {
  if (!surface || generation == 0 || tolerance <= 0.0 || line.nodes.size() < 2)
    return Core::Result<std::vector<FieldLine::FrontIntersection>>::Failure(
        Core::StatusCode::InvalidConfiguration,
        "front intersection requires line, surface, generation, and tolerance");
  std::vector<FieldLine::FrontIntersection> result;
  for (std::size_t i = 1; i < line.nodes.size(); ++i) {
    Vec3 left = FromNeutral(line.nodes[i - 1].positionM);
    Vec3 right = FromNeutral(line.nodes[i].positionM);
    double fLeft = surface(left), fRight = surface(right);
    if (!std::isfinite(fLeft) || !std::isfinite(fRight))
      return Core::Result<std::vector<FieldLine::FrontIntersection>>::Failure(
          Core::StatusCode::NumericalFailure,
          "front implicit function is nonfinite");
    double fraction = 0.5;
    if (std::abs(fLeft) <= 1.0e-14) {
      fraction = 0.0;
    } else if (std::abs(fRight) <= 1.0e-14) {
      fraction = 1.0;
    } else if (fLeft * fRight > 0.0) {
      // A tangency creates/annihilates two roots without a sign change.  The
      // midpoint check handles an exactly resolved tangent; production time
      // histories bracket it by adaptive subdivision before calling here.
      const double fMiddle = surface(0.5 * (left + right));
      if (std::abs(fMiddle) > 1.0e-12) continue;
    } else {
      double a = 0.0, b = 1.0;
      for (int iteration = 0; iteration < 100 &&
           Norm((b - a) * (right - left)) > tolerance; ++iteration) {
        const double middle = 0.5 * (a + b);
        const double fMiddle = surface(left + middle * (right - left));
        if (fLeft * fMiddle <= 0.0) { b = middle; fRight = fMiddle; }
        else { a = middle; fLeft = fMiddle; }
      }
      fraction = 0.5 * (a + b);
    }
    FieldLine::FrontIntersection intersection;
    std::ostringstream id;
    id << line.stableLineId << '-' << generation << '-' << (i - 1);
    intersection.stableId = id.str();
    intersection.frontGeneration = generation;
    intersection.timeS = time;
    intersection.arcLengthM = line.nodes[i - 1].arcLengthM + fraction *
        (line.nodes[i].arcLengthM - line.nodes[i - 1].arcLengthM);
    intersection.geometricIntersection = true;
    intersection.sourceEligible = eligible;
    if (!result.empty() && std::abs(result.back().arcLengthM -
                                    intersection.arcLengthM) <= tolerance)
      continue;
    result.push_back(intersection);
  }
  return Core::Result<std::vector<FieldLine::FrontIntersection>>::Success(result);
}

Core::Status AssignFiniteFluxTubeMeasure(FieldLine::LineRecord* line,
                                         double flux) {
  if (!line || !(flux > 0.0))
    return Core::Status::Failure(Core::StatusCode::InvalidConfiguration,
                                "finite tube requires line and positive flux");
  for (auto& node : line->nodes) {
    const double field = Norm(FromNeutral(node.magneticFieldT));
    if (!(field > 0.0))
      return Core::Status::Failure(Core::StatusCode::InvalidState,
                                  "finite tube encountered zero magnetic field");
    node.tubeAreaM2 = flux / field;
  }
  line->unsignedMagneticFluxWb = flux;
  line->quadratureWeight = flux;
  line->sourceMeasureStatus = FieldLine::SourceMeasureStatus::FiniteTube;
  return Core::Status::Success();
}

Core::Result<FieldLine::FieldLineSet> BuildFieldLineSet(
    std::string frame, double epoch, double begin, double end,
    std::string rotation, std::vector<FieldLine::LineRecord> lines) {
  std::sort(lines.begin(), lines.end(), [](const auto& a, const auto& b) {
    return a.stableLineId < b.stableLineId;
  });
  FieldLine::FieldLineSet set;
  set.frame = std::move(frame);
  set.epochS = epoch;
  set.historyBeginS = begin;
  set.historyEndS = end;
  set.rotationProvenance = std::move(rotation);
  set.lines = std::move(lines);
  set.bundleId = FieldLine::BundleIdentity(set);
  auto status = FieldLine::ValidateFieldLineSet(set);
  if (!status.ok())
    return Core::Result<FieldLine::FieldLineSet>::Failure(status.code,
                                                         status.message);
  return Core::Result<FieldLine::FieldLineSet>::Success(set);
}

} }  // namespace SEP::CoronalCME
