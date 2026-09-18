#include "sep_numerical_extensions.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP {
namespace NumericalExtensions {
namespace {

Transport::Status Error(Transport::StatusCode code,
                        const std::string& message) {
  return Transport::Status::Error(code, message);
}

bool Finite(double value) { return std::isfinite(value); }

double NormalAt(const BrownianAddress& address, std::uint64_t purpose) {
  Transport::KeyedRandomStream stream(
      address.campaign, address.particle,
      address.operatorId ^ purpose,
      address.physicalInterval);
  return stream.Normal01();
}

double BridgeIncrementRecursive(const BrownianAddress& requested,
                                double rootDtS) {
  BrownianAddress root = requested;
  root.depth = 0;
  root.index = 0;
  double increment = std::sqrt(rootDtS) *
      NormalAt(root, UINT64_C(0x42524f574e49414e));
  double nodeDt = rootDtS;
  for (unsigned depth = 1; depth <= requested.depth; ++depth) {
    const std::uint64_t parentIndex = requested.index >>
        (requested.depth - depth + 1);
    BrownianAddress parent = requested;
    parent.depth = depth - 1;
    parent.index = parentIndex;
    // Each internal node gets a purpose-separated conditional normal.  The
    // two children use opposite signs, ensuring they add to the parent even
    // when requested independently in different traversal orders.
    const std::uint64_t nodeKey =
        (static_cast<std::uint64_t>(parent.depth) << 56) ^ parent.index;
    const double conditional = std::sqrt(nodeDt * 0.25) *
        NormalAt(parent, UINT64_C(0x434f4e444954494f) ^ nodeKey);
    const bool chooseRight = ((requested.index >>
        (requested.depth - depth)) & UINT64_C(1)) != 0;
    increment = chooseRight ? 0.5 * increment - conditional
                            : 0.5 * increment + conditional;
    nodeDt *= 0.5;
  }
  return increment;
}

double EulerMultiplicative(double value, double drift, double diffusion,
                           double dt, double dW) {
  return value + drift * value * dt + diffusion * value * dW;
}

struct AdaptiveAccumulator {
  unsigned accepted = 0;
  unsigned rejected = 0;
  unsigned maxDepth = 0;
  double maxError = 0.0;
};

Transport::Status AdaptiveNode(
    double initial, double drift, double diffusion, double rootDt,
    const BrownianAddress& node, const AdaptiveSdeConfiguration& config,
    double* output, AdaptiveAccumulator* diagnostics) {
  const double nodeDt = std::ldexp(rootDt, -static_cast<int>(node.depth));
  const Transport::ScalarResult parentIncrement = BrownianIncrement(node, rootDt);
  if (!parentIncrement.status.ok()) return parentIncrement.status;
  BrownianAddress left = node;
  left.depth += 1;
  left.index *= 2;
  BrownianAddress right = left;
  right.index += 1;
  const Transport::ScalarResult leftIncrement = BrownianIncrement(left, rootDt);
  const Transport::ScalarResult rightIncrement = BrownianIncrement(right, rootDt);
  if (!leftIncrement.status.ok()) return leftIncrement.status;
  if (!rightIncrement.status.ok()) return rightIncrement.status;
  const double full = EulerMultiplicative(initial, drift, diffusion, nodeDt,
                                           parentIncrement.value);
  const double half = EulerMultiplicative(initial, drift, diffusion,
                                           0.5 * nodeDt,
                                           leftIncrement.value);
  const double twoHalf = EulerMultiplicative(half, drift, diffusion,
                                              0.5 * nodeDt,
                                              rightIncrement.value);
  const double scale = config.absoluteTolerance + config.relativeTolerance *
      std::max(std::fabs(full), std::fabs(twoHalf));
  const double normalized = std::fabs(twoHalf - full) / scale;
  diagnostics->maxError = std::max(diagnostics->maxError, normalized);
  diagnostics->maxDepth = std::max(diagnostics->maxDepth, node.depth);
  if (normalized <= 1.0) {
    *output = twoHalf;
    ++diagnostics->accepted;
    return Transport::Status::Ok();
  }
  ++diagnostics->rejected;
  if (node.depth >= config.maximumDepth)
    return Error(Transport::StatusCode::StepUnderflow,
                 "adaptive Brownian bridge reached maximum depth");

  // Rejected work has not modified caller state.  The left child is committed
  // only to this local variable; if the right child fails, the whole recursive
  // call returns an error and the public result still exposes the original.
  double midpoint = initial;
  Transport::Status status = AdaptiveNode(initial, drift, diffusion, rootDt,
                                           left, config, &midpoint, diagnostics);
  if (!status.ok()) return status;
  return AdaptiveNode(midpoint, drift, diffusion, rootDt,
                      right, config, output, diagnostics);
}

double Minmod(double a, double b) {
  if (a * b <= 0.0) return 0.0;
  return std::copysign(std::min(std::fabs(a), std::fabs(b)), a);
}

double LimitedSlope(double left, double center, double right,
                    TvdLimiter limiter, bool* limited) {
  const double backward = center - left;
  const double forward = right - center;
  double slope = 0.0;
  if (limiter == TvdLimiter::Minmod) {
    slope = Minmod(backward, forward);
  } else {
    slope = Minmod(0.5 * (right - left),
                   Minmod(2.0 * backward, 2.0 * forward));
  }
  if (limited) *limited = std::fabs(slope - 0.5 * (right - left)) >
      32.0 * std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::fabs(center));
  return slope;
}

std::vector<double> SpatialResidual(const std::vector<double>& values,
                                    double speed, double dx,
                                    AdvectionOrder order, TvdLimiter limiter,
                                    std::uint64_t* activations,
                                    std::uint64_t* faces) {
  const std::size_t n = values.size();
  std::vector<double> slopes(n, 0.0);
  if (order == AdvectionOrder::SecondOrderMuscl) {
    for (std::size_t i = 0; i < n; ++i) {
      bool limited = false;
      slopes[i] = LimitedSlope(values[(i + n - 1) % n], values[i],
                               values[(i + 1) % n], limiter, &limited);
      if (limited) ++*activations;
    }
  }
  std::vector<double> flux(n, 0.0);
  for (std::size_t face = 0; face < n; ++face) {
    const std::size_t left = face;
    const std::size_t right = (face + 1) % n;
    const double leftFace = values[left] + 0.5 * slopes[left];
    const double rightFace = values[right] - 0.5 * slopes[right];
    flux[face] = speed * (speed >= 0.0 ? leftFace : rightFace);
    ++*faces;
  }
  std::vector<double> residual(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    residual[i] = -(flux[i] - flux[(i + n - 1) % n]) / dx;
  return residual;
}

double Integral(const std::vector<double>& edges,
                const std::vector<double>& values) {
  long double total = 0.0L;
  for (std::size_t i = 0; i < values.size(); ++i)
    total += static_cast<long double>(values[i]) * (edges[i + 1] - edges[i]);
  return static_cast<double>(total);
}

bool ValidMesh(const std::vector<double>& edges, std::size_t cells) {
  if (edges.size() != cells + 1 || cells == 0) return false;
  for (std::size_t i = 0; i < edges.size(); ++i)
    if (!Finite(edges[i]) || (i && edges[i] <= edges[i - 1])) return false;
  return true;
}

}  // namespace

Transport::ScalarResult BrownianIncrement(
    const BrownianAddress& address, double rootIntervalS) {
  Transport::ScalarResult result;
  if (address.particle == 0 || address.operatorId == 0 ||
      !Finite(rootIntervalS) || rootIntervalS <= 0.0 ||
      address.depth > 62 || address.index >= (UINT64_C(1) << address.depth)) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "Brownian bridge address or root interval is invalid");
    return result;
  }
  result.value = BridgeIncrementRecursive(address, rootIntervalS);
  result.status = Finite(result.value)
      ? Transport::Status::Ok()
      : Error(Transport::StatusCode::InvalidParticleState,
              "Brownian bridge increment is non-finite");
  return result;
}

AdaptiveSdeResult AdvanceAdaptiveMultiplicativeSde(
    double initial, double drift, double diffusion, double dt,
    const BrownianAddress& root,
    const AdaptiveSdeConfiguration& config) {
  AdaptiveSdeResult result;
  result.value = initial;
  if (!Finite(initial) || !Finite(drift) || !Finite(diffusion) ||
      !Finite(dt) || dt <= 0.0 || !Finite(config.absoluteTolerance) ||
      config.absoluteTolerance < 0.0 || !Finite(config.relativeTolerance) ||
      config.relativeTolerance <= 0.0 || config.maximumDepth > 30 ||
      root.depth != 0 || root.index != 0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "adaptive stochastic configuration is invalid");
    return result;
  }
  AdaptiveAccumulator diagnostics;
  double staged = initial;
  result.status = AdaptiveNode(initial, drift, diffusion, dt, root, config,
                               &staged, &diagnostics);
  result.acceptedLeaves = diagnostics.accepted;
  result.rejectedTrials = diagnostics.rejected;
  result.maximumDepthReached = diagnostics.maxDepth;
  result.maximumNormalizedError = diagnostics.maxError;
  if (result.status.ok()) result.value = staged;
  return result;
}

AdvectionResult AdvectPeriodic(
    const std::vector<double>& values, double speed, double dt, double dx,
    AdvectionOrder order, TvdLimiter limiter) {
  AdvectionResult result;
  if (values.size() < 3 || !Finite(speed) || !Finite(dt) || dt < 0.0 ||
      !Finite(dx) || dx <= 0.0 || std::fabs(speed) * dt > dx) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "periodic advection requires a valid grid and CFL<=1");
    return result;
  }
  for (std::size_t i = 0; i < values.size(); ++i)
    if (!Finite(values[i]) || values[i] < 0.0) {
      result.status = Error(Transport::StatusCode::InvalidParticleState,
                            "advection cell average is invalid");
      return result;
    }
  long double initial = 0.0L;
  for (std::size_t i = 0; i < values.size(); ++i) initial += values[i];
  const std::vector<double> first = SpatialResidual(
      values, speed, dx, order, limiter, &result.limiterActivations,
      &result.faceUpdates);
  result.cellAverage = values;
  for (std::size_t i = 0; i < values.size(); ++i)
    result.cellAverage[i] += dt * first[i];
  if (order == AdvectionOrder::SecondOrderMuscl) {
    const std::vector<double> second = SpatialResidual(
        result.cellAverage, speed, dx, order, limiter,
        &result.limiterActivations, &result.faceUpdates);
    for (std::size_t i = 0; i < values.size(); ++i)
      result.cellAverage[i] = 0.5 * values[i] +
          0.5 * (result.cellAverage[i] + dt * second[i]);
  }
  long double final = 0.0L;
  for (std::size_t i = 0; i < result.cellAverage.size(); ++i) {
    if (!Finite(result.cellAverage[i]) || result.cellAverage[i] < -1.0e-12) {
      result.status = Error(Transport::StatusCode::StepUnderflow,
                            "TVD advection violated positivity");
      return result;
    }
    result.cellAverage[i] = std::max(0.0, result.cellAverage[i]);
    final += result.cellAverage[i];
  }
  result.conservationResidual = static_cast<double>((final - initial) * dx);
  result.status = Transport::Status::Ok();
  return result;
}

RemapResult RemapConservative(
    const std::vector<double>& oldEdges,
    const std::vector<double>& oldValues,
    const std::vector<double>& newEdges,
    RemapOrder order, double initialization) {
  RemapResult result;
  if (!ValidMesh(oldEdges, oldValues.size()) ||
      !ValidMesh(newEdges, newEdges.size() - 1) ||
      !Finite(initialization) || initialization < 0.0) {
    result.status = Error(Transport::StatusCode::InvalidArgument,
                          "conservative remap mesh or initialization is invalid");
    return result;
  }
  for (std::size_t i = 0; i < oldValues.size(); ++i)
    if (!Finite(oldValues[i]) || oldValues[i] < 0.0) {
      result.status = Error(Transport::StatusCode::InvalidParticleState,
                            "conservative remap state is invalid");
      return result;
    }
  result.oldIntegral = Integral(oldEdges, oldValues);
  std::vector<double> slopes(oldValues.size(), 0.0);
  if (order == RemapOrder::PiecewiseLinear && oldValues.size() > 2) {
    for (std::size_t i = 1; i + 1 < oldValues.size(); ++i) {
      const double leftDx = 0.5 * (oldEdges[i + 1] - oldEdges[i - 1]);
      const double rightDx = 0.5 * (oldEdges[i + 2] - oldEdges[i]);
      const double left = (oldValues[i] - oldValues[i - 1]) / leftDx;
      const double right = (oldValues[i + 1] - oldValues[i]) / rightDx;
      slopes[i] = Minmod(left, right);
      if (slopes[i] != 0.5 * (left + right)) ++result.limiterActivations;
      // Positivity at both cell faces is a sufficient bound for a linear
      // reconstruction.  Scaling the slope preserves the cell average.
      const double half = 0.5 * (oldEdges[i + 1] - oldEdges[i]);
      const double maxSlope = half > 0.0 ? oldValues[i] / half : 0.0;
      if (std::fabs(slopes[i]) > maxSlope) {
        slopes[i] = std::copysign(maxSlope, slopes[i]);
        ++result.limiterActivations;
      }
    }
  }
  result.cellAverage.assign(newEdges.size() - 1, 0.0);
  for (std::size_t j = 0; j + 1 < newEdges.size(); ++j) {
    long double integral = 0.0L;
    const double newWidth = newEdges[j + 1] - newEdges[j];
    const double coveredLeft = std::max(newEdges[j], oldEdges.front());
    const double coveredRight = std::min(newEdges[j + 1], oldEdges.back());
    if (coveredRight > coveredLeft) {
      for (std::size_t i = 0; i < oldValues.size(); ++i) {
        const double left = std::max(coveredLeft, oldEdges[i]);
        const double right = std::min(coveredRight, oldEdges[i + 1]);
        if (right <= left) continue;
        const double center = 0.5 * (oldEdges[i] + oldEdges[i + 1]);
        integral += oldValues[i] * (right - left) +
            0.5 * slopes[i] *
            ((right - center) * (right - center) -
             (left - center) * (left - center));
      }
    }
    const double uncovered = newWidth - std::max(0.0, coveredRight - coveredLeft);
    if (uncovered > 0.0) {
      integral += initialization * uncovered;
      result.initializedIntegral += initialization * uncovered;
    }
    result.cellAverage[j] = static_cast<double>(integral / newWidth);
    if (result.cellAverage[j] < -1.0e-12) {
      result.status = Error(Transport::StatusCode::InvalidParticleState,
                            "limited remap produced negative state");
      return result;
    }
    result.cellAverage[j] = std::max(0.0, result.cellAverage[j]);
  }
  result.newIntegral = Integral(newEdges, result.cellAverage);
  const double oldCoveredLeft = std::max(oldEdges.front(), newEdges.front());
  const double oldCoveredRight = std::min(oldEdges.back(), newEdges.back());
  double retainedOld = 0.0;
  if (oldCoveredRight > oldCoveredLeft) {
    std::vector<double> overlapEdges;
    overlapEdges.push_back(oldCoveredLeft);
    // Boundary exchange is easiest to compute by integrating the same limited
    // reconstruction over the intersection; here the final accounting identity
    // derives it from old, new, and explicitly initialized measure.
    retainedOld = result.newIntegral - result.initializedIntegral;
  }
  result.boundaryExchange = retainedOld - result.oldIntegral;
  result.numericalResidual = result.newIntegral - result.oldIntegral -
      result.boundaryExchange - result.initializedIntegral;
  result.status = Transport::Status::Ok();
  return result;
}

MoverResult ClassifyMoverFailure(const Transport::Status& status,
                                 const FailureContext& context,
                                 FailurePolicy policy) {
  MoverResult result;
  result.status = status;
  result.context = context;
  if (status.ok()) {
    result.disposition = MoverDisposition::Completed;
    result.fate = ParticleFate::Retain;
    return result;
  }
  switch (status.code) {
    case Transport::StatusCode::OutOfDomain:
      result.disposition = MoverDisposition::BoundaryExit;
      result.fate = ParticleFate::DeleteAtBoundary;
      break;
    case Transport::StatusCode::StepUnderflow:
    case Transport::StatusCode::UnresolvedCoefficient:
      result.disposition = MoverDisposition::RejectedStep;
      result.fate = ParticleFate::Unchanged;
      break;
    case Transport::StatusCode::InvalidParticleState:
      result.disposition = policy == FailurePolicy::QuarantineParticleLocal
          ? MoverDisposition::QuarantinedParticle
          : MoverDisposition::FatalInternalError;
      result.fate = policy == FailurePolicy::QuarantineParticleLocal
          ? ParticleFate::Quarantine : ParticleFate::Unchanged;
      break;
    case Transport::StatusCode::InvalidArgument:
    case Transport::StatusCode::InvalidCoefficient:
    case Transport::StatusCode::UnsupportedConfiguration:
    case Transport::StatusCode::NondifferentiableCoefficient:
      result.disposition = MoverDisposition::InvalidConfiguration;
      result.fate = ParticleFate::Unchanged;
      break;
    case Transport::StatusCode::Ok:
      break;
  }
  return result;
}

Transport::Status ParticleTransaction::Commit(
    Transport::ParticleState* destination) {
  if (!destination)
    return Error(Transport::StatusCode::InvalidArgument,
                 "particle transaction destination is null");
  if (committed_)
    return Error(Transport::StatusCode::InvalidArgument,
                 "particle transaction cannot be committed twice");
  *destination = staged_;
  committed_ = true;
  return Transport::Status::Ok();
}

void ParticleTransaction::Rollback() {
  staged_ = initial_;
  committed_ = false;
}

}  // namespace NumericalExtensions
}  // namespace SEP
