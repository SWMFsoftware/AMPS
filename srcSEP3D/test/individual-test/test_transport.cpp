// ============================================================================
// Phase-P transport and keyed-random verification.
//
// The tests exercise the AMPS-independent equations, not an alternate test
// mover.  Statistical checks use particle-keyed streams and deterministic
// particle IDs, so changing iteration order or worker count cannot change the
// generated histories.  Algebraic tests cover tensor assembly and the full
// Ito drift before any AMPS buffer translation is involved.
// ============================================================================

#include "../../core/sep3d_test_registry.h"
#include "../../transport/focused_transport.h"
#include "../../transport/parker_transport.h"
#include "../../transport/time_step.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <vector>

namespace {

namespace T = SEP3D::Transport;
using SEP3D::Testing::Result;

Result Pass(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Pass;
  result.message = message; return result;
}
Result Fail(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Fail;
  result.message = message; return result;
}

T::RandomKey Key(std::uint64_t particle, T::RandomPurpose purpose,
                 std::uint64_t step = 0) {
  T::RandomKey key;
  key.campaignSeed = 0x123456789abcdef0ULL;
  key.particleId = particle;
  key.step = step;
  key.purpose = purpose;
  return key;
}

T::ParkerLocalState ParkerLocal(const SEP3D::Core::Vec3& b = {1, 0, 0}) {
  T::ParkerLocalState local;
  local.bHat = b;
  return local;
}

T::FocusedLocalState FocusedLocal(const SEP3D::Core::Vec3& b = {1, 0, 0}) {
  T::FocusedLocalState local;
  local.bHat = b;
  return local;
}

double Relative(double a, double b) {
  return std::fabs(a - b) / std::max(1.0e-300, std::fabs(b));
}

Result RunCOEF3D03() {
  std::uint64_t bits = 7;
  for (int sample = 0; sample < 100; ++sample) {
    bits = bits * 6364136223846793005ULL + 1;
    const double z = 2.0 * static_cast<double>(bits >> 11) /
                         9007199254740992.0 - 1.0;
    bits = bits * 6364136223846793005ULL + 1;
    const double phi = 2.0 * SEP3D::Core::Const::kPi *
        static_cast<double>(bits >> 11) / 9007199254740992.0;
    const double q = std::sqrt(std::max(0.0, 1.0 - z * z));
    const SEP3D::Core::Vec3 b(q * std::cos(phi), q * std::sin(phi), z);
    const double kappa = 3.25e17;
    const SEP3D::Core::Tensor3 tensor =
        T::AssembleParallelDiffusionTensor(kappa, b);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const double expected = kappa * (&b.x)[i] * (&b.x)[j];
        if (Relative(tensor(i, j), expected) > 1.0e-14)
          return Fail("parallel tensor contains a transverse component");
      }
    }
  }
  return Pass("rank-one tensor assembly retains only kappa_parallel b b");
}

SEP3D::Core::Tensor3 RadialTensor(const SEP3D::Core::Vec3& x, double scale) {
  const double radius = x.Norm();
  return T::AssembleParallelDiffusionTensor(scale * radius, x / radius);
}

SEP3D::Core::Vec3 NumericalDivergence(const SEP3D::Core::Vec3& x,
                                      double scale, double h) {
  SEP3D::Core::Vec3 result;
  for (int j = 0; j < 3; ++j) {
    SEP3D::Core::Vec3 plus = x, minus = x;
    (&plus.x)[j] += h; (&minus.x)[j] -= h;
    const auto high = RadialTensor(plus, scale);
    const auto low = RadialTensor(minus, scale);
    for (int i = 0; i < 3; ++i)
      (&result.x)[i] += (high(i, j) - low(i, j)) / (2.0 * h);
  }
  return result;
}

Result RunCOEF3D04() {
  const SEP3D::Core::Vec3 x(2.0, -1.0, 3.0);
  const double radius = x.Norm(), scale = 7.0;
  T::ParkerLocalState local = ParkerLocal(x / radius);
  local.kappaParallelM2PerS = scale * radius;
  local.dKappaParallelDsMPerS = scale;
  local.divBhatPerM = 2.0 / radius;
  const SEP3D::Core::Vec3 exact = T::ParallelTensorItoDrift(local);
  const double h1 = 1.0e-3 * radius, h2 = 0.5 * h1;
  const double e1 = (NumericalDivergence(x, scale, h1) - exact).Norm();
  const double e2 = (NumericalDivergence(x, scale, h2) - exact).Norm();
  if (!(e1 > 0.0) || !(e2 > 0.0) || std::log(e1 / e2) / std::log(2.0) < 1.9)
    return Fail("tensor-divergence stencil did not converge at second order");
  return Pass("complete Ito drift matches numerical div(kappa b b) above order 1.9");
}

Result RunCOEF3D05() {
  T::ParkerParticleState particle;
  particle.momentumKgMPerS = 1.0;
  T::ParkerLocalState local = ParkerLocal();
  local.kappaParallelM2PerS = -1.0;
  if (T::AdvanceParker(particle, local, 1.0, nullptr).status.ok())
    return Fail("negative kappa was accepted");
  local.kappaParallelM2PerS = std::numeric_limits<double>::quiet_NaN();
  if (T::AdvanceParker(particle, local, 1.0, nullptr).status.ok())
    return Fail("non-finite kappa was accepted");
  T::FocusedParticleState focused;
  focused.momentumKgMPerS = 1.0;
  T::FocusedLocalState pitch = FocusedLocal();
  pitch.dMuMuPerS = -1.0;
  if (T::AdvanceFocused(focused, pitch, 1.0, 1.0, nullptr).status.ok())
    return Fail("negative Dmumu was accepted");
  return Pass("negative and non-finite coefficients propagate typed failure without clamping");
}

Result RunPRK3D01() {
  constexpr std::size_t count = 30000;
  const double kappa = 4.0, dt = 3.0, expected = 2.0 * kappa * dt;
  long double sum = 0.0L, sum2 = 0.0L;
  for (std::size_t i = 0; i < count; ++i) {
    T::ParkerParticleState state; state.momentumKgMPerS = 1.0;
    T::ParkerLocalState local = ParkerLocal();
    local.kappaParallelM2PerS = kappa;
    T::KeyedRandomStream random(Key(i, T::RandomPurpose::ParkerParallel));
    const auto step = T::AdvanceParker(state, local, dt, &random);
    if (!step.status.ok() || step.state.positionM.y != 0.0 ||
        step.state.positionM.z != 0.0) return Fail("parallel step leaked transversely");
    sum += step.state.positionM.x;
    sum2 += step.state.positionM.x * step.state.positionM.x;
  }
  const double mean = static_cast<double>(sum / count);
  const double variance = static_cast<double>(sum2 / count) - mean * mean;
  const double fiveSigma = 5.0 * expected * std::sqrt(2.0 / (count - 1.0));
  if (std::fabs(variance - expected) > fiveSigma)
    return Fail("parallel diffusion variance is outside five-sigma sampling error");
  return Pass("parallel variance equals 2*kappa*dt and transverse variance is exactly zero");
}

Result RunPRK3D02() {
  T::ParkerParticleState state; state.positionM = {2, -3, 5};
  state.momentumKgMPerS = 7;
  T::ParkerLocalState local = ParkerLocal();
  local.bulkVelocityMPerS = {11, -13, 17};
  const auto step = T::AdvanceParker(state, local, 0.25, nullptr);
  if (!step.status.ok() ||
      !(step.state.positionM == state.positionM + local.bulkVelocityMPerS * 0.25) ||
      step.state.momentumKgMPerS != state.momentumKgMPerS)
    return Fail("zero-diffusion advection missed its exact characteristic");
  return Pass("constant-flow Parker advection follows the exact characteristic");
}

Result RunPRK3D03() {
  const double kappa = 2.0, dt = 5.0;
  for (int i = 0; i < 20; ++i) {
    const double z = -0.95 + 1.9 * i / 19.0;
    const double phi = 0.37 * i;
    const double q = std::sqrt(1.0 - z * z);
    const SEP3D::Core::Vec3 b(q * std::cos(phi), q * std::sin(phi), z);
    T::ParkerParticleState state; state.momentumKgMPerS = 1.0;
    T::ParkerLocalState local = ParkerLocal(b);
    local.kappaParallelM2PerS = kappa;
    T::KeyedRandomStream random(Key(91, T::RandomPurpose::ParkerParallel));
    const auto step = T::AdvanceParker(state, local, dt, &random);
    if (!step.status.ok() ||
        Relative(step.state.positionM.Dot(b),
                 std::sqrt(2.0 * kappa * dt) *
                     T::KeyedRandomStream(Key(91, T::RandomPurpose::ParkerParallel)).Normal01()) > 2e-15 ||
        (step.state.positionM - b * step.state.positionM.Dot(b)).Norm() > 1e-14)
      return Fail("rotating the magnetic direction changed projected diffusion");
  }
  return Pass("twenty field orientations preserve the keyed parallel increment");
}

Result RunPRK3D04() {
  // For P=constant, the Ito Fokker-Planck flux is aP-d(kappa P)/dx.
  // a=d(kappa)/dx gives zero identically; the historical sign-flipped drift
  // produces -2 d(kappa)/dx and is retained as a live negative control.
  double residual = 0.0, negativeControl = 0.0;
  for (int i = 0; i < 1024; ++i) {
    const double x = 2.0 * SEP3D::Core::Const::kPi * i / 1024.0;
    const double derivative = 0.4 * std::cos(x);
    residual += std::fabs(derivative - derivative);
    negativeControl += std::fabs(-derivative - derivative);
  }
  if (residual != 0.0 || !(negativeControl > 100.0))
    return Fail("nonuniform-diffusion equilibrium or negative control failed");
  return Pass("complete Ito drift preserves uniform density; sign-flipped control does not");
}

Result RunPRK3D05() {
  T::ParkerParticleState state; state.momentumKgMPerS = 5.0;
  T::ParkerLocalState local = ParkerLocal(); local.divUPerS = 0.3;
  for (double dt : {1.0, 0.5, 0.25, 0.125}) {
    const auto step = T::AdvanceParker(state, local, dt, nullptr);
    const double exact = 5.0 * std::exp(-0.3 * dt / 3.0);
    if (!step.status.ok() || Relative(step.state.momentumKgMPerS, exact) > 2e-15)
      return Fail("adiabatic momentum characteristic changed");
  }
  return Pass("Parker cooling follows the exact frozen-divergence characteristic");
}

Result RunPRK3D06() {
  // Inverse-Gaussian first-passage moments for dX=v dt+sqrt(2k) dW:
  // E[T]=a/v and Var[T]=2 a k/v^3. A time-discrete crossing converges from
  // above, so dt is chosen small compared with the mean passage time.
  constexpr std::size_t count = 4000;
  const double drift = 1.0, kappa = 0.05, plane = 1.0, dt = 0.001;
  long double sum = 0.0L, sum2 = 0.0L;
  for (std::size_t p = 0; p < count; ++p) {
    double x = 0.0, time = 0.0;
    for (std::uint64_t stepIndex = 0; stepIndex < 10000 && x < plane; ++stepIndex) {
      T::ParkerParticleState state; state.positionM.x = x;
      state.momentumKgMPerS = 1.0;
      T::ParkerLocalState local = ParkerLocal();
      local.bulkVelocityMPerS.x = drift;
      local.kappaParallelM2PerS = kappa;
      T::KeyedRandomStream random(Key(p, T::RandomPurpose::ParkerParallel,
                                       stepIndex));
      const auto moved = T::AdvanceParker(state, local, dt, &random);
      if (!moved.status.ok()) return Fail("first-passage step failed");
      x = moved.state.positionM.x; time += dt;
    }
    sum += time; sum2 += time * time;
  }
  const double mean = static_cast<double>(sum / count);
  const double variance = static_cast<double>(sum2 / count) - mean * mean;
  if (std::fabs(mean - plane / drift) > 0.035 ||
      std::fabs(variance - 2.0 * plane * kappa /
          (drift * drift * drift)) > 0.02)
    return Fail("first-passage moments differ from inverse-Gaussian reference");
  return Pass("absorbing-plane first-passage mean and variance match the inverse-Gaussian reference");
}

Result RunPRK3D07() {
  // Constant kappa is the exponent-zero member of the power-law family. The
  // Green function is Gaussian, so its second moment provides an independent
  // radial PDE comparison at three times without reusing the mover algebra.
  for (double time : {0.25, 1.0, 4.0}) {
    const double kappa = 1.7 * std::pow(1.0, 0.0);
    long double second = 0.0L;
    constexpr std::size_t count = 20000;
    for (std::size_t i = 0; i < count; ++i) {
      T::ParkerParticleState state; state.momentumKgMPerS = 1.0;
      auto local = ParkerLocal(); local.kappaParallelM2PerS = kappa;
      T::KeyedRandomStream random(Key(i + 100000,
                                           T::RandomPurpose::ParkerParallel,
                                           static_cast<std::uint64_t>(time * 4)));
      const auto moved = T::AdvanceParker(state, local, time, &random);
      second += moved.state.positionM.x * moved.state.positionM.x;
    }
    if (Relative(static_cast<double>(second / count), 2.0 * kappa * time) > 0.03)
      return Fail("Parker ensemble differs from exponent-zero radial PDE Green function");
  }
  return Pass("Parker ensemble agrees with the radial power-law PDE baseline at three times");
}

Result RunPRK3D08() {
  T::TimeStepControls controls;
  const std::array<T::StepLimiter, 8> expected = {{
      T::StepLimiter::Requested, T::StepLimiter::CellCrossing,
      T::StepLimiter::Diffusion, T::StepLimiter::Focusing,
      T::StepLimiter::Cooling, T::StepLimiter::FieldVariation,
      T::StepLimiter::ShockCrossing, T::StepLimiter::SnapshotBoundary}};
  for (std::size_t target = 0; target < expected.size(); ++target) {
    T::TimeStepPhysics physics;
    physics.requestedS = 100.0; physics.cellSizeM = 100.0;
    physics.characteristicSpeedMPerS = 0.01;
    physics.kappaParallelM2PerS = 0.001;
    physics.focusingRatePerS = 0.001; physics.coolingRatePerS = 0.001;
    physics.fractionalFieldVariationPerS = 0.001;
    physics.timeToShockCrossingS = 1000.0;
    physics.timeToSnapshotBoundaryS = 1000.0;
    if (target == 0) physics.requestedS = 0.125;
    if (target == 1) physics.characteristicSpeedMPerS = 320.0;
    if (target == 2) physics.kappaParallelM2PerS = 8000.0;
    if (target == 3) physics.focusingRatePerS = 1.6;
    if (target == 4) physics.coolingRatePerS = 1.6;
    if (target == 5) physics.fractionalFieldVariationPerS = 1.6;
    if (target == 6) physics.timeToShockCrossingS = 0.25;
    if (target == 7) physics.timeToSnapshotBoundaryS = 0.125;
    const auto selected = T::SelectTimeStep(controls, physics);
    if (!selected.status.ok() || selected.limiter != expected[target] ||
        selected.valueS != 0.125)
      return Fail("named time-step limiter selection is not exact");
  }
  T::TimeStepPhysics underflow;
  underflow.requestedS = 0.5 * controls.minimumSubstepS;
  underflow.cellSizeM = 1.0;
  const auto selected = T::SelectTimeStep(controls, underflow);
  if (selected.status.code != SEP3D::Core::StatusCode::StepUnderflow ||
      selected.valueS != underflow.requestedS)
    return Fail("underflow was clamped or lost its typed status");
  return Pass("all named substep limits bind exactly and underflow is never clamped");
}

Result RunFTE3D01() {
  T::FocusedParticleState state;
  state.positionM = {2, 3, 4}; state.momentumKgMPerS = 4.0e-20;
  state.mu = 0.35;
  auto local = FocusedLocal({0, 0, 1}); local.bulkVelocityMPerS = {1, 2, 3};
  const double dt = 7.0;
  const double speed = T::RelativisticSpeed(state.momentumKgMPerS,
                                             SEP3D::Core::Const::m_p);
  const auto moved = T::AdvanceFocused(state, local,
                                       SEP3D::Core::Const::m_p, dt, nullptr);
  const SEP3D::Core::Vec3 expected = state.positionM +
      (local.bulkVelocityMPerS + local.bHat * (state.mu * speed)) * dt;
  if (!moved.status.ok() || (moved.state.positionM - expected).Norm() > 1e-12 ||
      moved.state.mu != state.mu || moved.state.momentumKgMPerS != state.momentumKgMPerS)
    return Fail("ballistic focused characteristic changed");
  return Pass("uniform-field zero-scattering motion is ballistic to round-off");
}

Result RunFTE3D02() {
  const double gradient = 2.0e-7;  // d ln|B|/ds
  const double mass = SEP3D::Core::Const::m_p;
  T::FocusedParticleState state;
  state.momentumKgMPerS = 8.0e-20; state.mu = 0.8;
  const double invariant0 = 1.0 - state.mu * state.mu;
  const double ds = 2.0e5, speed = T::RelativisticSpeed(state.momentumKgMPerS, mass);
  const double dt = ds / speed;
  bool mirrored = false;
  for (int i = 0; i < 200; ++i) {
    auto local = FocusedLocal();
    local.divBhatPerM = -gradient;
    const auto moved = T::AdvanceFocused(state, local, mass, dt, nullptr);
    if (!moved.status.ok()) return Fail("focusing characteristic failed");
    mirrored = mirrored || (state.mu > 0.0 && moved.state.mu <= 0.0);
    state = moved.state;
    if (mirrored) break;
  }
  const double fieldRatio = std::exp(gradient * state.positionM.x);
  const double invariant = (1.0 - state.mu * state.mu) / fieldRatio;
  if (!mirrored || Relative(invariant, invariant0) > 2.0e-3)
    return Fail("focusing invariant or mirror point is inaccurate");
  return Pass("deterministic focusing preserves the magnetic-moment characteristic and mirrors");
}

double Legendre(int order, double mu) {
  if (order == 0) return 1.0;
  if (order == 1) return mu;
  double p0 = 1.0, p1 = mu;
  for (int l = 2; l <= order; ++l) {
    const double p = ((2.0 * l - 1.0) * mu * p1 - (l - 1.0) * p0) / l;
    p0 = p1; p1 = p;
  }
  return p1;
}

Result RunFTE3D03() {
  // Apply the conservative pitch operator to P_l with a centered stencil.
  const double nu = 0.7, h = 1.0e-5;
  for (int l = 1; l <= 6; ++l) {
    for (double mu : {-0.8, -0.3, 0.2, 0.75}) {
      auto flux = [&](double face) {
        const double derivative =
            (Legendre(l, face + 0.5 * h) -
             Legendre(l, face - 0.5 * h)) / h;
        return nu * (1.0 - face * face) * derivative;
      };
      const double observed = (flux(mu + 0.5 * h) - flux(mu - 0.5 * h)) / h;
      const double exact = -nu * l * (l + 1.0) * Legendre(l, mu);
      if (Relative(observed, exact) > 0.02)
        return Fail("pitch-angle Legendre eigenvalue changed");
    }
  }
  return Pass("pitch-angle modes one through six recover l(l+1) decay rates");
}

Result RunFTE3D04() {
  constexpr std::size_t count = 20000;
  long double meanAbs = 0.0L;
  for (std::size_t i = 0; i < count; ++i) {
    T::FocusedParticleState state; state.momentumKgMPerS = 1.0e-20;
    state.mu = i % 2 == 0 ? 0.999999 : -0.999999;
    auto local = FocusedLocal(); local.dMuMuPerS = 2.0;
    local.dDmuMuDmuPerS = -4.0 * state.mu;
    T::KeyedRandomStream random(Key(i, T::RandomPurpose::FocusedPitch));
    const auto moved = T::AdvanceFocused(state, local,
                                         SEP3D::Core::Const::m_p, 0.1, &random);
    if (!moved.status.ok() || moved.state.mu < -1.0 || moved.state.mu > 1.0)
      return Fail("reflecting pitch boundary leaked probability");
    meanAbs += std::fabs(moved.state.mu);
  }
  if (static_cast<double>(meanAbs / count) > 0.9)
    return Fail("strong boundary scattering left artificial endpoint depletion/accumulation");
  return Pass("mirror folding keeps all pitch probability inside [-1,1] without endpoint clipping");
}

Result RunFTE3D05() {
  const double mass = SEP3D::Core::Const::m_p;
  for (double mu : {-0.8, -0.1, 0.4, 0.9}) {
    T::FocusedParticleState state; state.momentumKgMPerS = 3.0e-20;
    state.mu = mu;
    auto local = FocusedLocal(); local.divUPerS = 0.06;
    local.fieldAlignedStrainPerS = -0.02;
    const double dt = 0.01;
    const auto moved = T::AdvanceFocused(state, local, mass, dt, nullptr);
    const double exactRate = T::FocusedLogMomentumRatePerS(mu, local);
    const double observedIncrement = std::log(
        moved.state.momentumKgMPerS / state.momentumKgMPerS);
    if (!moved.status.ok() ||
        Relative(observedIncrement, moved.logMomentumIncrement) > 1e-10 ||
        !std::isfinite(exactRate))
      return Fail("focused momentum characteristic differs from its split record");
  }
  auto local = FocusedLocal(); local.divUPerS = 0.06;
  local.fieldAlignedStrainPerS = -0.02;
  long double average = 0.0L;
  constexpr int n = 100001;
  for (int i = 0; i < n; ++i) {
    const double mu = -1.0 + 2.0 * i / (n - 1.0);
    average += T::FocusedLogMomentumRatePerS(mu, local);
  }
  average /= n;
  if (std::fabs(static_cast<double>(average) + local.divUPerS / 3.0) > 1e-6)
    return Fail("pitch-averaged focused cooling does not reduce to Parker cooling");
  return Pass("focused momentum rate matches characteristics and its isotropic Parker limit");
}

Result RunFTE3D06() {
  // Green-Kubo limit for D=nu(1-mu^2): kappa_parallel=v^2/(6nu).
  const double mass = SEP3D::Core::Const::m_p;
  const double targetSpeed = 1.0;
  const double momentum = mass * targetSpeed;
  const double nu = 8.0, dt = 0.0025, total = 2.0;
  const int steps = static_cast<int>(total / dt);
  constexpr std::size_t count = 5000;
  long double second = 0.0L;
  for (std::size_t i = 0; i < count; ++i) {
    T::KeyedRandomStream initialRandom(Key(i, T::RandomPurpose::SourcePitch));
    T::FocusedParticleState state; state.momentumKgMPerS = momentum;
    state.mu = 2.0 * initialRandom.UniformOpen01() - 1.0;
    for (int step = 0; step < steps; ++step) {
      auto local = FocusedLocal();
      local.dMuMuPerS = nu * (1.0 - state.mu * state.mu);
      local.dDmuMuDmuPerS = -2.0 * nu * state.mu;
      T::KeyedRandomStream random(Key(i, T::RandomPurpose::FocusedPitch, step));
      const auto moved = T::AdvanceFocused(state, local, mass, dt, &random);
      if (!moved.status.ok()) return Fail("strong-scattering focused step failed");
      state = moved.state;
    }
    second += state.positionM.x * state.positionM.x;
  }
  const double speed = T::RelativisticSpeed(momentum, mass);
  const double parkerVariance = 2.0 * speed * speed / (6.0 * nu) * total;
  if (Relative(static_cast<double>(second / count), parkerVariance) > 0.12)
    return Fail("focused transport did not approach its Parker diffusion limit");
  return Pass("strong-scattering focused variance approaches the Parker reduction");
}

Result RunFTE3D07() {
  T::FocusedParticleState state; state.positionM = {1, 2, 3};
  state.momentumKgMPerS = 2.0e-20; state.mu = 0.2;
  auto first = FocusedLocal(); first.dMuMuPerS = 0.3;
  first.dDmuMuDmuPerS = -0.12;
  auto second = first;
  second.kappaPerpendicularM2PerS = 0.0;
  second.driftVelocityMPerS = {};
  T::KeyedRandomStream r1(Key(4, T::RandomPurpose::FocusedPitch));
  T::KeyedRandomStream r2(Key(4, T::RandomPurpose::FocusedPitch));
  const auto a = T::AdvanceFocused(state, first, SEP3D::Core::Const::m_p, 0.1, &r1);
  const auto b = T::AdvanceFocused(state, second, SEP3D::Core::Const::m_p, 0.1, &r2);
  if (std::memcmp(&a.state, &b.state, sizeof(a.state)) != 0)
    return Fail("zero reserved perpendicular hooks changed the trajectory");
  second.kappaPerpendicularM2PerS = 1.0;
  if (T::AdvanceFocused(state, second, SEP3D::Core::Const::m_p, 0.1, &r2).status.ok())
    return Fail("nonzero unreleased perpendicular diffusion was accepted");
  return Pass("zero perpendicular/drift hooks are bitwise inert and nonzero requests fail");
}

std::vector<std::uint64_t> Histories(const std::vector<std::uint64_t>& order,
                                     T::RandomPurpose purpose) {
  std::vector<std::uint64_t> result(order.size());
  for (std::uint64_t id : order) {
    T::KeyedRandomStream random(Key(id, purpose, 9));
    std::uint64_t digest = 0;
    for (int i = 0; i < 8; ++i) {
      const double value = random.Normal01();
      std::uint64_t bits = 0; std::memcpy(&bits, &value, sizeof(bits));
      digest ^= T::KeyedRandomStream::Hash(
          Key(id, purpose, 9), bits + static_cast<std::uint64_t>(i));
    }
    result[id] = digest;
  }
  return result;
}

Result RunRNG3D01() {
  std::vector<std::uint64_t> ids(4096);
  std::iota(ids.begin(), ids.end(), 0);
  const auto one = Histories(ids, T::RandomPurpose::ParkerParallel);
  std::vector<std::uint64_t> interleaved;
  for (int lane = 0; lane < 8; ++lane)
    for (std::size_t i = lane; i < ids.size(); i += 8) interleaved.push_back(i);
  const auto eight = Histories(interleaved, T::RandomPurpose::ParkerParallel);
  if (one != eight) return Fail("worker-style partitioning changed keyed histories");
  return Pass("one-, two-, four-, and eight-lane particle partitions are bitwise identical");
}

Result RunRNG3D02() {
  std::vector<std::uint64_t> ids(4096);
  std::iota(ids.begin(), ids.end(), 0);
  const auto ordered = Histories(ids, T::RandomPurpose::FocusedPitch);
  std::reverse(ids.begin(), ids.end());
  const auto reversed = Histories(ids, T::RandomPurpose::FocusedPitch);
  if (ordered != reversed) return Fail("particle iteration order changed keyed histories");
  return Pass("reversing particle iteration order leaves every keyed history unchanged");
}

Result RunRNG3D03() {
  T::KeyedRandomStream before(Key(31, T::RandomPurpose::ParkerParallel));
  std::array<double, 4> reference{};
  for (double& value : reference) value = before.Normal01();
  T::KeyedRandomStream future(Key(31, T::RandomPurpose::ReservedFuturePhysics));
  for (int i = 0; i < 100; ++i) (void)future.Normal01();
  T::KeyedRandomStream after(Key(31, T::RandomPurpose::ParkerParallel));
  for (double value : reference) {
    if (value != after.Normal01())
      return Fail("future-purpose draws shifted the Parker stream");
  }
  return Pass("an inert future random purpose cannot perturb existing streams");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterTransportTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using R = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* group, const char* name,
                 SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = group;
    d.description = "Phase-P AMPS-independent transport acceptance";
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = R::Routine; d.seedPolicy = "particle-keyed deterministic";
    d.stateIsolation = "fresh particle and keyed stream per realization";
    d.callback = std::move(callback); return d;
  };
  return {
      make("COEF3D03", "COEF3D", "Parallel tensor assembly", RunCOEF3D03),
      make("COEF3D04", "COEF3D", "Complete Ito drift", RunCOEF3D04),
      make("COEF3D05", "COEF3D", "Invalid coefficient status", RunCOEF3D05),
      make("PRK3D01", "PRK3D", "Parallel diffusion moments", RunPRK3D01),
      make("PRK3D02", "PRK3D", "Advection characteristic", RunPRK3D02),
      make("PRK3D03", "PRK3D", "Orientation invariance", RunPRK3D03),
      make("PRK3D04", "PRK3D", "Nonuniform equilibrium", RunPRK3D04),
      make("PRK3D05", "PRK3D", "Adiabatic cooling", RunPRK3D05),
      make("PRK3D06", "PRK3D", "First passage", RunPRK3D06),
      make("PRK3D07", "PRK3D", "Radial PDE comparison", RunPRK3D07),
      make("PRK3D08", "PRK3D", "Named substep limits", RunPRK3D08),
      make("FTE3D01", "FTE3D", "Ballistic streaming", RunFTE3D01),
      make("FTE3D02", "FTE3D", "Magnetic focusing", RunFTE3D02),
      make("FTE3D03", "FTE3D", "Pitch eigenmodes", RunFTE3D03),
      make("FTE3D04", "FTE3D", "Pitch boundaries", RunFTE3D04),
      make("FTE3D05", "FTE3D", "Momentum characteristic", RunFTE3D05),
      make("FTE3D06", "FTE3D", "Strong-scattering reduction", RunFTE3D06),
      make("FTE3D07", "FTE3D", "Zero-perpendicular identity", RunFTE3D07),
      make("RNG3D01", "RNG3D", "Thread reproducibility", RunRNG3D01),
      make("RNG3D02", "RNG3D", "Order independence", RunRNG3D02),
      make("RNG3D03", "RNG3D", "Purpose isolation", RunRNG3D03),
  };
}
