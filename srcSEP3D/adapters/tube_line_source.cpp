// ============================================================================
// srcSEP3D/adapters/tube_line_source.cpp
//
// Implementation of the 1-D-equivalent shock source on the observer-tube
// lines.  The physics, the weights and the frame are described in
// tube_line_source.h; this file follows the 1-D implementation
// (srcSEP/reduced_shock_particles.cpp, the reduced-shock branch of the
// source, and srcSEP/util/sep_reduced_source.cpp NetProtonRate) step by step,
// and the comments below name the 1-D counterpart of each step.
// ============================================================================
#include "tube_line_source.h"

#include "../transport/keyed_random.h"
#include "sep_coronal_cme/particle_source.h"
#include "sep_injection_spectrum.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace SEP3D {
namespace Adapters {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

// Relativistic momentum [kg m/s] of kinetic energy E [J]:
// (pc)^2 = E (E + 2 m c^2), the inverse of E = sqrt((pc)^2 + (mc^2)^2) - mc^2.
double MomentumFromEnergy(double energyJ, double massKg) {
  const double rest = massKg * Core::Const::c * Core::Const::c;
  return std::sqrt(energyJ * (energyJ + 2.0 * rest)) / Core::Const::c;
}

// The random stream of one (step, sample, species, purpose).  Each purpose
// (time, line, transverse position, spectrum, pitch, gyrophase, identity) has
// its own stream, so changing how one quantity is drawn never shifts the
// draws of another, and the candidates do not depend on which rank evaluates
// them.  The 1-D source keys its streams the same way (purpose-separated
// keys per event and macroparticle).
Transport::KeyedRandomStream Stream(std::uint64_t seed, std::uint64_t step,
                                    std::uint64_t sample, int species,
                                    Transport::RandomPurpose purpose) {
  Transport::RandomKey key;
  key.campaignSeed = seed;
  // particleId carries the sample and species: one stream per candidate.
  key.particleId = (sample << 8) ^ static_cast<std::uint64_t>(species & 0xff);
  key.step = step;
  key.substep = 0;
  key.purpose = purpose;
  return Transport::KeyedRandomStream(key);
}

// Crossing of a radius-monotone tube line with the front (release = false,
// level L = ShockSurface::Level) or with the release surface (release = true,
// level = ShockSurface::ReleaseLevel with f = fraction), at epochS.  The line
// runs from the inner sphere outward; the inner end is downstream (L < 0)
// once the front has passed it.
//   point  output: the crossing [m, frame], on the chord between two line
//          points (the polyline is the line representation everywhere)
//   chord  output: index i of the chord [p_i, p_{i+1}] holding the crossing;
//          p.size() (an out-of-range index) when the line does not cross,
//          i.e. unless level(first) < 0 < level(last)
// With a level that is monotone along the line (the heliocentric sphere and
// its release sphere on a radius-monotone line) the crossing is unique.
// For a non-monotone level the bisection below returns one sign change, not
// necessarily the first one along the line.
// Fails only on a line of fewer than two points or a provider error.
Core::Status LineCrossing(const Mesh::TracedTubeLine& line,
                          const ShockSurface& shock, double epochS,
                          bool release, double fraction, Core::Vec3* point,
                          std::size_t* chord) {
  const auto& p = line.pointsM;
  if (p.size() < 2) return Invalid("tube line has fewer than two points");
  auto level = [&](const Core::Vec3& x, double* value) {
    return release ? shock.ReleaseLevel(x, epochS, fraction, value)
                   : shock.Level(x, epochS, value);
  };
  // Binary search on the monotone radius: the level of a heliocentric
  // sphere family is monotone along a radius-monotone line.  The search keeps
  // level(p[low]) < 0 <= level(p[high]), so it costs O(log n) provider
  // queries instead of one per chord (the source evaluates every line at
  // every sample epoch).
  double l0 = 0.0, l1 = 0.0;
  Core::Status status = level(p.front(), &l0);
  if (!status.ok()) return status;
  status = level(p.back(), &l1);
  if (!status.ok()) return status;
  if (!(l0 < 0.0 && l1 > 0.0)) {
    *chord = p.size();  // no crossing on this line
    return Core::Status::OK();
  }
  std::size_t low = 0, high = p.size() - 1;
  while (high - low > 1) {
    const std::size_t middle = (low + high) / 2;
    double lm = 0.0;
    status = level(p[middle], &lm);
    if (!status.ok()) return status;
    if (lm < 0.0) low = middle; else high = middle;
  }
  double la = 0.0, lb = 0.0;
  status = level(p[low], &la);
  if (!status.ok()) return status;
  status = level(p[high], &lb);
  if (!status.ok()) return status;
  // la and lb are re-evaluated status checks of the bracketing chord ends;
  // their values are not used below (the bracket is already known).
  // Root of the level on the chord by bisection (the chord is the line).
  // 80 halvings shrink the bracket below 2^-53, the binary64 spacing of the
  // chord parameter near 1, so the root is located to roundoff of the chord.
  double a = 0.0, b = 1.0;
  for (int i = 0; i < 80; ++i) {
    const double m = 0.5 * (a + b);
    double lm = 0.0;
    status = level(p[low] + (p[high] - p[low]) * m, &lm);
    if (!status.ok()) return status;
    if (lm < 0.0) a = m; else b = m;
  }
  *point = p[low] + (p[high] - p[low]) * (0.5 * (a + b));
  *chord = low;
  return Core::Status::OK();
}

// Local line spacing s = sqrt(Phi/|B|) [m] at a point x of chord `chord`,
// linearly interpolated in arc fraction between the chord's two points
// (TracedTubeLine::localSpacingM is stored per point).  s is the side of the
// square cross-section (normal to B) that carries the element's flux Phi:
// s^2 |B| = Phi.  The fraction is clipped at 1 so a point slightly past the
// chord end does not extrapolate.
double LocalSpacing(const Mesh::TracedTubeLine& line, std::size_t chord,
                    const Core::Vec3& x) {
  const auto& p = line.pointsM;
  const double length = (p[chord + 1] - p[chord]).Norm();
  const double f = length > 0 ? (x - p[chord]).Norm() / length : 0.0;
  return line.localSpacingM[chord] +
         std::min(1.0, f) * (line.localSpacingM[chord + 1] - line.localSpacingM[chord]);
}

}  // namespace

Core::Status TubeLineNetRate(const Mesh::TracedTubeLine& line,
                             const ShockSurface& shock, double epochS,
                             double eta, double* ratePerS,
                             Core::Vec3* crossingM, ShockFrontState* front) {
  if (!ratePerS) return Invalid("tube line rate output is null");
  *ratePerS = 0.0;
  Core::Vec3 x;
  std::size_t chord = 0;
  // Crossing with the front itself (release = false).  The 1-D source finds
  // the same root on its moved material line (ReducedShock::ConnectLine).
  Core::Status status = LineCrossing(line, shock, epochS, false, 0.0, &x, &chord);
  if (!status.ok()) return status;
  if (chord >= line.pointsM.size()) return Core::Status::OK();
  // Local front state at the crossing.  x lies on the front to roundoff, so
  // the provider's projection (radial for the heliocentric sphere) returns
  // the state of the front point the line actually meets.  A point outside
  // the provider's support (e.g. beyond a finite cap) returns accepted = false.
  ShockFrontState state;
  status = shock.FrontAt(x, epochS, &state);
  if (!status.ok()) return status;
  if (crossingM) *crossingM = x;
  if (front) *front = state;
  // Only an admitted fast-mode shock with a valid jump is a source (1-D:
  // `if(!front.acceptedShock) return 0`).
  if (!state.accepted) return Core::Status::OK();
  const double bn = std::fabs(state.upstreamMagneticFieldT.Dot(state.outwardNormal));
  if (!(bn > 0.0) || !(state.upstreamProtonDensityM3 > 0.0) ||
      !(state.inflowMPerS > 0.0))
    return Invalid("accepted tube-line front lacks a physical upstream state "
                   "(tangent field or zero inflow); no obliquity floor is used");
  // eta n1 w1 Phi / |B1.n| (the 1-D NetProtonRate).  Units:
  // [m^-3][m/s][T m^2]/[T] = [1/s]; Phi/|B1.n| is the element's intercepted
  // area on the front, so the rate is eta times the upstream proton flux
  // through that area in the shock frame.
  *ratePerS = eta * state.upstreamProtonDensityM3 * state.inflowMPerS *
              line.fluxWb / bn;
  if (!std::isfinite(*ratePerS)) return Invalid("tube line rate is not finite");
  return Core::Status::OK();
}

// One step of the tube-line source (tube_line_source.h).  For each of the N
// stratified samples: evaluate every line's rate at the sample epoch, choose
// one line by rate, locate its release point, spread it over the element's
// footprint, draw the momentum and the pitch angle, and assign the weight
// w (sum R) dt / N.  Any failure aborts the whole batch, so a partial step is
// never injected (the 1-D source likewise rejects before any allocation).
TubeLineSourceResult GenerateTubeLineSource(const TubeLineSourceRequest& q) {
  TubeLineSourceResult result;
  if (!q.tube || !q.shock || !q.plasma || q.samplesPerStep == 0 ||
      !(q.dtS > 0.0) || !(q.massKg > 0.0) || !(q.injectionEfficiency > 0.0) ||
      !(q.releaseFraction > 0.0) || !(q.maximumEnergyJ > q.minimumEnergyJ) ||
      !(q.minimumEnergyJ > 0.0)) {
    result.status = Invalid("incomplete tube-line source request");
    return result;
  }
  const auto& lines = q.tube->lines;
  // Momentum limits of the source band, from its kinetic-energy limits.
  const double pmin = MomentumFromEnergy(q.minimumEnergyJ, q.massKg);
  const double pmax = MomentumFromEnergy(q.maximumEnergyJ, q.massKg);
  // Per-line rate and front state at the current sample epoch (reused
  // across samples; every entry is rewritten for every sample).
  std::vector<double> rates(lines.size());
  std::vector<ShockFrontState> fronts(lines.size());
  for (std::uint64_t j = 0; j < q.samplesPerStep; ++j) {
    // Stratified birth epoch: one draw per stratum of dt/N.
    //   phase = (j + u) / N in [j/N, (j+1)/N),   t_j = t0 + phase dt.
    // Stratification estimates the time integral of the (time-dependent)
    // source rate over the step without assuming it constant, as the 1-D
    // source's "stratified temporal Monte Carlo quadrature".
    auto time = Stream(q.campaignSeed, q.step, j, q.species,
                       Transport::RandomPurpose::TubeSourceTime);
    const double phase = (static_cast<double>(j) + time.UniformOpen01()) /
                         static_cast<double>(q.samplesPerStep);
    const double epoch = q.stepBeginS + phase * q.dtS;
    // R_l(t_j) for every line and their sum.  The lines are fixed curves in
    // the computation frame (no advection between samples); only the front
    // moves with the epoch.
    double total = 0.0;
    for (std::size_t l = 0; l < lines.size(); ++l) {
      const Core::Status status = TubeLineNetRate(
          lines[l], *q.shock, epoch, q.injectionEfficiency, &rates[l], nullptr,
          &fronts[l]);
      if (!status.ok()) { result.status = status; return result; }
      total += rates[l];
    }
    // Diagnostics of the first sample (its epoch lies in the first stratum).
    if (j == 0) {
      result.totalRateBeginPerS = total;
      result.acceptedLines = static_cast<std::uint64_t>(
          std::count_if(rates.begin(), rates.end(), [](double r) { return r > 0; }));
    }
    // No accepted shock on any line at this epoch: the stratum contributes
    // nothing, which is the correct value of its quadrature term (R = 0).
    if (!(total > 0.0)) { ++result.emptySamples; continue; }
    // Line by rate: inverse-CDF selection of l with P(l) = R_l / sum R.  The
    // last line takes the remainder, so the selection always returns a line.
    auto pick = Stream(q.campaignSeed, q.step, j, q.species,
                       Transport::RandomPurpose::TubeSourceLine);
    const double target = pick.UniformOpen01() * total;
    std::size_t line = 0;
    double cumulative = 0.0;
    for (; line + 1 < lines.size(); ++line) {
      cumulative += rates[line];
      if (target < cumulative) break;
    }
    // If roundoff in the running sum lets the loop fall through to a last
    // line of zero rate, step back to the nearest line that carries a source.
    // (A line where the loop breaks always has R_l > 0.)
    while (rates[line] == 0.0 && line > 0) --line;  // roundoff at the end
    // Front state at the selected line's shock crossing: it supplies the
    // compression ratio of the DSA law and the front speed of the pitch law
    // (the 1-D source also takes both from its shock connection; its normal
    // speed comes from the release point's normal foot, which for the
    // heliocentric sphere is the same uniform expansion speed).
    const ShockFrontState& front = fronts[line];
    // Release point: the line's crossing of the release surface
    // (1-D: ReferenceIntersections, the release sphere R + L_ref on the same
    // material line, upstream of the shock root).  On a radius-monotone line
    // the release crossing lies outward of the shock crossing.  A line whose
    // outer end lies inside the release sphere has no release point; that
    // fails the batch, as the 1-D source fails without an upstream L_ref
    // intersection (the tube must extend beyond (1+f) R_sh while the source
    // is active).
    Core::Vec3 release;
    std::size_t chord = 0;
    Core::Status status = LineCrossing(lines[line], *q.shock, epoch, true,
                                       q.releaseFraction, &release, &chord);
    if (!status.ok()) { result.status = status; return result; }
    if (chord >= lines[line].pointsM.size()) {
      result.status = Invalid("configured L_ref has no upstream intersection "
                              "on the same tube line");
      return result;
    }
    // Field and plasma at the release point (frame): b [T] and u = U' [m/s].
    // They are evaluated on the traced line, before the transverse spread
    // below; the pitch law uses these values for the whole footprint.
    Core::Vec3 b, u;
    status = q.plasma(release, epoch, &b, &u);
    if (!status.ok()) { result.status = status; return result; }
    const double bMagnitude = b.Norm();
    const Core::Vec3 bHat = b / bMagnitude;
    // Outward normal of the release surface.  The heliocentric release
    // sphere has the radial normal x/|x| (the coordinate origin is the Sun).
    const Core::Vec3 normal = release.Normalized();
    // Transverse position over the element footprint on the release sphere.
    //
    // The element's cross-section normal to B is a square of side s, with
    // s^2 |B| = Phi.  Project it along B onto the plane tangent to the
    // sphere: the tangent direction e2, which is perpendicular to both n and
    // B, keeps length s; the tangent direction e1 along the projection of b
    // onto the tangent plane is stretched by 1/cos(theta), cos(theta) =
    // |b.n|.  The footprint on the sphere is therefore an s/cos(theta) by s
    // rectangle of area s^2/cos(theta) = Phi/|B.n|, the same intercepted area
    // that sets R_l.  A uniform point of the rectangle is drawn and pushed
    // radially back onto the sphere of the release radius; for s << r this
    // is uniform over the footprint on the sphere.
    {
      const double spacing = LocalSpacing(lines[line], chord, release);
      // The 1e-3 floor only keeps the 1/cos(theta) stretch finite (at most
      // 1000 s) where B is nearly tangent to the release sphere; it is a
      // numerical guard, not a physical parameter.
      const double cosine = std::max(1e-3, std::fabs(bHat.Dot(normal)));
      // e1: the unit tangent along the projection of b; when b is radial
      // (no tangent projection) any tangent direction serves, chosen from
      // n x z, or n x x at the poles.
      Core::Vec3 e1 = bHat - normal * bHat.Dot(normal);
      if (e1.Norm() < 1e-12) e1 = normal.Cross(Core::Vec3(0, 0, 1));
      if (e1.Norm() < 1e-12) e1 = normal.Cross(Core::Vec3(1, 0, 0));
      e1 = e1.Normalized();
      const Core::Vec3 e2 = normal.Cross(e1);
      // Two independent uniform offsets, centred on the traced line:
      // a1 in [-s/(2 cos), s/(2 cos)) along e1, a2 in [-s/2, s/2) along e2.
      auto jitter = Stream(q.campaignSeed, q.step, j, q.species,
                           Transport::RandomPurpose::TubeSourceTransverse);
      const double a1 = (jitter.UniformOpen01() - 0.5) * spacing / cosine;
      const double a2 = (jitter.UniformOpen01() - 0.5) * spacing;
      const double radius = release.Norm();
      release = (release + e1 * a1 + e2 * a2).Normalized() * radius;
    }
    // Momentum from the local-compression DSA law of the selected line's
    // crossing: f(p) ~ p^-q, q = 3r/(r-1), so dN/dp ~ p^(2-q) on [pmin,pmax]
    // (MomentumSpectrum::LocalCompressionDsa, the same core law as 1-D).
    const auto dsa = SEP::CoronalCME::MomentumSpectrum::LocalCompressionDsa(
        pmin, pmax, front.compressionRatio);
    if (!dsa.ok()) { result.status = Invalid(dsa.status.message); return result; }
    auto spectral = Stream(q.campaignSeed, q.step, j, q.species,
                           Transport::RandomPurpose::TubeSourceSpectrum);
    double p = 0.0, weight = 1.0;
    if (q.logUniformMomentum) {
      // The 1-D proposal g(p) = 1/(p L) with the exact weight of the same
      // DSA law, evaluated in x = p/pmin (sep_injection_spectrum).
      //   x = exp(u ln(pmax/pmin)) is log-uniform on [1, pmax/pmin];
      //   w(x) = f(x) / g(x) = f(x) x ln(pmax/pmin), f the normalized DSA
      //   density in x (ImportanceWeightForLogUniformProposal).
      // Working in x = p/pmin keeps the power law representable: in SI
      // units p^(2-q) overflows binary64 for weak shocks with large q
      // (the 1-D source gives the same reason).  The clamp catches the one
      // ulp exp() may add at u -> 1.  The proposal puts equal numbers of
      // macroparticles in equal logarithmic momentum intervals, so the steep
      // high-energy tail is sampled far more often than by direct sampling;
      // w restores the DSA law, and E[w] = 1 keeps the expected physical
      // number unchanged.
      SEP::Injection::Spectrum scaled;
      scaled.measure = SEP::Injection::Measure::Momentum;
      scaled.minimum = 1.0;
      scaled.maximum = pmax / pmin;
      scaled.powerIndex = -dsa.value.NumberExponent();
      double x = std::exp(spectral.UniformOpen01() * std::log(scaled.maximum));
      if (x > scaled.maximum) x = scaled.maximum;
      const auto w = SEP::Injection::ImportanceWeightForLogUniformProposal(scaled, x);
      if (!w.status.ok() || !std::isfinite(w.value) || !(w.value > 0.0)) {
        result.status = Invalid("log-uniform DSA importance weight is not finite positive");
        return result;
      }
      p = pmin * x;
      weight = w.value;
    } else {
      // Direct sampling (constant-statistical-weight): inverse CDF of the
      // DSA law itself, weight 1.
      p = dsa.value.InverseCdf(spectral.UniformOpen01());
    }
    // Relativistic speed v = p c^2 / E_total = c p / sqrt(p^2 + m^2 c^2).
    const double speed = Core::Const::c * p /
        std::sqrt(p * p + q.massKg * q.massKg * Core::Const::c * Core::Const::c);
    // Advection of the release surface frame: U.n minus the front speed at
    // the foot (the 1-D -normalFoot.normalSpeed + U.n).
    // A particle moving with plasma-frame speed v and cosine mu along B has
    // the normal velocity
    //   v_n(mu) = advection + v mu (b.n),   advection = U'.n - V_n,
    // relative to a surface that moves with the front's normal speed V_n.
    // For the radial normal U'.n = U.n (Omega x x is tangential), so the value
    // is the same in the corotating and the inertial frame.
    const double advection = u.Dot(normal) - front.normalSpeedMPerS;
    const double bDotN = bHat.Dot(normal);
    auto pitch = Stream(q.campaignSeed, q.step, j, q.species,
                        Transport::RandomPurpose::TubeSourcePitch);
    const double pitchDraw = pitch.UniformOpen01();
    // Parker: isotropic mu, the angular marginal the Parker equation evolves
    // (it has no pitch-angle state; 1-D stores the same projection).
    double mu = 2.0 * pitchDraw - 1.0;
    if (q.focusedPitchLaw) {
      // Focused movers: the positive-normal-flux law of the 1-D source,
      // density ~ max(0, v_n(mu)) on [-1, 1]: an isotropic population counted
      // by the flux it carries outward through the release surface
      // (BuildFocusedEscapeDistribution, 64 trapezoid intervals plus the
      // kink at v_n = 0, so the integration is exact).
      // mu is measured along B (polarity +1 with the signed b.n).  The whole
      // band must escape: test its minimum speed first, as in 1-D.
      // The largest normal velocity, advection + v |b.n|, grows with v, so
      // escape at pmin implies escape at every momentum of the band.  The
      // source is a prescribed NET escaping spectrum; dropping or redrawing
      // slow particles would change its number and its DSA marginal, so an
      // incompatible band fails instead (the 1-D source does the same).
      const double minimumSpeed = Core::Const::c * pmin /
          std::sqrt(pmin * pmin + q.massKg * q.massKg * Core::Const::c * Core::Const::c);
      const auto preflight = SEP::CoronalCME::BuildFocusedEscapeDistribution(
          minimumSpeed, 1, bDotN, advection, 64, 0);
      if (!preflight.ok()) { result.status = Invalid(preflight.status.message); return result; }
      if (preflight.value.status == SEP::CoronalCME::FocusedEscapeStatus::NoFocusedEscape) {
        result.status = Invalid("the net DSA source band includes momenta with no "
                                "outward release-surface escape");
        return result;
      }
      // The law at the realized speed.  SampleFocusedMu accepts only an
      // Admissible law, so a realized momentum without escape (excluded by
      // the preflight above) would still fail here rather than be injected.
      const auto law = SEP::CoronalCME::BuildFocusedEscapeDistribution(
          speed, 1, bDotN, advection, 64, 0);
      if (!law.ok()) { result.status = Invalid(law.status.message); return result; }
      const auto sampled = SEP::CoronalCME::SampleFocusedMu(law.value, pitchDraw);
      if (!sampled.ok()) { result.status = Invalid(sampled.status.message); return result; }
      mu = sampled.value;
    }
    auto gyro = Stream(q.campaignSeed, q.step, j, q.species,
                       Transport::RandomPurpose::TubeSourceGyrophase);
    InjectedParticle injected;
    injected.status = Core::Status::OK();
    // A particle born at phase tau of the step moves only the remaining
    // (1 - tau) dt in its first mover call (AMPS MoveParticle), so the birth
    // times keep their stratified positions inside the step.
    injected.remainingFirstStepFraction = 1.0 - phase;
    // Stable identity from its own keyed stream: a deterministic function of
    // (seed, step, sample, species), identical on every rank.  0 is the
    // "no particle" identity the movers reject, so it is mapped to 1.
    injected.particle.stableId = Transport::KeyedRandomStream::Hash(
        Stream(q.campaignSeed, q.step, j, q.species,
               Transport::RandomPurpose::TubeSourceIdentity).state().key, 0);
    if (injected.particle.stableId == 0) injected.particle.stableId = 1;
    injected.particle.species = q.species;
    injected.particle.positionM = release;
    injected.particle.momentumKgMPerS = p;
    injected.particle.mu = mu;
    injected.particle.gyrophaseRad = 2.0 * Core::Const::kPi * gyro.UniformOpen01();
    // W = w (sum_l R_l(t_j)) dt / N physical protons (see the header): the
    // line choice by R_l / sum R and this weight together give each line the
    // 1-D weight R_l dt / N per stratum (1-D: weight * rate * dt / N).
    injected.particle.statisticalWeight =
        weight * total * q.dtS / static_cast<double>(q.samplesPerStep);
    // Keys of the particle's own mover streams start at this step.
    injected.particle.completedStep = q.step;
    if (!std::isfinite(injected.particle.statisticalWeight) ||
        !(injected.particle.statisticalWeight > 0.0)) {
      result.status = Invalid("tube-line source weight is not finite positive");
      return result;
    }
    result.representedParticles += injected.particle.statisticalWeight;
    result.particles.push_back(std::move(injected));
  }
  result.status = Core::Status::OK();
  return result;
}

}  // namespace Adapters
}  // namespace SEP3D
