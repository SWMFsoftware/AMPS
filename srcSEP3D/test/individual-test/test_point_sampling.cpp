// PS3D01-PS3D04: AMPS-independent tests of the point-sampling kernels
// (output/point_sampling.h).  Every expected value below is computed from an
// independent closed form or by hand from the fixture, never by calling the
// routine under test a second time.
//
//   PS3D01  energy channel grid: exact linear/logarithmic edges, closed-range
//           channel lookup, centres, and rejected grids.
//   PS3D02  relativistic kinematics: E(p) and v(p) against E = sqrt(p^2c^2 +
//           m^2c^4) - mc^2, the relativistic velocity addition in the
//           collinear and perpendicular cases, and the two-frame evaluation.
//   PS3D03  accumulation, normalization and the Tecplot row: hand-summed
//           density, flux, field-aligned flux, intensity and fluence of a
//           two-call fixture in cm units, out-of-range tails, and the
//           cumulative fluence over two intervals.
//   PS3D04  point-file and product-list grammar, including file:line errors.

#include "sep3d_test_registry.h"

#include "point_sampling.h"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace {

using SEP3D::Testing::Result;
namespace C = SEP3D::Core;
namespace PS = SEP3D::Output::PointSampling;

Result Pass(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Pass;
  result.message = message; return result;
}
Result Fail(const std::string& message) {
  Result result; result.status = SEP3D::Testing::Status::Fail;
  result.message = message; return result;
}

bool Close(double actual, double expected, double relative = 1.0e-12) {
  return std::fabs(actual - expected) <=
         relative * std::max(std::fabs(expected), 1.0e-300);
}

std::string Num(double value) {
  std::ostringstream out; out.precision(17); out << value; return out.str();
}

// ---- PS3D01 ----------------------------------------------------------------
Result RunPS3D01() {
  std::vector<double> edges;
  if (!PS::BuildEnergyEdges(0.0, 10.0, 5, PS::EnergyScale::Linear, &edges).ok())
    return Fail("linear 0-10 MeV / 5 grid rejected");
  const double linear[] = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0};
  for (int i = 0; i < 6; ++i)
    if (!Close(edges[i], linear[i]) && !(linear[i] == 0.0 && edges[i] == 0.0))
      return Fail("linear edge " + std::to_string(i) + " = " + Num(edges[i]));

  if (!PS::BuildEnergyEdges(1.0, 1000.0, 3, PS::EnergyScale::Logarithmic,
                            &edges).ok())
    return Fail("log 1-1000 MeV / 3 grid rejected");
  const double logarithmic[] = {1.0, 10.0, 100.0, 1000.0};
  for (int i = 0; i < 4; ++i)
    if (!Close(edges[i], logarithmic[i], 1.0e-13))
      return Fail("log edge " + std::to_string(i) + " = " + Num(edges[i]));
  if (edges.front() != 1.0 || edges.back() != 1000.0)
    return Fail("log end points are not exact");

  // Closed range: the upper channel owns an interior edge, E_max is in the
  // last channel, values outside are typed tails, NaN is never in range.
  struct Probe { double e; int channel; };
  const Probe probes[] = {{1.0, 0},     {9.99, 0},  {edges[1], 1},
                          {999.0, 2},   {1000.0, 2}, {1000.5, PS::kAboveRange},
                          {0.999, PS::kBelowRange},
                          {std::nan(""), PS::kBelowRange}};
  for (const Probe& probe : probes)
    if (PS::ChannelIndex(edges, probe.e) != probe.channel)
      return Fail("channel of " + Num(probe.e) + " is " +
                  std::to_string(PS::ChannelIndex(edges, probe.e)) +
                  ", expected " + std::to_string(probe.channel));

  PS::SpeciesChannels species;
  species.scale = PS::EnergyScale::Logarithmic;
  species.edgesMeV = edges;
  if (!Close(PS::ChannelCenterMeV(species, 0), std::sqrt(10.0), 1.0e-13))
    return Fail("logarithmic centre is not the geometric mean");
  species.scale = PS::EnergyScale::Linear;
  if (!Close(PS::ChannelCenterMeV(species, 0), 5.5, 1.0e-13))
    return Fail("linear centre is not the arithmetic mean");

  if (PS::BuildEnergyEdges(5.0, 5.0, 2, PS::EnergyScale::Linear, &edges).ok() ||
      PS::BuildEnergyEdges(1.0, 2.0, 0, PS::EnergyScale::Linear, &edges).ok() ||
      PS::BuildEnergyEdges(0.0, 2.0, 2, PS::EnergyScale::Logarithmic, &edges).ok() ||
      PS::BuildEnergyEdges(-1.0, 2.0, 2, PS::EnergyScale::Linear, &edges).ok())
    return Fail("an empty, zero-interval, log-from-zero or negative grid was accepted");
  return Pass("exact linear/log edges, closed-range lookup, centres and grid "
              "rejections");
}

// ---- PS3D02 ----------------------------------------------------------------
Result RunPS3D02() {
  const double c = C::Const::c;
  const double m = C::Const::m_p;
  const double restMeV = m * c * c / PS::kJoulePerMeV;
  // Independent inverse: for kinetic energy T, p c = sqrt(T^2 + 2 T m c^2).
  for (double kineticMeV : {1.0e-3, 1.0, 100.0, 1.0e4}) {
    const double pcMeV = std::sqrt(kineticMeV * kineticMeV +
                                   2.0 * kineticMeV * restMeV);
    const double p = pcMeV * PS::kJoulePerMeV / c;
    const double energy = PS::KineticEnergyFromMomentumJ(p, m) / PS::kJoulePerMeV;
    if (!Close(energy, kineticMeV, 1.0e-12))
      return Fail("E(p) at " + Num(kineticMeV) + " MeV gives " + Num(energy));
    const double speed = PS::SpeedFromMomentumMPerS(p, m);
    const double expectedSpeed = pcMeV / (kineticMeV + restMeV) * c;
    if (!Close(speed, expectedSpeed, 1.0e-12))
      return Fail("v(p) at " + Num(kineticMeV) + " MeV");
    const double fromSpeed = PS::KineticEnergyFromSpeedJ(speed, m) / PS::kJoulePerMeV;
    if (!Close(fromSpeed, kineticMeV, 1.0e-9))
      return Fail("E(v) at " + Num(kineticMeV) + " MeV gives " + Num(fromSpeed));
  }

  // Velocity addition.  Rest particle -> u = U.  Collinear 0.5c (+) 0.5c =
  // 0.8c.  Perpendicular: u_par = U, u_perp = v'/gamma_U.
  const C::Vec3 U(0.5 * c, 0.0, 0.0);
  if (!Close(PS::BoostVelocity(C::Vec3(), U).x, 0.5 * c))
    return Fail("a particle at rest in the moving frame does not move with U");
  if (!Close(PS::BoostVelocity(C::Vec3(0.5 * c, 0.0, 0.0), U).x, 0.8 * c))
    return Fail("collinear 0.5c + 0.5c is not 0.8c");
  const C::Vec3 perpendicular = PS::BoostVelocity(C::Vec3(0.0, 0.3 * c, 0.0), U);
  if (!Close(perpendicular.x, 0.5 * c) ||
      !Close(perpendicular.y, 0.3 * c * std::sqrt(1.0 - 0.25)))
    return Fail("perpendicular velocity addition is wrong");

  // Two frames.  With U = 0 both frames coincide; v_par = mu v.
  const double p = std::sqrt(1.0 + 2.0 * restMeV) * PS::kJoulePerMeV / c;  // 1 MeV
  PS::ParticleInput particle;
  particle.weight = 1.0; particle.massKg = m; particle.momentumKgMPerS = p;
  particle.mu = 0.6;
  particle.plasmaVelocityDirection = C::Vec3(0.6, 0.8, 0.0);  // mu along +x
  PS::ParticleKinematics k;
  if (!PS::EvaluateKinematics(particle, C::Vec3(), C::Vec3(1, 0, 0), &k).ok())
    return Fail("two-frame evaluation rejected a valid particle");
  if (!Close(k.energyMeV[PS::kPlasmaFrame], 1.0, 1.0e-12) ||
      !Close(k.energyMeV[PS::kInertialFrame], 1.0, 1.0e-9) ||
      !Close(k.parallelSpeedMPerS[PS::kPlasmaFrame],
             0.6 * k.speedMPerS[PS::kPlasmaFrame]) ||
      !Close(k.parallelSpeedMPerS[PS::kInertialFrame],
             0.6 * k.speedMPerS[PS::kPlasmaFrame], 1.0e-12))
    return Fail("with U = 0 the plasma and inertial frames differ");
  // A 400 km/s wind along +x adds to the field-aligned (x) velocity:
  // inertial v_par = (v_x + U)/(1 + v_x U/c^2) by collinear addition of x.
  const double wind = 4.0e5;
  if (!PS::EvaluateKinematics(particle, C::Vec3(wind, 0, 0), C::Vec3(1, 0, 0),
                              &k).ok())
    return Fail("two-frame evaluation rejected a wind");
  const double vx = 0.6 * k.speedMPerS[PS::kPlasmaFrame];
  const double gammaU = 1.0 / std::sqrt(1.0 - wind * wind / (c * c));
  const double expectedX = (vx + wind) / (1.0 + vx * wind / (c * c));
  const double vy = 0.8 * k.speedMPerS[PS::kPlasmaFrame];
  const double expectedY = vy / gammaU / (1.0 + vx * wind / (c * c));
  if (!Close(k.parallelSpeedMPerS[PS::kInertialFrame], expectedX, 1.0e-12) ||
      !Close(k.speedMPerS[PS::kInertialFrame],
             std::sqrt(expectedX * expectedX + expectedY * expectedY), 1.0e-12))
    return Fail("inertial-frame velocity with a wind is wrong");
  if (!(k.energyMeV[PS::kInertialFrame] > k.energyMeV[PS::kPlasmaFrame]))
    return Fail("a particle streaming with the wind did not gain inertial energy");

  particle.mu = 1.5;
  if (PS::EvaluateKinematics(particle, C::Vec3(), C::Vec3(1, 0, 0), &k).ok())
    return Fail("|mu| > 1 was accepted");
  particle.mu = 0.0;
  if (PS::EvaluateKinematics(particle, C::Vec3(), C::Vec3(2, 0, 0), &k).ok())
    return Fail("a non-unit magnetic direction was accepted");
  return Pass("E(p), v(p), E(v), velocity addition and two-frame kinematics "
              "match independent closed forms");
}

// ---- PS3D03 ----------------------------------------------------------------
//
// Fixture: one point, one species, linear channels [0,5) and [5,10] MeV,
// cell volume V = 4 m^3, time step dt = 10 s, U = 0 (both frames equal).
// Call 1: A (w=2, 1 MeV, mu=0.5), B (w=3, 7 MeV, mu=-1), C (w=1, 20 MeV).
// Call 2: A only.  N = 2 calls.
//   channel 0: n = (2/4 + 2/4)/2 = 0.5 m^-3          -> 5e-7 cm^-3
//              F = (2 vA/4 + 2 vA/4)/2 = 0.5 vA      -> x 1e-4 cm^-2 s^-1
//              F| = 0.5 * 0.5 vA                       (mu = 0.5)
//              j = F / (4 pi 5 MeV)
//              fluence = 2 (2 vA/4) 10 = 10 vA m^-2  -> x 1e-4 cm^-2
//              samples = 2
//   channel 1: n = (3/4)/2, F = (3 vB/4)/2, F| = -F, fluence = 3 vB/4 * 10,
//              samples = 1
//   above range: density (1/4)/2, samples 1; below range: none.
Result RunPS3D03() {
  const double c = C::Const::c;
  const double m = C::Const::m_p;
  const double restMeV = m * c * c / PS::kJoulePerMeV;
  auto momentum = [&](double kineticMeV) {
    return std::sqrt(kineticMeV * kineticMeV + 2.0 * kineticMeV * restMeV) *
           PS::kJoulePerMeV / c;
  };
  auto speed = [&](double kineticMeV) {
    return std::sqrt(kineticMeV * kineticMeV + 2.0 * kineticMeV * restMeV) /
           (kineticMeV + restMeV) * c;
  };
  PS::SpeciesChannels species;
  species.ampsIndex = 0; species.symbol = "H_PLUS"; species.massKg = m;
  species.scale = PS::EnergyScale::Linear;
  if (!PS::BuildEnergyEdges(0.0, 10.0, 2, species.scale, &species.edgesMeV).ok())
    return Fail("fixture grid rejected");
  const std::vector<PS::SpeciesChannels> table = {species};
  const PS::Layout layout(1, table);
  std::vector<double> sums(layout.size(), 0.0);
  const double volume = 4.0, dt = 10.0;

  struct Particle { double w, e, mu; };
  auto add = [&](const Particle& item) {
    PS::ParticleInput input;
    input.weight = item.w; input.massKg = m;
    input.momentumKgMPerS = momentum(item.e); input.mu = item.mu;
    input.plasmaVelocityDirection =
        C::Vec3(item.mu, std::sqrt(1.0 - item.mu * item.mu), 0.0);
    PS::ParticleKinematics k;
    const C::Status status =
        PS::EvaluateKinematics(input, C::Vec3(), C::Vec3(1, 0, 0), &k);
    if (status.ok()) PS::Accumulate(layout, species, 0, 0, item.w, volume, dt, k, &sums);
    return status.ok();
  };
  const Particle a{2.0, 1.0, 0.5}, b{3.0, 7.0, -1.0}, above{1.0, 20.0, 0.0};
  if (!add(a) || !add(b) || !add(above) || !add(a))
    return Fail("fixture particle rejected");

  const double vA = speed(1.0), vB = speed(7.0);
  const double n = 2.0;
  const double fourPi = 4.0 * C::Const::kPi;
  for (int frame = 0; frame < PS::kFrames; ++frame) {
    auto at = [&](std::size_t channel, std::size_t q) {
      return sums[layout.Channel(0, 0, frame, channel, q)];
    };
    const double tol = frame == PS::kPlasmaFrame ? 1.0e-12 : 1.0e-9;
    if (!Close(at(0, PS::kDensity) / n, 0.5) ||
        !Close(at(0, PS::kFlux) / n, 0.5 * vA, tol) ||
        !Close(at(0, PS::kFieldAlignedFlux) / n, 0.25 * vA, tol) ||
        !Close(at(0, PS::kFluence), 10.0 * vA, tol) ||
        at(0, PS::kSamples) != 2.0 ||
        !Close(at(1, PS::kDensity) / n, 0.375) ||
        !Close(at(1, PS::kFlux) / n, 0.375 * vB, tol) ||
        !Close(at(1, PS::kFieldAlignedFlux) / n, -0.375 * vB, tol) ||
        !Close(at(1, PS::kFluence), 7.5 * vB, tol) ||
        at(1, PS::kSamples) != 1.0)
      return Fail(std::string("hand-summed channel sums differ in the ") +
                  (frame == 0 ? "plasma" : "inertial") + " frame");
    if (!Close(sums[layout.Tail(0, 0, frame, PS::kAboveDensity)] / n, 0.125) ||
        sums[layout.Tail(0, 0, frame, PS::kAboveSamples)] != 1.0 ||
        sums[layout.Tail(0, 0, frame, PS::kBelowSamples)] != 0.0)
      return Fail("out-of-range tails are wrong");
  }

  // Tecplot product: the first data row of the zone (channel 0) in cm units,
  // after folding this interval into an empty cumulative fluence twice (two
  // identical intervals -> cumulative = 2 x interval fluence).
  std::vector<double> cumulative(layout.size(), 0.0);
  PS::FoldFluence(layout, table, sums, &cumulative);
  PS::FoldFluence(layout, table, sums, &cumulative);
  PS::ProductSelection products;
  if (!PS::ParseProducts("density, flux, spectrum, fluence", &products).ok())
    return Fail("product list rejected");
  PS::IntervalInfo info;
  info.outputFileNumber = 7; info.sampleCalls = 2;
  info.firstSampleTimeS = 100.0; info.lastSampleTimeS = 110.0;
  PS::PointDefinition point;
  point.name = "earth"; point.positionM = C::Vec3(C::Const::AU, 0, 0);
  const std::string text = PS::FormatTecplot(layout, {point}, table, products,
                                             info, sums, cumulative);
  if (text.find("ZONE T=\"earth H_PLUS\", I=2") == std::string::npos ||
      text.find("DATASETAUXDATA amps_output_file_number=\"7\"") == std::string::npos ||
      text.find("DATASETAUXDATA sample_calls=\"2\"") == std::string::npos)
    return Fail("Tecplot zone header or dataset metadata missing");
  std::istringstream lines(text);
  std::string line;
  std::vector<double> row;
  bool inZone = false;
  while (std::getline(lines, line)) {
    if (line.rfind("ZONE", 0) == 0) { inZone = true; continue; }
    if (!inZone || line.rfind("AUXDATA", 0) == 0) continue;
    std::istringstream fields(line);
    double value;
    while (fields >> value) row.push_back(value);
    break;
  }
  // E_min E_max E_center | n F F| j fl_int fl fl_I samples | (inertial ...)
  if (row.size() != 3 + 2 * 8) return Fail("unexpected Tecplot row width " +
                                           std::to_string(row.size()));
  const double flux = 0.5 * vA * 1.0e-4;
  const double expected[] = {0.0, 5.0, 2.5, 0.5e-6, flux, 0.5 * flux,
                             flux / (fourPi * 5.0), 10.0 * vA * 1.0e-4,
                             20.0 * vA * 1.0e-4,
                             20.0 * vA * 1.0e-4 / (fourPi * 5.0), 2.0};
  for (int i = 0; i < 11; ++i) {
    const bool match = expected[i] == 0.0 ? row[i] == 0.0
                                          : Close(row[i], expected[i], 1.0e-9);
    if (!match)
      return Fail("Tecplot column " + std::to_string(i) + " = " + Num(row[i]) +
                  ", expected " + Num(expected[i]));
  }
  return Pass("hand-summed density, flux, field-aligned flux, intensity, "
              "interval/cumulative fluence and tails in cm units");
}

// ---- PS3D04 ----------------------------------------------------------------
Result RunPS3D04() {
  std::vector<PS::PointDefinition> points;
  const std::string file =
      "! sampling points\n"
      "1.5e11 0 0 earth   # Earth\n"
      "\n"
      "7.5e10 1e9 -2e9\n"
      "  3e10 0 0 inner\n";
  if (!PS::ParsePointFile(file, "points.txt", &points).ok() || points.size() != 3)
    return Fail("valid point file rejected");
  if (points[0].name != "earth" || points[1].name != "1" ||
      points[2].name != "inner" || points[1].positionM.z != -2.0e9 ||
      points[1].origin != "points.txt:4")
    return Fail("names, unnamed index or origin are wrong");
  const struct { const char* text; const char* where; } bad[] = {
      {"1 2\n", "points.txt:1"},
      {"0 0 0\n1 2 x\n", "points.txt:2"},
      {"1 2 3 a b\n", "points.txt:1"},
      {"! nothing\n", "points.txt"}};
  for (const auto& item : bad) {
    const C::Status status = PS::ParsePointFile(item.text, "points.txt", &points);
    if (status.ok() || status.message.find(item.where) != 0)
      return Fail(std::string("bad point file not rejected at ") + item.where +
                  ": " + status.message);
  }

  PS::ProductSelection products;
  if (!PS::ParseProducts("fluence,density", &products).ok() ||
      !products.fluence || !products.density || products.flux ||
      products.spectrum)
    return Fail("product subset parsed incorrectly");
  if (PS::ParseProducts("flux, flux", &products).ok() ||
      PS::ParseProducts("flux, intensity", &products).ok() ||
      PS::ParseProducts("", &products).ok())
    return Fail("a repeated, unknown or empty product list was accepted");
  return Pass("point-file rows, default index names, file:line errors and "
              "product lists");
}

}  // namespace

std::vector<SEP3D::Testing::Descriptor> RegisterPointSamplingTests() {
  using D = SEP3D::Testing::Descriptor;
  using I = SEP3D::Testing::InitializationLevel;
  using RC = SEP3D::Testing::RuntimeClass;
  auto make = [](const char* id, const char* name, const char* description,
                 SEP3D::Testing::TestCallback callback) {
    D d; d.id = id; d.name = name; d.group = "PS3D";
    d.description = description;
    d.initialization = I::None; d.supportedBuildModes = "standalone-no-AMPS";
    d.runtime = RC::Routine; d.seedPolicy = "deterministic";
    d.stateIsolation = "pure functions on local fixtures";
    d.callback = std::move(callback); return d;
  };
  return {
      make("PS3D01", "Point-sampling energy channels",
           "Exact linear/log edges, closed-range lookup, rejected grids.", RunPS3D01),
      make("PS3D02", "Point-sampling kinematics",
           "Relativistic E(p), v(p), velocity addition, plasma/inertial frames.", RunPS3D02),
      make("PS3D03", "Point-sampling accumulation and product",
           "Hand-summed density, flux, intensity, fluence and Tecplot row.", RunPS3D03),
      make("PS3D04", "Point-sampling input grammar",
           "Point-file rows, default names, errors, product lists.", RunPS3D04),
  };
}
