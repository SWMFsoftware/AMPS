// srcSEP3D point sampling: numerical kernels and the Tecplot product.
// Physics, units and the data layout are documented in point_sampling.h.

#include "point_sampling.h"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>

namespace SEP3D {
namespace Output {
namespace PointSampling {
namespace {

Core::Status Invalid(const std::string& message) {
  return Core::Status(Core::StatusCode::InvalidInput, message);
}

bool FiniteVector(const Core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

// Strict finite-number parse of a whole token (no trailing characters).
bool ParseNumber(const std::string& token, double* value) {
  if (token.empty()) return false;
  errno = 0;
  char* end = nullptr;
  const double parsed = std::strtod(token.c_str(), &end);
  if (errno == ERANGE || end == token.c_str() || *end != '\0' ||
      !std::isfinite(parsed)) return false;
  *value = parsed;
  return true;
}

// SI -> product units.  1 m^-3 = 1e-6 cm^-3; 1 m^-2 = 1e-4 cm^-2.
constexpr double kPerM3ToPerCm3 = 1.0e-6;
constexpr double kPerM2ToPerCm2 = 1.0e-4;
constexpr double kFourPi = 4.0 * Core::Const::kPi;
using C3 = Core::Vec3;

}  // namespace

Core::Status BuildEnergyEdges(double minimumMeV, double maximumMeV,
                              unsigned channels, EnergyScale scale,
                              std::vector<double>* edgesMeV) {
  if (edgesMeV == nullptr) return Invalid("energy-edge output is null");
  if (!std::isfinite(minimumMeV) || !std::isfinite(maximumMeV) ||
      channels == 0 || maximumMeV <= minimumMeV || minimumMeV < 0.0)
    return Invalid("energy channels need finite 0 <= E_min < E_max and at "
                   "least one interval");
  if (scale == EnergyScale::Logarithmic && minimumMeV <= 0.0)
    return Invalid("a logarithmic energy scale needs E_min > 0");
  std::vector<double> edges(static_cast<std::size_t>(channels) + 1);
  const double n = static_cast<double>(channels);
  for (unsigned i = 0; i <= channels; ++i) {
    const double fraction = static_cast<double>(i) / n;
    edges[i] = scale == EnergyScale::Linear
        ? minimumMeV + fraction * (maximumMeV - minimumMeV)
        : minimumMeV * std::pow(maximumMeV / minimumMeV, fraction);
  }
  // Exact end points: the configured range is closed, [E_min, E_max].
  edges.front() = minimumMeV;
  edges.back() = maximumMeV;
  for (std::size_t i = 1; i < edges.size(); ++i)
    if (!(edges[i] > edges[i - 1]))
      return Invalid("energy channels are not strictly increasing (range too "
                     "narrow for the requested number of intervals)");
  *edgesMeV = std::move(edges);
  return Core::Status::OK();
}

int ChannelIndex(const std::vector<double>& edgesMeV, double energyMeV) {
  if (edgesMeV.size() < 2 || !(energyMeV >= edgesMeV.front()))
    return kBelowRange;  // also catches NaN
  if (energyMeV > edgesMeV.back()) return kAboveRange;
  if (energyMeV == edgesMeV.back())
    return static_cast<int>(edgesMeV.size()) - 2;
  // First edge strictly greater than E; the channel is the one before it.
  const auto upper =
      std::upper_bound(edgesMeV.begin(), edgesMeV.end(), energyMeV);
  return static_cast<int>(upper - edgesMeV.begin()) - 1;
}

double ChannelCenterMeV(const SpeciesChannels& species, std::size_t channel) {
  const double low = species.edgesMeV[channel];
  const double high = species.edgesMeV[channel + 1];
  return species.scale == EnergyScale::Logarithmic ? std::sqrt(low * high)
                                                   : 0.5 * (low + high);
}

double KineticEnergyFromMomentumJ(double momentumKgMPerS, double massKg) {
  const double c = Core::Const::c;
  const double pc = momentumKgMPerS * c;
  const double restEnergy = massKg * c * c;
  return pc * pc / (std::sqrt(pc * pc + restEnergy * restEnergy) + restEnergy);
}

double SpeedFromMomentumMPerS(double momentumKgMPerS, double massKg) {
  const double c = Core::Const::c;
  const double pc = momentumKgMPerS * c;
  const double restEnergy = massKg * c * c;
  return pc * c / std::sqrt(pc * pc + restEnergy * restEnergy);
}

double KineticEnergyFromSpeedJ(double speedMPerS, double massKg) {
  const double c = Core::Const::c;
  const double beta2 = (speedMPerS / c) * (speedMPerS / c);
  const double root = std::sqrt(1.0 - beta2);
  return massKg * c * c * beta2 / (root * (1.0 + root));
}

Core::Vec3 BoostVelocity(const Core::Vec3& velocityInMovingFrame,
                         const Core::Vec3& frameVelocity) {
  const double c2 = Core::Const::c * Core::Const::c;
  const double frameGamma = 1.0 / std::sqrt(1.0 - frameVelocity.NormSq() / c2);
  const double dot = frameVelocity.Dot(velocityInMovingFrame) / c2;
  const Core::Vec3 numerator = velocityInMovingFrame / frameGamma +
      frameVelocity + frameVelocity * (frameGamma / (1.0 + frameGamma) * dot);
  return numerator / (1.0 + dot);
}

Core::Status EvaluateKinematics(const ParticleInput& particle,
                                const Core::Vec3& bulkVelocityMPerS,
                                const Core::Vec3& magneticDirection,
                                ParticleKinematics* kinematics) {
  if (kinematics == nullptr) return Invalid("kinematics output is null");
  if (!std::isfinite(particle.weight) || particle.weight <= 0.0 ||
      !std::isfinite(particle.massKg) || particle.massKg <= 0.0 ||
      !std::isfinite(particle.momentumKgMPerS) ||
      particle.momentumKgMPerS < 0.0 || !std::isfinite(particle.mu) ||
      std::fabs(particle.mu) > 1.0 ||
      !FiniteVector(particle.plasmaVelocityDirection))
    return Invalid("sampled particle has a non-finite or unphysical state");
  if (!FiniteVector(bulkVelocityMPerS) ||
      bulkVelocityMPerS.Norm() >= Core::Const::c)
    return Invalid("sampled cell has a non-finite or superluminal bulk velocity");
  if (!FiniteVector(magneticDirection) ||
      std::fabs(magneticDirection.Norm() - 1.0) > 1.0e-12)
    return Invalid("sampled cell has no unit magnetic direction");

  // Plasma frame: the stored transport state.
  const double speed =
      SpeedFromMomentumMPerS(particle.momentumKgMPerS, particle.massKg);
  ParticleKinematics out;
  out.energyMeV[kPlasmaFrame] =
      KineticEnergyFromMomentumJ(particle.momentumKgMPerS, particle.massKg) /
      kJoulePerMeV;
  out.speedMPerS[kPlasmaFrame] = speed;
  out.parallelSpeedMPerS[kPlasmaFrame] = particle.mu * speed;

  // Inertial frame: boost the plasma-frame velocity vector with U.
  Core::Vec3 plasmaVelocity;
  if (speed > 0.0) {
    const double length = particle.plasmaVelocityDirection.Norm();
    if (!(length > 0.0))
      return Invalid("sampled particle has momentum but no velocity direction");
    plasmaVelocity = particle.plasmaVelocityDirection * (speed / length);
  }
  const Core::Vec3 inertial = BoostVelocity(plasmaVelocity, bulkVelocityMPerS);
  const double inertialSpeed = inertial.Norm();
  if (!std::isfinite(inertialSpeed) || inertialSpeed >= Core::Const::c)
    return Invalid("boosted particle velocity is not below c");
  out.energyMeV[kInertialFrame] =
      KineticEnergyFromSpeedJ(inertialSpeed, particle.massKg) / kJoulePerMeV;
  out.speedMPerS[kInertialFrame] = inertialSpeed;
  out.parallelSpeedMPerS[kInertialFrame] = inertial.Dot(magneticDirection);
  *kinematics = out;
  return Core::Status::OK();
}

Layout::Layout(std::size_t points, const std::vector<SpeciesChannels>& species)
    : points_(points) {
  std::size_t offset = 0;
  for (const SpeciesChannels& item : species) {
    const std::size_t frameStride =
        item.channels() * kChannelQuantities + kTailQuantities;
    channels_.push_back(item.channels());
    speciesOffset_.push_back(offset);
    speciesStride_.push_back(kFrames * frameStride);
    offset += kFrames * frameStride;
  }
  pointStride_ = offset;
}

std::size_t Layout::FrameBase(std::size_t point, std::size_t slot,
                              int frame) const {
  const std::size_t frameStride =
      channels_[slot] * kChannelQuantities + kTailQuantities;
  return point * pointStride_ + speciesOffset_[slot] +
         static_cast<std::size_t>(frame) * frameStride;
}

std::size_t Layout::Channel(std::size_t point, std::size_t slot, int frame,
                            std::size_t channel, std::size_t q) const {
  return FrameBase(point, slot, frame) + channel * kChannelQuantities + q;
}

std::size_t Layout::Tail(std::size_t point, std::size_t slot, int frame,
                         std::size_t q) const {
  return FrameBase(point, slot, frame) +
         channels_[slot] * kChannelQuantities + q;
}

void Accumulate(const Layout& layout, const SpeciesChannels& species,
                std::size_t point, std::size_t slot,
                double weight, double cellVolumeM3, double timeStepS,
                const ParticleKinematics& kinematics,
                std::vector<double>* sums) {
  std::vector<double>& s = *sums;
  const double perVolume = weight / cellVolumeM3;
  for (int frame = 0; frame < kFrames; ++frame) {
    const double speed = kinematics.speedMPerS[frame];
    const int channel = ChannelIndex(species.edgesMeV, kinematics.energyMeV[frame]);
    if (channel >= 0) {
      const std::size_t k = static_cast<std::size_t>(channel);
      s[layout.Channel(point, slot, frame, k, kDensity)] += perVolume;
      s[layout.Channel(point, slot, frame, k, kFlux)] += perVolume * speed;
      s[layout.Channel(point, slot, frame, k, kFieldAlignedFlux)] +=
          perVolume * kinematics.parallelSpeedMPerS[frame];
      s[layout.Channel(point, slot, frame, k, kFluence)] +=
          perVolume * speed * timeStepS;
      s[layout.Channel(point, slot, frame, k, kSamples)] += 1.0;
    } else {
      const bool below = channel == kBelowRange;
      s[layout.Tail(point, slot, frame, below ? kBelowDensity : kAboveDensity)] +=
          perVolume;
      s[layout.Tail(point, slot, frame, below ? kBelowFlux : kAboveFlux)] +=
          perVolume * speed;
      s[layout.Tail(point, slot, frame, below ? kBelowSamples : kAboveSamples)] +=
          1.0;
    }
  }
}

void FoldFluence(const Layout& layout,
                 const std::vector<SpeciesChannels>& species,
                 const std::vector<double>& sums,
                 std::vector<double>* cumulative) {
  for (std::size_t p = 0; p < layout.points(); ++p)
    for (std::size_t slot = 0; slot < species.size(); ++slot)
      for (int frame = 0; frame < kFrames; ++frame)
        for (std::size_t k = 0; k < species[slot].channels(); ++k) {
          const std::size_t i = layout.Channel(p, slot, frame, k, kFluence);
          (*cumulative)[i] += sums[i];
        }
}

Core::Status ParseProducts(const std::string& text, ProductSelection* products) {
  if (products == nullptr) return Invalid("product output is null");
  ProductSelection parsed;
  std::set<std::string> seen;
  std::stringstream input(text);
  std::string token;
  while (std::getline(input, token, ',')) {
    const std::size_t first = token.find_first_not_of(" \t");
    const std::size_t last = token.find_last_not_of(" \t");
    token = first == std::string::npos ? std::string()
                                       : token.substr(first, last - first + 1);
    if (!seen.insert(token).second)
      return Invalid("product '" + token + "' is listed twice");
    if (token == "density") parsed.density = true;
    else if (token == "flux") parsed.flux = true;
    else if (token == "spectrum") parsed.spectrum = true;
    else if (token == "fluence") parsed.fluence = true;
    else
      return Invalid("unknown product '" + token +
                     "' (use density, flux, spectrum, fluence)");
  }
  if (seen.empty()) return Invalid("at least one product is required");
  *products = parsed;
  return Core::Status::OK();
}

std::string FormatTecplot(const Layout& layout,
                          const std::vector<PointDefinition>& points,
                          const std::vector<SpeciesChannels>& species,
                          const ProductSelection& products,
                          const IntervalInfo& info,
                          const std::vector<double>& sums,
                          const std::vector<double>& cumulative) {
  const double n = static_cast<double>(info.sampleCalls);
  const char* frameName[kFrames] = {"plasma", "inertial"};
  std::ostringstream out;
  out << std::scientific << std::setprecision(10);
  out << "TITLE=\"srcSEP3D point sampling, AMPS output " << info.outputFileNumber
      << "\"\n";
  // Variables: channel geometry, then each selected product in both frames.
  out << "VARIABLES=\"E_min [MeV]\", \"E_max [MeV]\", \"E_center [MeV]\"";
  for (int frame = 0; frame < kFrames; ++frame) {
    const std::string f = frameName[frame];
    if (products.density) out << ", \"n_" << f << " [cm^-3]\"";
    if (products.flux)
      out << ", \"F_" << f << " [cm^-2 s^-1]\", \"F_parallel_" << f
          << " [cm^-2 s^-1]\"";
    if (products.spectrum)
      out << ", \"j_" << f << " [cm^-2 s^-1 sr^-1 MeV^-1]\"";
    if (products.fluence)
      out << ", \"fluence_interval_" << f << " [cm^-2]\", \"fluence_" << f
          << " [cm^-2]\", \"fluence_intensity_" << f
          << " [cm^-2 sr^-1 MeV^-1]\"";
    out << ", \"samples_" << f << "\"";
  }
  out << '\n';
  out << "DATASETAUXDATA amps_output_file_number=\"" << info.outputFileNumber
      << "\"\n";
  out << "DATASETAUXDATA sample_calls=\"" << info.sampleCalls << "\"\n";
  out << "DATASETAUXDATA first_sample_time_s=\"" << info.firstSampleTimeS
      << "\"\n";
  out << "DATASETAUXDATA last_sample_time_s=\"" << info.lastSampleTimeS
      << "\"\n";

  for (std::size_t p = 0; p < points.size(); ++p) {
    for (std::size_t slot = 0; slot < species.size(); ++slot) {
      const SpeciesChannels& sp = species[slot];
      out << "ZONE T=\"" << points[p].name << " " << sp.symbol << "\", I="
          << sp.channels() << ", DATAPACKING=POINT\n";
      out << "AUXDATA point=\"" << points[p].name << "\"\n";
      out << "AUXDATA species=\"" << sp.symbol << "\"\n";
      out << "AUXDATA position_x_m=\"" << points[p].positionM.x << "\"\n";
      out << "AUXDATA position_y_m=\"" << points[p].positionM.y << "\"\n";
      out << "AUXDATA position_z_m=\"" << points[p].positionM.z << "\"\n";
      if (points[p].cellKnown) {
        // The sampled AMPS cell: centre and edge lengths [m], written with
        // 17 significant digits so the doubles round-trip exactly (a reader
        // can then identify the cell bit-for-bit, e.g. by srcSEP3D's
        // position-hash cell id in the cells.csv publication).
        out << std::setprecision(17);
        const C3& center = points[p].cellCenterM;
        const C3& size = points[p].cellSizeM;
        out << "AUXDATA cell_center_x_m=\"" << center.x << "\"\n";
        out << "AUXDATA cell_center_y_m=\"" << center.y << "\"\n";
        out << "AUXDATA cell_center_z_m=\"" << center.z << "\"\n";
        out << "AUXDATA cell_size_x_m=\"" << size.x << "\"\n";
        out << "AUXDATA cell_size_y_m=\"" << size.y << "\"\n";
        out << "AUXDATA cell_size_z_m=\"" << size.z << "\"\n";
        out << std::setprecision(10);
      }
      // Particles outside [E_min, E_max]: interval-mean density [cm^-3],
      // flux [cm^-2 s^-1] and the number of macroparticle samples.
      for (int frame = 0; frame < kFrames; ++frame) {
        const std::string f = frameName[frame];
        out << "AUXDATA below_range_density_" << f << "_cm3=\""
            << sums[layout.Tail(p, slot, frame, kBelowDensity)] / n *
                   kPerM3ToPerCm3 << "\"\n";
        out << "AUXDATA below_range_flux_" << f << "_cm2_s=\""
            << sums[layout.Tail(p, slot, frame, kBelowFlux)] / n *
                   kPerM2ToPerCm2 << "\"\n";
        out << "AUXDATA below_range_samples_" << f << "=\""
            << static_cast<std::uint64_t>(
                   sums[layout.Tail(p, slot, frame, kBelowSamples)]) << "\"\n";
        out << "AUXDATA above_range_density_" << f << "_cm3=\""
            << sums[layout.Tail(p, slot, frame, kAboveDensity)] / n *
                   kPerM3ToPerCm3 << "\"\n";
        out << "AUXDATA above_range_flux_" << f << "_cm2_s=\""
            << sums[layout.Tail(p, slot, frame, kAboveFlux)] / n *
                   kPerM2ToPerCm2 << "\"\n";
        out << "AUXDATA above_range_samples_" << f << "=\""
            << static_cast<std::uint64_t>(
                   sums[layout.Tail(p, slot, frame, kAboveSamples)]) << "\"\n";
      }
      for (std::size_t k = 0; k < sp.channels(); ++k) {
        const double widthMeV = sp.edgesMeV[k + 1] - sp.edgesMeV[k];
        out << sp.edgesMeV[k] << ' ' << sp.edgesMeV[k + 1] << ' '
            << ChannelCenterMeV(sp, k);
        for (int frame = 0; frame < kFrames; ++frame) {
          auto at = [&](std::size_t q) {
            return sums[layout.Channel(p, slot, frame, k, q)];
          };
          const double flux = at(kFlux) / n * kPerM2ToPerCm2;
          if (products.density) out << ' ' << at(kDensity) / n * kPerM3ToPerCm3;
          if (products.flux)
            out << ' ' << flux << ' '
                << at(kFieldAlignedFlux) / n * kPerM2ToPerCm2;
          if (products.spectrum) out << ' ' << flux / (kFourPi * widthMeV);
          if (products.fluence) {
            const double total =
                cumulative[layout.Channel(p, slot, frame, k, kFluence)] *
                kPerM2ToPerCm2;
            out << ' ' << at(kFluence) * kPerM2ToPerCm2 << ' ' << total << ' '
                << total / (kFourPi * widthMeV);
          }
          out << ' ' << static_cast<std::uint64_t>(at(kSamples));
        }
        out << '\n';
      }
    }
  }
  return out.str();
}

Core::Status ParsePointFile(const std::string& text,
                            const std::string& displayName,
                            std::vector<PointDefinition>* points) {
  if (points == nullptr) return Invalid("point output is null");
  std::vector<PointDefinition> parsed;
  std::istringstream input(text);
  std::string line;
  std::size_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::size_t comment = line.find_first_of("!#");
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream fields(line);
    std::vector<std::string> tokens;
    std::string token;
    while (fields >> token) tokens.push_back(token);
    if (tokens.empty()) continue;
    const std::string where = displayName + ":" + std::to_string(lineNumber);
    if (tokens.size() != 3 && tokens.size() != 4)
      return Invalid(where + ": expected 'x y z [name]' (metres), found " +
                     std::to_string(tokens.size()) + " fields");
    double xyz[3];
    for (int i = 0; i < 3; ++i)
      if (!ParseNumber(tokens[i], &xyz[i]))
        return Invalid(where + ": '" + tokens[i] + "' is not a finite number");
    PointDefinition point;
    point.positionM = Core::Vec3(xyz[0], xyz[1], xyz[2]);
    point.name = tokens.size() == 4 ? tokens[3] : std::to_string(parsed.size());
    point.origin = where;
    parsed.push_back(point);
  }
  if (parsed.empty())
    return Invalid(displayName + ": the point file contains no points");
  *points = std::move(parsed);
  return Core::Status::OK();
}

}  // namespace PointSampling
}  // namespace Output
}  // namespace SEP3D
