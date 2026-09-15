#ifndef _SEP_PHYSICAL_UNITS_H_
#define _SEP_PHYSICAL_UNITS_H_

#include <cmath>
#include <stdexcept>

namespace SEP {
namespace Units {

// These deliberately small wrappers prevent dimensionally different values
// from being interchanged at new geometry/source interfaces.  Legacy AMPS APIs
// still consume doubles, so callers cross that boundary explicitly with
// Value(); the unit-bearing type remains visible everywhere calculations are
// assembled and reviewed.
template <class Tag>
class Quantity {
 public:
  explicit Quantity(double value) : value_(value) {}
  double Value() const { return value_; }

 private:
  double value_;
};

struct EnergyJTag {};
struct MomentumKgMPerSTag {};
struct NumberDensityPerM3Tag {};
struct SpeedMPerSTag {};
struct AreaM2Tag {};
struct VolumeM3Tag {};
struct TimeSTag {};
struct LengthMTag {};
struct MagneticFieldTTag {};

typedef Quantity<EnergyJTag> EnergyJ;
typedef Quantity<MomentumKgMPerSTag> MomentumKgMPerS;
typedef Quantity<NumberDensityPerM3Tag> NumberDensityPerM3;
typedef Quantity<SpeedMPerSTag> SpeedMPerS;
typedef Quantity<AreaM2Tag> AreaM2;
typedef Quantity<VolumeM3Tag> VolumeM3;
typedef Quantity<TimeSTag> TimeS;
typedef Quantity<LengthMTag> LengthM;
typedef Quantity<MagneticFieldTTag> MagneticFieldT;

// One electron volt is exact in SI because the elementary charge is exact.
// Keeping this conversion here prevents MeV values from reaching SI momentum
// routines through an unlabelled double.
inline EnergyJ EnergyFromMeV(double energy_MeV) {
  const double joule_per_MeV = 1.602176634e-13;
  if (!std::isfinite(energy_MeV) || energy_MeV < 0.0) {
    throw std::invalid_argument("energy_MeV must be finite and non-negative");
  }
  return EnergyJ(energy_MeV * joule_per_MeV);
}

}  // namespace Units
}  // namespace SEP

#endif
