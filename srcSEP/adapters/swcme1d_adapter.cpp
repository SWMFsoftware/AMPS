#include "swcme1d_adapter.h"

#include "../util/sep_background_snapshot.h"

// The makefile's adapter-only rule supplies AMPS/src/models/swcme as an
// include root.  This is the sole srcSEP production translation unit that
// consumes the canonical 1-D SWCME C++ interface.
#include "swcme1d.hpp"

#include <algorithm>
#include <cmath>

namespace {

swcme1d::Model model;
swcme1d::StepState state{};
bool state_prepared = false;
bool clamp_sheath = true;

void AssertSwcmeMayWrite() {
  SEP::Background::SnapshotStore::Instance().AssertProviderMayWrite(
      SEP::Background::Provider::Swcme);
}

}  // namespace

namespace SEP {
namespace SW1DAdapter {

void Configure(Scenario scenario) {
  AssertSwcmeMayWrite();

  if (scenario == Scenario::Fast) {
    model.SetAmbient(400.0, 6.0, 5.0, 1.2e5)
        .SetCME(1.05, 1900.0, 8.0e-8)
        .SetGeometry(0.12, 0.22)
        .SetSmoothing(0.010, 0.020, 0.030)
        .SetSheathEjecta(1.25, 2.0, 1.12, 0.50, 0.80);
  } else {
    model.SetAmbient(380.0, 5.0, 4.5, 1.0e5)
        .SetCME(1.05, 950.0, 3.0e-8)
        .SetGeometry(0.08, 0.18)
        .SetSmoothing(0.015, 0.030, 0.050)
        .SetSheathEjecta(1.15, 1.5, 1.08, 0.60, 0.90);
  }

  clamp_sheath = true;
  state_prepared = false;
}

void PrepareState(double epoch_seconds) {
  AssertSwcmeMayWrite();
  state = model.prepare_step(epoch_seconds);
  state_prepared = true;
}

void EnableSheathClamp(bool enabled) { clamp_sheath = enabled; }

double ShockRadiusM() {
  return state_prepared ? state.r_sh_m : 0.0;
}

double ShockSpeedMPerS() {
  return state_prepared ? state.V_sh_ms : 0.0;
}

double CompressionRatio() {
  return state_prepared ? state.rc : 1.0;
}

double DlnB_Dr_at_r(double radius_m) {
  if (!state_prepared || !std::isfinite(radius_m) || radius_m <= 0.0) {
    return 0.0;
  }

  const double radius_au = radius_m / swcme1d::AU;
  const double k = state.k_AU;
  const double k2r2 = (k * radius_au) * (k * radius_au);
  return (-2.0 / radius_au +
          (k * k * radius_au) / (1.0 + k2r2)) /
         swcme1d::AU;
}

bool QueryAtRadius(double radius_m,
                   double& number_density_m3,
                   double& speed_m_s,
                   double& divergence_s_inv,
                   bool apply_sheath_clamp) {
  if (!state_prepared || !std::isfinite(radius_m)) return false;

  double radius = std::max(radius_m, 1.05 * swcme1d::Rs);
  double radial_field_t = 0.0;
  double azimuthal_field_t = 0.0;
  double field_magnitude_t = 0.0;
  double divergence = 0.0;
  model.evaluate_radii_with_B_div(
      state, &radius, &number_density_m3, &speed_m_s, &radial_field_t,
      &azimuthal_field_t, &field_magnitude_t, &divergence, 1);

  if (!std::isfinite(number_density_m3) || number_density_m3 < 0.0) {
    number_density_m3 = 0.0;
  }
  if (!std::isfinite(speed_m_s)) speed_m_s = 0.0;
  if (!std::isfinite(divergence)) divergence = 0.0;

  if (clamp_sheath && apply_sheath_clamp && state.rc > 1.0 &&
      radius <= state.r_sh_m && radius >= state.r_le_m) {
    const double upstream_density =
        swcme1d::Model::density_upstream(state, radius);
    if (std::isfinite(upstream_density)) {
      number_density_m3 = std::max(number_density_m3, upstream_density);
    }

    const double minimum_speed = state.V_up_ms;
    const double maximum_speed = std::max(state.V_up_ms, state.V_sh_ms);
    speed_m_s = std::min(std::max(speed_m_s, minimum_speed), maximum_speed);
  }

  divergence_s_inv = divergence;
  return true;
}

void WriteRadialProfileFromR(const double* radii_m,
                             int count,
                             const char* file_name,
                             double epoch_seconds) {
  if (!state_prepared || radii_m == nullptr || count <= 0 ||
      file_name == nullptr) {
    return;
  }
  model.write_tecplot_radial_profile_from_r(
      state, radii_m, count, file_name, epoch_seconds);
}

void WriteShockVsTime(double duration_seconds,
                      int sample_count,
                      const char* file_name) {
  if (sample_count <= 0 || file_name == nullptr) return;
  model.write_tecplot_shock_vs_time(
      duration_seconds, sample_count, file_name);
}

}  // namespace SW1DAdapter
}  // namespace SEP
