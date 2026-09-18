// ============================================================================
// srcSEP -> canonical SWCME adapter.
//
// This header deliberately contains no SWCME include and exposes no SWCME
// type.  AMPS translation units may include it without adding
// AMPS/src/models/swcme to their include path.  The provider-specific header is
// confined to swcme1d_adapter.cpp, matching srcSEP3D's adapter boundary.
// ============================================================================

#ifndef SEP_ADAPTERS_SWCME1D_ADAPTER_H
#define SEP_ADAPTERS_SWCME1D_ADAPTER_H

namespace SEP {
namespace SW1DAdapter {

enum class Scenario {
  Fast,
  Slow
};

// Configure the process-wide 1-D SW+CME model.  Call PrepareState() after
// configuration and once for each physical background epoch.
void Configure(Scenario scenario);
void PrepareState(double epoch_seconds);

// Optional monotonic stabilization applied only inside the modeled sheath.
void EnableSheathClamp(bool enabled = true);

// Provider-neutral scalar views used by the rest of srcSEP.  Keeping these
// accessors here prevents the provider's prepared-state type from leaking
// through sep.h.
double ShockRadiusM();
double ShockSpeedMPerS();
double CompressionRatio();
double DlnB_Dr_at_r(double radius_m);

// Query density [m^-3], radial speed [m/s], and divergence [s^-1] at radius_m.
bool QueryAtRadius(double radius_m,
                   double& number_density_m3,
                   double& speed_m_s,
                   double& divergence_s_inv,
                   bool apply_sheath_clamp = true);

// Diagnostics remain owned by the adapter so their SWCME StepState argument
// never becomes part of the srcSEP application interface.
void WriteRadialProfileFromR(const double* radii_m,
                             int count,
                             const char* file_name,
                             double epoch_seconds);
void WriteShockVsTime(double duration_seconds,
                      int sample_count,
                      const char* file_name);

}  // namespace SW1DAdapter
}  // namespace SEP

#endif  // SEP_ADAPTERS_SWCME1D_ADAPTER_H
