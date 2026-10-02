#ifndef SEP3D_OUTPUT_SHOCK_HISTORY_H
#define SEP3D_OUTPUT_SHOCK_HISTORY_H

// Portable publication boundary for the native propagation telemetry. MPI and
// AMPS ownership belong to main_lib.cpp; this writer receives already reduced
// facts, never computes a second CME trajectory and never creates particles.
#include "../core/sep3d_types.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

namespace SEP3D { namespace Output {
struct ShockHistorySample {
  // Sphere radius/speed, or SSE apex distance/speed. Off-axis arrival must
  // use the finite directional front; this CSV is not an observer-hit test.
  double timeS=0.0, radiusM=0.0, speedMPerS=0.0;
  double mpiRadiusSpreadM=0.0, mpiClockSpreadS=0.0;
  std::uint64_t tick=0, generation=0, particles=0, injections=0;
  bool active=false;
  std::string providerIdentity, configurationFingerprint;
};

inline Core::Status CheckShockHistorySample(const ShockHistorySample& sample, double dt) {
  if (!std::isfinite(dt) || dt<=0 || !std::isfinite(sample.timeS) || sample.timeS<0 ||
      !std::isfinite(sample.radiusM) || sample.radiusM<=0 ||
      !std::isfinite(sample.speedMPerS) || sample.speedMPerS<=0 ||
      !std::isfinite(sample.mpiRadiusSpreadM) || sample.mpiRadiusSpreadM<0 || sample.mpiRadiusSpreadM>1.0 ||
      !std::isfinite(sample.mpiClockSpreadS) || sample.mpiClockSpreadS<0 || sample.mpiClockSpreadS>1e-9 ||
      std::fabs(sample.timeS-static_cast<double>(sample.tick)*dt)>1e-8*std::max(1.0,sample.timeS) ||
      sample.generation==0 || sample.providerIdentity.empty() || sample.configurationFingerprint.empty())
    return Core::Status(Core::StatusCode::InvalidInput,"invalid native shock history clock/geometry/identity/MPI agreement");
  if (sample.particles!=0 || sample.injections!=0)
    return Core::Status(Core::StatusCode::ConfigurationConflict,"shock-propagation requires zero particles and cumulative injections");
  return Core::Status::OK();
}

class ShockHistoryWriter {
 public:
  Core::Status Open(const std::filesystem::path& file, double dt) {
    if (stream_.is_open() || hasRow_)
      return Core::Status(Core::StatusCode::InvalidTransition,"native writer can only open once");
    dt_=dt;
    std::error_code error;
    if (!std::isfinite(dt) || dt<=0 || std::filesystem::exists(file,error) || error)
      return Core::Status(Core::StatusCode::InvalidInput,"native history path already exists; choose a fresh --output-dir");
    if (!file.parent_path().empty()) std::filesystem::create_directories(file.parent_path(),error);
    if (error) return Core::Status(Core::StatusCode::Error,"cannot create native history directory: "+error.message());
    stream_.open(file);
    if (!stream_) return Core::Status(Core::StatusCode::Error,"cannot open native history");
    stream_ << "time_s,tick,shock_radius_m,shock_speed_m_s,shock_active,generation,particle_count,injected_particle_count,mpi_radius_spread_m,mpi_clock_spread_s\n";
    stream_ << std::setprecision(17);
    return Core::Status::OK();
  }
  Core::Status Append(const ShockHistorySample& sample) {
    if (!stream_.is_open())
      return Core::Status(Core::StatusCode::InvalidTransition,"native history writer is not open");
    Core::Status status=CheckShockHistorySample(sample,dt_);
    if (!status.ok()) return status;
    if ((!hasRow_ && (sample.tick!=0 || sample.timeS!=0)) ||
        (hasRow_ && (sample.tick!=last_.tick+1 || sample.generation<=last_.generation ||
                    sample.radiusM<=last_.radiusM || sample.providerIdentity!=last_.providerIdentity ||
                    sample.configurationFingerprint!=last_.configurationFingerprint)))
      return Core::Status(Core::StatusCode::InvalidInput,"native history must start at zero and retain contiguous ticks, increasing generations/radii and provider identity");
    stream_ << sample.timeS << ',' << sample.tick << ',' << sample.radiusM << ',' << sample.speedMPerS << ','
            << (sample.active?1:0) << ',' << sample.generation << ',' << sample.particles << ',' << sample.injections << ','
            << sample.mpiRadiusSpreadM << ',' << sample.mpiClockSpreadS << '\n';
    // Close ownership is explicit. Flushing each boundary leaves useful partial
    // telemetry after a failed run without qualifying it as complete evidence.
    stream_.flush();
    if (!stream_) return Core::Status(Core::StatusCode::Error,"native history write failed");
    last_=sample; hasRow_=true;
    return Core::Status::OK();
  }
  Core::Status Close() {
    if (!stream_.is_open() || !hasRow_)
      return Core::Status(Core::StatusCode::InvalidTransition,"native history cannot close without observed rows");
    stream_.flush(); const bool good=static_cast<bool>(stream_); stream_.close();
    if (!good || stream_.fail()) return Core::Status(Core::StatusCode::Error,"native history close failed");
    return Core::Status::OK();
  }
 private:
  std::ofstream stream_;
  double dt_=0.0;
  bool hasRow_=false;
  ShockHistorySample last_;
};

}}  // namespace SEP3D::Output
#endif
