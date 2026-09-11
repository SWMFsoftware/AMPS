#ifndef SWCME_SEP_INTERFACE_HPP
#define SWCME_SEP_INTERFACE_HPP

// ============================================================================
// swcme_sep_interface.hpp
// ----------------------------------------------------------------------------
// Thin, stable SWCME-to-SEP/AMPS adapter.
//
// The physics remains owned by swcme1d::Model / swcme3d::Model.  This adapter
// does not duplicate Parker, shock, connectivity, or DSA equations.  Its job is
// to package those validated results into a small SI-unit contract convenient
// for a particle/focused-transport caller.
//
// Hot-loop contract
// -----------------
//   * prepare(t) is called once per SWCME field-update cadence;
//   * PreparedStep is immutable/read-only while particle threads query it;
//   * evaluate_background() performs no dynamic allocation for an N=1 query;
//   * source-surface construction/connectivity are update-cadence operations,
//     not per-particle operations, and may allocate diagnostic/source vectors.
//
// This separation is intentional: AMPS can cache BackgroundState queries as it
// already does for fields, while shock/source state is refreshed only when the
// background step changes.
// ============================================================================

#include "swcme1d.hpp"
#include "swcme3d.hpp"
#include "swcme_sep_source.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace swcme {
namespace sep {

struct BackgroundState {
  ModelStatus status;
  std::array<double,3> position_m{{0.0,0.0,0.0}};
  double density_m3 = 0.0;
  std::array<double,3> velocity_m_s{{0.0,0.0,0.0}};
  std::array<double,3> magnetic_T{{0.0,0.0,0.0}};
  double magnetic_magnitude_T = 0.0;
  double div_velocity_s_inv = 0.0;
};

struct SourceSurface {
  ModelStatus status;
  double time_s = 0.0;
  double total_surface_area_m2 = 0.0;
  double active_surface_area_m2 = 0.0;
  std::size_t patch_count = 0;
  std::size_t active_patch_count = 0;
  std::vector<SEPSourceState> patches;
};

inline std::string spectrum_configuration_manifest(const SpectrumConfig& c) {
  std::ostringstream out;
  out.setf(std::ios::scientific);
  out << std::setprecision(17)
      << "sep_source_contract_version=1\n"
      << "sep_particle_mass_kg=" << c.particle_mass_kg << '\n'
      << "sep_charge_number=" << c.charge_number << '\n'
      << "sep_energy_min_MeV=" << c.kinetic_energy_min_MeV << '\n'
      << "sep_energy_max_MeV=" << c.kinetic_energy_max_MeV << '\n'
      << "sep_reference_energy_MeV=" << c.reference_energy_MeV << '\n'
      << "sep_normalization=" << normalization_mode_name(c.normalization) << '\n'
      << "sep_reference_differential_intensity_SI=";
  if (std::isfinite(c.reference_differential_intensity_SI))
    out << c.reference_differential_intensity_SI;
  else
    out << "NA";
  out << '\n';
  return out.str();
}

class Interface1D {
public:
  using PreparedStep = swcme1d::StepState;

  explicit Interface1D(const swcme1d::Params& params = swcme1d::Params{},
                       const SpectrumConfig& spectrum = SpectrumConfig{})
      : model_(params), spectrum_(spectrum) {
    const ModelStatus spectrum_status=validate_spectrum_config(spectrum_);
    if (spectrum_status.failure()) throw std::invalid_argument(spectrum_status.summary());
    const swcme::config::ValidationResult model_status=model_.validate();
    if (!model_status.ok()) throw std::invalid_argument(model_status.summary("swcme SEP 1D"));
  }

  const swcme1d::Model& model() const noexcept { return model_; }
  const SpectrumConfig& spectrum_config() const noexcept { return spectrum_; }

  PreparedStep prepare(double time_s) const { return model_.prepare_step(time_s); }

  ModelStatus evaluate_background(const PreparedStep& step,
                                  double radius_m,
                                  BackgroundState& out) const {
    // Reject a foreign prepared state before clearing the caller's output.
    // This ordering is part of PST02: an ownership failure is transactional
    // and cannot erase a previously valid AMPS background record.
    const ModelStatus ownership=model_.validate_prepared_state(
        step,"SEP Interface1D background");
    if (!ownership.ok()) return ownership;
    out=BackgroundState{};
    out.position_m={{radius_m,0.0,0.0}};
    double n=0.0,V=0.0,Br=0.0,Bphi=0.0,Bmag=0.0,divV=0.0;
    const ModelStatus status=model_.evaluate_radii_with_B_div_checked(
        step,&radius_m,&n,&V,&Br,&Bphi,&Bmag,&divV,1);
    out.status=status;
    if (status.failure()) return status;
    out.density_m3=n;
    out.velocity_m_s={{V,0.0,0.0}};
    // In the +X equatorial reference used by the common 1-D/3-D regression,
    // the Parker azimuthal basis is +Y; signed Bphi therefore maps directly to
    // the Y component and preserves the production Parker polarity.
    out.magnetic_T={{Br,Bphi,0.0}};
    out.magnetic_magnitude_T=Bmag;
    out.div_velocity_s_inv=divV;
    return status;
  }

  ModelStatus source_at_shock(const PreparedStep& step,
                              SEPSourceState& out) const {
    // The checked 1-D source path preserves `out` on ownership failure; only a
    // state prepared by this adapter's Model can reach spectrum conversion.
    acceleration::ShockAccelerationState a;
    const ModelStatus ownership=
        model_.shock_acceleration_state_checked(step,a);
    if (!ownership.ok()) return ownership;
    return make_source_state(a,spectrum_,false,false,0,0.0,out);
  }

  std::string resolved_manifest() const {
    return swcme1d::resolved_configuration_manifest(model_.GetParams())+
           spectrum_configuration_manifest(spectrum_);
  }

private:
  swcme1d::Model model_;
  SpectrumConfig spectrum_;
};

class Interface3D {
public:
  using PreparedStep = swcme3d::StepState;

  explicit Interface3D(const swcme3d::Params& params,
                       const SpectrumConfig& spectrum = SpectrumConfig{})
      : model_(params), params_(params), spectrum_(spectrum) {
    const ModelStatus spectrum_status=validate_spectrum_config(spectrum_);
    if (spectrum_status.failure()) throw std::invalid_argument(spectrum_status.summary());
    const swcme::config::ValidationResult model_status=model_.validate();
    if (!model_status.ok()) throw std::invalid_argument(model_status.summary("swcme SEP 3D"));
  }

  const swcme3d::Model& model() const noexcept { return model_; }
  const SpectrumConfig& spectrum_config() const noexcept { return spectrum_; }

  PreparedStep prepare(double time_s) const { return model_.prepare_step(time_s); }

  // Single-point AMPS hot-loop query.  All outputs are SI.  No std::vector or
  // heap object is created: the existing checked batch API is called with N=1
  // and stack scalars, preserving the explicit status propagation from Fix 11.
  ModelStatus evaluate_background(const PreparedStep& step,
                                  const std::array<double,3>& position_m,
                                  BackgroundState& out) const {
    // Perform ownership validation before resetting `out`, preserving the
    // previous record exactly when an AMPS caller accidentally mixes adapters.
    const ModelStatus ownership=model_.validate_prepared_state(
        step,"SEP Interface3D background");
    if (!ownership.ok()) return ownership;
    out=BackgroundState{};
    out.position_m=position_m;
    const double x=position_m[0],y=position_m[1],z=position_m[2];
    double n=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0,divV=0.0;
    const ModelStatus status=model_.evaluate_cartesian_with_B_div_checked(
        step,&x,&y,&z,&n,&vx,&vy,&vz,&bx,&by,&bz,&divV,1);
    out.status=status;
    if (status.failure()) return status;
    out.density_m3=n;
    out.velocity_m_s={{vx,vy,vz}};
    out.magnetic_T={{bx,by,bz}};
    out.magnetic_magnitude_T=std::hypot(bx,std::hypot(by,bz));
    out.div_velocity_s_inv=divV;
    if (!std::isfinite(out.magnetic_magnitude_T)) {
      out.status=ModelStatus::make(StatusCode::NonFiniteResult,
                                   "SEP Interface3D background |B|");
      return out.status;
    }
    return status;
  }

  ModelStatus source_at_direction(const PreparedStep& step,
                                  const std::array<double,3>& direction,
                                  SEPSourceState& out,
                                  std::size_t source_id=0,
                                  double patch_area_m2=0.0) const {
    // Guard the outer interface so the caller's source record remains
    // byte-for-byte unchanged on mismatch.  The model repeats the check at its
    // checked physics boundary as defense in depth for direct users.
    const ModelStatus ownership=model_.validate_prepared_state(
        step,"SEP Interface3D directional source");
    if (!ownership.ok()) return ownership;
    acceleration::ShockAccelerationState acceleration_state;
    const double u[3]={direction[0],direction[1],direction[2]};
    const ModelStatus acceleration_status=
        model_.shock_acceleration_state_checked(step,u,acceleration_state);
    if (!acceleration_status.ok()) {
      out=SEPSourceState{};
      out.status=acceleration_status;
      return out.status;
    }
    return make_source_state(acceleration_state,spectrum_,false,false,
                             source_id,patch_area_m2,out);
  }

  // Build one source record per shock-mesh cell.  Physics is sampled at the
  // direction of the physical triangle centroid using the canonical local
  // shock-state path; the triangle area comes only from the validated mesh.
  // Area fractions are normalized over cells that actually host a physical
  // SOURCE-mode fast shock.  Sub-fast portions remain present/auditable with
  // active=false and exactly zero source weight.
  ModelStatus build_shock_surface_source(const PreparedStep& step,
                                         std::size_t n_theta,
                                         std::size_t n_phi,
                                         SourceSurface& out) const {
    // Mesh construction is an unchecked value-returning API, so validate at
    // the adapter boundary before clearing `out` or allocating mesh storage.
    const ModelStatus ownership=model_.validate_prepared_state(
        step,"SEP Interface3D surface source");
    if (!ownership.ok()) return ownership;
    out=SourceSurface{};
    out.time_s=step.time_s;
    try {
      const swcme3d::ShockMesh mesh=model_.build_shock_mesh(step,n_theta,n_phi);
      swcme3d::TriMetrics metrics;
      model_.compute_triangle_metrics(mesh,metrics);
      out.patch_count=metrics.area.size();
      out.patches.reserve(out.patch_count);

      long double total_area=0.0L;
      long double active_area=0.0L;
      for (std::size_t i=0;i<metrics.area.size();++i) {
        total_area+=static_cast<long double>(metrics.area[i]);
        std::array<double,3> u{{metrics.cx[i],metrics.cy[i],metrics.cz[i]}};
        const double norm=std::hypot(u[0],std::hypot(u[1],u[2]));
        if (!std::isfinite(norm) || !(norm>0.0)) {
          out.status=ModelStatus::make(StatusCode::InvalidMesh,
                                       "SEP source centroid direction",i);
          return out.status;
        }
        u[0]/=norm; u[1]/=norm; u[2]/=norm;
        SEPSourceState source;
        const ModelStatus status=source_at_direction(step,u,source,i,metrics.area[i]);
        if (status.failure()) { out.status=status; return status; }
        if (source.active) {
          active_area+=static_cast<long double>(metrics.area[i]);
          ++out.active_patch_count;
        }
        out.patches.push_back(source);
      }

      out.total_surface_area_m2=static_cast<double>(total_area);
      out.active_surface_area_m2=static_cast<double>(active_area);
      if (!std::isfinite(out.total_surface_area_m2) ||
          !(out.total_surface_area_m2>0.0) ||
          !std::isfinite(out.active_surface_area_m2) ||
          out.active_surface_area_m2<0.0) {
        out.status=ModelStatus::make(StatusCode::NonFiniteResult,
                                     "SEP source surface area");
        return out.status;
      }

      for (SEPSourceState& source : out.patches) {
        source.active_surface_area_m2=out.active_surface_area_m2;
        if (source.active && out.active_surface_area_m2>0.0) {
          source.area_fraction=source.patch_area_m2/out.active_surface_area_m2;
          source.relative_patch_weight=
              source.relative_source_weight_per_area*source.area_fraction;
        } else {
          source.area_fraction=0.0;
          source.relative_patch_weight=0.0;
        }
      }
      out.status=ModelStatus::success();
      return out.status;
    } catch (const std::exception&) {
      out.status=ModelStatus::make(StatusCode::InvalidMesh,
                                   "SEP build_shock_surface_source");
      return out.status;
    }
  }

  // Convert the selected production Parker-line/shock intersection directly to
  // an SEP source record.  No second shock equation is evaluated here; the
  // direction is taken from the already selected cobpoint and routed through
  // the same source_at_direction() path used by surface patches.
  ModelStatus source_at_observer_cobpoint(
      const PreparedStep& step,
      const std::array<double,3>& observer_m,
      SEPSourceState& out,
      swcme3d::ConnectivityState* connectivity_out=nullptr,
      const swcme3d::ConnectivityOptions& options=swcme3d::ConnectivityOptions{}) const {
    // Connectivity is evaluated only after ownership succeeds.  This leaves
    // both optional connectivity output and the SEP source output untouched
    // when the prepared state belongs to another model.
    const ModelStatus ownership=model_.validate_prepared_state(
        step,"SEP Interface3D observer source");
    if (!ownership.ok()) return ownership;
    const double observer[3]={observer_m[0],observer_m[1],observer_m[2]};
    const swcme3d::ConnectivityState connectivity=
        model_.observer_connectivity(step,observer,options);
    if (connectivity_out) *connectivity_out=connectivity;

    if (!connectivity.connected || connectivity.roots.empty()) {
      out=SEPSourceState{};
      out.time_s=step.time_s;
      out.spectrum=spectrum_;
      out.connection_evaluated=true;
      out.connected=false;
      out.status=ModelStatus::make(StatusCode::NoConnection,
                                   "SEP observer cobpoint");
      return out.status;
    }
    if (connectivity.selected_root>=connectivity.roots.size()) {
      out=SEPSourceState{};
      out.status=ModelStatus::make(StatusCode::GeometryFailure,
                                   "SEP selected cobpoint index");
      return out.status;
    }

    const swcme3d::ConnectivityRoot& root=
        connectivity.roots[connectivity.selected_root];
    std::array<double,3> u{{root.position_m[0],root.position_m[1],root.position_m[2]}};
    const double norm=std::hypot(u[0],std::hypot(u[1],u[2]));
    if (!std::isfinite(norm) || !(norm>0.0)) {
      out=SEPSourceState{};
      out.status=ModelStatus::make(StatusCode::GeometryFailure,
                                   "SEP cobpoint direction");
      return out.status;
    }
    u[0]/=norm; u[1]/=norm; u[2]/=norm;

    acceleration::ShockAccelerationState acceleration_state;
    const double ud[3]={u[0],u[1],u[2]};
    const ModelStatus acceleration_status=
        model_.shock_acceleration_state_checked(step,ud,acceleration_state);
    if (!acceleration_status.ok()) {
      out=SEPSourceState{};
      out.status=acceleration_status;
      return out.status;
    }
    return make_source_state(acceleration_state,spectrum_,true,true,0,0.0,out);
  }

  std::string resolved_manifest() const {
    return swcme3d::resolved_configuration_manifest(params_)+
           spectrum_configuration_manifest(spectrum_);
  }

private:
  swcme3d::Model model_;
  swcme3d::Params params_;
  SpectrumConfig spectrum_;
};

}  // namespace sep
}  // namespace swcme

#endif
