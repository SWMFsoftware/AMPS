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
//
// Transport-context contract
// --------------------------
// Pressure, Parker path length, and magnetic focusing are obtained from the
// same swcme_solarwind.hpp helpers used by direct model/connectivity queries.
// The adapter only packages those SI values and prepared-state provenance;
// there is no AMPS-specific physical formula or unit convention to drift.
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

  // Preserve the provenance of the prepared state used for this record.  AMPS
  // can archive or assert these tokens at its field-update boundary; the
  // adapter copies them verbatim and never derives a new identity.
  ModelIdentity owner_model_identity = 0;
  ConfigurationDigest configuration_digest = 0;

  std::array<double,3> position_m{{0.0,0.0,0.0}};
  double density_m3 = 0.0;
  double pressure_Pa = 0.0;
  std::array<double,3> velocity_m_s{{0.0,0.0,0.0}};
  std::array<double,3> magnetic_T{{0.0,0.0,0.0}};
  double magnetic_magnitude_T = 0.0;
  double div_velocity_s_inv = 0.0;
  double focusing_length_m = 0.0;
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
  // Version 2 adds pressure, focusing, and optional observer path length to
  // the serialized source contract.  Recording the schema change prevents an
  // older AMPS campaign reader from silently assigning the shifted columns to
  // the version-1 meanings.
  out << std::setprecision(17)
      << "sep_source_contract_version=2\n"
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

// Attach transport quantities that are not part of the common shock-source
// algebra.  Both adapters call the same production pressure and Parker
// focusing functions; this helper performs only record assembly and finite
// checks, so no physical equation is duplicated in the AMPS-facing layer.
inline ModelStatus attach_source_transport_context(
    const solarwind::PreparedState& solar_wind,
    double radius_m,
    double sin_theta,
    SEPSourceState& out) {
  out.upstream_pressure_Pa=solarwind::proton_pressure_Pa(
      solar_wind,out.upstream_density_m3);
  out.focusing_length_m=solarwind::parker_focusing_length_m(
      solar_wind,radius_m,sin_theta);
  if (!std::isfinite(out.upstream_pressure_Pa) ||
      !(out.upstream_pressure_Pa>0.0) ||
      !std::isfinite(out.focusing_length_m) ||
      !(out.focusing_length_m>0.0)) {
    out.status=ModelStatus::make(StatusCode::NonFiniteResult,
                                 "SEP source transport context");
  }
  return out.status;
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
    out.owner_model_identity=step.owner_model_identity;
    out.configuration_digest=step.configuration_digest;
    out.position_m={{radius_m,0.0,0.0}};
    double n=0.0,V=0.0,Br=0.0,Bphi=0.0,Bmag=0.0,divV=0.0;
    const ModelStatus status=model_.evaluate_radii_with_B_div_checked(
        step,&radius_m,&n,&V,&Br,&Bphi,&Bmag,&divV,1);
    out.status=status;
    if (status.failure()) return status;
    out.density_m3=n;
    out.pressure_Pa=solarwind::proton_pressure_Pa(
        step.common.solar_wind,n);
    out.velocity_m_s={{V,0.0,0.0}};
    // In the +X equatorial reference used by the common 1-D/3-D regression,
    // the Parker azimuthal basis is +Y; signed Bphi therefore maps directly to
    // the Y component and preserves the production Parker polarity.
    out.magnetic_T={{Br,Bphi,0.0}};
    out.magnetic_magnitude_T=Bmag;
    out.div_velocity_s_inv=divV;
    out.focusing_length_m=solarwind::parker_focusing_length_m(
        step.common.solar_wind,radius_m,
        step.common.solar_wind.reference_sin_theta);
    if (!std::isfinite(out.pressure_Pa) || !(out.pressure_Pa>0.0) ||
        !std::isfinite(out.focusing_length_m) ||
        !(out.focusing_length_m>0.0)) {
      out.status=ModelStatus::make(StatusCode::NonFiniteResult,
                                   "SEP Interface1D transport context");
    }
    return out.status;
  }

  ModelStatus source_at_shock(const PreparedStep& step,
                              SEPSourceState& out) const {
    // The checked 1-D source path preserves `out` on ownership failure; only a
    // state prepared by this adapter's Model can reach spectrum conversion.
    acceleration::ShockAccelerationState a;
    const ModelStatus ownership=
        model_.shock_acceleration_state_checked(step,a);
    if (!ownership.ok()) return ownership;
    const ModelStatus source_status=
        make_source_state(a,spectrum_,false,false,0,0.0,out);
    if (!source_status.ok()) return source_status;
    return attach_source_transport_context(
        step.common.solar_wind,a.radius_m,
        step.common.solar_wind.reference_sin_theta,out);
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
    out.owner_model_identity=step.owner_model_identity;
    out.configuration_digest=step.configuration_digest;
    out.position_m=position_m;
    const double x=position_m[0],y=position_m[1],z=position_m[2];
    double n=0.0,vx=0.0,vy=0.0,vz=0.0,bx=0.0,by=0.0,bz=0.0,divV=0.0;
    const ModelStatus status=model_.evaluate_cartesian_with_B_div_checked(
        step,&x,&y,&z,&n,&vx,&vy,&vz,&bx,&by,&bz,&divV,1);
    out.status=status;
    if (status.failure()) return status;
    out.density_m3=n;
    out.pressure_Pa=solarwind::proton_pressure_Pa(
        step.common.solar_wind,n);
    out.velocity_m_s={{vx,vy,vz}};
    out.magnetic_T={{bx,by,bz}};
    out.magnetic_magnitude_T=std::hypot(bx,std::hypot(by,bz));
    out.div_velocity_s_inv=divV;
    const double radius_m=std::hypot(x,std::hypot(y,z));
    const std::array<double,3> radial_hat={{x/radius_m,y/radius_m,z/radius_m}};
    const std::array<double,3> solar_axis_hat={{
        step.solar_axis_hat[0],step.solar_axis_hat[1],step.solar_axis_hat[2]}};
    const double sin_theta=solarwind::parker_sin_colatitude(
        solar_axis_hat,radial_hat);
    out.focusing_length_m=solarwind::parker_focusing_length_m(
        step.common.solar_wind,radius_m,sin_theta);
    if (!std::isfinite(out.magnetic_magnitude_T) ||
        !std::isfinite(out.pressure_Pa) || !(out.pressure_Pa>0.0) ||
        !std::isfinite(out.focusing_length_m) ||
        !(out.focusing_length_m>0.0)) {
      out.status=ModelStatus::make(StatusCode::NonFiniteResult,
                                   "SEP Interface3D transport context");
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
    const ModelStatus source_status=make_source_state(
        acceleration_state,spectrum_,false,false,source_id,patch_area_m2,out);
    if (!source_status.ok()) return source_status;

    // The production acceleration record supplies the physical surface point.
    // Derive only its normalized direction here, then delegate pressure and
    // focusing to the common solar-wind core used by direct SWCME queries.
    const double radius_m=std::hypot(
        out.position_m[0],std::hypot(out.position_m[1],out.position_m[2]));
    if (!std::isfinite(radius_m) || !(radius_m>0.0)) {
      out.status=ModelStatus::make(StatusCode::GeometryFailure,
                                   "SEP source radius");
      return out.status;
    }
    const std::array<double,3> radial_hat={{
        out.position_m[0]/radius_m,out.position_m[1]/radius_m,
        out.position_m[2]/radius_m}};
    const std::array<double,3> solar_axis_hat={{
        step.solar_axis_hat[0],step.solar_axis_hat[1],step.solar_axis_hat[2]}};
    return attach_source_transport_context(
        step.common.solar_wind,radius_m,
        solarwind::parker_sin_colatitude(solar_axis_hat,radial_hat),out);
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

    // Reuse the public directional adapter so pressure, focusing, spectrum,
    // and shock scalars have one mapping path.  Connectivity contributes only
    // its observer-specific flags and the already computed Parker arc length.
    const ModelStatus source_status=source_at_direction(step,u,out);
    if (!source_status.ok()) return source_status;
    out.connection_evaluated=true;
    out.connected=true;
    out.field_line_path_length_m=root.path_length_m;
    if (!std::isfinite(out.field_line_path_length_m) ||
        out.field_line_path_length_m<0.0) {
      out.status=ModelStatus::make(StatusCode::NonFiniteResult,
                                   "SEP observer field-line path length");
    }
    return out.status;
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
