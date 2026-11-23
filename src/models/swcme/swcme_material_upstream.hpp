#pragma once

#include "swcme_constants.hpp"
#include "swcme_kinematics.hpp"
#include "swcme_shock.hpp"
#include "swcme_status.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace swcme {
namespace material1d {

// Versioned prescribed SHOCK_ONLY ambient. Positions and vectors are HCI/SI,
// centered on the Sun. A radial solution is exposed in Cartesian components
// so a curved field line can use the actual strain instead of inferring it
// from a radial derivative or a density history. This is a source-free mass
// and ideal-induction profile with a declared external momentum force.
inline const char* profile_key() { return "swcme1d-adiabatic-outer-wind-v1"; }
using Vector = std::array<double,3>;
using Tensor = std::array<double,9>;  // Entry (i,j) is d(vector_i)/d(x_j).

struct Config {
  // Zero/default quantities require explicit ownership and validation. The
  // 0.30–1.20 AU test deck is not a default for another event.
  double inner_radius_m = 0.0, outer_radius_m = 0.0, reference_radius_m = 0.0;
  double wind_m_s = 0.0, reference_electron_density_m3 = 0.0;
  double reference_electron_temperature_k = 0.0, reference_proton_temperature_k = 0.0;
  double reference_field_t = 0.0, reference_sin_colatitude = -1.0;
  int polarity = 0;
  Vector rotation_axis{{0.0,0.0,0.0}};
  double rotation_rad_s = std::numeric_limits<double>::quiet_NaN();
  double solar_gravity_m3_s2 = 0.0, gamma = 0.0;
};

struct Sample {
  Vector position_m{}, velocity_m_s{}, magnetic_field_t{};
  Tensor velocity_gradient_s_inv{}, magnetic_gradient_t_m_inv{};
  Vector pressure_gradient_pa_m_inv{}, mass_density_gradient_kg_m4{};
  Vector electron_temperature_gradient_k_m{}, proton_temperature_gradient_k_m{};
  double electron_density_m3 = 0.0, proton_density_m3 = 0.0;
  double electron_temperature_k = 0.0, proton_temperature_k = 0.0;
  double mass_density_kg_m3 = 0.0, pressure_pa = 0.0, sound_speed_m_s = 0.0;
  Vector external_force_n_m3{};
  double external_work_w_m3 = 0.0, heat_w_m3 = 0.0;
};

struct Motion {
  Vector position_m{}, velocity_m_s{};
  Tensor deformation{};  // dX_i/d(initial_position_j), dimensionless.
  double jacobian = 0.0;
};

namespace detail {
inline Vector scale(const Vector& v,double a) { return {{a*v[0],a*v[1],a*v[2]}}; }
inline double norm(const Vector& v) { return std::hypot(std::hypot(v[0],v[1]),v[2]); }
inline Vector cross(const Vector& a,const Vector& b) {
  return {{a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}};
}
inline double dot(const Vector& a,const Vector& b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
template<std::size_t N> inline bool finite(const std::array<double,N>& data) {
  for(double x:data) if(!std::isfinite(x)) return false;
  return true;
}
inline Config checked(Config c) {
  for(double positive:{c.inner_radius_m,c.outer_radius_m,c.reference_radius_m,c.wind_m_s,
       c.reference_electron_density_m3,c.reference_electron_temperature_k,
       c.reference_proton_temperature_k,c.reference_field_t,c.solar_gravity_m3_s2}) {
    if(!std::isfinite(positive)||positive<=0.0)
      throw std::invalid_argument("material upstream: explicit finite positive SI anchors are required");
  }
  if(c.inner_radius_m<solarwind::MIN_RADIUS_M||c.outer_radius_m<=c.inner_radius_m||
     c.reference_radius_m<c.inner_radius_m||c.reference_radius_m>c.outer_radius_m||
     !std::isfinite(c.gamma)||c.gamma<=1.0||!std::isfinite(c.reference_sin_colatitude)||
     c.reference_sin_colatitude<0.0||c.reference_sin_colatitude>1.0||
     (c.polarity!=1&&c.polarity!=-1)||!std::isfinite(c.rotation_rad_s)||c.rotation_rad_s<0.0)
    throw std::invalid_argument("material upstream: invalid radial coverage, EOS, polarity, latitude or rotation");
  const double axis_norm=norm(c.rotation_axis);
  if(!std::isfinite(axis_norm)||axis_norm<=0.0)
    throw std::invalid_argument("material upstream: a nondegenerate HCI rotation axis is required");
  for(double& component:c.rotation_axis) component/=axis_norm;
  return c;
}
}  // namespace detail

class Profile {
 public:
  explicit Profile(const Config& config) : config_(detail::checked(config)) {
    const double winding=config_.rotation_rad_s*config_.reference_radius_m*
        config_.reference_sin_colatitude/config_.wind_m_s;
    br_ref_=config_.polarity*config_.reference_field_t/std::hypot(1.0,winding);
    // Verify the physical profiles at BOTH endpoints before exposing this
    // owner. Interior operations retain their independent finite checks.
    Sample endpoint;
    for(double radius:{config_.inner_radius_m,config_.outer_radius_m})
      throw_if_error(sample({{radius,0.0,0.0}},endpoint));
  }
  const Config& config() const { return config_; }

  // Evaluate the prescribed undisturbed ambient. A local RH surface query
  // is separate; this function never returns a downstream plasma volume.
  ModelStatus sample(const Vector& position,Sample& out) const {
    const double r=detail::norm(position);
    const auto coverage=check_radius(r);
    if(!coverage.ok()) return coverage;
    return sample_on_ray(position,r,detail::scale(position,1.0/r),out);
  }

 private:
  friend class Event;
  // A spherical front has an authoritative kinematic radius. Its Cartesian
  // projection can change the recomputed norm by a rounding unit. Surface
  // queries use the checked canonical radius rather than enlarging coverage.
  ModelStatus sample_on_ray(const Vector& position,double r,const Vector& er,Sample& out) const {
    const double ratio=config_.reference_radius_m/r;
    const double density_ratio=ratio*ratio;
    const double thermal_ratio=std::pow(density_ratio,config_.gamma-1.0);
    Sample value;
    value.position_m=position;
    value.velocity_m_s=detail::scale(er,config_.wind_m_s);
    value.electron_density_m3=value.proton_density_m3=
        config_.reference_electron_density_m3*density_ratio;
    value.electron_temperature_k=config_.reference_electron_temperature_k*thermal_ratio;
    value.proton_temperature_k=config_.reference_proton_temperature_k*thermal_ratio;
    // Quasineutral e/p, no alpha component, electron inertia neglected as in
    // the existing SWCME fluid mass convention. Both pressures contribute
    // to the scalar EOS and therefore to the same canonical RH fast speed.
    value.mass_density_kg_m3=constants::PROTON_MASS_KG*value.proton_density_m3;
    value.pressure_pa=constants::BOLTZMANN_J_K*value.proton_density_m3*
        (value.electron_temperature_k+value.proton_temperature_k);
    value.sound_speed_m_s=std::sqrt(config_.gamma*value.pressure_pa/value.mass_density_kg_m3);
    value.pressure_gradient_pa_m_inv=detail::scale(er,-2.0*config_.gamma*value.pressure_pa/r);
    value.mass_density_gradient_kg_m4=detail::scale(er,-2.0*value.mass_density_kg_m3/r);
    value.electron_temperature_gradient_k_m=detail::scale(er,-2.0*(config_.gamma-1.0)*value.electron_temperature_k/r);
    value.proton_temperature_gradient_k_m=detail::scale(er,-2.0*(config_.gamma-1.0)*value.proton_temperature_k/r);

    const double br=br_ref_*density_ratio;
    const double k=config_.rotation_rad_s/config_.wind_m_s;
    const Vector winding=detail::cross(config_.rotation_axis,position);
    // B = Br (e_r - k axis cross x). This produces B_phi proportional to
    // 1/r on each ray, preserves steady ideal induction and has a continuous
    // polar limit without constructing a singular azimuthal unit vector.
    for(std::size_t i=0;i<3;++i) value.magnetic_field_t[i]=br*(er[i]-k*winding[i]);
    for(std::size_t j=0;j<3;++j) {
      Vector basis{};
      basis[j]=1.0;
      const Vector wind_derivative=detail::cross(config_.rotation_axis,basis);
      for(std::size_t i=0;i<3;++i) {
        const double radial_derivative=((i==j?1.0:0.0)-er[i]*er[j])/r;
        value.velocity_gradient_s_inv[3*i+j]=config_.wind_m_s*radial_derivative;
        value.magnetic_gradient_t_m_inv[3*i+j]=br*(radial_derivative-k*wind_derivative[i])-
            (2.0*er[j]/r)*value.magnetic_field_t[i];
      }
    }
    // U is constant on a radial parcel, so the material acceleration is
    // zero. Balance pressure, central gravity and the FULL vector Lorentz
    // force with f_ext [N/m^3]. The magnetic force can be transverse away
    // from the equator. Its magnitude is not a scalar heating term.
    const Tensor& db=value.magnetic_gradient_t_m_inv;
    const Vector curl{{db[7]-db[5],db[2]-db[6],db[3]-db[1]}};
    const Vector lorentz=detail::scale(detail::cross(curl,value.magnetic_field_t),
        1.0/constants::VACUUM_PERMEABILITY_N_A2);
    const double gravity_density=value.mass_density_kg_m3*(config_.solar_gravity_m3_s2/r)/r;
    for(std::size_t i=0;i<3;++i)
      value.external_force_n_m3[i]=value.pressure_gradient_pa_m_inv[i]+gravity_density*er[i]-lorentz[i];
    value.external_work_w_m3=detail::dot(value.velocity_m_s,value.external_force_n_m3);
    // Adiabatic internal energy evolves with p div(U); the separate force
    // work enters total energy, in W/m^3, without changing the parcel entropy.
    value.heat_w_m3=0.0;
    for(double primitive:{value.electron_density_m3,value.electron_temperature_k,
        value.proton_temperature_k,value.mass_density_kg_m3,value.pressure_pa,
        value.sound_speed_m_s,std::abs(br)}) {
      if(!std::isfinite(primitive)||primitive<=0.0)
        return ModelStatus::make(StatusCode::NonFiniteResult,"material upstream primitive");
    }
    if(!detail::finite(value.velocity_m_s)||!detail::finite(value.magnetic_field_t)||
       !detail::finite(value.velocity_gradient_s_inv)||!detail::finite(value.magnetic_gradient_t_m_inv)||
       !detail::finite(value.pressure_gradient_pa_m_inv)||!detail::finite(value.mass_density_gradient_kg_m4)||
       !detail::finite(value.electron_temperature_gradient_k_m)||!detail::finite(value.proton_temperature_gradient_k_m)||
       !detail::finite(value.external_force_n_m3)||!std::isfinite(value.external_work_w_m3))
      return ModelStatus::make(StatusCode::NonFiniteResult,"material upstream derivatives or work");
    out=value;
    return ModelStatus::success();
  }

 public:
  ModelStatus move(const Vector& initial,double elapsed_s,Motion& out) const {
    const double a=detail::norm(initial);
    const auto initial_status=check_radius(a);
    if(!initial_status.ok()) return initial_status;
    if(!std::isfinite(elapsed_s))
      return ModelStatus::make(StatusCode::NonFiniteInput,"material trajectory elapsed time");
    const double r=a+config_.wind_m_s*elapsed_s;
    const auto final_status=check_radius(r);
    if(!final_status.ok()) return final_status;
    // Along this exact radial map, r varies monotonically with elapsed time
    // in both directions. Both endpoints bound the whole path. There is no
    // clipped endpoint or fractional sampling step when coverage is missing.
    const Vector er=detail::scale(initial,1.0/a);
    const double stretch=r/a;
    Motion value;
    value.position_m=elapsed_s==0.0?initial:detail::scale(initial,stretch);
    value.velocity_m_s=detail::scale(er,config_.wind_m_s);
    value.jacobian=stretch*stretch;
    for(std::size_t i=0;i<3;++i) for(std::size_t j=0;j<3;++j) {
      const double radial=er[i]*er[j];
      value.deformation[3*i+j]=radial+stretch*((i==j?1.0:0.0)-radial);
    }
    if(!detail::finite(value.position_m)||!detail::finite(value.velocity_m_s)||
       !detail::finite(value.deformation)||!std::isfinite(value.jacobian)||value.jacobian<=0.0)
      return ModelStatus::make(StatusCode::NonFiniteResult,"material trajectory deformation");
    const auto rounded_status=check_radius(detail::norm(value.position_m));
    if(!rounded_status.ok()) return rounded_status;
    out=value;
    return ModelStatus::success();
  }

  bool covers_radius(double r) const {
    return std::isfinite(r)&&r>=config_.inner_radius_m&&r<=config_.outer_radius_m;
  }

 private:
  ModelStatus check_radius(double radius) const {
    if(!std::isfinite(radius))
      return ModelStatus::make(StatusCode::NonFiniteInput,"material upstream radius");
    if(!covers_radius(radius))
      return ModelStatus::make_value(StatusCode::OutsideModelDomain,"material upstream radius",radius);
    return ModelStatus::success();
  }
  const Config config_;
  double br_ref_;
};

struct EventConfig {
  // Time is elapsed seconds from the event launch, including both ends.
  std::string identity;
  double begin_s = std::numeric_limits<double>::quiet_NaN();
  double end_s = std::numeric_limits<double>::quiet_NaN();
  kinematics::Config front;
};

struct Shock {
  Vector position_m{}, outward_normal{};
  double epoch_s = 0.0;
  Sample ambient;
  shock::JumpResult jump;
};

class Event {
 public:
  Event(const Config& wind,const EventConfig& event)
      : profile_(wind),event_(event),fingerprint_(digest(profile_.config(),event_)) {
    const auto& front=event_.front;
    const bool mode=front.mode==kinematics::Mode::Ballistic||front.mode==kinematics::Mode::DBM||
        front.mode==kinematics::Mode::DataDriven;
    if(event_.identity.empty()||!std::isfinite(event_.begin_s)||!std::isfinite(event_.end_s)||
       event_.begin_s<0.0||event_.end_s<event_.begin_s||!mode||
       front.extrapolation!=kinematics::ExtrapolationPolicy::OutsideTime||
       front.Vsw_m_s!=profile_.config().wind_m_s||kinematics::validate(front)!=kinematics::Status::Ok)
      throw std::invalid_argument("material upstream event: invalid identity, time coverage or canonical kinematics");
    for(double epoch:{event_.begin_s,event_.end_s})
      if(kinematics::evaluate(front,epoch).status!=kinematics::Status::Ok)
        throw std::invalid_argument("material upstream event: kinematics does not cover the declared interval");
  }

  const Profile& profile() const { return profile_; }
  const EventConfig& config() const { return event_; }
  ConfigurationDigest fingerprint() const { return fingerprint_; }

  bool covers_time(double epoch) const {
    return std::isfinite(epoch)&&epoch>=event_.begin_s&&epoch<=event_.end_s;
  }
  ModelStatus sample(const Vector& position,double epoch,Sample& out) const {
    if(!covers_time(epoch)) return outside_time();
    return profile_.sample(position,out);
  }
  ModelStatus move(const Vector& position,double begin,double end,Motion& out) const {
    if(!covers_time(begin)||!covers_time(end)) return outside_time();
    return profile_.move(position,end-begin,out);
  }

  ModelStatus advance_fixed(const Vector& position,double begin,double end,
      double fixed_dt_s,Motion& out) const {
    const auto span=check_fixed_span(begin,end,fixed_dt_s);
    if(!span.ok()) return span;
    // The input duration owns a physical step. The difference of two large
    // represented epochs can round to a different duration even when end is
    // exactly the correctly rounded begin+dt. Preserve dt for motion too.
    return profile_.move(position,fixed_dt_s,out);
  }

  ModelStatus require_step(double begin,double end,double fixed_dt_s,
      const std::vector<Vector>& material_vertices,bool require_front) const {
    const auto span=check_fixed_span(begin,end,fixed_dt_s);
    if(!span.ok()) return span;
    for(const auto& vertex:material_vertices) {
      Motion candidate;
      const auto status=profile_.move(vertex,fixed_dt_s,candidate);
      if(!status.ok()) return status;
    }
    // Canonical 1-D ballistic, DBM and monotone PCHIP front histories do not
    // turn back through the interval. Bound front coverage by their endpoints
    // if the caller needs shock queries for this step. Insertions and release
    // positions remain separate caller-owned queries through checked APIs.
    if(require_front) for(double epoch:{begin,end}) {
      const auto state=kinematics::evaluate(event_.front,epoch);
      if(state.status!=kinematics::Status::Ok||!profile_.covers_radius(state.radius_m))
        return ModelStatus::make(StatusCode::OutsideModelDomain,"material upstream front coverage");
    }
    return ModelStatus::success();
  }

  ModelStatus shock_at(const Vector& ray,double epoch,Shock& out) const {
    if(!covers_time(epoch)) return outside_time();
    const double ray_norm=detail::norm(ray);
    if(!std::isfinite(ray_norm)||ray_norm<=0.0)
      return ModelStatus::make(StatusCode::DegenerateVector,"material upstream shock ray");
    const auto state=kinematics::evaluate(event_.front,epoch);
    if(state.status!=kinematics::Status::Ok||!profile_.covers_radius(state.radius_m))
      return ModelStatus::make(StatusCode::OutsideModelDomain,"material upstream front coverage");
    Shock candidate;
    candidate.epoch_s=epoch;
    for(std::size_t i=0;i<3;++i) candidate.outward_normal[i]=ray[i]/ray_norm;
    candidate.position_m=detail::scale(candidate.outward_normal,state.radius_m);
    const auto sampled=profile_.sample_on_ray(candidate.position_m,state.radius_m,
        candidate.outward_normal,candidate.ambient);
    if(!sampled.ok()) return sampled;
    shock::PrimitiveState upstream;
    upstream.rho_kg_m3=candidate.ambient.mass_density_kg_m3;
    upstream.pressure_Pa=candidate.ambient.pressure_pa;
    upstream.velocity_m_s=candidate.ambient.velocity_m_s;
    upstream.magnetic_T=candidate.ambient.magnetic_field_t;
    // The prescribed front and normal are inertial. The existing solver
    // transforms into the shock-normal relative flow, classifies a genuine
    // fast front, and supplies one-sided bulk states in HCI. Species heating
    // downstream is unavailable. No new source or sheath law is inferred.
    candidate.jump=shock::solve_ideal_mhd_fast_shock(upstream,candidate.outward_normal,
        state.speed_m_s,profile_.config().gamma);
    if(!candidate.jump.solver_converged)
      return ModelStatus::make(StatusCode::ShockSolverFailure,shock::solve_status_name(candidate.jump.status));
    out=candidate;
    return ModelStatus::success();
  }

 private:
  ModelStatus check_fixed_span(double begin,double end,double dt) const {
    if(!covers_time(begin)||!covers_time(end)) return outside_time();
    if(!std::isfinite(dt)||dt<=0.0||end<=begin||end!=begin+dt)
      return ModelStatus::make(StatusCode::InvalidNumericalStep,"complete fixed AMPS step is required");
    return ModelStatus::success();
  }
  static ModelStatus outside_time() {
    return ModelStatus::make(StatusCode::OutsideModelDomain,"material upstream event time coverage");
  }
  static ConfigurationDigest digest(const Config& c,const EventConfig& e) {
    ConfigurationDigestBuilder d;
    d.add_string(profile_key());
    // The typed event API owns a length-delimited string. Hash every byte
    // rather than treating an embedded zero as the end of the event identity.
    d.add_uint64(e.identity.size());
    for(char byte:e.identity) d.add_byte(static_cast<unsigned char>(byte));
    for(double x:{c.inner_radius_m,c.outer_radius_m,c.reference_radius_m,c.wind_m_s,
        c.reference_electron_density_m3,c.reference_electron_temperature_k,c.reference_proton_temperature_k,
        c.reference_field_t,c.reference_sin_colatitude,static_cast<double>(c.polarity),c.rotation_rad_s,
        c.solar_gravity_m3_s2,c.gamma,e.begin_s,e.end_s,e.front.r0_m,e.front.V0_m_s,e.front.Vsw_m_s,e.front.Gamma_m_inv})
      d.add_double(x);
    for(double x:c.rotation_axis) d.add_double(x);
    d.add_uint64(static_cast<unsigned>(e.front.mode));
    d.add_uint64(static_cast<unsigned>(e.front.extrapolation));
    for(const auto* knots:{&e.front.data_time_s,&e.front.data_radius_m}) {
      d.add_uint64(knots->size());
      for(double x:*knots) d.add_double(x);
    }
    return d.value();
  }
  const Profile profile_;
  const EventConfig event_;
  const ConfigurationDigest fingerprint_;
};

}  // namespace material1d
}  // namespace swcme
