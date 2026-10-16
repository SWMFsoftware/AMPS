#pragma once
// C++17 declarations for specification v1.2 Sections 16.12 and 17.
// This is an interface contract, not an implemented AMPS physics library.
#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace amps_scattering {
using Vec3 = std::array<double,3>;
using Id = std::uint64_t;
enum class OperatorConvention { Standard, RuffoloPhi, EpremHalfD };
enum class MuReference { AlongB, OutwardField };
enum class FrameKind { HostInertial, Plasma, Wave, Corotating, Shock, SourceMixed };
enum class ClockKind { HostTime, PlasmaTime, WaveTime, SourceMixedTime };
enum class Side { Minus, Plus };
enum class Kernel { EulerMaruyama, Milstein, RandomRotation, HemisphericRotation,
                    FiniteVolumeChain, TensorIto, ElasticWaveSplit, PoissonRedraw };
enum class TimeDirection { Forward, Backward };
enum class Admissibility { Admissible, AdmissibleOnDomain, InadmissibleAsPrinted, Unestablished };
enum class Provides { ScalarDmumu, MuMomentumTensor, Spectrum, SpatialClosure,
                      TransportSet, GrowthRate, Kernel, Estimator, ReferenceRecord };
enum class RuntimeState { ReadyExplicitInputs, ReadyPublishedVariant, RequiresExternalInput,
                          RequiresUserDecision, RequiresSourceOrCodeAudit, ReferenceDataOnly,
                          BlockedDependency };
enum class Error { Ok, InvalidInput, InvalidField, OutsideDomain, MissingCoverage,
                   UnresolvedSource, Inadmissible, IncompatibleKernel, Nonconverged,
                   Nonprogress, UnsupportedFrameMapping, VersionMismatch };
enum class ZeroClass { Regular, DriftJump, UnboundedDrift, InternalBarrier };
enum class TailPolicy { Unavailable, DeclaredAnalyticTail, PhysicalTruncation };
enum class EnergyKind { MagneticOnly, TotalAlfvenEquipartition, SourceDefined };
enum class LambdaStatus { Finite, Divergent, Unavailable, Nonconverged };
enum class SpatialMeasure { Volume, FluxTube };
enum class MomentumMeasure { GyrotropicP2, FlatDpDmu, SourceDefined };

struct Diagnostic { Error code=Error::Ok; std::string message;
                    std::vector<std::string> blocking_items; };
template<class T> struct Result { std::optional<T> value; Diagnostic diagnostic; };
struct Species { Id id; double mass_kg, charge_C; };
// Authoritative p and mu in the named local frame. v and gamma are derived once.
struct GyroState { Species species; double p_kg_m_s, mu; Id frame_id; };
struct FourMomentum { Species species; double total_energy_J; Vec3 p_kg_m_s; Id frame_id; };
struct Kinematics { double speed_m_s, gamma, total_energy_J, kinetic_energy_J; };
struct FrameContract { FrameKind kind; Id frame_id; MuReference mu_reference;
                       ClockKind clock; int retained_order_u_over_c; };
struct SampledMeasure { SpatialMeasure spatial; MomentumMeasure momentum;
                       std::string jacobian_equation; };
struct Provenance { std::string model_id, source_key, equation, source_revision;
                    std::vector<std::string> register_items; bool derived=false; };
struct DomainContract { double mu_min,mu_max; bool includes_min,includes_max;
                        ZeroClass zero; std::string predicate_id;
                        std::string boundary_extension; };
struct SnapshotKey { Id background, turbulence, topology; double time_s; };
struct SpectrumView;
struct LocalContext { double B_T; Vec3 b_hat; double inv_focusing_length_per_m;
                      Vec3 plasma_velocity_m_s; double VA_m_s;
                      SnapshotKey snapshot; const SpectrumView* spectrum=nullptr; };
struct Readiness { RuntimeState state; Admissibility admissibility;
                   bool inputs_ok, domain_ok, capability_ok, pair_certified;
                   Diagnostic diagnostic; };
struct ScalarValue { double D_per_s;
                     std::optional<double> derivative_minus_per_s, derivative_plus_per_s;
                     double quadrature_absolute_error_per_s; };
struct TensorValue {
  double D_mumu_per_s, D_mup_kg_m_s2, D_pp_kg2_m2_s3;
  // All six derivatives: order D_mumu, D_mup, D_pp. Coordinate units follow SI.
  std::array<double,3> d_dmu, d_dp;
  double quadrature_error_scaled; bool derivatives_available;
};
struct LambdaValue { LambdaStatus status; double lambda_parallel_unfocused_m;
                     double absolute_error_m; Provenance provenance; };

class ScalarCoefficient {
 public:
  virtual ~ScalarCoefficient() = default;
  virtual Provenance provenance() const = 0;
  virtual FrameContract frame() const = 0;
  virtual DomainContract domain() const = 0;
  virtual Readiness readiness(Kernel, const GyroState&, const LocalContext&) const = 0;
  virtual Result<ScalarValue> evaluate(const GyroState&, const LocalContext&) const = 0;
};
class TensorCoefficient {
 public:
  virtual ~TensorCoefficient() = default;
  virtual Provenance provenance() const = 0;
  virtual FrameContract frame() const = 0;
  virtual DomainContract domain() const = 0;
  virtual Readiness readiness(Kernel, const GyroState&, const LocalContext&) const = 0;
  virtual Result<TensorValue> evaluate(const GyroState&, const LocalContext&) const = 0;
};
// The caller cannot pass a TensorCoefficient to a scalar kernel by accident.
Result<Kinematics> kinematics(const GyroState&);
Result<LambdaValue> lambda_parallel(const ScalarCoefficient&, const GyroState&, const LocalContext&);

struct WavePopulation { Id id; int direction_along_B; int polarization_source;
                        std::string polarization_convention; Id frame_id;
                        double phase_speed_along_B_m_s; };
struct SpectrumMetadata { std::string canonical_quantity, conversion_id, source_quantity;
                          std::string input_units,input_sidedness,input_component_content;
                          double k_min_per_m,k_max_per_m, transverse_variance_T2;
                          TailPolicy tail; EnergyKind energy_kind; SnapshotKey snapshot; };
// Positive |k| inside each explicitly identified wave population.
struct SpectrumView { SpectrumMetadata metadata; std::vector<WavePopulation> populations;
                      const void* immutable_data=nullptr; };
Result<double> spectral_power_T2_m(const SpectrumView&,Id population,double k_abs_per_m);
// Per real particle in the plasma frame; host multiplies by the statistical weight.
// Binned entries sum to the particle change; waves receive its exact negative.
struct Exchange { Id event_id,population_id,spectral_bin_id;
                  SnapshotKey snapshot;
                  double delta_particle_energy_J,delta_particle_parallel_momentum_kg_m_s; };
struct RandomState { Id particle_id,stream_id,event_counter;
                     // Independent Wiener drivers parameterized in physical host time.
                     std::vector<double> pending_bridge_dt_host_s,pending_bridge_dW_host_sqrt_s; };
struct StepRequest { double dt_host_s; TimeDirection time_direction;
                     Kernel kernel; RandomState random; };
struct StepResult { GyroState state; RandomState random;
                    double dt_host_used_s,dt_coefficient_used_s;
                    double log_weight_increment; unsigned reflections;
                    std::vector<Exchange> exchange; Diagnostic diagnostic; };
// Pure scattering. Focusing, streaming, adiabatic work and shock handling are
// supplied by the audited host composer, not silently added by these functions.
Result<StepResult> scatter_scalar(const ScalarCoefficient&,const GyroState&,
                                 const LocalContext&,const StepRequest&);
Result<StepResult> scatter_tensor(const TensorCoefficient&,const GyroState&,
                                 const LocalContext&,const StepRequest&);
struct StepLimit { double max_dt_host_s; std::string limiting_process;
                   Diagnostic diagnostic; };
Result<StepLimit> scattering_limit(const ScalarCoefficient&,const GyroState&,
                                  const LocalContext&,Kernel);
// Rows of Q are origin states: Q_ij=T_ij / stationary_mass_i.
struct FvGenerator { std::vector<double> edges,centres,stationary_mass,T_face;
                     std::vector<double> rate_left,rate_right,rate_diagonal;
                     double quadrature_error; std::string internal_extension; };
Result<FvGenerator> fv_reference(const ScalarCoefficient&,const GyroState&,
                                const LocalContext&,const std::vector<double>& edges);
Result<GyroState> boost_along_field(const GyroState&,Id target_frame,double u_along_B_m_s);
// Exact nonparallel boosts require a full momentum vector, including gyrophase.
Result<FourMomentum> boost_vector(const FourMomentum&,Id target_frame,Vec3 frame_velocity_m_s);
struct RunRecord { std::string specification_digest, data_digest, code_revision;
                   SampledMeasure measure; FrameContract frame; Kernel kernel;
                   TimeDirection time_direction; SnapshotKey snapshot;
                   std::vector<Provenance> components;
                   std::vector<std::string> decisions, released_pairs, observables;
                   std::string error_budget_id, rng_contract, process_ownership; };
} // namespace amps_scattering
