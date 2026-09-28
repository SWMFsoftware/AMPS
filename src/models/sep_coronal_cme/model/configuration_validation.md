## 14. Proposed configuration contract

Outside the separately owned raw `[swcme]` section, every dimensional value is
a bare number in the fixed unit encoded by the key suffix. SI suffixes are the
default (`_m`, `_m_per_s`, `_tesla`, `_rad`, `_j`); an explicitly named
domain unit such as `_ev` is legal only where the schema defines it. For
example, the production parser accepts `interface_radius_m`,
`initial_rate_m_per_s`, and `lateral_axis_tilt_rad`; values such as `2.5 Rs`, `500 km/s`,
`3 deg`, and `5/3` are invalid.
The blocks below specify the schema shape, not a runnable deck: `REQUIRED`,
`REQUIRED_OR_ZERO`, and alternatives separated by `|` are documentation
metasyntax and must be replaced by one legal literal value.

The raw `[swcme]` block belongs only to the existing SWCME authority and is
absent when `background.authority=analytic-coronal-composite` and
`shock.authority=standalone-mhd-ellipsoid`; mixing it into this stand-alone
profile is a configuration error.

Schemas 1--4 retain their exact parser behavior, enum meanings, and resolved
fingerprints. Schema-5 fields are rejected in older schemas. Parsing is
two-pass: first lex sections and assignments and obtain `run.schema_version`;
then invoke the version-specific section, key, inactive-sentinel, and
required-field tables. Human-authored schema-5 decks supply every physical
choice explicitly. Inactive numeric branches use zero and inactive string/file
branches use `none`, unless the selector forbids the section entirely.

```ini
[run]
schema_version = 5
intent = production-shock-injection | analytic-verification
transport = ballistic-verification | parker | focused-pitch-angle-diffusion | focused-discrete-scattering
transport_frame = inertial | rigid-corotating
start_time_s = REQUIRED
end_time_s = REQUIRED
maximum_steps = REQUIRED_OR_ZERO
campaign_seed_u64 = REQUIRED
random_stream_layout = keyed-v1

[domain]
solar_radius_m = REQUIRED
qualified_source_inner_radius_m = REQUIRED
outer_radius_m = REQUIRED
inner_boundary = solar-sphere-absorb
outer_boundary = escape
outer_boundary_geometry = sun-centered-sphere
coordinate_frame = REQUIRED
epoch_utc = REQUIRED

[mesh]
model = amps-cartesian-amr
domain_min_x_m = REQUIRED
domain_min_y_m = REQUIRED
domain_min_z_m = REQUIRED
domain_max_x_m = REQUIRED
domain_max_y_m = REQUIRED
domain_max_z_m = REQUIRED
global_target_cell_size_m = REQUIRED
solar_surface_target_cell_size_m = REQUIRED
solar_refinement_decay_length_m = REQUIRED
tube_centerline_source = none | field-line-requests
tube_target_cell_size_m = REQUIRED_OR_ZERO
tube_transverse_refinement_decay_length_m = REQUIRED_OR_ZERO
maximum_refinement_level = REQUIRED

[particle_numerics]
time_step_model = fixed-upper-bound | adaptive-local
maximum_time_step_s = REQUIRED
spatial_accuracy_factor = REQUIRED
pitch_angle_accuracy_factor = REQUIRED_OR_ZERO
pitch_endpoint_regularization = REQUIRED_OR_ZERO
shock_motion_accuracy_factor = REQUIRED
species_binding = stable-id-plus-compiled-slot-verified
species_numerics = explicit-all-compiled

[species.ID]
compiled_slot = REQUIRED
chemical_symbol = REQUIRED
expected_mass_kg = REQUIRED
expected_charge_c = REQUIRED
mass_number = REQUIRED_OR_ZERO
transport_role = charged-sep | initialization-only
base_macroparticle_weight = REQUIRED
minimum_population = REQUIRED_OR_ZERO
target_population = REQUIRED_OR_ZERO
maximum_population = REQUIRED_OR_ZERO

[population_control]
model = disabled | amps-conservative-split-merge
check_cadence_steps = REQUIRED_OR_ZERO
minimum_total_particles = REQUIRED_OR_ZERO
target_total_particles = REQUIRED_OR_ZERO
maximum_total_particles = REQUIRED_OR_ZERO
minimum_particles_per_active_block = REQUIRED_OR_ZERO
maximum_particles_per_active_block = REQUIRED_OR_ZERO
maximum_weight_ratio_per_merge_group = REQUIRED_OR_ZERO
conservation = species-weight-momentum-relativistic-energy

[active_corridor]
model = full-domain | swept-field-line-tube
line_ids = none | REQUIRED
distance_metric = none | minimum-euclidean-to-polyline
width_profile = none | constant | tabulated-by-line-coordinate
constant_half_width_m = REQUIRED_OR_ZERO
width_profile_file = none | REQUIRED
geometric_buffer_m = REQUIRED_OR_ZERO
block_rule = full-domain | intersects-swept-tube-bounds
connectivity_closure = none | sun-to-source-to-observer-face-connected
amr_closure = none | ancestors-children-and-ghost-neighbors
lateral_particle_policy = not-applicable | escape

[background]
authority = analytic-coronal-composite

[solar_rotation]
model = rigid | latitude-dependent-verification
rigid_input_rotation_rate_rad_per_s = REQUIRED_OR_ZERO
input_rate_convention = sidereal | synodic
synodic_conversion_ephemeris_file = none | REQUIRED
differential_rotation_coefficients_file = none | REQUIRED

[pfss]
lower_boundary = analytic-harmonics | magnetogram-coefficients
outer_boundary_radius_m = REQUIRED
maximum_degree = REQUIRED
harmonic_normalization = orthonormal-complex
coefficients_file = REQUIRED
flux_balance = reject | remove-monopole
spectral_filter = none | heat-kernel
apodization_degree = REQUIRED_OR_ZERO
input_grid = coefficient-file | sine-latitude | uniform-latitude
quadrature = coefficient-file | cell-area | gauss-legendre
cache_radial_points = REQUIRED
cache_colatitude_points = REQUIRED
cache_longitude_points = REQUIRED
cache_interpolation_order = REQUIRED

[open_flux_calibration]
model = compare-only | magnetogram-scale | not-applicable-verification
reference_definition = rotation-mean-unsigned-radial-field | none
construction_reference_asset_file = none | REQUIRED
construction_data_use_role = none | construction
construction_folding_correction_asset_file = none | REQUIRED
construction_reference_radius_handling = none | sample-position | pre-normalized-to-comparison-radius
construction_comparison_radius_m = REQUIRED_OR_ZERO
qualification_reference_asset_file = none | REQUIRED
qualification_data_use_role = none | qualification
qualification_folding_correction_asset_file = none | REQUIRED
qualification_reference_radius_handling = none | sample-position | pre-normalized-to-comparison-radius
qualification_comparison_radius_m = REQUIRED_OR_ZERO
nested_sphere_radii_m = REQUIRED | none
maximum_signed_to_unsigned_flux_fraction = REQUIRED_OR_ZERO
maximum_nested_unsigned_flux_spread_fraction = REQUIRED_OR_ZERO
maximum_angular_quadrature_relative_error = REQUIRED_OR_ZERO
maximum_qualification_relative_mismatch = REQUIRED_OR_ZERO
scale_factor_minimum = REQUIRED_OR_ZERO
scale_factor_maximum = REQUIRED_OR_ZERO
apply_stage = photospheric-coefficients-before-pfss | none
radius_ensemble_file = none | REQUIRED
radius_ensemble_member_id = none | REQUIRED
radius_ensemble_selection_role = none | preregistered-topology-and-open-flux
radius_ensemble_prior_weight = REQUIRED_OR_ZERO

[formation_height_validation]
model = none | preregistered-cartesian-candidate-product
candidate_product_asset_file = none | REQUIRED
product_axes = none | magnetogram-scale-radii-wind-front
event_formation_constraint_asset_file = none | REQUIRED
event_formation_constraint_representation = none | frequency-time-likelihood | preinferred-height-time-likelihood
event_formation_data_use_role = none | qualification | withheld-validation
density_conditioning = none | per-wind-density-member
joint_comparison = none | d6-topology-coronal-holes-d1-d2-typeii-euv
candidate_weighting = none | preregistered-joint-likelihood
sep_output_selection = forbidden

[current_sheet]
model = finite-shell-schatten | none
interface_coupling = direct-pfss-scs | overlap-minimized | not-applicable
interface_radius_m = REQUIRED_OR_ZERO
outer_radial_radius_m = REQUIRED_OR_ZERO
sector_mapping = field-line-traced | none
maximum_degree = REQUIRED_OR_ZERO
unsigned_boundary_fit = nonnegative-constrained | none
fit_quadrature = cell-area | gauss-legendre | none
fit_angular_oversampling = REQUIRED_OR_ZERO
maximum_normal_field_error_tesla = REQUIRED_OR_ZERO
maximum_unsigned_flux_relative_error = REQUIRED_OR_ZERO
maximum_negative_area_fraction = 0
radialization_gate = outer-zonal-power-and-latitude-flatness | diagnostic-only | not-applicable
maximum_outer_zonal_nonmonopole_power_fraction = REQUIRED_OR_ZERO
latitude_diagnostic_radii_m = REQUIRED | none
latitude_minimum_unmasked_longitude_fraction = REQUIRED_OR_ZERO
maximum_unsigned_radial_flux_rms_fraction = REQUIRED_OR_ZERO
maximum_unsigned_radial_flux_p95_to_p05_ratio = REQUIRED_OR_ZERO
interface_particle_rule = not-applicable | reject-unverified | resolved-transition
resolved_transition_width_m = REQUIRED_OR_ZERO
resolved_transition_profile = none | quintic-vector-potential
transition_vector_potential = none | signed-mie
transition_vector_potential_gauge = none | zero-mean-mie-v1
transition_hcs_policy = not-applicable | exclude-clearance
transition_hcs_clearance_m = REQUIRED_OR_ZERO
maximum_transition_normal_trace_jump_fraction = REQUIRED_OR_ZERO
transition_minimum_field_tesla = REQUIRED_OR_ZERO
transition_jump_support_fraction = REQUIRED_OR_ZERO
transition_clearance_convergence_asset_file = none | REQUIRED
maximum_verified_kink_angle_rad = REQUIRED_OR_ZERO
weak_field_policy = not-applicable | diagnostic | exclude-source
weak_field_reference_tesla = REQUIRED_OR_ZERO
weak_field_relative_threshold = REQUIRED_OR_ZERO

[current_sheet_transport]
model = not-applicable | sector-confined | ideal-coordinate-crossing
pure_hcs_minimum_clearance_m = REQUIRED_OR_ZERO

[closed_field_plasma]
model = isothermal-hydrostatic | polytropic-hydrostatic
reference_radius_m = REQUIRED
base_state = species-density-temperature
base_proton_number_density_m3 = REQUIRED
base_proton_temperature_k = REQUIRED
base_electron_temperature_k = REQUIRED
base_alpha_temperature_k = REQUIRED_OR_ZERO
closed_polytropic_index = REQUIRED_OR_ZERO
base_normalization = prescribed | global-separatrix-pressure-scale-verification | footpoint-separatrix-constrained
base_normalization_asset_file = none | REQUIRED
composition_source = solar-wind
velocity_frame = corotating
force_model = gravity-plus-centrifugal | gravity-only-approximation
maximum_centrifugal_to_gravity_ratio = REQUIRED

[open_closed_interface]
representation = sharp-one-sided | finite-width-volume
policy = diagnostic-kinematic | bounded-approximation | stationary-td-equilibrium
state_origin = analytic-composite | equilibrium-solver | imported-equilibrium
transition_width_m = REQUIRED_OR_ZERO
thickness_ensemble_m = none | REQUIRED
interface_velocity_model = corotating-topology-surface | versioned-asset
interface_velocity_asset_file = none | REQUIRED
uncertainty_asset_file = none | REQUIRED
maximum_absolute_traction_jump_pa = REQUIRED_OR_ZERO
maximum_relative_traction_jump = REQUIRED_OR_ZERO
traction_jump_quantile = REQUIRED_OR_ZERO
maximum_absolute_mass_flux_jump_kg_per_m2_per_s = REQUIRED_OR_ZERO
maximum_relative_mass_flux_jump = REQUIRED_OR_ZERO
mass_flux_jump_quantile = REQUIRED_OR_ZERO
maximum_absolute_volume_momentum_residual_n_per_m3 = REQUIRED_OR_ZERO
maximum_relative_volume_momentum_residual = REQUIRED_OR_ZERO
volume_residual_quantile = REQUIRED_OR_ZERO
minimum_convergence_levels = REQUIRED
maximum_residual_refinement_change = REQUIRED
null_neighborhood_field_threshold_tesla = REQUIRED
maximum_absolute_normal_field_tesla = REQUIRED
maximum_relative_normal_field = REQUIRED
maximum_absolute_relative_normal_velocity_m_per_s = REQUIRED
maximum_relative_normal_velocity_fast_mach = REQUIRED

[open_open_interface]
model = none-continuous-single-family | versioned-interface-catalog
representation = not-applicable | sharp-one-sided | finite-width-volume
policy = not-applicable | diagnostic-kinematic | bounded-approximation
interface_catalog_file = none | REQUIRED
uncertainty_asset_file = none | REQUIRED
thickness_ensemble_m = none | REQUIRED
maximum_absolute_traction_jump_pa = REQUIRED_OR_ZERO
maximum_relative_traction_jump = REQUIRED_OR_ZERO
traction_jump_quantile = REQUIRED_OR_ZERO
maximum_absolute_mass_flux_jump_kg_per_m2_per_s = REQUIRED_OR_ZERO
maximum_relative_mass_flux_jump = REQUIRED_OR_ZERO
mass_flux_jump_quantile = REQUIRED_OR_ZERO
maximum_absolute_volume_momentum_residual_n_per_m3 = REQUIRED_OR_ZERO
maximum_relative_volume_momentum_residual = REQUIRED_OR_ZERO
volume_residual_quantile = REQUIRED_OR_ZERO
minimum_convergence_levels = REQUIRED_OR_ZERO
maximum_residual_refinement_change = REQUIRED_OR_ZERO

[plasma_sheet]
model = none | gaussian-tube-contrast
density_contrast = REQUIRED_OR_ZERO
angular_half_width_rad = REQUIRED_OR_ZERO
thermodynamic_rule = none | fixed-temperature
balance_policy = not-applicable | diagnostic-kinematic | bounded-approximation
uncertainty_asset_file = none | REQUIRED
angular_half_width_ensemble_rad = none | REQUIRED
maximum_absolute_volume_momentum_residual_n_per_m3 = REQUIRED_OR_ZERO
maximum_relative_volume_momentum_residual = REQUIRED_OR_ZERO
volume_residual_quantile = REQUIRED_OR_ZERO
minimum_convergence_levels = REQUIRED_OR_ZERO
maximum_residual_refinement_change = REQUIRED_OR_ZERO

[solar_wind]
model = flux-tube-polytropic | empirical-kinematic
scientific_role = analytic-verification | sensitivity | event-nominal
base_reference_radius_m = REQUIRED
energy_closure = isothermal | polytropic | empirical-profile
polytropic_index = REQUIRED_OR_ZERO
temperature_model = uniform | tube-from-target-speed | versioned-profile
uniform_base_temperature_k = REQUIRED_OR_ZERO
target_speed_definition = none | asymptotic | finite-radius
target_speed_relation = none | wsa-versioned
target_speed_coefficients_file = none | REQUIRED
target_speed_radius_m = REQUIRED_OR_ZERO
kinematic_profile_manifest_file = none | REQUIRED
kinematic_profile_manifest_schema = none | sep-kinematic-wind-profile-v1
kinematic_interpolation = none | quintic-hermite-c2-certified
inner_outer_blend_inner_radius_m = REQUIRED_OR_ZERO
inner_outer_blend_outer_radius_m = REQUIRED_OR_ZERO
mass_flux_model = uniform-base-density | mass-per-magnetic-flux | radial-flux-density-with-mapped-field | colocated-density-velocity-field
mass_flux_coefficients_file = none | REQUIRED
mass_flux_data_use_role = none | construction
base_number_density_m3 = REQUIRED_OR_ZERO
mass_per_magnetic_flux_kg_per_s_per_wb = REQUIRED_OR_ZERO
outer_mass_flux_density_kg_per_m2_per_s = REQUIRED_OR_ZERO
mass_flux_reference_radius_m = REQUIRED_OR_ZERO
composition = proton-electron | proton-electron-alpha | versioned-mixture
alpha_to_proton_ratio = REQUIRED
composition_asset_file = none | REQUIRED
composition_uncertainty_model = none | covariance | named-ensemble
composition_covariance_or_ensemble_file = none | REQUIRED
electron_density_conversion = quasineutral-composition-and-charge-state
momentum_residual_policy = transonic-solve | report-and-gate
momentum_residual_acceleration_floor_m_per_s2 = REQUIRED_OR_ZERO
maximum_normalized_momentum_residual = REQUIRED
maximum_absolute_momentum_residual_m_per_s2 = REQUIRED
d7_qualification_asset_file = none | REQUIRED
d7_qualification_data_use_role = none | qualification
d7_withheld_validation_asset_file = none | REQUIRED
d7_withheld_validation_data_use_role = none | withheld-validation
d7_comparison_radii_m = none | REQUIRED
d7_maximum_relative_density_mismatch = REQUIRED_OR_ZERO
d7_maximum_relative_speed_mismatch = REQUIRED_OR_ZERO
d7_maximum_relative_mass_flux_mismatch = REQUIRED_OR_ZERO
d7_maximum_overlap_absolute_log_density_mismatch = REQUIRED_OR_ZERO
d7_maximum_overlap_covariance_normalized_mismatch = REQUIRED_OR_ZERO
required_consumer_coverage_policy = fail-required-event-support | diagnostic-mask-sensitivity
d7_minimum_open_magnetic_flux_coverage_fraction = REQUIRED_OR_ZERO
d7_minimum_open_area_coverage_fraction = REQUIRED_OR_ZERO
d7_minimum_source_flux_coverage_fraction = REQUIRED_OR_ZERO
d7_minimum_observer_exposure_coverage_fraction = REQUIRED_OR_ZERO
d7_minimum_export_support_coverage_fraction = REQUIRED_OR_ZERO

[plasma_eos]
model = ideal-single-fluid
adiabatic_index = REQUIRED
electron_mass_in_density = include | neglect-recorded
near_parallel_threshold_rad = REQUIRED

[source_surface_coupling]
model = rotating-footpoint-conservative
winding_construction = radial-from-scs-boundary
longitude_mapping = full-jacobian | uniform-speed-special-case
minimum_forward_mapping_jacobian = REQUIRED_OR_ZERO
fold_action = fail | mark-invalid-diagnostic
azimuthal_flow = angular-momentum-conserving
coupled_wind_mapping_tolerance = REQUIRED

[turbulence]
open_field_model = prescribed | wkb-outward
closed_field_model = invalid | prescribed-bidirectional | direct-mean-free-path
direction_basis = geometric-outward-inward
reference_radius_m = REQUIRED
wave_energy_at_reference_j_m3 = REQUIRED_OR_ZERO
wave_energy_radial_exponent = REQUIRED_OR_ZERO
outward_cross_helicity = REQUIRED_OR_ZERO
correlation_length_at_reference_m = REQUIRED_OR_ZERO
correlation_length_radial_exponent = REQUIRED_OR_ZERO
spectral_model = kolmogorov | kraichnan | power-law
custom_spectral_index = REQUIRED_OR_ZERO
closed_turbulence_asset_file = none | REQUIRED
maximum_delta_b_over_b = REQUIRED
delta_b_over_b_action = fatal-small-amplitude-closure | diagnostic-empirical
maximum_wave_pressure_to_thermal_pressure_diagnostic = REQUIRED
maximum_wave_pressure_to_total_momentum_pressure_diagnostic = REQUIRED
wave_force_absolute_tolerance_asset_file = REQUIRED
maximum_absolute_wave_acceleration_m_per_s2 = REQUIRED
maximum_wave_force_fraction = REQUIRED
wave_force_quantile = REQUIRED
maximum_wave_force_quantile_fraction = REQUIRED

[transport_coefficients]
parallel_mfp_model = none | single-power-law | smooth-broken-power-law
parker_parallel_diffusion = none | from-mean-free-path
perpendicular_diffusion_model = none | fixed-ratio | versioned-asset
perpendicular_to_parallel_ratio = REQUIRED_OR_ZERO
perpendicular_diffusion_coefficients_file = none | REQUIRED
drift_model = none
drift_coefficients_file = none
focused_collision_model = none | pitch-angle-diffusion | discrete-scattering
pitch_angle_shape = none | isotropic | quasilinear-versioned
pitch_angle_coefficients_file = none | REQUIRED
discrete_scattering_kernel = none | isotropic-poisson
discrete_scattering_kernel_file = none
reference_radius_m = REQUIRED_OR_ZERO
reference_rigidity_v = REQUIRED_OR_ZERO
reference_mfp_m = REQUIRED_OR_ZERO
rigidity_exponent = REQUIRED_OR_ZERO
single_radial_exponent = REQUIRED_OR_ZERO
break_radius_m = REQUIRED_OR_ZERO
inner_radial_exponent = REQUIRED_OR_ZERO
outer_radial_exponent = REQUIRED_OR_ZERO
transition_sharpness = REQUIRED_OR_ZERO
out_of_domain_policy = not-applicable | fail | ballistic-verification-only

[shock]
authority = standalone-mhd-ellipsoid
surface_mode = solar-anchored-dome
geometry_record = shock_geometry
jump_model = ideal-mhd-oblique
activation = fast-magnetosonic
critical_mach_model = edmiston-kennel-1984 | none
critical_mach_table_file = REQUIRED | none
critical_mach_convention = fast | total-alfven | normal-alfven | none
critical_mach_out_of_domain_policy = fail-preflight | diagnostic-inapplicable | exclude-source-budgeted | not-applicable
maximum_criticality_excluded_area_fraction = REQUIRED_OR_ZERO
maximum_criticality_excluded_incident_number_fraction = REQUIRED_OR_ZERO
maximum_criticality_excluded_incident_kinetic_energy_fraction = REQUIRED_OR_ZERO
initial_fast_shock_requirement = none | any-fast-patch | minimum-fast-area-fraction
initial_fast_area_region = whole-clipped-dome | apex-cone
initial_fast_cone_half_angle_rad = REQUIRED_OR_ZERO
minimum_initial_fast_area_fraction = REQUIRED_OR_ZERO
activation_event_time_tolerance_s = REQUIRED
provider_update_cadence_s = REQUIRED

[shock_geometry]
model = triaxial-ellipsoid
axis_convention = radial-principal-axis
initial_parameterization = center-and-principal-axes | apex-and-principal-axes
reference_time_s = REQUIRED
direction_frame = REQUIRED
direction_longitude_rad = REQUIRED
direction_latitude_rad = REQUIRED
lateral_axis_tilt_rad = REQUIRED
orientation_evolution = fixed
evolution = independent-component-laws | tabulated-snapshots
initial_apex_radius_m = REQUIRED_OR_ZERO
snapshots_file = none | REQUIRED
snapshot_interpolation = none | monotone-cubic-c1
snapshot_extrapolation = reject

[shock_geometry.center_distance]
initial_value_m = REQUIRED_OR_ZERO
kinematics = smooth-rate-transition | inactive
initial_rate_m_per_s = REQUIRED_OR_ZERO
final_rate_m_per_s = REQUIRED_OR_ZERO
transition_start_time_s = REQUIRED_OR_ZERO
transition_duration_s = REQUIRED_OR_ZERO
after_transition = constant-final-rate | none

[shock_geometry.radial_semiaxis]
# Same fields as shock_geometry.center_distance.

[shock_geometry.lateral_semiaxis_1]
# Same fields as shock_geometry.center_distance.

[shock_geometry.lateral_semiaxis_2]
# Same fields as shock_geometry.center_distance.

[cme_piston]
enabled = false | true
geometry_source = none | independent-record
geometry_record = none | piston_geometry

[piston_geometry]
model = none | triaxial-ellipsoid-versioned-file
record_file = none | REQUIRED
minimum_front_separation_m = REQUIRED_OR_ZERO

[source]
enabled = false | true
surface_requirement = fast | supercritical
scientific_role = event-nominal | sensitivity | verification
population_semantics = net-first-passage-upstream-released
release_model = empirical-net-first-passage-release
release_boundary = upstream-reference-first-passage | front-conormal-flux-no-through-flow-sensitivity | shock-adjacent-absorbing-verification
shock_reentry_policy = absorb-post-reference-return-loss-ledger | conormal-no-through-flow | absorb-shock-adjacent-verification
reference_surface_distance_m = REQUIRED_OR_ZERO
focused_release_phase_space = outward-first-passage-flux-weighted | unrestricted-isotropic-verification | not-applicable
focused_release_mu_c_minimum_abs_bdotn = REQUIRED_OR_ZERO
focused_release_pitch_root_absolute_tolerance = REQUIRED_OR_ZERO
focused_release_support_root_relative_momentum_tolerance = REQUIRED_OR_ZERO
focused_release_support_root_time_tolerance_s = REQUIRED_OR_ZERO
focused_release_normal_flux_quadrature_absolute_tolerance_m_per_s = REQUIRED_OR_ZERO
focused_release_normal_flux_quadrature_relative_tolerance = REQUIRED_OR_ZERO
no_focused_escape_policy = fail-preflight | exclude-budgeted | not-applicable
maximum_no_escape_number_fraction = REQUIRED_OR_ZERO
maximum_no_escape_energy_fraction = REQUIRED_OR_ZERO
closed_field_policy = diagnose-only
pfss_scs_transition_policy = exclude | convergence-qualified
spectrum_model = local-compression-dsa | fixed-phase-space-power-law
normalization_model = physical-rate | flux-fraction
physical_rate_patch_distribution = none | incoming-species-flux | area | versioned-asset
physical_rate_patch_asset_file = none | REQUIRED
energy_budget_model = shock-frame-kinetic-conversion | versioned-available-flux
energy_budget_asset_file = none | REQUIRED
test_particle_energy_fraction_limit = REQUIRED
maximum_nonthermal_energy_fraction = REQUIRED
placement_kernel = none | upstream-top-hat | upstream-quintic
placement_kernel_support = none | upstream-one-sided
placement_thickness_cell_fraction = REQUIRED_OR_ZERO
maximum_placement_thickness_m = REQUIRED_OR_ZERO
minimum_conormal_diffusivity_m2_per_s = REQUIRED_OR_ZERO
radial_envelope = hard-cutoff | quintic-taper
taper_start_apex_radius_m = REQUIRED_OR_ZERO
zero_source_apex_radius_m = REQUIRED

[consumer_acceptance_budgets]
model = disabled-verification | preregistered-event-grade
transition_consumer_budget_asset_file = none | REQUIRED
front_return_budget_asset_file = none | REQUIRED
birth_energy_edges_file = none | REQUIRED
birth_time_edges_file = none | REQUIRED
cohort_partition = none | stable-species-birth-energy-patch-lineage-time
finite_footprint_measure = none | unsigned-magnetic-flux
maximum_transition_source_area_fraction = REQUIRED_OR_ZERO
maximum_transition_counterfactual_number_rate_fraction = REQUIRED_OR_ZERO
maximum_transition_counterfactual_energy_rate_fraction = REQUIRED_OR_ZERO
maximum_transition_finite_footprint_flux_fraction = REQUIRED_OR_ZERO
maximum_transition_runtime_represented_number_loss_fraction = REQUIRED_OR_ZERO
maximum_transition_runtime_represented_birth_energy_loss_fraction = REQUIRED_OR_ZERO
maximum_front_return_represented_number_loss_fraction = REQUIRED_OR_ZERO
maximum_front_return_represented_birth_energy_loss_fraction = REQUIRED_OR_ZERO
exceedance_action = diagnostic-only-verification | mark-not-event-grade

[source_species.ID]
upstream_abundance_model = background-density | relative-to-species | versioned-asset
relative_to_species_id = none | REQUIRED
relative_abundance = REQUIRED_OR_ZERO
abundance_asset_file = none | REQUIRED
physical_rate_authority = inactive | scalar | versioned-asset
physical_particle_rate_per_s = REQUIRED_OR_ZERO
physical_rate_asset_file = none | REQUIRED
upstream_release_fraction = REQUIRED_OR_ZERO
energy_coordinate = kinetic-per-particle | kinetic-per-nucleon
minimum_energy_ev = REQUIRED
maximum_energy_ev = REQUIRED
fixed_phase_space_power_index = REQUIRED_OR_ZERO
samples_per_step = REQUIRED
momentum_direction = isotropic | pitch-angle-law
source_frame = upstream-plasma | outward-wave | inward-wave
pitch_angle_asset_file = none | REQUIRED

[observer.ID]
position_mode = fixed-cartesian | ephemeris-file
coordinate_frame = REQUIRED
ephemeris_file = none | REQUIRED
velocity_source = coordinate-frame-stationary | ephemeris-file
position_x_m = REQUIRED_OR_ZERO
position_y_m = REQUIRED_OR_ZERO
position_z_m = REQUIRED_OR_ZERO
collection_geometry = volume-sphere | surface-disk | surface-sphere
estimator = volume-residence | surface-crossing
geometry_radius_m = REQUIRED
crossing_sense = not-applicable | inward | outward | both
surface_normal_x = REQUIRED_OR_ZERO
surface_normal_y = REQUIRED_OR_ZERO
surface_normal_z = REQUIRED_OR_ZERO
angular_acceptance = omnidirectional | fixed-look-cone | pitch-angle-bins
look_axis_x = REQUIRED_OR_ZERO
look_axis_y = REQUIRED_OR_ZERO
look_axis_z = REQUIRED_OR_ZERO
look_cone_half_angle_rad = REQUIRED_OR_ZERO
sampling_cadence_s = REQUIRED
accumulation_interval_s = REQUIRED
species = REQUIRED
energy_coordinate = kinetic-per-particle | kinetic-per-nucleon
energy_grid = linear | logarithmic | explicit-edges
minimum_energy_ev = REQUIRED_OR_ZERO
maximum_energy_ev = REQUIRED_OR_ZERO
energy_channel_count = REQUIRED_OR_ZERO
energy_edges_file = none | REQUIRED
pitch_angle_grid = none | linear | explicit-edges
pitch_angle_channel_count = REQUIRED_OR_ZERO
pitch_angle_edges_file = none | REQUIRED
reported_intensity = directional-per-sr | accepted-solid-angle-integrated | omnidirectional
empty_bin_policy = physical-zero-with-validity
instrument_response = none | online-versioned
instrument_response_file = none | REQUIRED

[output]
directory = REQUIRED
initialization_mesh_file = sep3d-initialization-mesh.dat
initialization_data_file = sep3d-initialization-data.dat
time_dependent_cadence_s = REQUIRED
write_tecplot = true
invalid_value_sentinel = -1.7976931348623157e308
require_validity_columns = true

[field_line_export]
enabled = false | true
bundle_path = none | REQUIRED
bundle_schema = 3
sampling = adaptive
geometry_time_model = stationary-in-transport-frame | rigid-rotation-from-trace-time
maximum_segment_length_m = REQUIRED
maximum_relative_field_change = REQUIRED
maximum_relative_plasma_change = REQUIRED
time_sampling = static-background
shock_intersections = false | true
shock_history_start_s = REQUIRED_OR_ZERO
shock_history_end_s = REQUIRED_OR_ZERO
shock_history_cadence_s = REQUIRED_OR_ZERO
write_text_diagnostics = false | true

[field_line.ID]
seed_mode = observer-connected | photospheric-footpoint | cartesian
observer_id = none | REQUIRED
observer_mapping = none | static-volume-intersection | time-resolved-volume-intersection
observer_mapping_cadence_s = REQUIRED_OR_ZERO
observer_mapping_event_time_tolerance_s = REQUIRED_OR_ZERO
observer_mapping_out_of_tolerance = not-applicable | fail
seed_x_m = REQUIRED_OR_ZERO
seed_y_m = REQUIRED_OR_ZERO
seed_z_m = REQUIRED_OR_ZERO
photospheric_longitude_rad = REQUIRED_OR_ZERO
photospheric_latitude_rad = REQUIRED_OR_ZERO
arc_length_orientation = geometric-outward
trace_branches = both-from-seed | outward-from-photosphere
trace_time_s = REQUIRED
outer_radius_m = REQUIRED
connection_tolerance_m = REQUIRED_OR_ZERO
measure = characteristic | flux-tube | quadrature
measure_group_id = none | REQUIRED
represented_magnetic_flux_wb = REQUIRED_OR_ZERO
cross_section_model = none | traced-reference-footprint
cross_section_asset_file = none | REQUIRED
cross_section_boundary_samples = REQUIRED_OR_ZERO
maximum_cross_section_relative_error = REQUIRED_OR_ZERO
```

`sep-kinematic-wind-profile-v1` is a checksummed manifest, not a monolithic
numeric table with one global coordinate declaration. It contains a stable
manifest ID and one or more channel records. Each channel record contains all
of the following fields; omission is a parse error rather than permission to
infer a value from another channel:

| Channel field | Required meaning |
|---|---|
| `channel_id`, `channel_role`, `species_id` | Stable identity; role is `inner-mass-density`, `inner-electron-density`, `outer-velocity`, or `species-temperature`. `species_id` is required only where the quantity is species resolved. |
| `value_file`, `content_checksum`, `value_units` | Immutable numeric payload, its content identity, and the exact SI unit. |
| `velocity_component`, `velocity_reference_frame` | For an outer-velocity channel, the only legal pairs are `radial,inertial`, `field-aligned,inertial`, and `field-aligned,corotating`; both fields are `not-applicable` for density or temperature channels. The pair states exactly which scalar is stored instead of treating “field aligned” as an implicit frame conversion. |
| `minimum_radial_projection` | A channel-local dimensionless guard in `(0,1]` for `radial,inertial`; exactly zero and inactive for both field-aligned pairs and all nonvelocity channels. |
| `consumer_selector_kind`, `consumer_selector_members`, `excluded_line_or_tube_ids` | Deterministic routing by `stable-topology-class` or `stable-line-or-tube-ids`. A topology selector may subtract explicit stable IDs so a broad radial population and a special field-aligned subset can coexist. Expanded selector domains are pairwise disjoint and cover every required consumer exactly once; order, priority, and fallback matching are forbidden. |
| `abscissa`, `abscissa_units` | `heliocentric-radius` in meters or `oriented-field-line-arclength` in meters. This choice belongs to the channel, so a radial density and an arc-length velocity may coexist. |
| `coordinate_frame`, `trace_epoch_utc` | Frame and epoch of the geometry to which the channel applies. |
| `line_ids`, `topology_class`, `background_fingerprint` | Exact line/tube support. Stable line IDs and the complete immutable background fingerprint are mandatory for field-aligned velocity or arc-length abscissa; population-level radial products still declare their topology class. |
| `arc_length_origin` | `not-applicable` for radius; otherwise a named geometric origin and orientation consistent with the exported line. A bare numeric offset is not an identity. |
| `support_segments` | Closed, nonoverlapping support intervals carrying stable segment IDs and topology identity. Gaps and turning points are explicit boundaries. |
| `data_source`, `instrument`, `processing_version`, `source_epoch_utc`, `data_use_role` | Complete provenance and one role: `construction`, `qualification`, or `withheld-validation`. Active background channels must be `construction`; qualification and withheld records cannot be consumed as construction. |
| `uncertainty_model`, `covariance_asset_file`, `covariance_checksum` | `covariance` or `named-ensemble` plus an immutable payload that includes cross-channel correlations used in the overlap test. |
| `node_value`, `node_first_derivative`, `node_second_derivative` | The values and derivatives needed by the declared quintic-Hermite interpolant at every node, in units implied by the value and abscissa. |

`kinematic_interpolation=quintic-hermite-c2-certified` has one exact meaning.
On each support segment, adjacent node values and their common first and second
derivatives define the quintic Hermite polynomial. Preparation evaluates its
value, first derivative, second derivative, and all real interior extrema. The
channel passes only if adjacent pieces are `C2`, every required physical value
is strictly positive, and the polynomial does not overshoot the closed range of
its two endpoint values. A genuine resolved extremum therefore appears as an
explicit node. Certification uses the stored coefficients directly and is
repeated after unit/frame conversion; a generic cubic spline or a library
"shape-preserving" option is not an equivalent implementation.

For a channel tabulated against heliocentric radius on a particular line,
`r(s)` must be strictly monotone on each support segment. A zero derivative,
turning point, or reversal terminates the segment; the manifest must split the
line there and may not select a branch by nearest radius. Arc-length channels
instead use their recorded origin, outward orientation, and trace epoch.
Neither coordinate may be silently reinterpreted as the other.

The schema-5 configuration families and option ownership are:

| Family | Application-facing selection | Shared-model contract | Qualification status |
|---|---|---|---|
| Analytic verification | `run.intent=analytic-verification`; analytic harmonic/PFSS inputs | Uniform/transonic verification wind, diagnostic interfaces, manufactured source or disabled source | Verification only; application owns mesh, time step, output, and AMPS bindings. |
| Stand-alone event candidate | `run.intent=production-shock-injection`; `background.authority=analytic-coronal-composite` | Empirical manifest wind, independent D6/D7 qualification, bounded or qualified equilibrium interface, finite-reference source | Event-nominal only after all Stage-12/13 data and release gates pass. |
| Declared sensitivity | Production-capable executable with component `scientific_role=sensitivity` | Target-speed, diagnostic interface, conormal source, or budgeted exclusions as explicitly selected | Never relabeled event-nominal; all changed authorities are fingerprinted. |
| Field-aligned consumer | `srcSEP` selects `field_line_input.provider=sep-field-line-bundle` | Reads the immutable shared bundle; no local background, source, or observer redefinition | Inherits bundle qualification and adds independent 1-D numerical/parity gates. |
| SWCME alternative backend | Application selects the sibling SWCME authority and its raw `[swcme]` block | No `sep_coronal_cme` background/shock sections may be mixed into that family | Qualified by the SWCME/application contract, not by silently falling back to this model. |

Shared physical options and validation belong to `sep_coronal_cme`; neutral
snapshot and bundle types belong to `sep_common`; mesh, MPI, particle
allocation, output, and backend selection remain application-adapter options.

Schema 5 exposes only capabilities whose governing operator is specified and
whose release tests are named. The following remain explicit roadmap items and
must return typed `NotImplemented` if encountered through a migrated or
programmatic configuration; they never fall back to a nearby enum:

| Reserved capability | First enabling gate |
|---|---|
| Finite-thickness HCS and cross-sector/HCS drift | Stage 11A |
| Common-flux-surface PFSS/SCS transition | future schema plus signed-gauge reconstruction, free-boundary sheet, `CPL3D11--12`, and D9 zero-crossing-flux gate |
| General nonradial winding with `R_w<R_scs` | generalized Piola pushforward, rederived wind/map coupling, `CPL3D10`, and a future schema |
| Physical free-escape-boundary SEP source | normalized diffusion-advection escape flux, pitch-angle/placement contract, and dedicated source tests in a future schema |
| Downstream shock/sheath crossing and transport | Stage 11B |
| Closed-loop shock injection | closed-loop loss-cone, precipitation, source-budget, and two-footpoint tests |
| General shock-frame isotropic/gyrotropic source | full-velocity or qualified gyro-averaged source operator |
| Non-rigid or topology-changing field-line geometry | time-dependent bundle schema and 3-D/1-D parity gate |
| Momentum/species-dependent integrated-Peclet release surfaces | future source-surface discriminator, per-species/per-momentum/per-patch/per-time geometry, offset caps and all topology/clearance gates described below |
| Return-to-release renewal after a front encounter | future source return-operator discriminator, normalized dwell/transition kernel, downstream or sheath authority, and closed number/momentum/energy/shock-work ledgers |
| Foreshock-modified scattering | distinct future empirical-MFP-proxy and self-generated-wave branches, each with the provenance, coverage, and coefficient gates described below |
| Gradient/curvature drift | future `relativistic-guiding-center` transport discriminator with a jointly derived spatial/momentum operator, reduction, weak-field, HCS/separatrix, and interface test suite |
| Deflecting and rotating CME geometry | future vector-center and SO(3)-attitude history discriminators with analytic derivatives, covariance, and swept-geometry tests |
| Optional impulsive/flare-associated source | future additional-source-component discriminator with independent phase-space law, normalization, provenance, and ledgers |
| Continuous wind plausibility envelope | future wind-plausibility discriminator with a versioned quantity/frame/support/uncertainty asset and certified full-support extrema tests |
| Arbitrary versioned discrete-scattering kernel | normalized transition law, detailed-balance/invariant measure where claimed, autocorrelation/MFP relation, endpoint behavior, and dedicated regression tests in a future schema |

The names in the following list are **reserved future-schema discriminants**,
not accepted schema-5 input.  A schema-5 parser must reject every one before
allocation with a typed `NotImplemented` result.  They extend the existing
`[source]`, `[transport_coefficients]`, `[shock_geometry]`, and `[solar_wind]`
authorities; they do not create parallel authorities whose values could
silently disagree with those sections.

- A future `[source]` value
  `release_boundary=upstream-integrated-peclet-first-passage` would extend the
  existing release-boundary discriminator and replace, rather than
  supplement, the schema-5 fixed-distance construction.  For stable species
  \(a\), momentum coordinate \(p\), front patch \(\sigma\), and event time
  \(t\), its upstream normal offset \(L_a(p,\sigma,t)\) must be the
  sign-certified root of

  \[
  \int_0^{L_a(p,\sigma,t)}
    \frac{u^{\rm in}_{1n}(d,\sigma,t)}
         {\kappa_{nn,a}(d,p,\sigma,t)}\,{\rm d}d
    = \mathcal P_{\rm target},
  \]

  using the positive front-frame upstream inflow and coefficients derived from
  the same transport authority as the mover. The Parker branch uses its
  diffusion tensor directly; a focused branch uses only the declared
  diffusion-approximation diagnostic
  \(\kappa_\parallel=v\lambda_\parallel/3\), together with its independently
  declared \(\kappa_\perp\), and never a second surface-only coefficient.
  This fixed-integrated-P discriminator owns exactly one finite positive scalar
  `target_integrated_peclet`, applied uniformly to every supported stratum. A
  species/momentum/patch/time-dependent target would be a different calibrated
  source model with a separately named discriminator, inference asset,
  normalization, and tests; it cannot be smuggled into this branch through a
  table. A paired ordered
  `minimum_reference_surface_distance_m`/`maximum_reference_surface_distance_m`
  bounds the admissible root, and `offset_bound_action` is either fatal or a
  preregistered non-renormalizing sensitivity exclusion;
  it may not also accept `reference_surface_distance_m` as an independent
  physical value.  Every generated surface has a stable identity containing
  species, momentum interval, patch, time interval, coefficient/background
  generation, and root tolerance.  Registered minimum/maximum offsets are
  admissibility bounds, not clipping instructions: reaching either cap is a
  typed failure or a preregistered, non-renormalized sensitivity exclusion.
  The resolved offset must also pass front contact, normal-reach, fold,
  self-overlap, other-surface overlap, transition/HCS clearance, active-mask,
  placement-support, and outer-boundary-clearance tests for the entire time
  interval. Existing transition, HCS, active-domain, and placement clearances
  remain their single authorities; the future source branch cannot copy them
  into inconsistent local fields. The source density consequently becomes
  momentum dependent in its spatial factor; a code path that retains a momentum-independent
  `phi_ref(x,t)` is invalid. The branch remains sensitivity-only until its
  target-depth inference, uncertainty, and cross-event transfer are qualified.
- A future `[source]` value
  `shock_reentry_policy=return-to-release-renewal-kernel` would extend the
  existing re-entry-policy discriminator and replace the schema-5 absorbing return
  policy. Its `return_renewal_kernel_file` is the sole checksummed conditional
  kernel authority and must include terminal absorption
  probability, downstream/sheath state, residence-time distribution, exit
  patch, pitch/gyrophase convention, momentum/energy transition, and the
  kernel's normalization measure.  It must close represented number,
  momentum, energy, and shock-work ledgers and distinguish repeated cycles by
  lineage.  Immediate re-emission with unchanged momentum, or a fresh draw
  from the original source spectrum, is not a physical renewal kernel and is
  forbidden because it double counts acceleration already folded into the
  empirical first-passage source.
- Future foreshock scattering uses a discriminated tuple in the existing
  `[turbulence]` and `[transport_coefficients]` authorities; it never modifies
  the schema-5 mean free path behind the parser. A future
  `parallel_mfp_model=foreshock-distance-proxy` consumes one checksummed
  `foreshock_mfp_proxy_coefficients_file` and is a sensitivity-only positive,
  smooth reduction factor `F_foreshock` with `0<F_foreshock<=1`, applied as
  `lambda_parallel=F_foreshock*lambda_parallel,ambient`. It depends on upstream
  distance, rigidity/species, shock patch, obliquity/Mach state, and time; its
  support begins at or outside the finite reference surface, it must recover
  `F_foreshock=1` at the outer edge, and it may not be labeled self-consistent
  wave growth.
  A future `turbulence.open_field_model=offline-self-generated-wave-asset`
  consumes a checksummed, versioned asset
  containing \(W^\pm(k,\mathbf x,t)\), units, frame, wave-number convention,
  grid/interpolation, coverage mask, resonance mapping, solver/version,
  numerical tolerances, conservation residuals, and immutable background and
  shock fingerprints. A product claimed to replay a self-consistent coupled
  solution additionally binds the source spectrum/normalization,
  species/weights, represented particle distribution, transport and return
  policies, coupling iteration, and coupling cadence; without those bindings
  it is only a prescribed external wave field. It requires the paired future
  `parallel_mfp_model=from-wave-spectrum` and derives \(D_{\mu\mu}\) and/or
  \(\lambda_\parallel\) once; those coefficients cannot also be independently
  prescribed. A truly
  `turbulence.open_field_model=runtime-coupled-self-generated-waves` branch
  requires a separately
  implemented wave-growth/damping/transport solver and coupled conservation
  tests.  In particular, this plan makes **no assumption that `srcSEP`
  currently contains a self-excited-wave solver**; a separately qualified
  solver may produce the offline asset, but its output is never inferred from
  the present application name or directory.
- A future
  `transport_coefficients.drift_model=relativistic-guiding-center` must select
  one jointly derived phase-space operator.  It must define frame, charge
  sign, relativistic magnetic gradient/curvature coefficients, the associated
  energy/momentum evolution, and validity guards for magnetization, field
  smoothness, weak field, transition layers, separatrices, and the HCS.
  Appending a drift velocity and an independent `q E dot v_d` energy term to
  the existing focused mover is forbidden unless the derivation proves that
  the existing plasma-frame momentum operator contains no overlapping work.
  Invalid regions require a typed alternate operator or termination policy,
  never coefficient clipping.
- Future deflection/rotation extends `[shock_geometry]` with two independent
  but synchronized records. A future
  `evolution=vector-trajectory-and-shape-history` replaces the fixed direction
  plus scalar center-distance authority with a Sun-centered vector trajectory;
  `orientation_evolution=versioned-so3-attitude-history` replaces the fixed
  attitude authority. The former owns
  \(\mathbf c(t)\), \(\dot{\mathbf c}(t)\), and any nonradial deflection; the
  latter owns a proper attitude rotation \(R(t)\in SO(3)\) and angular
  velocity. The symbol `Q` remains reserved for the ellipsoid quadratic-form
  tensor of Section 8 and is not itself an orthogonal attitude matrix. A
  changing attitude alone cannot stand in for center deflection.  Interpolation
  must retain the complete derivative: if
  \(\mathbf c=d_c\hat{\mathbf e}_r\), then
  \(\dot{\mathbf c}=\dot d_c\hat{\mathbf e}_r+
  d_c\dot{\hat{\mathbf e}}_r\). Omitting the directional term is invalid. It
  must remain on SO(3), provide consistent analytic derivatives, preserve the
  ellipsoid nesting/solar-anchoring contracts, and carry a checksummed
  event-specific covariance or ensemble; componentwise quaternion splines and
  a universal undocumented angular range are invalid.
- A future `[source]` value
  `release_model=empirical-net-first-passage-plus-impulsive-components` must
  extend the existing source as a
  separately identified component rather than overwrite the shock
  first-passage normalization. The selector composes the existing shock
  planner with a separately versioned impulsive-source provider; it does not
  duplicate that provider's physical fields in `[source]`. The component
  requires a physical release region and
  magnetic-connectivity rule, event-time profile, species-resolved momentum
  and pitch law, frame/Jacobian, number/energy normalization and budget,
  provenance/data-use role, and component-specific ledgers.  Overlap in
  space/time/species with the shock component must be explicitly partitioned
  or modeled jointly so the total source is not double counted.  This is an
  optional attribution/sensitivity capability, not a prerequisite for the
  schema-5 shock-only baseline. Every impulsive birth must also be classified
  relative to the moving front when a front generation exists; a legal earlier
  birth carries a typed no-front-yet state rather than a guessed side. Birth
  on the unmodeled downstream side is
  rejected unless an independently validated downstream provider owns that
  state. A later front encounter must select a typed policy: terminal
  absorb-and-ledger in the upstream-only approximation, or a separately
  validated downstream-transfer/renewal operator. A terminal interaction must
  not overwrite the immutable impulsive source origin. Reacceleration must
  retain impulsive ancestry and shock-work ledgers and may not be counted again in
  the calibrated shock first-passage source.
- A future `[solar_wind]` value
  `plausibility_model=continuous-versioned-envelope` must consume one
  checksummed `plausibility_asset_file` with an explicit qualification or
  withheld-validation data-use role. The single asset contains separately
  checksummed, single-quantity channels and explicit gate records. Each
  resolved gate pairs exactly one selected velocity channel (for example
  inertial radial speed or corotating field-aligned speed) with exactly one
  dynamically compatible acceleration channel. A quasi-steady profile uses
  the field-aligned advective acceleration; a time-dependent profile instead
  uses the full consistently projected material derivative. Each channel owns
  its quantity, units, reference frame, radial/temporal/topological support,
  uncertainty envelope, interpolant, extrema certificate, and tolerance;
  asset-level covariance may bind the speed and acceleration observations.
  A speed bound is never reused as an acceleration bound. Validation uses
  exact extrema only where they are proved for each selected interpolant;
  otherwise certified interval arithmetic or interval branch-and-bound
  encloses every global extremum over the full required support. Passing a
  finite node/radius list or an uncertified optimizer is insufficient.
  The asset may impose event- and topology-specific acceleration/deceleration
  bounds, but the implementation must not invent a universal monotonic-wind
  rule or hidden maximum speed.

`[species.ID]`, `[source_species.ID]`, `[observer.ID]`, and
`[field_line.ID]` are repeatable;
their string IDs are unique and stable across output, restart, and bundle
exchange. For `source.enabled=false`, the `[source]` section contains only that
selector and all `[source_species.ID]` sections are absent. For an enabled
production source, there is exactly one source-species record for every
compiled AMPS species and every species must have `transport_role=charged-sep`;
branch-inactive rate, abundance, spectrum, and angular fields use zero/`none`.

A `traced-reference-footprint` asset is a versioned SI geometry record, not an
image or an area scalar. It contains frame, epoch, an oriented triangulated
reference patch transverse to the local field, ordered boundary connectivity,
interior quadrature nodes/weights, and a content checksum. The named line seed
lies inside its own patch. Validation rejects self-intersection, repeated or
inverted triangles, a field-tangent reference patch, nonpositive signed flux,
insufficient boundary sampling, tube-map folds, and overlap with another
footprint in the same measure group. Refining both triangulation and boundary
traces must reduce flux/volume/overlap residuals below the configured error.

The shared, immutable `ObserverOptions` record is the single physical
authority; the `srcSEP3D` parser and adapter create and bind it but do not own
a second definition. A field-line request with
`seed_mode=observer-connected` must name that record and uses its position as
the seed. A photospheric or Cartesian seed may either set `observer_id=none`
or name an observer solely for an independent physical-overlap mapping; the
observer never moves that seed. During export, the mapper appends the accepted
time-resolved line/volume interval(s), closest coordinate, perpendicular
separation, represented tube volume, tolerance, and validity to that
observer's bundle record. `srcSEP` consumes this mapping and does not ask for
or infer a second observer position. Linear, logarithmic, and explicit-edge
energy branches are mutually exclusive. Energy per nucleon requires the
validated integer `mass_number` of every selected species. `volume-sphere`
requires the residence estimator; either surface geometry requires the crossing
estimator and a complete geometric acceptance. `observer.species` is an
ordered comma-separated list of stable `[species.ID]` suffixes, not chemical
symbols or mutable numerical indices. The global finite invalid sentinel is interpreted only when its
companion validity flag is zero; physical zero and an empty-but-valid bin
remain numeric zero.

Cross-option validation is transactional. In particular:

- all `_time_s` values are elapsed SI seconds from `domain.epoch_utc` (from
  the bundle epoch in `srcSEP`), `start_time_s<end_time_s`, and a nonzero
  `maximum_steps` is only a safety cap. Background, shock, source, ephemeris,
  observer, and imported histories cover the closed run interval. The
  unsigned campaign seed and `keyed-v1` stream layout are fingerprinted and
  restart-identical;
- `run.intent=production-shock-injection` requires `source.enabled=true` and
  the complete source/source-species contract. `source.enabled=false` is legal
  only for `analytic-verification`; a production deck may still be inspected
  with the `--initialization-only` CLI while retaining its enabled source,
  because that mode stops before particle injection/iteration;
- `domain.solar_radius_m < domain.qualified_source_inner_radius_m <
  domain.outer_radius_m`; `solar_wind.base_reference_radius_m` lies from the
  solar radius through the qualified source radius, and all source-supporting
  patches lie at or above the latter. A source-enabled production composite
  additionally requires
  `qualified_source_inner_radius_m<R_i<R_scs<domain.outer_radius_m` and
  `R_b<domain.outer_radius_m`, so every magnetic-authority partition is
  nonempty and the outer Parker region lies inside the physical domain. A
  source-active interval clipped by its configured zero-source radius may be
  empty and is reported as such rather than written with inverted bounds;
- each Cartesian mesh minimum is below its corresponding maximum, the solar
  sphere, the complete exact outer sphere, and every observer collection
  region lie inside those bounds. Global/solar target sizes and decay length
  are positive. Tube values are both zero when their source is `none` and both
  positive otherwise. `maximum_refinement_level` can realize the requested
  finest target.
  The tube centerline is the same stable field-line set used by the selected
  field-line requests; a second geometrically inferred
  Parker spiral is forbidden. Refinement studies vary every active target and
  decay length and must leave the physical domain fixed;
- `finite-shell-schatten` requires `R_sun < R_i <= R_b`, `R_i < R_scs`,
  `sector_mapping=field-line-traced`, positive SCS resolution/tolerances, and
  a current-sheet transport policy other than `not-applicable`; there is no
  required ordering between `R_b` and `R_scs`. `direct-pfss-scs` requires
  `R_i=R_b`, zero transition width, transition profile `none`, and
  analytic-verification intent. Production requires
  `radialization_gate=outer-zonal-power-and-latitude-flatness`, positive
  flatness bounds, `0<=maximum_outer_zonal_nonmonopole_power_fraction<=1`,
  `0<latitude_minimum_unmasked_longitude_fraction<=1`, an RMS bound greater
  than zero, a percentile-ratio bound at least one, and a strictly increasing
  list of unique diagnostic radii from `R_scs` through the exact spherical
  domain. The list grammar is one or more bare SI decimal scalars separated by
  commas; surrounding ASCII whitespace is ignored, empty/duplicate tokens and
  more than 64 entries are rejected, and the fingerprint serializes the
  sorted IEEE-754 values in canonical hexadecimal form. The user must include
  10 and 20 solar radii when each is exterior and in-domain. No radius, including
  `R_scs=2.5 R_sun`, bypasses those outcome gates;
- `radialization_gate=diagnostic-only` is legal only with
  `run.intent=analytic-verification` and finite-shell SCS. It requires the same
  valid exterior radius list and positive mask-coverage fraction, including
  applicable 10/20-solar-radius entries, but sets the zonal, RMS, and
  percentile acceptance thresholds to zero/inactive. All metrics are emitted
  without an outcome comparison. Malformed/nonfinite fields or metrics,
  insufficient mask coverage, failed constrained fit, signed/unsigned-flux or
  divergence failure, unintended nulls, and failed connectivity remain fatal;
  `diagnostic-only` suppresses only the production radialization thresholds,
  not numerical or topological validity;
- `current_sheet.model=none` requires `run.intent=analytic-verification`: all SCS fields are
  zero/`none`, `radialization_gate=not-applicable`,
  `current_sheet_transport.model=not-applicable`, and the Parker
  start radius is derived as `R_b` (`R_i=R_scs=R_b` conceptually);
- production requires an independently checksummed D6 qualification reference
  for the same rotation/time interval as the magnetic map. `compare-only`
  applies no scale, requires `apply_stage=none`, zero scale bounds, no
  construction reference or construction folding-correction asset, and
  `construction_data_use_role=none`; its qualification reference has
  `qualification_data_use_role=qualification` and must pass D6.
  `magnetogram-scale` requires a checksummed construction reference with
  `construction_data_use_role=construction`,
  `apply_stage=photospheric-coefficients-before-pfss`, ordered positive scale
  bounds, and the two-pass lifecycle in Section 13.1. It derives one positive
  factor within its declared bounds and applies it exactly once to the
  flux-balanced photospheric coefficients before PFSS filtering. The required
  D6 qualification reference has `qualification_data_use_role=qualification`
  and is independent of the construction reference: the two assets have
  disjoint observation identities/data intervals (or an explicitly
  preregistered, nonoverlapping partition), provenance, and covariance. One
  asset, observation, or covariance realization may not both determine the
  scale and qualify or weight that candidate. Scaling any later field,
  reusing a scaled state as Pass A, qualifying against the construction
  residual, or tuning the factor to SEP output is forbidden. At least two
  unique, increasing nested-sphere radii lie at or outside `R_scs` and inside
  the exact spherical domain. Signed/unsigned, nested-spread, and independent
  angular-quadrature tolerances are positive and all must pass. The radii list
  uses the same bounded bare-SI list grammar and canonical serialization as
  the latitude diagnostic list. The radius convention is validated separately
  for construction and qualification: `sample-position` evaluates
  `mean(r(t)^2*abs(B_r(t)))` and requires its comparison radius to be zero;
  `pre-normalized-to-comparison-radius` requires a positive comparison radius
  inside the exact domain, exact membership in the nested-sphere list, and
  matching asset metadata. Pass A always evaluates the deterministic
  `r_cmp` rule in Section 13.1. `maximum_qualification_relative_mismatch`
  gates only the independent qualification product; a construction residual
  can never be used as a candidate-product weight. Analytic verification with
  `model=not-applicable-verification` makes all construction and qualification
  reference/correction/radius fields `none` or zero. Radius-ensemble fields are either all
  inactive or name one checksummed ensemble, stable member ID, positive prior
  weight, and `preregistered-topology-and-open-flux` role. The selected member's
  stored `(R_b,R_i,R_scs,resolved_transition_width_m)` tuple must equal the
  active PFSS/current-sheet configuration bit-for-bit after SI normalization,
  and the entered prior weight must equal that member's immutable stored prior;
  the runtime cannot relabel an arbitrary radius choice. A member is selected
  before SEP comparison and its identity is restart/fingerprint state; radius
  scans remain topology/uncertainty ensembles, not amplitude calibration;
- an enabled `[formation_height_validation]` names two immutable assets. The
  candidate-product asset enumerates the complete Cartesian product of
  magnetogram realization, field-scale mode/value, `(R_b,R_i,R_scs)`,
  density/wind member, and front-geometry/kinematics member, including each
  checksum and preregistered prior weight. Omitting an inconvenient tuple or
  adding one after viewing D1/D2 is invalid. The event-constraint asset carries
  the event-specific type-II/EUV likelihood in height and time, including
  fundamental/harmonic, projection, density-model, covariance, and data-use
  metadata. `density_conditioning=per-wind-density-member` is mandatory.
  `frequency-time-likelihood` is preferred: for every density/wind member the
  implementation computes the local electron plasma frequency, evaluates the
  declared fundamental/harmonic radio lane directly in frequency--time space,
  and carries the observational and density-member covariance. A
  `preinferred-height-time-likelihood` asset must declare the density model,
  checksum, fundamental/harmonic choice, projection treatment, and full
  covariance used by that inference. It is not statistically independent
  across candidate members that share or depend on that density inference;
  those members require a joint covariance or explicitly conditional
  likelihood. The constraint used to qualify or weight candidates has
  `event_formation_data_use_role=qualification`. A
  `withheld-validation` constraint is reported only after the candidate set
  and weights are frozen and cannot gate, select, or reweight it. Every tuple
  is rebuilt completely and compared jointly with the independent D6
  qualification product, coronal-hole/topology observations, D1, D2, and the
  formation likelihood. A D6 construction asset or any observation used to
  determine a magnetic scale cannot also supply candidate weight. Candidate
  weights and rejection rules are frozen in the assets; SEP intensities,
  spectra, onset, or fluence are forbidden selectors. `model=none` makes all
  other fields inactive and cannot support a claim that campaign C calibrated
  or validated shock-formation height;
- production transport across the PFSS/SCS join requires a resolved transition
  with `R_i<R_b` and positive width contained in the PFSS/SCS analytic overlap, endpoint
  gauge/normal-flux matching, the `quintic-vector-potential` profile,
  `transition_vector_potential=signed-mie`, the
  `zero-mean-mie-v1` gauge, and convergence. The signed flux must close before
  either potential is built; the unsigned `B_tilde` cannot be supplied to this
  reconstruction. Production therefore uses
  `overlap-minimized`; `direct-pfss-scs` and an unverified zero-thickness
  operator are verification-only;
- the first-release production transition requires
  `transition_hcs_policy=exclude-clearance`, positive clearance, a checksummed
  clearance-convergence asset with the exact refinement/metric contents in
  Section 6.4, a positive transition minimum field,
  `0<transition_jump_support_fraction<=1`, and a nonnegative normal-trace-jump
  tolerance. The configured clearance must equal the named production member
  of that asset. The transition clearance is the
  single authority for mesh validity, source, movers, observers, and exports
  near `S_tr`. Future `common-flux-surface-qualified` and
  `finite-thickness-qualified` spellings are not schema-5 grammar and return a
  precise `NotImplemented` migration diagnostic if encountered.
  `current_sheet_transport.pure_hcs_minimum_clearance_m` applies only to the
  pure finite-SCS/Parker HCS outside the transition: it is positive for
  `sector-confined`, while `ideal-coordinate-crossing` requires zero and the
  separately qualified ideal-HCS crossing operator. `not-applicable` is legal
  only for the no-SCS branch. The source transition selector cannot override
  either topology policy;
- `consumer_acceptance_budgets.model=preregistered-event-grade` requires the
  transition budget and both cohort-edge assets, conditionally requires the
  front-return asset as stated below, and uses the exact cohort partition and unsigned-flux
  footprint measure defined in Sections 6.4 and 10.3, and
  `exceedance_action=mark-not-event-grade`. Every maximum is finite in
  `[0,1]`; asset strata may tighten but never loosen the scalar envelope.
  Surface area and counterfactual number/energy rates use all source gates
  except the clearance mask. Finite observer/export footprints use physical
  unsigned magnetic flux. Characteristic lines and point observers instead
  receive a typed valid/rejected state and cannot enter a flux-fraction test.
  Runtime transition and front-return number fractions use represented
  physical weight in stable birth cohorts, never macroparticle count. Each
  energy-loss numerator is the sum of the immutable **birth kinetic energy**
  of particles lost through that channel divided by the corresponding cohort
  birth-energy denominator. Event-time kinetic-energy and momentum loss remain
  dimensional, unbounded diagnostics and are never substituted into a bounded
  fraction. Zero denominators produce the typed inapplicable states specified
  above. A removing front policy requires the front-return asset and finite
  registered limits in `[0,1]`; that asset and both `maximum_front_return_*`
  limits consume only `DelayedFrontReturn` numerators. Substituting
  `ImmediateShockAdjacentReturn`, a verification-only term, or an alias is a
  schema/integrity failure. Zero is legal and means that any measured loss
  fails event-grade. A conormal/no-through-flow or other
  nonabsorbing policy requires that asset `none` and both front-return maxima
  zero while retaining nonloss contact diagnostics. The disabled budget model
  is verification-only, makes every asset/maximum inactive, and cannot be
  labeled event-grade;
- every structured wind requires the full longitude Jacobian and D7.
  `flux-tube-polytropic` requires an isothermal/polytropic energy closure,
  `transonic-solve`, `kinematic_profile_manifest_file=none`,
  `kinematic_profile_manifest_schema=none`, `kinematic_interpolation=none`,
  a zero radial-projection guard, and either `uniform` or
  `tube-from-target-speed` temperature. The
  uniform branch is `analytic-verification`; target speed is `sensitivity`,
  uses the coupled transonic solve, and cannot be labeled `event-nominal`.
  `finite-radius` requires a positive asset-matched target radius, while
  `asymptotic` requires that radius to be zero. Target-speed plus
  `uniform-base-density` is forbidden for production/sensitivity reporting.
  `empirical-kinematic` requires `empirical-profile`, `versioned-profile`, a
  checksummed `sep-kinematic-wind-profile-v1` manifest whose independently
  checksummed channel records provide the required inner-density,
  outer-velocity, and species-temperature products,
  `quintic-hermite-c2-certified` interpolation, an explicit
  `R_w0<=r_a<r_b` overlap wholly supported by both products, exactly one
  positive mass-per-flux authority, no target-speed relation/coefficients, no
  extrapolation, and `report-and-gate`. Its construction records all have the
  `construction` role and carry uncertainties plus a covariance or named
  ensemble. `mass-per-magnetic-flux` requires a direct positive `eta_m`;
  `colocated-density-velocity-field` requires co-located positive density and
  velocity with nonzero mapped `B`; and
  `radial-flux-density-with-mapped-field` requires positive `F_m` and nonzero
  mapped `B_r` at the same position/epoch. Speed alone, or radial mass flux
  without that field, is underdetermined and fails. Exactly one normalization
  branch is active, and its represented flux/solid-angle measure is present;
- empirical profile discriminants are never inferred from file contents.
  Each channel independently declares its physical role/species, checksum,
  units, velocity component/reference frame and channel-local radial-projection
  guard (when applicable), deterministic consumer selector, abscissa, spatial
  coordinate frame, trace epoch,
  line catalogue and stable line IDs, immutable-background fingerprint,
  arc-length origin/orientation, support segments, provenance/data-use role,
  and uncertainty/covariance or ensemble identity. Consequently an inner
  density channel may use heliocentric radius while an outer velocity channel
  uses field-line arc length; no global coordinate label may silently coerce
  either channel. An outer-velocity channel with
  `velocity_component=radial, velocity_reference_frame=inertial` requires its
  own `0<minimum_radial_projection<=1` and applies the guard at every
  interpolation/evaluation point before `u_s=u_r/(t_hat dot e_r)`. A
  `field-aligned,inertial` channel stores `u_inertial dot t_hat` and subtracts
  `(Omega_F cross x) dot t_hat` exactly once; a `field-aligned,corotating`
  channel stores `u_s` directly. Both require positive resulting `u_s`, a zero
  inactive radial-projection guard, stable line IDs, and an exact immutable-
  background fingerprint. A `radial,corotating` scalar is noncanonical and
  unsupported by schema 5 (a heliocentric rigid rotation does not change the
  radial component), so that pair is rejected rather than assigned a second
  meaning. Here field aligned always means the outward geometric tangent
  `t_hat`, never the signed magnetic direction `b_hat`.
  Before preparation, topology selectors are expanded against the frozen line/
  tube catalogue and their explicit exclusions. More than one match is typed
  `AmbiguousKinematicRoute`; zero matches follow
  `required_consumer_coverage_policy` (fatal for event-nominal, explicit masked
  census only for sensitivity). Routing is independent of record/support order
  and component type. Failure of the selected channel's support, frame,
  projection, or certificate never falls back to another channel.
  `field-line-arclength` likewise requires the line
  catalogue/fingerprint plus its stored origin, orientation, and trace epoch;
  `heliocentric-radius` does not authorize arc-length evaluation. A radial
  asset can never satisfy the field-aligned branch by relabelling metadata.
  A radius-tabulated channel used on a traced line requires monotone `r(s)` on
  every support segment; every turning point, zero derivative, or reversal
  splits the channel into stable, nonoverlapping segments. Every required
  consumer interval must remain inside one supported topology segment without
  crossing a data discontinuity or line-ID event;
- every empirical channel supplies finite node value, first derivative, and
  second derivative. `quintic-hermite-c2-certified` means the unique interval
  quintic determined by those six endpoint data, not a generic spline alias.
  Validation checks value, first-derivative, and second-derivative continuity
  at every shared node, evaluates all real in-interval extrema, proves the
  declared positive quantities remain positive and within their adjacent-node
  no-overshoot envelope, and rejects uncertified intervals. A physical
  extremum is represented by an explicit node; changing interpolator or
  silently clipping a negative/overshooting polynomial is forbidden;
- empirical electron density requires the declared composition and charge
  states. `versioned-mixture` requires its checksummed composition and
  uncertainty/covariance asset; fixed proton/electron(/alpha) branches require
  the compatible zero/nonzero alpha field and reject conflicting composition
  files. The conversion to mass density follows quasineutrality and the
  `electron_mass_in_density` policy; no hidden mean-molecular-weight constant
  is permitted;
- an empirical wind may be `sensitivity` or `event-nominal`, not analytic
  verification. Its acceleration floor, absolute/normalized momentum limits,
  and overlap absolute-log/covariance-normalized mismatch limits are strictly
  positive. Event-nominal additionally requires
  `required_consumer_coverage_policy=fail-required-event-support`, a
  checksummed D7 qualification asset with
  `d7_qualification_data_use_role=qualification`, independent of every
  construction channel, and unique in-support
  comparison radii, positive preregistered density/speed/mass-flux bounds, and
  coverage fractions in `(0,1]`. Missing support for any required source,
  finite observer footprint, or export interval is fatal even if aggregate
  magnetic-flux coverage passes. `diagnostic-mask-sensitivity` is legal only
  for sensitivity and retains all uncovered physical measures. The D7 census
  is weighted by open magnetic flux/area, source incident-number/energy flux,
  observer exposure, and export support rather than trace count. It also
  serializes every stable rejected tube/line ID with a typed reason; those IDs
  and reasons are part of the D7 fingerprint and are never replaced by an
  aggregate trace count. A separate checksummed
  `d7_withheld_validation_asset_file` may have only
  `d7_withheld_validation_data_use_role=withheld-validation`; it is evaluated
  after the event-nominal decision is frozen and can never gate, select, or
  reweight a candidate. Construction, qualification, and withheld observation
  identities/data intervals and covariance partitions must be disjoint or
  explicitly modeled jointly. The radius
  list uses the canonical bare-SI list grammar and contains `1.1`, `2`, and `5`
  times the configured solar radius whenever each value lies in the physical,
  construction, and qualification support. D7 mass flux is always evaluated
  at `mass_flux_reference_radius_m`; a target-radius value additionally exists
  only for finite-radius target-speed mode. Every inactive branch field is
  exactly zero/`none`. These are the complete schema-5 D7 gates: required
  consumer support and weighted coverage, independent finite-radius
  density/speed/mass-flux comparisons, inner/outer overlap agreement,
  certified positive/no-overshoot `C2` interpolation, and the pointwise plus
  quantile momentum-residual gates. They do **not** constitute a continuous
  observational acceleration envelope, do not assert monotonic speed, and do
  not impose an undocumented universal speed ceiling. A full-support
  observational plausibility envelope is the reserved future branch specified
  above, not a hidden strengthening of D7;
- production requires `fold_action=fail`; `mark-invalid-diagnostic` is
  analytic-verification only, deactivates every affected cell/source/observer,
  and never licenses particle transport across a folded map;
- schema 5 requires
  `source_surface_coupling.winding_construction=radial-from-scs-boundary` and
  derives the manifest value `R_w=R_scs`; there is no second equality-only
  input authority. A distinct `R_w<R_scs` would require the
  future generalized solenoidal Piola pushforward, rederived nonradial
  wind/map coupling, and its own schema/tests; it cannot be emulated by
  changing the initial radius of the radial formulas;
- an isothermal open wind sets `solar_wind.polytropic_index=0` as an inactive
  sentinel and uses its explicit temperature; a polytropic open wind requires
  `1<gamma_w<3/2`. An isothermal closed plasma likewise sets
  `closed_polytropic_index=0`; a polytropic closed plasma requires
  `gamma_c>1` and positive enthalpy over every closed-loop point. Neither wind
  index is passed to the characteristic-speed or shock-EOS kernels;
- `open_closed_interface.policy=diagnostic-kinematic` is the default for the
  independently assembled analytic PFSS/SCS, open-wind, and closed-plasma
  background and is restricted to analytic verification or declared
  sensitivity work. A production event selects `bounded-approximation` or
  `stationary-td-equilibrium`; neither spelling converts a diagnostic state
  into an equilibrium. All three policies enforce the pointwise one-sided
  `B_n` and interface-relative `w_n=(u-v_I) dot n` bounds. The interface
  velocity asset is required exactly for `versioned-asset`; otherwise it is
  `none` and the corotating topology-surface velocity is derived once.
  `sharp-one-sided` requires zero transition width, no thickness ensemble,
  zero volume-residual limits and quantile, and emits every component and norm
  of the full vector traction jump plus the signed, absolute, and normalized
  mass-flux jump. Its absolute/relative mass-flux bounds are finite and
  nonnegative, `0<mass_flux_jump_quantile<1`, and the pointwise maximum is
  always gated. `finite-width-volume` requires positive width, a
  positive ordered thickness ensemble containing the selected width, zero
  traction and mass-flux-jump limits/quantiles, an EOS-consistent layer, and the full volume momentum
  residual. A finite layer may not select `stationary-td-equilibrium`.
  `diagnostic-kinematic` sets the active **force-balance** residual bounds to
  zero as non-gating sentinels but still supplies the active representation's
  force-residual quantile and at least two refinement levels so the reported
  distribution is reproducible. A sharp diagnostic retains its separately
  gated mass-flux-jump bounds and quantile because mass continuity is not a
  force-balance claim.
  `bounded-approximation` requires a checksummed uncertainty asset, at least
  three mesh/thickness convergence levels, strictly positive active absolute
  and relative residual bounds, and pointwise compliance across the
  preregistered uncertainty and convergence ensemble.
  `stationary-td-equilibrium` requires a sharp interface and
  `state_origin=equilibrium-solver|imported-equilibrium`; it rejects
  `analytic-composite`, requires nonnegative solver-tolerance traction bounds,
  nonnegative solver-tolerance mass-flux-jump bounds, at least three
  convergence levels, and gates the full vector traction and mass-flux jump at
  every point. For every active metric, `0<quantile<1`, the refinement-change bound
  is nonnegative, and a quantile or area mean cannot hide a failing pointwise
  maximum. Near registered null/cusp neighborhoods a relative diagnostic may
  be inapplicable, but the absolute gate remains finite and normative.
  `global-separatrix-pressure-scale-verification` is sensitivity-only. A
  `footpoint-separatrix-constrained` asset must preserve two-footpoint
  hydrostatic consistency but is not an equilibrium proof; fitting one scalar
  scale or copying temperature across the interface cannot satisfy the vector
  balance by construction;
- `open_open_interface.model=none-continuous-single-family` requires
  `representation=not-applicable`, `policy=not-applicable`, no interface or
  uncertainty assets, no thickness ensemble, and every residual limit,
  quantile, and convergence count set to its inactive zero value. This branch
  is legal only when every open-field channel is continuous and belongs to one
  declared plasma family; a fast/slow label, discontinuous manifest channel,
  or finite open--open transition activates
  `model=versioned-interface-catalog`. The catalog has stable interface and
  adjacent-family IDs, frame/epoch/background fingerprints, geometry and
  velocity provenance, and a checksum. `sharp-one-sided` has zero thickness
  and volume-residual controls and emits the full vector MHD traction jump in
  a deterministic surface basis plus the signed, absolute, and normalized
  mass-flux jump. The sharp branch requires finite nonnegative absolute and
  relative mass-flux-jump bounds, `0<mass_flux_jump_quantile<1`, and pointwise
  maximum compliance for both diagnostic and bounded policies.
  `finite-width-volume` requires a positive
  ordered thickness ensemble containing the selected thickness, zero traction
  and mass-flux-jump controls, an explicit force-term inventory, and the full volume momentum
  residual. `diagnostic-kinematic` leaves the active **force-balance**
  residual non-gating but still emits its pointwise distribution and a
  reproducible refinement sequence; a sharp diagnostic continues to gate the
  separately configured mass-flux-continuity bounds;
  `bounded-approximation` requires a checksummed uncertainty asset, at least
  three refinement/thickness levels, positive active absolute/relative bounds,
  and pointwise compliance across the preregistered ensemble. No open--open
  branch may claim stationary equilibrium. These rules apply to fast/slow and
  all other open-plasma class boundaries. Plasma-sheet edges retain their
  `[plasma_sheet]` controls, but their result enters the same D8 interface
  product with a stable interface-class tag. A global or signed mean never
  establishes balance, a quantile cannot hide a failed maximum, and a sharp
  traction jump cannot replace a smooth-layer volume residual;
- exactly one resolution-independent mass-flux normalization is active. A
  per-numerical-tube absolute `kg/s` input is not legal. Every structured
  production branch passes the D7 density/mass-flux comparison; replacing a
  uniform base density by a uniform mass-per-flux input is not an exemption;
- `composition=proton-electron` requires zero `alpha_to_proton_ratio` and zero
  closed-plasma alpha base temperature; `proton-electron-alpha` requires both
  to be positive. Production recommends
  `electron_mass_in_density=include`; the recorded neglect branch changes the
  density, characteristic speeds, and fingerprint and is verification-only;
- every time-step accuracy factor lies strictly in `(0,1)`, the fixed upper
  bound is positive, pitch-angle fields are active only for the diffusion
  mover, and exactly one stable `[species.ID]` record with a unique in-range
  compiled slot and positive base macroparticle weight is installed for every
  compiled species before allocation. Its expected symbol, mass, and charge
  must match the AMPS molecular table; duplicate symbols are legal when slots,
  stable IDs, and roles differ. A positive integer `mass_number` is mandatory
  for every selected per-nucleon ion and zero for electrons or species for
  which per-nucleon output is forbidden;
- disabled population control requires all global, per-block, and per-species
  population bounds and cadences to be zero; enabled population control requires
  `minimum_total_particles <= target_total_particles <= maximum_total_particles`,
  ordered per-block bounds, positive cadence, ordered per-species
  minimum/target/maximum values whose sums are compatible with the global
  bounds, and a merge/split implementation
  that closes species-resolved represented weight, momentum, and relativistic
  energy. It never merges different stable species IDs or compiled slots, even
  when their chemical symbols match;
- a swept active corridor names existing stable line IDs, has a positive
  width plus buffer, uses conservative block/bounding-volume intersection,
  closes AMR ancestors/children/ghost neighbors, and passes the face-connected
  Sun--source--observer test. It requires enabled field-line export and its
  `line_ids` are exactly a subset of the existing `[field_line.ID]` requests.
  `full-domain` makes every geometric corridor field inactive;
- both open-field turbulence branches require a positive reference radius,
  reference wave energy, reference correlation length, and the selected
  spectral model. The prescribed branch additionally requires its wave-energy
  radial exponent and `-1<=outward_cross_helicity<=1`; the WKB branch derives
  radial wave energy from wave action, sets the energy exponent to zero as an
  inactive sentinel, and requires `outward_cross_helicity=1` so `w_in=0`.
  The correlation-length radial exponent remains active in both branches;
  `power-law` requires its custom spectral index. A zero/missing WKB reference
  amplitude is invalid rather than an implicit no-turbulence fallback.
  `prescribed-bidirectional` closed turbulence requires
  a checksummed two-footpoint asset; `direct-mean-free-path` requires no wave
  asset and marks closed-field wave quantities inapplicable rather than zero.
  Every wave-applicable state emits `delta B/B`, `p_w/p`,
  `p_w/(p+rho*u_s^2)`, signed wave acceleration, and the normalized
  wave-force residual. The force residual is the universal uncoupled-wind
  gate. `delta_b_over_b_action=fatal-small-amplitude-closure` is required for
  WKB and any QLT/small-amplitude coefficient; only an empirical amplitude
  branch with an independently prescribed coefficient may select the
  diagnostic action. `p_w/p` is diagnostic rather than a universal fatal
  threshold. The wave-force absolute tolerance comes from a checksummed
  derivative/refinement asset, is nonnegative, and is never inserted into a
  denominator. Both force-fraction bounds are nonnegative,
  `0<wave_force_quantile<1`, and the pointwise combined maximum plus the
  separately bounded quantile must pass before publication.
  `closed_field_model=invalid`
  additionally requires a topology/reachability proof that no source, active
  mover cell, observer support, or exported line enters the closed region;
- for `single-power-law`, `reference_mfp_m` is `lambda_0` at
  `reference_radius_m`; for `smooth-broken-power-law`, the reference radius
  must equal `break_radius_m` and `reference_mfp_m` is `lambda_b`. A ballistic
  out-of-domain policy is legal only for an active finite coefficient model in
  analytic verification and is encoded as a typed infinite-MFP state, never a
  large silent number. `parallel_mfp_model=none` instead requires
  `out_of_domain_policy=not-applicable`;
- `parker` requires a non-`none` mean-free-path model and the spatial-diffusion branch and derives
  `kappa_parallel=v*lambda_parallel/3`; focused pitch-angle diffusion requires
  a non-`none` mean-free-path model and
  `focused_collision_model=pitch-angle-diffusion`; focused discrete scattering
  requires a non-`none` mean-free-path model and the schema-5
  `isotropic-poisson` kernel with no file. A `versioned-kernel` spelling is a
  reserved future capability and fails before allocation. Inactive collision
  branches are `none`, and neither perpendicular
  diffusion nor drift is inferred from the parallel mean free path;
- `ballistic-verification` is legal only with analytic-verification intent (or
  an imported bundle explicitly qualified for it); it requires
  `parallel_mfp_model=none`, zero spatial diffusion, no focused collision,
  perpendicular diffusion, or drift, and all coefficient scalars/assets
  inactive. It is a manufactured-test mover, not a production SEP closure;
- schema-5 transport requires `drift_model=none` for both Parker and focused
  movers. A nonzero perpendicular coefficient is rejected for production until a conservative
  open/closed-separatrix crossing operator has passed its event and flux tests;
  it remains legal only in analytic domains proven to contain no reachable
  separatrix. Strict 3-D/1-D parity requires both perpendicular diffusion and
  drift to be `none` and either an infinitesimal characteristic comparison, a
  manufactured transversely uniform finite tube, or demonstrated convergence
  of a traced-footprint quadrature subdivision;
- a nonzero plasma-sheet model requires contrast `>=1`, positive width, and
  `fixed-temperature`; it modifies the selected base or outer normalization,
  never both. It selects `diagnostic-kinematic` for analytic/sensitivity use
  or `bounded-approximation` for event-grade use and is always evaluated with
  the smooth-layer volume momentum residual, never a sharp traction jump. The
  bounded policy requires a checksummed uncertainty asset, an ordered positive
  half-width ensemble containing the selected width, at least three
  radial/angular convergence levels, positive dimensional and normalized
  pointwise bounds, `0<volume_residual_quantile<1`, and a nonnegative
  refinement-change bound. Model `none` requires `not-applicable` and zero or
  `none` for every residual-control field;
- production uses `plasma_eos.adiabatic_index=1.6666666666666667` (the
  floating-point representation of the physical value `5/3`); another positive EOS index
  is verification-only. The sole critical-Mach selector and table asset live
  in `[shock]`, and their declared EOS/Mach convention must match the plasma
  EOS. A `none` model requires table and convention `none`, policy
  `not-applicable`, and all three exclusion budgets zero; an enabled model
  requires exactly one matching fast/total-Alfvén/normal-Alfvén convention and
  a table covering the entire obliquity interval `[0,pi/2]`.
  EOS/convention/obliquity-coverage mismatch is always fatal.
  `surface_requirement=fast` may use `diagnostic-inapplicable` for a
  beta-domain or verified exact-normal-Alfvén-asymptotic miss, with all three
  exclusion budgets zero/inactive and no source change.
  `surface_requirement=supercritical`
  requires a non-`none`, successfully checksummed table plus either
  `fail-preflight` or `exclude-source-budgeted`; the latter requires positive
  preregistered area, incident-number, and incident-kinetic-energy limits,
  zero source on missed patches,
  no redistribution to neighboring patches, and complete-history coverage
  below all three Section 9.4 measures. Silent clamping/extrapolation is
  forbidden;
- `weak_field_policy=exclude-source` excludes a patch exactly when
  `abs(B)<weak_field_relative_threshold*weak_field_reference_tesla`; both
  factors are positive, the excluded area/rate are reported, and `diagnostic`
  never changes source eligibility;
- a focused run must satisfy the selected PFSS/SCS kink, pure-HCS, and
  transition-sheet topology policies. Under `exclude-clearance`, neither a
  particle-support region nor a production line/footprint may touch the mask;
- closed-field injection requires an explicit closed turbulence/mover model;
- production requires exactly one rigid sidereal rotation authority in
  `[solar_rotation]`; the static closed provider, open-tube frame, Parker map,
  observers, and exporter consume its resolved value read-only. Provider-local
  input copies are unknown keys. `model=rigid` requires a positive
  `rigid_input_rotation_rate_rad_per_s`, requires the differential asset to be
  `none`, and requires the conversion ephemeris exactly when the input
  convention is synodic. `SepCoronalCmeConfiguration::Create()` performs that
  conversion before provider validation, publishes the derived immutable
  `resolved_sidereal_rotation_rate_rad_per_s`, and fingerprints the raw value,
  convention, resolved value, rotation axis/frame, and asset checksum.
  `latitude-dependent-verification` requires zero rigid input, convention
  `sidereal`, no synodic asset, and one checksummed differential-rotation asset;
  it is analytic-verification only until a selectable time-dependent closed
  provider and 3-D mapping exist;
- exactly one ellipsoid radial parameterization and one evolution branch is
  active. Component-law fields are required only for
  `independent-component-laws`; for `tabulated-snapshots` their kinematics are
  `inactive`, all numbers are zero, and the snapshot file is required. Schema
  5 orientation is fixed: tabulated snapshots may vary center and semiaxes but
  not the radial principal basis or tilt;
- the derived candidate-front apex is nondecreasing over the production source
  interval; reaching the zero-source radius latches source termination;
- an enabled piston requires the complete independent `piston_geometry` asset;
  the center-versus-apex exclusivity rule applies separately inside the front
  and piston records, not between those two surfaces;
- an enabled source has exactly one stable-ID record per compiled species and
  schema 5 unconditionally rejects a neutral compiled slot before allocation.
  It also requires
  `population_semantics=net-first-passage-upstream-released` and
  `release_model=empirical-net-first-passage-release`; no other interpretation
  is inferred from the DSA slope. The preferred event branch requires
  `release_boundary=upstream-reference-first-passage`, a positive finite
  `reference_surface_distance_m`,
  `shock_reentry_policy=absorb-post-reference-return-loss-ledger`, and a
  non-`none` normalized placement kernel with
  `placement_kernel_support=upstream-one-sided`; its full support remains
  strictly between the reference surface's upstream side and the active-domain
  outer boundary. Symmetric or downstream-overlapping kernels are invalid.
  The distance is smaller than every sampled local normal reach and is
  invariant under mesh/placement refinement. An offset-surface fold,
  overlap, front contact, active-mask contact, or incomplete history is fatal.
  An event-nominal source must select the preferred reference-surface branch.
  In schema 5 this distance is one fixed, positive SI length shared by species,
  momentum, patch, and event time. It is neither recomputed from the selected
  transport coefficient nor silently moved when that coefficient changes.
  `absorb-post-reference-return-loss-ledger` is likewise literal: after a
  committed positive-distance first passage, the first subsequent front hit is
  terminal and is recorded as `DelayedFrontReturn`; schema 5 performs no
  delayed re-emission or renewal. These fixed-length/absorbing semantics are
  part of the resolved source fingerprint and must remain invariant across
  restart.

  D10 additionally reports, without changing either surface geometry or
  release normalization, the dimensionless integrated transport-depth
  diagnostic

  \[
  \mathcal P_a(p,\sigma,t)=
  \int_0^{L_{\rm ref}}
    \frac{u^{\rm in}_{1n}(d,\sigma,t)}
         {\kappa_{nn,a}(d,p,\sigma,t)}\,{\rm d}d .
  \]

  Here \(u^{\rm in}_{1n}>0\) is the front-frame upstream inflow, and
  \(\kappa_{nn}=\mathbf n\mathbin{\cdot}\boldsymbol\kappa
  \mathbin{\cdot}\mathbf n\) is evaluated from the same coefficient authority
  and on the same normal ray as transport. Parker uses its transport tensor;
  focused transport reports the explicitly labeled diffusion-approximation
  proxy formed from \(\kappa_\parallel=v\lambda_\parallel/3\) and the
  independently declared perpendicular coefficient. This proxy never changes
  the focused mover. For constant quantities this reduces
  to \(\mathcal P=L_{\rm ref}u^{\rm in}_{1n}/\kappa_{nn}\). Therefore
  \(\mathcal P<1\) means that the fixed surface lies *within* one conventional
  diffusion length; it does not prove that front returns are rare. A
  branch with no declared spatial-diffusion interpretation reports a typed
  inapplicable state. Where the diagnostic is applicable, its sign convention
  and finite positive inflow and \(\kappa_{nn}\) must be certified over the
  complete ray; failure is a fatal background/coefficient inconsistency. D10
  also publishes the measured
  finite-horizon no-front-return fraction from the cohort ledgers, with the
  exact evaluation horizon. It must not call this a time-asymptotic survival or
  return probability, because it also depends on pitch-angle memory, moving
  geometry, boundaries, and the finite horizon;

  Schema 5 deliberately exposes no input key that can tune the D10
  interpolation or transport-depth quadrature into a different physical
  result. Production uses the versioned implementation-owned
  `D10TransportDepthQuadratureV1` policy. Its ray interpolation rule, absolute
  and relative error bounds, subdivision limit, and finite-value checks are
  fixed in the checksummed release-evidence tolerance profile and emitted in
  the resolved manifest. Missing or mismatched policy evidence is fatal. A
  verification harness may tighten those values to demonstrate convergence,
  but an event deck cannot override them. Making the policy runtime-selectable
  would require a later schema and explicit keys. The resolved algorithm ID
  and tolerance-profile checksum, rather than unowned hidden constants, enter
  `release_calibration_fingerprint`;
  `front-conormal-flux-no-through-flow-sensitivity` is Parker-only, requires
  `scientific_role=sensitivity`,
  `shock_reentry_policy=conormal-no-through-flow`, zero reference distance,
  `placement_kernel=none`, `placement_kernel_support=none`, zero placement
  controls, and a finite strictly positive
  `minimum_conormal_diffusivity_m2_per_s`. Before forming its direction it
  evaluates `n dot kappa dot n`; a nonfinite value or one not greater than the
  configured minimum is a typed rejection, so the code never divides by a
  vanishing conormal diffusivity. The direction is then exactly
  `kappa*n/(n dot kappa dot n)`. It must pass a conormal-flux manufactured
  solution and may not claim physical reflection. Every nonconormal branch
  sets `minimum_conormal_diffusivity_m2_per_s=0`.
  `shock-adjacent-absorbing-verification` requires
  `scientific_role=verification`, its matching re-entry policy, and a
  non-`none` kernel with `placement_kernel_support=upstream-one-sided`; it is legal only with
  `run.intent=analytic-verification` and is rejected by event-nominal and
  observation-comparison products;
  `flux-fraction` activates
  `upstream_release_fraction` and an
  upstream abundance; `physical-rate` activates exactly one scalar or
  checksummed rate asset plus a nonnegative patch partition, never both.
  Scalar physical-rate authority requires a positive
  `physical_particle_rate_per_s` and no rate asset; asset authority requires a
  checksummed rate asset and a zero scalar. The flux-fraction branch requires
  `physical_rate_authority=inactive`, zero scalar, no rate asset,
  `physical_rate_patch_distribution=none`, and no patch asset. An
  `area` or `incoming-species-flux` patch distribution forbids the patch asset;
  `versioned-asset` requires it;
  Per-nucleon energy requires a valid
  mass number, a fixed spectrum requires a positive phase-space index, and the
  source frame/angular/placement branches are complete. Every production
  source-species record has positive `samples_per_step` and a normalization
  that is positive on at least one eligible source interval; a species may be
  locally zero because no patch is active, but an identically zero record may
  not silently exclude a compiled species. Number- and
  nonthermal-energy-budget gates must pass before weights are initialized.
  For focused transport the preferred branch requires
  `focused_release_phase_space=outward-first-passage-flux-weighted`; the
  parser rejects `unrestricted-isotropic-verification` outside the absorbing
  verification branch. The authoritative normalization is
  `Z=integral a(mu)*[w_n(mu)]_+ dmu`, evaluated in the declared transport/source
  frame with the surface velocity transformed into that same frame; a derived
  pitch cutoff `mu_c` is diagnostic only.
  `outward-first-passage-flux-weighted` requires all six
  `focused_release_*tolerance*`/`focused_release_mu_c_minimum_abs_bdotn`
  controls above to be finite and strictly positive; both relative tolerances
  are less than one. Parker, conormal, and unrestricted-isotropic branches set
  all six to their inactive zero value. The `mu_c` threshold controls only
  whether that diagnostic is reported: crossing it never selects physical
  pitch support or changes `Z`. The root and quadrature controls bound a
  sign-certified calculation; they are never used as `Z<=epsilon` escape
  criteria, normalization floors, or permission to clamp a negative estimate.
  When `abs(b dot n)` exceeds the registered diagnostic threshold, the
  ordinary signed-pitch diagnostic may be reported.
  In the tangent-field limit it never divides by `b dot n`: if
  `(u-V_ref) dot n>0`, the entire declared pitch support remains eligible and
  is flux weighted by the positive normal speed; otherwise there is no focused
  escape. Tests exercise positive, zero, and negative tangent advection and a
  nonparallel transport-frame transformation. Pitch integration is partitioned
  at seed-support boundaries, pitch-law knots, and every isolated root of
  unclipped `w_n=0`; momentum and event-time intervals are split at every
  sign-certified admissibility boundary within their configured root
  tolerances. A root-isolation or quadrature failure is
  `UnresolvedPositiveFluxNumerics` and is fatal before allocation, never
  `NoFocusedEscape`. `NoFocusedEscape` is evaluated per species/momentum
  interval/patch/time interval with the signed magnetic orientation. A
  `fail-preflight` policy rejects any positive candidate rate on unavailable
  support and is mandatory for `scientific_role=event-nominal`.
  `exclude-budgeted` requires `scientific_role=sensitivity`, both no-escape
  fraction bounds in `[0,1)`, debits number and energy before sampling without
  renormalization, and fails when either bound is exceeded. Parker and the conormal branch
  require `no_focused_escape_policy=not-applicable` and zero no-escape bounds;
- source accounting uses the following exact, schema-stable ledger names:
  `CandidateReferenceRelease`, `NoFocusedEscapeExcluded`,
  `CommittedFirstPassageRelease`, `ImmediateShockAdjacentReturn`,
  `DelayedFrontReturn`, `NetFirstPassageRelease`, and the time-dependent
  `SurvivingUpstreamInventory`. Their cohort identity,
  represented number, immutable birth kinetic energy, event-time kinetic
  energy, and momentum columns are not aliases and cannot be renamed or
  collapsed. The input spelling
  `shock_reentry_policy=absorb-post-reference-return-loss-ledger` is a policy
  selector, not a second ledger name: every front hit after a committed
  positive-distance flight debits exactly one `DelayedFrontReturn` record.
  `PostReferenceFrontReturnLoss` is not a schema-5 ledger term and is rejected
  if supplied by an asset or restart. The source identities
  `CandidateReferenceRelease = NoFocusedEscapeExcluded +
  CommittedFirstPassageRelease` and
  `NetFirstPassageRelease = CommittedFirstPassageRelease` are checked per
  stable species/birth-energy/patch/lineage/time cohort. At ledger time `t`,
  `SurvivingUpstreamInventory(t)` equals committed represented content through
  `t` minus every terminal physical loss through `t`; it is never relabeled as
  a first-passage release rate. Only the shock-adjacent absorbing verification
  branch may publish `GrossShockAdjacentEmission` and
  `HitEscapeSurfaceBeforeShock`; those diagnostic-horizon terms cannot alter
  any preferred-branch release identity. All source, surface,
  coordinate/storage, and mover-frame transformations are explicit and are
  fingerprinted with the source record;
- the schema-5 resolved source record contains one
  `release_calibration_fingerprint`. This is a derived manifest field, not a
  second parser authority. Its canonical hash covers the raw fixed
  `reference_surface_distance_m`, unit/frame and normal convention, source and
  placement semantics, the selected `run.transport` mover and
  `run.transport_frame`, shock/background/coefficient generations, stable
  species and momentum support, D10 diagnostic algorithm ID and
  release-evidence tolerance-profile checksum, and the exact calibration
  horizon. In schema 5 that horizon binds the raw closed
  `[run.start_time_s, run.end_time_s]` interval together with the declared
  cohort follow-up, right-censoring, or survival/competing-risk interpretation;
  generic time support is not a substitute for it. The hash also covers
  focused-release root/quadrature tolerances, and the checksums/data-use roles
  of any preregistered campaign
  evidence used to choose the distance. If no external calibration evidence is
  used, that state is explicit. The scalar input remains the sole geometry
  authority; descriptive provenance cannot override it. Changing any covered
  item changes the fingerprint and invalidates restart compatibility and any
  previously qualified D10/front-return product;
- the analytic kinetic-conversion budget requires a complete admissible
  downstream state and `energy_budget_asset_file=none`; the asset branch
  requires a positive checksummed patch/time-dependent available-flux record.
  Production enforces
  `0<maximum_nonthermal_energy_fraction<=test_particle_energy_fraction_limit<1`;
  both limits are preregistered release inputs, not hidden defaults;
- Parker transport requires the upstream-plasma frame and uses the
  `1/(4*pi*p^2)` phase-space Jacobian; isotropy is a property of its transported
  distribution, not a license to create a shock-adjacent half-population.
  Focused transport uses the normalized joint
  `H(p,mu,t)/(2*pi*p^2)` first-passage law. Its seed angular density may be
  isotropic or versioned, but the realized law is restricted and weighted by
  positive outward normal flux. A source-frame boost carries its complete
  momentum-space Jacobian into the local-plasma mover frame and cannot preserve
  only energy-bin centers.
  `outward-wave`/`inward-wave` are focused-only parallel boosts and require a
  positive corresponding directional wave population at every injecting
  patch; a general shock-frame source is unavailable in schema 5;
- `pfss_scs_transition_policy=convergence-qualified` requires a registered
  width/profile ensemble; otherwise transition patches are diagnostic-only
  with zero injection and a separate ledger reason;
- observer IDs are unique; each application accepts only its compatible
  position/estimator branch; energy edges are strictly increasing and match
  the declared species/energy coordinate. Residence estimators require a
  positive volume and accumulation interval; crossing estimators require a
  positive derived area, angular acceptance, crossing sense, and geometric
  factor. A disk alone requires one normalized fixed normal; a sphere uses its
  local outward normal and forbids fixed-normal fields. The radius uniquely
  determines volume/disk area/sphere area. `position_mode=fixed-cartesian`
  requires `velocity_source=coordinate-frame-stationary`, no ephemeris file,
  and the three explicit coordinates. `position_mode=ephemeris-file` requires
  the ephemeris velocity source, a checksummed file covering the complete
  sampling horizon, a strictly monotone detector-acquisition-clock to run-time
  mapping, and zero fixed-coordinate fields. A fixed-look cone requires a normalized look axis
  and half-angle in `(0,pi]`; omnidirectional sampling makes those fields zero.
  Pitch edges are strictly increasing within `[-1,1]` and are focused-only.
  `reported_intensity=omnidirectional` requires either the explicitly isotropic
  Parker distribution or focused pitch bins that cover all of `[-1,1]`
  without gaps and full gyrophase acceptance. Partial pitch/look coverage may
  report only `directional-per-sr` or
  `accepted-solid-angle-integrated`; no undocumented angular inversion is
  allowed.
  Binning uses detector-rest-frame momentum and electromagnetic field. Online
  response requires a checksummed response asset;
- the output cadence is positive, both initialization filenames remain below
  the resolved output directory, Tecplot initialization output is enabled for
  schema 5, and every potentially inapplicable physical column has a typed
  validity companion. The one manifest-owned finite sentinel is never a
  physical default and is rejected if it is nonfinite;
- if `field_line_export.enabled=false`, that section contains only the selector
  and no `[field_line.ID]` sections exist. In this branch, mesh tube refinement
  is `none` with zero tube size/decay, and a swept corridor is illegal. Enabled
  export requires at least one stable string line ID; tube refinement may use
  those same authoritative requests. Schema 5 exports geometry that is static
  in its transport frame or follows the exact declared rigid rotation; it does
  not evolve topology. `stationary-in-transport-frame` means node coordinates
  are constant only in `run.transport_frame`, with exact stored transformations
  to requested epochs. `rigid-rotation-from-trace-time` requires the rigid-
  corotating authority and its single validated rotation rate; an inertial
  frame receives the resulting time-dependent rotation rather than ignoring
  `timeS`. Shock and observer-mapping histories may vary in time.
  `observer_id=none` requires `observer_mapping=none`; a named observer
  requires the appropriate static or time-resolved mapping branch;
- a bundle intended for `srcSEP` moving-source transport requires
  `shock_intersections=true`, an ordered history interval covering every
  possible source-active time within the run, positive history cadence, and
  event-refined intersections independent of that output cadence. A
  background-only analytic bundle may disable intersections and sets all three
  history scalars to zero;
- an observer-connected schema-5 line requires a volume-sphere observer and a
  finite `flux-tube` or `quadrature` measure; characteristic lines cannot own
  count-based observer volume. `static-volume-intersection` is legal only when
  observer and line geometry are stationary in the same frame and sets mapping
  cadence/event-time tolerance to zero. Any ephemeris,
  a fixed inertial observer against a rotating line, or other relative motion
  requires `time-resolved-volume-intersection`, positive mapping cadence and
  event-time tolerance, stable overlap-component IDs, and complete topology-
  segmented coverage of the run horizon. Interpolation across an overlap
  creation/annihilation event is forbidden and out-of-tolerance action is
  `fail`;
- every production-exported line has one unique composite-footpoint sector,
  positive transition-sheet clearance, and agreement with the SCS construction
  sector outside the overlap. A line or finite footprint touching the
  exclusion, transition null/kink, or ambiguous/mismatched sector is rejected
  as `topology-inapplicable`, serialized in the connection summary, and never
  shifted to a nearby line;
- a `characteristic` line has `measure_group_id=none`, zero represented flux,
  and no cross-section asset. A `flux-tube` line has no group ID but has a
  positive flux and a valid traced-reference-footprint. Every `quadrature`
  line names a measure group, has a positive flux and nonoverlapping traced
  footprint, and derives its weight as its flux divided by the group's summed
  flux; no independent weight is accepted. Observer volume and shock-source
  allocation use this identical footprint/flux measure and must close within
  `maximum_cross_section_relative_error`;
- `seed_mode=observer-connected` or `cartesian` requires
  `trace_branches=both-from-seed`: the tracer integrates both field signs,
  identifies the Sun-connected branch, joins it through the seed, and only
  then orients arc length geometrically outward. A photospheric seed requires
  `outward-from-photosphere`. Observer-connected seeds require a positive
  `connection_tolerance_m`; other seed modes set it to zero and do not acquire
  an observer by proximity. Every `field_line.outer_radius_m` lies above all
  of its mapped observer support and at or below `domain.outer_radius_m`, and
  the selected open trace must actually reach it;
- disabled source syntax consists only of `source.enabled=false` and has no
  source-species records; it is legal only for analytic verification and does
  not invent placeholder physical values;
- inactive fields are zero/`none`, not ignored; and
- `taper_start_apex_radius_m < zero_source_apex_radius_m` for a taper.

The resolved manifest records all raw and derived inputs: critical points,
PFSS/SCS and sector-map hashes, coefficient-file content hashes, topology and
longitude-map hashes, minimum `A_phi`, first fold radius if diagnostic,
structured-wind residuals, D6--D10 products and tolerance assets, the single
resolved solar-rotation authority, source-population semantics, geometry and
activation histories, whole-dome and gate-region areas, and source-envelope/
loss accounting, including the derived `release_calibration_fingerprint`,
integrated transport-depth strata, and finite-horizon return census. File identity is based on
canonical role and content checksum, not path alone.

The `srcSEP` consumer uses a separate runtime input and never parses the
`srcSEP3D` input indirectly. It reuses the exact schema-5 definitions of
`[particle_numerics]`, repeated `[species.ID]`, and
`[population_control]`. Numerical weights/population bounds are local controls,
while each species record binds one stable bundle identity to the consumer's
own verified compiled slot. All background, source-intersection, physical
coefficient, source-species, and observer definitions come from the immutable
bundle. These bindings and controls participate in the `srcSEP` fingerprint.
Its complete application-level delta is:

```ini
[run]
schema_version = 5
intent = imported-field-line-transport
transport = ballistic-verification | parker | focused-pitch-angle-diffusion | focused-discrete-scattering
transport_frame = bundle
start_time_s = REQUIRED
end_time_s = REQUIRED
maximum_steps = REQUIRED_OR_ZERO
campaign_seed_u64 = REQUIRED
random_stream_layout = keyed-v1

[particle_numerics]
time_step_model = fixed-upper-bound | adaptive-local
maximum_time_step_s = REQUIRED
spatial_accuracy_factor = REQUIRED
pitch_angle_accuracy_factor = REQUIRED_OR_ZERO
pitch_endpoint_regularization = REQUIRED_OR_ZERO
shock_motion_accuracy_factor = REQUIRED
species_binding = stable-id-plus-compiled-slot-verified
species_numerics = explicit-all-compiled

[species.ID]
compiled_slot = REQUIRED
chemical_symbol = REQUIRED
expected_mass_kg = REQUIRED
expected_charge_c = REQUIRED
mass_number = REQUIRED_OR_ZERO
transport_role = charged-sep | initialization-only
base_macroparticle_weight = REQUIRED
minimum_population = REQUIRED_OR_ZERO
target_population = REQUIRED_OR_ZERO
maximum_population = REQUIRED_OR_ZERO

[population_control]
model = disabled | amps-conservative-split-merge
check_cadence_steps = REQUIRED_OR_ZERO
minimum_total_particles = REQUIRED_OR_ZERO
target_total_particles = REQUIRED_OR_ZERO
maximum_total_particles = REQUIRED_OR_ZERO
minimum_particles_per_active_block = REQUIRED_OR_ZERO
maximum_particles_per_active_block = REQUIRED_OR_ZERO
maximum_weight_ratio_per_merge_group = REQUIRED_OR_ZERO
conservation = species-weight-momentum-relativistic-energy

[field_line_input]
provider = sep-field-line-bundle
bundle_path = REQUIRED
line_ids = REQUIRED
required_bundle_schema = 3
spatial_extrapolation = reject
temporal_extrapolation = reject
shock_source = none | imported-intersections
coefficient_authority = bundle-manifest
observer_definitions = bundle
verify_checksums = true
require_background_fingerprint = REQUIRED
unresolved_discontinuity = reject

[line_mesh.ID]
line_id = REQUIRED
start_s_m = REQUIRED
end_s_m = REQUIRED
point_count = REQUIRED
resampling = conservative-positive-one-sided

[output]
directory = REQUIRED
cadence_s = REQUIRED
write_line_diagnostics = false | true
invalid_value_sentinel = -1.7976931348623157e308
require_validity_columns = true
```

`line_ids` preserves user order for reporting but does not change stable IDs,
random streams, or physical results. `srcSEP` validates the bundle before
allocating particles: schema, checksums, units, coordinate frame, epoch,
species compatibility, line openness, monotone arc length, observer bounds,
source normalization, and coefficient-model compatibility are all release
gates. A mismatch is a configuration error, never a warning followed by a
best-effort run. `srcSEP` rejects local `[observer.ID]`, `[source_species.ID]`,
or background-physics sections because they would create a second authority.
`field_line_input.shock_source=imported-intersections` requires the bundle
manifest flag `shock_intersections=true` and complete source-active history;
`shock_source=none` requires `shock_intersections=false` and creates no source
from the bundle. Either mismatch is fatal. A background-only bundle therefore
remains usable for source-disabled transport/verification without inventing an
empty intersection authority.
It may select only a mover that the bundle's stored coefficients and
discontinuity policy support. Each selected line has exactly one
`[line_mesh.ID]` record. Its interval lies inside the bundle's validated arc-
length domain, contains every mapped observer assigned to that line, and has
at least two points. Conservative positive-state resampling is performed
separately on every one-sided smooth segment; no mesh interval spans an HCS,
separatrix, radial interface, invalid state, or source event. Thus the 1-D
length, starting point, and total point count are explicit numerical choices
without becoming a second authority for the physical field-line geometry.

---
