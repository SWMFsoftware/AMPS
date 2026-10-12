#ifndef SEP_COMMON_SEP_INPUT_LAYOUT_H
#define SEP_COMMON_SEP_INPUT_LAYOUT_H

// ===========================================================================
// Shared SEP runtime-input vocabulary (content schema 6)
// ===========================================================================
//
// srcSEP (field-line transport) and srcSEP3D (three-dimensional transport)
// read one runtime-input vocabulary written in the generic AMPS block format
// (src/INPUT_FILE.md, reader src/general/input/amps_input.h).  This header is
// the single table of that vocabulary: every block and every
// application-owned key that either application understands, each marked with
// the applications that implement it.  Both parsers walk their documents
// against this table, so the two applications cannot drift into different
// spellings of the same setting.
//
// Design decisions (user, 2026-10-11; src/INPUT_FILE.md section 5):
//
//  D1  Both applications use the block '#begin application: sep' and the
//      content schema 6.  Every earlier schema (srcSEP 1-4, srcSEP3D 1-5) is
//      retired; a file that names one stops with a "retired schema" message.
//  D2  The turbulence description lives in '#begin turbulence' inside
//      'application: sep'.  It is the only block whose CONTENT differs
//      between the applications: srcSEP3D describes a prescribed spectrum,
//      srcSEP its self-consistent Alfven-wave model.  srcSEP's turbulence
//      command-line options remain and override the file.
//  D3  The vocabulary is identical; an application rejects any setting it
//      cannot honour with a typed "not supported by <application>" error and
//      never ignores it.  This applies equally to whole blocks (for example
//      mesh.active_region in srcSEP) and to individual values (for example a
//      srcSEP observer of kind fixed-cartesian).
//  D4  Every physical quantity is written in exactly one place.  When
//      'physics: swcme' is present it owns the ambient solar wind (speed,
//      rotation, density, temperature, field, polarity), so those keys are
//      illegal in 'physics: background' / 'parker'.  'source: shock' owns the
//      particle source (energy range, reference spectrum, injection
//      efficiency), so the SWCME 'source.*' keys are illegal in
//      'physics: swcme'; each application synthesizes them for the SWCME
//      resolver from 'source: shock' (SynthesizeSwcmeSource()).
//  D5  Observers share one vocabulary built on 'kind'.  srcSEP implements
//      kind = field-line-radius (heliocentric_radius_m on the field line);
//      srcSEP3D implements the spatial kinds.
//
// Units.  Application keys carry their SI unit in the key name (src/INPUT_FILE.md
// section 4); model blocks (physics: swcme, physics: parallel_diffusion,
// physics: shock_front) keep the vocabularies of their libraries.
//
// This header depends only on the C++ standard library so that the parsers
// can validate a complete input before AMPS, MPI or a mesh exists.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace SEP {
namespace InputLayout {

// The one content schema shared by both applications (decision D1).
constexpr unsigned kSchemaVersion = 6;

// Applications that implement a block or key.  The values are bit masks so a
// table row can name either application or both.
enum Support : unsigned {
  kSrcSep = 1u,    // srcSEP, field-line transport
  kSrcSep3d = 2u,  // srcSEP3D, three-dimensional transport
  kBoth = 3u
};

inline const char* ApplicationName(unsigned application) {
  return application == kSrcSep ? "srcSEP"
       : application == kSrcSep3d ? "srcSEP3D" : "srcSEP and srcSEP3D";
}

// One block of the shared layout.  'parent' is the qualified path of the
// enclosing block ("" = document root), 'name' is "" for an unnamed block and
// "*" for any one instance name ('#begin observer: earth').  'model' marks a
// block whose keys belong to a shared library rather than to this table.
struct BlockRule {
  const char* parent;
  const char* kind;
  const char* name;
  unsigned support;
  bool model;
};

constexpr BlockRule kBlocks[] = {
    // Run control: time step, stop condition, seed, restart files.
    {"", "run", "", kBoth, false},
    {"run", "restart", "", kSrcSep3d, false},
    // Application section: content schema, transport choice and numerics.
    {"", "application", "sep", kBoth, false},
    {"application.sep", "transport", "", kSrcSep3d, false},
    {"application.sep", "population_control", "", kSrcSep3d, false},
    {"application.sep", "storage", "", kSrcSep3d, false},
    // The only block whose content differs between the applications (D2).
    {"application.sep", "turbulence", "", kBoth, false},
    // AMR support mesh and the embedded field line.
    {"", "mesh", "", kBoth, false},
    {"mesh", "domain", "", kBoth, false},
    {"mesh", "field_line", "", kBoth, false},
    {"mesh", "refinement", "solar", kBoth, false},
    {"mesh", "refinement", "tube", kBoth, false},
    {"mesh", "active_region", "", kSrcSep3d, false},
    {"mesh", "memory", "", kSrcSep3d, false},
    // Physics.  physics: background selects the plasma/IMF provider; the
    // nested parker block is the analytic Parker record of srcSEP3D.
    {"", "physics", "background", kBoth, false},
    {"physics.background", "parker", "", kSrcSep3d, false},
    {"", "physics", "shock", kBoth, false},
    {"", "physics", "swcme", kBoth, true},
    {"", "physics", "parallel_diffusion", kBoth, true},
    {"", "physics", "shock_front", kSrcSep3d, true},
    // Particle source, species weights, observers, output products.
    {"", "source", "shock", kBoth, false},
    {"", "species", "", kBoth, false},
    {"", "observer", "*", kBoth, false},
    {"", "output", "", kBoth, false},
    // Fixed point sampling in the AMPS output cycle (srcSEP3D; user decisions
    // 2026-10-11).  One 'species: SYMBOL' block per sampled AMPS species holds
    // its energy channels; 'point: NAME' blocks and/or a point file hold the
    // locations.  srcSEP samples along field lines and does not implement it.
    {"", "sampling", "points", kSrcSep3d, false},
    {"sampling.points", "species", "*", kSrcSep3d, false},
    {"sampling.points", "point", "*", kSrcSep3d, false},
};

// One application-owned key.  'block' is the qualified path of the block that
// holds it; "observer.*" stands for every '#begin observer: ID' block.
struct KeyRule {
  const char* block;
  const char* key;
  unsigned support;
};

constexpr KeyRule kKeys[] = {
    // ---- run ---------------------------------------------------------------
    {"run", "time_step_s", kBoth},
    {"run", "time_step_mode", kSrcSep3d},
    {"run", "maximum_particle_speed_m_s", kSrcSep3d},
    {"run", "time_step_margin_factor", kSrcSep3d},
    {"run", "maximum_time_steps", kBoth},
    {"run", "stop_shock_radius_m", kSrcSep3d},
    {"run", "campaign_seed", kSrcSep3d},
    {"run.restart", "input_path", kSrcSep3d},
    {"run.restart", "output_path", kSrcSep3d},

    // ---- application: sep ----------------------------------------------------
    // transport = parker | focused-diffusion | focused-scattering selects the
    // particle-transport formulation in both applications.  srcSEP maps it to
    // its production movers parker, fte-dmumu and fte-mfp (the identical
    // Parker, continuous D_mumu and discrete nu = v/lambda_parallel
    // formulations; srcSEP/PRODUCTION_MOVER_API.md).
    {"application.sep", "schema_version", kBoth},
    {"application.sep", "transport", kBoth},
    {"application.sep", "intent", kSrcSep3d},
    {"application.sep", "background_cadence_steps", kSrcSep3d},
    {"application.sep", "injection_cadence_steps", kSrcSep3d},
    {"application.sep.transport", "cell_crossing_fraction", kSrcSep3d},
    {"application.sep.transport", "diffusion_fraction", kSrcSep3d},
    {"application.sep.transport", "focusing_fraction", kSrcSep3d},
    {"application.sep.transport", "cooling_fraction", kSrcSep3d},
    {"application.sep.transport", "field_variation_fraction", kSrcSep3d},
    {"application.sep.transport", "shock_crossing_fraction", kSrcSep3d},
    {"application.sep.transport", "minimum_substep_s", kSrcSep3d},
    {"application.sep.transport", "maximum_substeps", kSrcSep3d},
    {"application.sep.transport", "pitch_angle_scheme", kSrcSep3d},
    {"application.sep.transport", "spatial_diffusion_model", kSrcSep3d},
    {"application.sep.transport", "pitch_angle_diffusion_model", kSrcSep3d},
    {"application.sep.transport", "mean_free_path_model", kSrcSep3d},
    {"application.sep.transport", "constant_dmumu_per_s", kSrcSep3d},
    {"application.sep.transport", "constant_mean_free_path_m", kSrcSep3d},
    {"application.sep.transport", "mean_free_path_reference_m", kSrcSep3d},
    {"application.sep.transport", "mean_free_path_reference_radius_m", kSrcSep3d},
    {"application.sep.transport", "mean_free_path_reference_rigidity_v", kSrcSep3d},
    {"application.sep.transport", "mean_free_path_radial_exponent", kSrcSep3d},
    {"application.sep.transport", "mean_free_path_rigidity_exponent", kSrcSep3d},
    {"application.sep.transport",
     "spatial_quadrature_absolute_tolerance_m2_per_s", kSrcSep3d},
    {"application.sep.transport", "spatial_quadrature_relative_tolerance", kSrcSep3d},
    {"application.sep.transport", "spatial_quadrature_maximum_recursion", kSrcSep3d},
    {"application.sep.transport", "focused_scattering_frame", kSrcSep3d},
    {"application.sep.transport", "maximum_scattering_events_per_substep", kSrcSep3d},
    {"application.sep.transport", "perpendicular_diffusion", kSrcSep3d},
    {"application.sep.transport", "constant_kappa_perpendicular_m2_per_s", kSrcSep3d},
    {"application.sep.transport", "kappa_perpendicular_to_parallel_ratio", kSrcSep3d},
    {"application.sep.transport", "drifts", kSrcSep3d},
    {"application.sep.population_control", "mode", kSrcSep3d},
    {"application.sep.population_control",
     "minimum_particles_per_cell_per_species", kSrcSep3d},
    {"application.sep.population_control",
     "target_particles_per_cell_per_species", kSrcSep3d},
    {"application.sep.population_control",
     "maximum_particles_per_cell_per_species", kSrcSep3d},
    {"application.sep.population_control", "cadence_steps", kSrcSep3d},
    {"application.sep.storage", "magnetic_gradient", kSrcSep3d},
    {"application.sep.storage", "velocity_gradient", kSrcSep3d},
    {"application.sep.storage", "sampling_bytes_per_cell", kSrcSep3d},

    // ---- application: sep / turbulence (decision D2) -----------------------
    // srcSEP3D: prescribed (or SWMF) turbulence spectrum and amplitude law.
    {"application.sep.turbulence", "authority", kSrcSep3d},
    {"application.sep.turbulence", "model", kSrcSep3d},
    {"application.sep.turbulence", "amplitude_model", kSrcSep3d},
    {"application.sep.turbulence", "delta_b_over_b", kSrcSep3d},
    {"application.sep.turbulence",
     "wave_energy_density_at_reference_j_per_m3", kSrcSep3d},
    {"application.sep.turbulence", "wave_energy_density_radial_exponent", kSrcSep3d},
    {"application.sep.turbulence", "normalized_cross_helicity", kSrcSep3d},
    {"application.sep.turbulence", "reference_radius_m", kSrcSep3d},
    {"application.sep.turbulence", "k_min_per_m", kSrcSep3d},
    {"application.sep.turbulence", "k_max_per_m", kSrcSep3d},
    {"application.sep.turbulence", "k_min_radial_exponent", kSrcSep3d},
    {"application.sep.turbulence", "k_max_radial_exponent", kSrcSep3d},
    {"application.sep.turbulence", "spectral_index", kSrcSep3d},
    {"application.sep.turbulence", "correlation_length_m", kSrcSep3d},
    {"application.sep.turbulence", "correlation_length_radial_exponent", kSrcSep3d},
    {"application.sep.turbulence", "validity_cadence_s", kSrcSep3d},
    {"application.sep.turbulence", "missing_data", kSrcSep3d},
    {"application.sep.turbulence", "resonance_range", kSrcSep3d},
    {"application.sep.turbulence", "self_consistent_3d", kSrcSep3d},
    // srcSEP: self-consistent Alfven-wave turbulence on the field line
    // (SEP::Turbulence::Configuration, srcSEP/TURBULENCE_MODEL.md).  Each key
    // has a command-line option that overrides it (srcSEP/README.md).
    {"application.sep.turbulence", "source", kSrcSep},
    {"application.sep.turbulence", "representation", kSrcSep},
    {"application.sep.turbulence", "coupling", kSrcSep},
    {"application.sep.turbulence", "advection", kSrcSep},
    {"application.sep.turbulence", "reflection", kSrcSep},
    {"application.sep.turbulence", "cascade", kSrcSep},
    {"application.sep.turbulence", "shock_injection", kSrcSep},
    {"application.sep.turbulence", "inner_boundary", kSrcSep},
    {"application.sep.turbulence", "inner_boundary_value_si", kSrcSep},
    {"application.sep.turbulence", "outer_boundary", kSrcSep},
    {"application.sep.turbulence", "outer_boundary_value_si", kSrcSep},
    {"application.sep.turbulence", "reflection_coefficient", kSrcSep},
    {"application.sep.turbulence", "cascade_coefficient", kSrcSep},
    {"application.sep.turbulence", "perpendicular_correlation_length_m", kSrcSep},
    {"application.sep.turbulence", "electron_heating_fraction", kSrcSep},
    {"application.sep.turbulence", "spectral_k_min_per_m", kSrcSep},
    {"application.sep.turbulence", "spectral_k_max_per_m", kSrcSep},
    {"application.sep.turbulence", "spectral_bins", kSrcSep},
    {"application.sep.turbulence", "cfl_safety", kSrcSep},
    {"application.sep.turbulence", "operator_accuracy_safety", kSrcSep},
    {"application.sep.turbulence", "maximum_source_fraction", kSrcSep},
    {"application.sep.turbulence", "maximum_cascade_fraction", kSrcSep},
    {"application.sep.turbulence", "minimum_substep_s", kSrcSep},
    {"application.sep.turbulence", "maximum_substeps", kSrcSep},
    {"application.sep.turbulence", "conservation_relative_tolerance", kSrcSep},
    {"application.sep.turbulence", "shock_source_efficiency", kSrcSep},
    {"application.sep.turbulence", "shock_source_plus_fraction", kSrcSep},
    {"application.sep.turbulence", "diagnostic_cadence_steps", kSrcSep},

    // ---- mesh ----------------------------------------------------------------
    {"mesh", "global_cell_size_m", kBoth},
    {"mesh", "minimum_cell_size_m", kBoth},
    {"mesh", "maximum_level", kBoth},
    {"mesh", "cells_per_block_edge", kSrcSep3d},
    {"mesh", "memory_budget_bytes", kSrcSep3d},
    {"mesh", "block_overhead_bytes", kSrcSep3d},
    {"mesh.domain", "inner_radius_m", kBoth},
    {"mesh.domain", "outer_radius_m", kBoth},
    {"mesh.domain", "preset", kSrcSep3d},
    {"mesh.domain", "inner_boundary", kSrcSep3d},
    {"mesh.domain", "outer_radius_mode", kSrcSep3d},
    {"mesh.domain", "box_geometry", kSrcSep3d},
    {"mesh.domain", "corner_direction_x", kSrcSep3d},
    {"mesh.domain", "corner_direction_y", kSrcSep3d},
    {"mesh.domain", "corner_direction_z", kSrcSep3d},
    {"mesh.domain", "corner_margin_m", kSrcSep3d},
    {"mesh.domain", "outer_boundary", kSrcSep3d},
    {"mesh.domain", "coordinate_frame", kSrcSep3d},
    {"mesh.domain", "origin_x_m", kSrcSep3d},
    {"mesh.domain", "origin_y_m", kSrcSep3d},
    {"mesh.domain", "origin_z_m", kSrcSep3d},
    {"mesh.field_line", "origin_x_m", kBoth},
    {"mesh.field_line", "origin_y_m", kBoth},
    {"mesh.field_line", "origin_z_m", kBoth},
    {"mesh.field_line", "initial_x_m", kBoth},
    {"mesh.field_line", "initial_y_m", kBoth},
    {"mesh.field_line", "initial_z_m", kBoth},
    {"mesh.field_line", "length_m", kBoth},
    {"mesh.field_line", "point_count", kBoth},
    {"mesh.field_line", "start_mode", kSrcSep3d},
    {"mesh.field_line", "end_radius_m", kSrcSep3d},
    {"mesh.refinement.solar", "enabled", kBoth},
    {"mesh.refinement.solar", "surface_cell_size_m", kBoth},
    {"mesh.refinement.solar", "transition_outer_radius_m", kBoth},
    {"mesh.refinement.solar", "profile", kBoth},
    {"mesh.refinement.solar", "exponent", kBoth},
    {"mesh.refinement.solar", "anchor", kSrcSep3d},
    {"mesh.refinement.tube", "enabled", kBoth},
    {"mesh.refinement.tube", "reference_radius_m", kBoth},
    {"mesh.refinement.tube", "radius_at_reference_m", kBoth},
    {"mesh.refinement.tube", "radius_mode", kBoth},
    {"mesh.refinement.tube", "center_cell_size_m", kBoth},
    {"mesh.refinement.tube", "transverse_profile", kBoth},
    {"mesh.refinement.tube", "transverse_exponent", kBoth},
    {"mesh.refinement.tube", "source_longitude_rad", kSrcSep3d},
    {"mesh.refinement.tube", "source_colatitude_rad", kSrcSep3d},
    {"mesh.active_region", "mode", kSrcSep3d},
    {"mesh.active_region", "reference_radius_m", kSrcSep3d},
    {"mesh.active_region", "radius_at_reference_m", kSrcSep3d},
    {"mesh.active_region", "radius_mode", kSrcSep3d},
    {"mesh.active_region", "buffer_blocks", kSrcSep3d},
    {"mesh.active_region", "solar_sphere_radius_m", kSrcSep3d},
    {"mesh.memory", "base_cell_bytes", kSrcSep3d},
    {"mesh.memory", "base_node_bytes", kSrcSep3d},
    {"mesh.memory", "block_structure_bytes", kSrcSep3d},
    {"mesh.memory", "communication_bytes_per_block", kSrcSep3d},
    {"mesh.memory", "particle_bytes", kSrcSep3d},
    {"mesh.memory", "particles_per_cell", kSrcSep3d},
    {"mesh.memory", "halo_fraction", kSrcSep3d},
    {"mesh.memory", "safety_margin_fraction", kSrcSep3d},

    // ---- physics ---------------------------------------------------------------
    // provider: srcSEP implements swcme (the standalone SWCME StepState along
    // the Parker field line); srcSEP3D implements all values.
    {"physics.background", "provider", kBoth},
    {"physics.background", "model_id", kSrcSep3d},
    // Analytic Parker record.  The first seven keys are owned by
    // physics: swcme when that block is present (D4, SwcmeOwnedParkerKey()).
    {"physics.background.parker", "reference_radius_m", kSrcSep3d},
    {"physics.background.parker", "radial_field_at_reference_t", kSrcSep3d},
    {"physics.background.parker", "solar_rotation_rate_rad_per_s", kSrcSep3d},
    {"physics.background.parker", "solar_wind_speed_m_per_s", kSrcSep3d},
    {"physics.background.parker", "magnetic_polarity", kSrcSep3d},
    {"physics.background.parker", "number_density_at_one_au_m3", kSrcSep3d},
    {"physics.background.parker", "temperature_k", kSrcSep3d},
    {"physics.background.parker", "validity_cadence_s", kSrcSep3d},
    // authority: srcSEP implements swcme; srcSEP3D none | swcme.
    {"physics.shock", "authority", kBoth},

    // ---- source: shock (owns the particle source, D4) --------------------------
    {"source.shock", "enabled", kBoth},
    {"source.shock", "samples_per_step", kBoth},
    {"source.shock", "minimum_energy_j", kBoth},
    {"source.shock", "maximum_energy_j", kBoth},
    {"source.shock", "injection_efficiency", kBoth},
    {"source.shock", "spectrum_reference_energy_j", kBoth},
    {"source.shock", "spectrum_particle_mass_kg", kBoth},
    {"source.shock", "spectrum_charge_number", kBoth},
    {"source.shock", "spectrum_normalization", kBoth},
    {"source.shock", "spectrum_reference_intensity_per_m2_s_sr_j", kBoth},
    {"source.shock", "physical_particle_rate_per_s", kSrcSep3d},
    {"source.shock", "spectrum_model", kSrcSep3d},
    {"source.shock", "phase_space_power_index", kSrcSep3d},
    {"source.shock", "normalization", kSrcSep3d},
    {"source.shock", "normalization_radius_m", kSrcSep3d},
    {"source.shock", "weight_model", kSrcSep3d},
    {"source.shock", "maximum_events_per_species_per_step", kSrcSep3d},
    // Explicit particle-source model (srcSEP3D): reduced-shock-front is the
    // only model; swcme is rejected as not ready (user decision 2026-10-12).
    {"source.shock", "model", kSrcSep3d},

    // ---- species ----------------------------------------------------------------
    // The common statistical weight of every compiled AMPS species.
    {"species", "macroparticle_weight", kBoth},

    // ---- observer: ID (decision D5) ----------------------------------------------
    {"observer.*", "kind", kBoth},
    {"observer.*", "minimum_energy_j", kBoth},
    {"observer.*", "maximum_energy_j", kBoth},
    {"observer.*", "energy_bins", kBoth},
    {"observer.*", "energy_spacing", kBoth},
    {"observer.*", "pitch_angle_bins", kBoth},
    {"observer.*", "heliocentric_radius_m", kSrcSep},
    {"observer.*", "normalization", kSrcSep3d},
    {"observer.*", "position_x_m", kSrcSep3d},
    {"observer.*", "position_y_m", kSrcSep3d},
    {"observer.*", "position_z_m", kSrcSep3d},
    {"observer.*", "follows_trajectory", kSrcSep3d},
    {"observer.*", "velocity_x_m_per_s", kSrcSep3d},
    {"observer.*", "velocity_y_m_per_s", kSrcSep3d},
    {"observer.*", "velocity_z_m_per_s", kSrcSep3d},
    {"observer.*", "collection_radius_m", kSrcSep3d},
    {"observer.*", "shell_radius_m", kSrcSep3d},
    {"observer.*", "cadence_s", kSrcSep3d},
    {"observer.*", "minimum_mu", kSrcSep3d},
    {"observer.*", "maximum_mu", kSrcSep3d},
    {"observer.*", "species", kSrcSep3d},
    {"observer.*", "products", kSrcSep3d},

    // ---- output -----------------------------------------------------------------
    {"output", "initialization_mesh_tecplot_file", kBoth},
    {"output", "initialization_field_line_tecplot_file", kBoth},
    {"output", "initialization_data_tecplot_file", kBoth},
    {"output", "cadence_steps", kSrcSep3d},
    {"output", "checkpoint_cadence_steps", kSrcSep3d},
    {"output", "directory", kSrcSep3d},
    {"output", "prefix", kSrcSep3d},

    // ---- sampling: points (srcSEP3D) ---------------------------------------------
    // products: comma list of density, flux, spectrum, fluence.  point_file:
    // rows "x y z [name]" [m], relative to the file holding the block.
    {"sampling.points", "products", kSrcSep3d},
    {"sampling.points", "point_file", kSrcSep3d},
    // Channel energies are kinetic energies per particle in MeV: the one
    // deliberate exception to the SI key suffixes (user decision 2026-10-11).
    {"sampling.points.species.*", "minimum_energy_mev", kSrcSep3d},
    {"sampling.points.species.*", "maximum_energy_mev", kSrcSep3d},
    {"sampling.points.species.*", "energy_intervals", kSrcSep3d},
    {"sampling.points.species.*", "energy_scale", kSrcSep3d},
    // Point position [m], AMPS heliocentric Cartesian frame.
    {"sampling.points.point.*", "position_x_m", kSrcSep3d},
    {"sampling.points.point.*", "position_y_m", kSrcSep3d},
    {"sampling.points.point.*", "position_z_m", kSrcSep3d},
};

// ---------------------------------------------------------------------------
// Lookup helpers
// ---------------------------------------------------------------------------

// Qualified path of a block: the parent path joined with "kind" or
// "kind.name" ('#begin refinement: solar' inside mesh -> mesh.refinement.solar).
inline std::string QualifiedPath(const std::string& parent,
                                 const std::string& kind,
                                 const std::string& name) {
  std::string own = kind;
  if (!name.empty()) own += "." + name;
  return parent.empty() ? own : parent + "." + own;
}

// Table row of a block, or nullptr when the block is not part of the shared
// layout at all (an error in both applications).
inline const BlockRule* FindBlock(const std::string& parent,
                                  const std::string& kind,
                                  const std::string& name) {
  for (const BlockRule& rule : kBlocks) {
    if (parent != rule.parent || kind != rule.kind) continue;
    const std::string ruleName = rule.name;
    if (ruleName == "*" ? !name.empty() : ruleName == name) return &rule;
  }
  return nullptr;
}

// The table's block path for a qualified block path: every instance of a
// table block named "*" maps to its wildcard row -- 'observer.earth' to
// "observer.*", 'sampling.points.species.H_PLUS' to
// "sampling.points.species.*", 'sampling.points.point.earth' to
// "sampling.points.point.*".
inline std::string TableBlockPath(const std::string& blockPath) {
  static const char* const kInstanceParents[] = {
      "observer.", "sampling.points.species.", "sampling.points.point."};
  for (const char* parent : kInstanceParents) {
    const std::string prefix = parent;
    if (blockPath.size() > prefix.size() && blockPath.rfind(prefix, 0) == 0 &&
        blockPath.find('.', prefix.size()) == std::string::npos)
      return prefix + "*";
  }
  return blockPath;
}

// Table row of an application-owned key, or nullptr for an unknown key.
inline const KeyRule* FindKey(const std::string& blockPath,
                              const std::string& key) {
  const std::string table = TableBlockPath(blockPath);
  for (const KeyRule& rule : kKeys)
    if (table == rule.block && key == rule.key) return &rule;
  return nullptr;
}

// Human-readable list of the blocks allowed under 'parent', used in the
// diagnostic for a misplaced or unknown block.
inline std::string AllowedChildren(const std::string& parent) {
  std::string out;
  for (const BlockRule& rule : kBlocks) {
    if (parent != rule.parent) continue;
    std::string label = rule.kind;
    const std::string name = rule.name;
    if (name == "*") label += ": NAME";
    else if (!name.empty()) label += ": " + name;
    out += (out.empty() ? "" : ", ") + label;
  }
  return out.empty() ? std::string("none") : out;
}

// Typed D3 diagnostic for a setting that is part of the shared vocabulary but
// is not implemented by the application reading the file.
inline std::string NotSupportedMessage(const std::string& what,
                                       unsigned reader) {
  const unsigned other = reader == kSrcSep ? kSrcSep3d : kSrcSep;
  return what + " is a " + ApplicationName(other) +
         " setting that is not supported by " + ApplicationName(reader) +
         " (shared SEP input schema 6, src/INPUT_FILE.md section 5)";
}

// ---------------------------------------------------------------------------
// Decision D4: single source of truth
// ---------------------------------------------------------------------------

// Keys of physics.background.parker that physics: swcme owns when present.
// Returns the owning SWCME key for the diagnostic, or nullptr when 'key' is
// not SWCME-owned.  The analytic Parker provider then receives the canonical
// SWCME values: wind speed, rotation rate, the radial field derived from the
// total |B| at 1 AU (expressed at 1 AU), density at 1 AU, proton temperature
// and the radial polarity.
inline const char* SwcmeOwnedParkerKey(const std::string& key) {
  struct Owned { const char* parker; const char* swcme; };
  static const Owned kOwned[] = {
      {"reference_radius_m", "ambient.magnetic_field_1au (the 1-AU reference)"},
      {"radial_field_at_reference_t", "ambient.magnetic_field_1au"},
      {"solar_rotation_rate_rad_per_s", "the SWCME solar rotation rate"},
      {"solar_wind_speed_m_per_s", "ambient.wind_speed"},
      {"magnetic_polarity", "parker.radial_polarity"},
      {"number_density_at_one_au_m3", "ambient.density_1au"},
      {"temperature_k", "ambient.proton_temperature"}};
  for (const Owned& owned : kOwned)
    if (key == owned.parker) return owned.swcme;
  return nullptr;
}

// SWCME 'source.*' keys that source: shock owns (D4).  Returns the source:
// shock key that replaces 'swcmeKey', or nullptr when it is not a source key.
inline const char* SourceShockKeyForSwcmeKey(const std::string& swcmeKey) {
  struct Owned { const char* swcme; const char* source; };
  static const Owned kOwned[] = {
      {"source.particle_mass", "spectrum_particle_mass_kg"},
      {"source.charge_number", "spectrum_charge_number"},
      {"source.energy_min", "minimum_energy_j"},
      {"source.energy_max", "maximum_energy_j"},
      {"source.reference_energy", "spectrum_reference_energy_j"},
      {"source.injection_efficiency", "injection_efficiency"},
      {"source.normalization", "spectrum_normalization"},
      {"source.reference_intensity_si",
       "spectrum_reference_intensity_per_m2_s_sr_j"}};
  for (const Owned& owned : kOwned)
    if (swcmeKey == owned.swcme) return owned.source;
  return nullptr;
}

// SWCME keys that only the three-dimensional resolver (swcme3d_input.hpp,
// geometry_keys()) defines: CME/shock shape and direction, the rotation axis
// and rate, and the source-surface discretization.  srcSEP's SWCME1D model is
// spherically symmetric about the field line and has a fixed Carrington
// rotation rate, so it rejects these keys with a D3 diagnostic (the rotation
// rate is accepted when it equals the fixed SWCME1D value).  The list is
// cross-checked against swcme::input3d::detail::geometry_keys() by srcSEP3D
// CFG3D18.
inline const std::vector<std::string>& SwcmeThreeDimensionalKeys() {
  static const std::vector<std::string> keys = {
      "geometry.shape", "geometry.axis_ratio_y", "geometry.axis_ratio_z",
      "geometry.half_width_rad", "geometry.cme_direction_x",
      "geometry.cme_direction_y", "geometry.cme_direction_z",
      "geometry.solar_rotation_axis_x", "geometry.solar_rotation_axis_y",
      "geometry.solar_rotation_axis_z",
      "parker.solar_rotation_rate_rad_per_s",
      "surface.theta_intervals", "surface.phi_points"};
  return keys;
}

inline bool IsSwcmeThreeDimensionalKey(const std::string& key) {
  for (const std::string& candidate : SwcmeThreeDimensionalKeys())
    if (candidate == key) return true;
  return false;
}

// The physical source description written in source: shock.  Both
// applications fill this record while walking the block and hand it to
// SynthesizeSwcmeSource() when physics: swcme is present.
struct SwcmeSourceSettings {
  double minimumEnergyJ = 0.0;              // kinetic energy [J]
  double maximumEnergyJ = 0.0;              // kinetic energy [J]
  double referenceEnergyJ = 0.0;            // spectrum normalization energy [J]
  double particleMassKg = 0.0;              // reference-spectrum mass [kg]
  long chargeNumber = 0;                    // reference-spectrum charge number
  double injectionEfficiency = 0.0;         // dimensionless, [0, 1]
  std::string normalization;                // relative-only | reference-differential-intensity
  double referenceIntensitySi = 0.0;        // [m^-2 s^-1 sr^-1 J^-1]
};

// One synthesized SWCME assignment ("key = value").  The caller attaches the
// input line of the source: shock key that supplies it, so the SWCME
// resolver's diagnostics still point at the line a user has to edit.
struct SynthesizedAssignment {
  std::string key;
  std::string value;
};

// Shortest decimal text that round-trips a double exactly: the smallest
// precision p in [1, 17] whose "%.*g" text parses back to the identical
// binary value (17 digits always suffice for IEEE binary64).  A synthesized
// SWCME value is therefore bit-identical to the SI value the application
// itself uses, and a written "1.602176634e-14" stays readable as such.
inline std::string ExactDecimal(double value) {
  char buffer[64];
  for (int precision = 1; precision <= 17; ++precision) {
    std::snprintf(buffer, sizeof(buffer), "%.*g", precision, value);
    if (std::strtod(buffer, nullptr) == value) break;
  }
  return buffer;
}

// Translate the SI source: shock description into the canonical SWCME
// 'source.*' assignments.  Units are named in each value ("J", "kg", "si"),
// which the SWCME resolvers convert with their own unit tables
// (swcme1d_input.hpp: J -> MeV through MEV_TO_J).  The normalization
// spellings map relative-only -> relative_only and
// reference-differential-intensity -> reference_differential_intensity;
// an unknown spelling is passed through unchanged so that the SWCME
// resolver reports it with its own vocabulary.
inline std::vector<SynthesizedAssignment> SynthesizeSwcmeSource(
    const SwcmeSourceSettings& source) {
  std::string normalization = source.normalization;
  if (normalization == "relative-only") normalization = "relative_only";
  else if (normalization == "reference-differential-intensity")
    normalization = "reference_differential_intensity";
  return {
      {"source.particle_mass", ExactDecimal(source.particleMassKg) + " kg"},
      {"source.charge_number", std::to_string(source.chargeNumber)},
      {"source.energy_min", ExactDecimal(source.minimumEnergyJ) + " J"},
      {"source.energy_max", ExactDecimal(source.maximumEnergyJ) + " J"},
      {"source.reference_energy", ExactDecimal(source.referenceEnergyJ) + " J"},
      {"source.injection_efficiency", ExactDecimal(source.injectionEfficiency)},
      {"source.normalization", normalization},
      {"source.reference_intensity_si",
       ExactDecimal(source.referenceIntensitySi) + " si"}};
}

// The source: shock keys that feed SynthesizeSwcmeSource(): required when
// physics: swcme is present, rejected otherwise (a reduced shock-front run
// has no SWCME source and must not carry inert spectrum settings).
inline const std::vector<std::string>& SwcmeSourceKeys() {
  static const std::vector<std::string> keys = {
      "spectrum_reference_energy_j", "spectrum_particle_mass_kg",
      "spectrum_charge_number", "spectrum_normalization",
      "spectrum_reference_intensity_per_m2_s_sr_j"};
  return keys;
}

}  // namespace InputLayout
}  // namespace SEP

#endif  // SEP_COMMON_SEP_INPUT_LAYOUT_H
