#include "test_framework.hpp"

#include <swcme_sep_interface.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// V6 validates the exact serialized boundary available to an external AMPS
// build.  The repository intentionally does not contain AMPS itself, so the
// transport calculation below is a validation-owned consumer smoke test, not
// a claim of coupled-model flux skill.  EVT01 requires real coupled artifacts
// and provenance before a RELEASE_VALIDATION campaign can report V6 PASS.

struct ParsedAMPSSource {
  bool active=false;
  double time_s=0.0;
  double source_weight=0.0;
  double compression=1.0;
  double fast_mach=0.0;
  double q_phase_space=std::numeric_limits<double>::quiet_NaN();
};

struct TransportControls {
  std::uint64_t random_seed=0x56365f414d50535fULL;
  double time_step_s=600.0;
  std::size_t time_steps=288;
  std::size_t macro_particles_per_source=4096;
};

struct AttributionRow {
  const char* name;
  bool source_dependent;
  bool three_dimensional;
  bool perpendicular_diffusion;
  std::vector<double> profile;
};

std::vector<std::string> split_csv(const std::string& record) {
  std::vector<std::string> fields;
  std::size_t begin=0;
  while (true) {
    const std::size_t comma=record.find(',',begin);
    fields.push_back(record.substr(begin,comma-begin));
    if (comma==std::string::npos) break;
    begin=comma+1;
  }
  return fields;
}

// This parser is intentionally independent of SEPSourceState and consumes
// only the public CSV schema an AMPS integration sees.  Column positions are
// checked against the complete header so a reordered/omitted field cannot be
// hidden by compiling producer and consumer against one C++ record type.
ParsedAMPSSource parse_amps_source(const std::string& header,
                                   const std::string& record) {
  const auto names=split_csv(header);
  const auto values=split_csv(record);
  if (names.size()!=37 || values.size()!=names.size() ||
      names[0]!="status" || names[1]!="active" ||
      names[5]!="time_s" || names[16]!="relative_source_weight_per_area" ||
      names[18]!="compression" || names[20]!="fast_mach" ||
      names[27]!="q_phase_space") {
    throw std::runtime_error("unsupported SWCME/AMPS source schema");
  }
  if (values[0]!="OK")
    throw std::runtime_error("AMPS source record carries failure status");
  ParsedAMPSSource out;
  out.active=values[1]=="1";
  out.time_s=std::stod(values[5]);
  out.source_weight=std::stod(values[16]);
  out.compression=std::stod(values[18]);
  out.fast_mach=std::stod(values[20]);
  out.q_phase_space=values[27]=="NA"
      ? std::numeric_limits<double>::quiet_NaN() : std::stod(values[27]);
  return out;
}

// FNV-1a is used as a compact regression fingerprint, not as a security hash.
// The comparison also retains the complete bytes, so a collision could not
// make unequal source records pass the identity assertion.
std::uint64_t record_hash(const std::string& bytes) {
  std::uint64_t hash=1469598103934665603ULL;
  for (unsigned char byte:bytes) {
    hash^=byte;
    hash*=1099511628211ULL;
  }
  return hash;
}

std::uint64_t next_random(std::uint64_t& state) {
  state^=state>>12;
  state^=state<<25;
  state^=state>>27;
  return state*2685821657736338717ULL;
}

double uniform01(std::uint64_t& state) {
  return static_cast<double>(next_random(state)>>11)*0x1.0p-53;
}

// A small deterministic transport consumer verifies that serialized source
// records can drive both field-aligned and 3-D control paths.  It deliberately
// contains no SWCME physics and is not an AMPS accuracy reference.  Identical
// seeds make 1-D and 3-D/no-perpendicular-diffusion paths bitwise comparable;
// enabling perpendicular diffusion consumes the same controls but introduces
// a documented cross-field survival factor.
std::vector<double> transport_smoke(
    const std::vector<ParsedAMPSSource>& sources,
    const TransportControls& controls,
    bool source_dependent,
    bool three_dimensional,
    bool perpendicular_diffusion) {
  std::vector<double> profile(controls.time_steps,0.0);
  std::uint64_t random_state=controls.random_seed;
  for (const ParsedAMPSSource& source:sources) {
    if (!source.active) continue;
    const double source_scale=source_dependent
        ? source.source_weight*(source.compression-1.0)*source.fast_mach
        : 1.0;
    for (std::size_t particle=0;
         particle<controls.macro_particles_per_source;++particle) {
      // One field-aligned random variate is consumed in every case. The
      // three-dimensional flag alone cannot alter a no-perpendicular run.
      const double parallel_scatter=0.85+0.30*uniform01(random_state);
      double cross_field_factor=1.0;
      if (three_dimensional && perpendicular_diffusion) {
        const double cross_field_displacement=2.0*uniform01(random_state)-1.0;
        cross_field_factor=std::exp(-1.5*cross_field_displacement*
                                    cross_field_displacement);
      }
      const std::size_t injection=static_cast<std::size_t>(
          std::max(0.0,std::floor(source.time_s/controls.time_step_s)));
      for (std::size_t step=injection;step<profile.size();++step) {
        const double age=(step-injection)*controls.time_step_s;
        profile[step]+=source_scale*parallel_scatter*cross_field_factor*
                       std::exp(-age/(12.0*3600.0));
      }
    }
  }
  const double normalization=static_cast<double>(
      controls.macro_particles_per_source);
  for (double& value:profile) value/=normalization;
  return profile;
}

bool exact_profile(const std::vector<double>& left,
                   const std::vector<double>& right) {
  return left==right;
}

bool different_profile(const std::vector<double>& left,
                       const std::vector<double>& right) {
  if (left.size()!=right.size()) return true;
  for (std::size_t i=0;i<left.size();++i)
    if (left[i]!=right[i]) return true;
  return false;
}

bool physically_plausible_profile(const std::vector<double>& profile) {
  return !profile.empty() &&
      std::all_of(profile.begin(),profile.end(),[](double value) {
        return std::isfinite(value) && value>=0.0;
      }) && *std::max_element(profile.begin(),profile.end())>0.0;
}

void make_equivalent_params(swcme1d::Params& one,swcme3d::Params& three) {
  one.V_sw_kms=410.0;
  one.n1AU_cm3=5.5;
  one.B1AU_nT=5.2;
  one.T_K=1.25e5;
  // Spell out the common closure as part of the coupling fixture.  The 1-D
  // and 3-D parameter types deliberately have independent defaults, so a V6
  // equivalence case must copy every source-affecting control explicitly.
  one.gamma_ad=5.0/3.0;
  one.sin_theta=1.0;
  one.kinematics_mode=swcme::kinematics::Mode::DBM;
  one.r0_Rs=20.0;
  one.V0_sh_kms=1450.0;
  one.Gamma_kmInv=5.0e-8;
  one.region_mode=swcme::regions::Mode::ShockOnly;
  one.shock_acceleration_mode=swcme::acceleration::Mode::Source;
  one.relative_source_weight_per_area=1.0;

  three.V_sw_kms=one.V_sw_kms;
  three.n1AU_cm3=one.n1AU_cm3;
  three.B1AU_nT=one.B1AU_nT;
  three.T_K=one.T_K;
  three.gamma_ad=one.gamma_ad;
  three.sin_theta=one.sin_theta;
  three.kinematics_mode=one.kinematics_mode;
  three.r0_Rs=one.r0_Rs;
  three.V0_sh_kms=one.V0_sh_kms;
  three.Gamma_kmInv=one.Gamma_kmInv;
  three.region_mode=one.region_mode;
  three.shock_acceleration_mode=one.shock_acceleration_mode;
  three.relative_source_weight_per_area=one.relative_source_weight_per_area;
  three.shape=swcme3d::ShockShape::Sphere;
  three.cme_dir[0]=1.0; three.cme_dir[1]=0.0; three.cme_dir[2]=0.0;
  three.solar_rotation_axis[0]=0.0;
  three.solar_rotation_axis[1]=0.0;
  three.solar_rotation_axis[2]=1.0;
}

}  // namespace

void test_v6(swcme_test::Context& context) {
  std::cout << "V6 SEP-facing AMPS integration\n"
            << "  fixture_kind=COUPLING_SMOKE_TEST "
               "coupled_flux_skill=NOT_EVALUATED\n";

  swcme1d::Params one_params;
  swcme3d::Params three_params;
  make_equivalent_params(one_params,three_params);
  swcme::sep::SpectrumConfig spectrum;
  spectrum.kinetic_energy_min_MeV=2.0;
  spectrum.kinetic_energy_max_MeV=500.0;
  spectrum.reference_energy_MeV=20.0;

  const swcme::sep::Interface1D one(one_params,spectrum);
  const swcme::sep::Interface3D three(three_params,spectrum);
  const std::array<double,5> times_s{{2.0*3600.0,6.0*3600.0,12.0*3600.0,
                                      24.0*3600.0,36.0*3600.0}};
  std::string standalone_one_history_bytes;
  std::string adapter_one_history_bytes;
  std::string standalone_three_history_bytes;
  std::string adapter_three_history_bytes;
  std::vector<ParsedAMPSSource> parsed_one_sources;
  std::vector<ParsedAMPSSource> parsed_three_sources;
  double previous_radius=0.0;

  for (double time_s:times_s) {
    const auto one_step=one.prepare(time_s);
    const auto three_step=three.prepare(time_s);
    swcme::sep::SEPSourceState one_source,three_source;
    context.expect_true(one.source_at_shock(one_step,one_source).ok(),
                        "1-D AMPS-facing source record succeeds");
    context.expect_true(three.source_at_direction(
        three_step,{{1.0,0.0,0.0}},three_source).ok(),
        "3-D AMPS-facing source record succeeds");

    // Build the standalone record directly from the public model result. This
    // proves that crossing Interface1D adds no hidden normalization or units.
    swcme::acceleration::ShockAccelerationState acceleration;
    context.expect_true(one.model().shock_acceleration_state_checked(
        one_step,acceleration).ok(),"standalone SWCME source state succeeds");
    swcme::sep::SEPSourceState standalone;
    context.expect_true(swcme::sep::make_source_state(
        acceleration,spectrum,false,false,0,0.0,standalone).ok(),
        "standalone state converts to public SEP record");
    context.expect_true(swcme::sep::attach_source_transport_context(
        one_step.common.solar_wind,acceleration.radius_m,
        one_step.common.solar_wind.reference_sin_theta,standalone).ok(),
        "standalone source transport context succeeds");

    // Repeat the standalone construction through the 3-D production model.
    // This makes record identity a true producer/boundary assertion for both
    // dimensional adapters instead of assuming their independent floating-
    // point paths must always round every diagnostic to the same final bit.
    swcme::acceleration::ShockAccelerationState acceleration_three;
    const double direction[3]={1.0,0.0,0.0};
    context.expect_true(three.model().shock_acceleration_state_checked(
        three_step,direction,acceleration_three).ok(),
        "standalone 3-D SWCME source state succeeds");
    swcme::sep::SEPSourceState standalone_three;
    context.expect_true(swcme::sep::make_source_state(
        acceleration_three,spectrum,false,false,0,0.0,standalone_three).ok(),
        "standalone 3-D state converts to public SEP record");
    const double acceleration_three_radius=std::hypot(
        acceleration_three.position_m[0],
        std::hypot(acceleration_three.position_m[1],
                   acceleration_three.position_m[2]));
    context.expect_true(swcme::sep::attach_source_transport_context(
        three_step.common.solar_wind,acceleration_three_radius,1.0,
        standalone_three).ok(),
        "standalone 3-D source transport context succeeds");

    const std::string standalone_record=
        swcme::sep::serialize_source_csv(standalone);
    const std::string standalone_three_record=
        swcme::sep::serialize_source_csv(standalone_three);
    const std::string one_record=swcme::sep::serialize_source_csv(one_source);
    const std::string three_record=swcme::sep::serialize_source_csv(three_source);
    context.expect_true(standalone_record==one_record,
                        "standalone 1-D and AMPS-boundary records are byte-identical");
    context.expect_true(standalone_three_record==three_record,
                        "standalone 3-D and AMPS-boundary records are byte-identical");
    context.expect_true(record_hash(standalone_record)==record_hash(one_record),
                        "standalone/adapter 1-D source-record hashes are identical");
    context.expect_true(record_hash(standalone_three_record)==record_hash(three_record),
                        "standalone/adapter 3-D source-record hashes are identical");

    const ParsedAMPSSource parsed=parse_amps_source(
        swcme::sep::source_csv_header(),one_record);
    const ParsedAMPSSource parsed_three=parse_amps_source(
        swcme::sep::source_csv_header(),three_record);
    context.expect_true(parsed.active,"source history remains physically active");
    context.expect_true(parsed.time_s==time_s,"serialized source time is exact");
    context.expect_true(parsed.compression>1.0 && parsed.compression<4.0,
                        "source compression is on physical fast-shock branch");
    context.expect_true(parsed.fast_mach>1.0,
                        "source history remains super-fast");
    context.expect_true(std::isfinite(parsed.q_phase_space) &&
                        parsed.q_phase_space>2.0,
                        "source history carries a physical DSA slope");
    // The dimensional models are independent implementations, so their full
    // records may differ by a final rounding bit in diagnostics that AMPS does
    // not use.  The transport inputs themselves must nevertheless agree to a
    // tight numerical tolerance before identical-control comparisons begin.
    const double source_scale=std::max(
        {std::abs(parsed.source_weight),std::abs(parsed_three.source_weight),1.0});
    const double compression_scale=std::max(
        {std::abs(parsed.compression),std::abs(parsed_three.compression),1.0});
    const double mach_scale=std::max(
        {std::abs(parsed.fast_mach),std::abs(parsed_three.fast_mach),1.0});
    const double slope_scale=std::max(
        {std::abs(parsed.q_phase_space),std::abs(parsed_three.q_phase_space),1.0});
    context.expect_true(parsed.active==parsed_three.active &&
                        parsed.time_s==parsed_three.time_s,
                        "1-D/3-D source activation and time controls agree");
    context.expect_near(parsed.source_weight,parsed_three.source_weight,
                        2.0e-14*source_scale,
                        "1-D/3-D source weights agree");
    context.expect_near(parsed.compression,parsed_three.compression,
                        2.0e-14*compression_scale,
                        "1-D/3-D compression ratios agree");
    context.expect_near(parsed.fast_mach,parsed_three.fast_mach,
                        2.0e-14*mach_scale,
                        "1-D/3-D fast Mach numbers agree");
    context.expect_near(parsed.q_phase_space,parsed_three.q_phase_space,
                        2.0e-14*slope_scale,
                        "1-D/3-D DSA slopes agree");
    context.expect_true(one_source.position_m[0]>previous_radius,
                        "source radius increases monotonically");
    previous_radius=one_source.position_m[0];
    parsed_one_sources.push_back(parsed);
    parsed_three_sources.push_back(parsed_three);
    standalone_one_history_bytes+=standalone_record+'\n';
    adapter_one_history_bytes+=one_record+'\n';
    standalone_three_history_bytes+=standalone_three_record+'\n';
    adapter_three_history_bytes+=three_record+'\n';
  }
  context.expect_true(standalone_one_history_bytes==adapter_one_history_bytes,
                      "complete standalone/adapter 1-D histories are byte-identical");
  context.expect_true(standalone_three_history_bytes==adapter_three_history_bytes,
                      "complete standalone/adapter 3-D histories are byte-identical");
  context.expect_true(record_hash(standalone_one_history_bytes)==
                      record_hash(adapter_one_history_bytes) &&
                      record_hash(standalone_three_history_bytes)==
                      record_hash(adapter_three_history_bytes),
                      "complete producer/boundary history hashes are identical");
  // Emit both pre-transport fingerprints as campaign evidence.  They are not
  // expected to equal one another because 1-D and 3-D retain independent
  // roundoff paths; each printed hash has already been proven identical to
  // its corresponding adapter-side history above.
  std::cout << "  boundary_history_hash_1d=0x" << std::hex
            << record_hash(standalone_one_history_bytes)
            << " boundary_history_hash_3d=0x"
            << record_hash(standalone_three_history_bytes) << std::dec << '\n';

  const TransportControls controls;
  std::vector<AttributionRow> matrix;
  matrix.push_back({"fixed_1d_no_perp",false,false,false,
      transport_smoke(parsed_one_sources,controls,false,false,false)});
  matrix.push_back({"fixed_3d_no_perp",false,true,false,
      transport_smoke(parsed_three_sources,controls,false,true,false)});
  matrix.push_back({"source_1d_no_perp",true,false,false,
      transport_smoke(parsed_one_sources,controls,true,false,false)});
  matrix.push_back({"source_3d_no_perp",true,true,false,
      transport_smoke(parsed_three_sources,controls,true,true,false)});
  matrix.push_back({"source_3d_perp",true,true,true,
      transport_smoke(parsed_three_sources,controls,true,true,true)});

  context.expect_true(exact_profile(matrix[0].profile,matrix[1].profile),
                      "fixed-source 1-D and 3-D/no-perp transport are exact");
  context.expect_true(exact_profile(matrix[2].profile,matrix[3].profile),
                      "source-dependent 1-D and 3-D/no-perp transport are exact");
  context.expect_true(different_profile(matrix[0].profile,matrix[2].profile),
                      "activating SWCME source dependence changes transport");
  context.expect_true(different_profile(matrix[3].profile,matrix[4].profile),
                      "perpendicular diffusion produces attributable 3-D change");
  for (const AttributionRow& row:matrix) {
    context.expect_true(physically_plausible_profile(row.profile),
                        std::string(row.name)+" profile is finite/nonnegative");
    std::cout << "  attribution=" << std::left << std::setw(22) << row.name
              << " final=" << std::scientific << std::setprecision(8)
              << row.profile.back() << " peak="
              << *std::max_element(row.profile.begin(),row.profile.end())
              << " seed=0x" << std::hex << controls.random_seed << std::dec
              << " dt_s=" << controls.time_step_s
              << " steps=" << controls.time_steps
              << " particles_per_source="
              << controls.macro_particles_per_source << '\n';
  }
}
