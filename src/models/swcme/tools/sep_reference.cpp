// Reference SWCME -> SEP source exporter used by the validation campaign.
//
// This is intentionally a small consumer of the public swcme_sep_interface.hpp
// contract.  It is not a second physics implementation.  Its purpose is to
// generate a human-readable time history that AMPS developers can compare with
// their direct call path before particle transport is enabled.

#include <swcme_sep_interface.hpp>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
  std::string output_path;
  bool print_manifest=false;
  double start_hours=0.0;
  double end_hours=48.0;
  double step_hours=1.0;
  double probe_time_hours=std::numeric_limits<double>::quiet_NaN();
  double v0_kms=std::numeric_limits<double>::quiet_NaN();
  double gamma_km_inv=std::numeric_limits<double>::quiet_NaN();
  double half_width_deg=std::numeric_limits<double>::quiet_NaN();
  double observer_lon_deg=0.0;
  double observer_lat_deg=0.0;
  std::vector<double> energies_MeV{1.0,10.0,100.0};
};

std::vector<double> parse_energies(const std::string& text) {
  std::vector<double> values;
  std::stringstream input(text);
  std::string token;
  while (std::getline(input,token,',')) {
    if (token.empty()) continue;
    values.push_back(std::stod(token));
  }
  if (values.empty()) throw std::invalid_argument("--energies requires at least one value");
  return values;
}

Options parse(int argc,char** argv) {
  Options o;
  for (int i=1;i<argc;++i) {
    const std::string a=argv[i];
    auto need=[&](const char* name)->const char* {
      if (i+1>=argc) throw std::invalid_argument(std::string(name)+" requires a value");
      return argv[++i];
    };
    if (a=="--output") o.output_path=need("--output");
    else if (a=="--start-hours") o.start_hours=std::stod(need("--start-hours"));
    else if (a=="--end-hours") o.end_hours=std::stod(need("--end-hours"));
    else if (a=="--step-hours") o.step_hours=std::stod(need("--step-hours"));
    else if (a=="--probe-time-hours") o.probe_time_hours=std::stod(need("--probe-time-hours"));
    else if (a=="--v0-kms") o.v0_kms=std::stod(need("--v0-kms"));
    else if (a=="--gamma-km-inv") o.gamma_km_inv=std::stod(need("--gamma-km-inv"));
    else if (a=="--half-width-deg") o.half_width_deg=std::stod(need("--half-width-deg"));
    else if (a=="--observer-lon-deg") o.observer_lon_deg=std::stod(need("--observer-lon-deg"));
    else if (a=="--observer-lat-deg") o.observer_lat_deg=std::stod(need("--observer-lat-deg"));
    else if (a=="--energies") o.energies_MeV=parse_energies(need("--energies"));
    else if (a=="--print-manifest") o.print_manifest=true;
    else if (a=="--help") {
      std::cout
        << "usage: sep_reference [options]\n"
        << "  --output FILE       write CSV to FILE (default stdout)\n"
        << "  --start-hours H     first launch-relative time (default 0)\n"
        << "  --end-hours H       last launch-relative time (default 48)\n"
        << "  --step-hours H      cadence (default 1)\n"
        << "  --probe-time-hours H  print one machine-readable source/cobpoint probe\n"
        << "  --v0-kms V          override CME apex reference speed\n"
        << "  --gamma-km-inv G    override DBM drag coefficient\n"
        << "  --half-width-deg D  override finite-SSE half width\n"
        << "  --observer-lon-deg D  1-AU observer longitude (default 0)\n"
        << "  --observer-lat-deg D  1-AU observer latitude (default 0)\n"
        << "  --energies E1,E2    relative spectrum sample energies in MeV\n"
        << "  --print-manifest    print resolved SWCME+SEP configuration and exit\n";
      std::exit(0);
    } else {
      throw std::invalid_argument("unknown option: "+a);
    }
  }
  if (!(o.step_hours>0.0) || o.end_hours<o.start_hours)
    throw std::invalid_argument("invalid time interval/cadence");
  return o;
}

} // namespace

int main(int argc,char** argv) {
  try {
    const Options options=parse(argc,argv);

    // Use the standardized science defaults: finite SSE, SHOCK_ONLY + SOURCE.
    // A 1-AU +X observer is only a reference consumer; event configurations
    // should construct the same interface with their own validated Params.
    swcme3d::Params params;
    if (std::isfinite(options.v0_kms)) params.V0_sh_kms=options.v0_kms;
    if (std::isfinite(options.gamma_km_inv)) params.Gamma_kmInv=options.gamma_km_inv;
    if (std::isfinite(options.half_width_deg))
      params.half_width_rad=options.half_width_deg*swcme::constants::PI/180.0;
    swcme::sep::SpectrumConfig spectrum;
    spectrum.kinetic_energy_min_MeV=1.0;
    spectrum.kinetic_energy_max_MeV=1000.0;
    spectrum.reference_energy_MeV=10.0;
    spectrum.normalization=swcme::sep::NormalizationMode::RelativeOnly;
    const swcme::sep::Interface3D interface(params,spectrum);

    if (options.print_manifest) {
      std::cout << interface.resolved_manifest();
      return 0;
    }

    const double lon=options.observer_lon_deg*swcme::constants::PI/180.0;
    const double lat=options.observer_lat_deg*swcme::constants::PI/180.0;
    const std::array<double,3> observer{{
        swcme::constants::AU_M*std::cos(lat)*std::cos(lon),
        swcme::constants::AU_M*std::cos(lat)*std::sin(lon),
        swcme::constants::AU_M*std::sin(lat)}};

    // Probe mode is used by campaign parameter sweeps.  It emits simple
    // KEY=value lines so the Python runner can capture a requested metric with
    // a regular expression without parsing the full reference CSV.
    if (std::isfinite(options.probe_time_hours)) {
      const double time_s=options.probe_time_hours*3600.0;
      const auto step=interface.prepare(time_s);
      swcme::sep::SEPSourceState source;
      swcme3d::ConnectivityState connectivity;
      const swcme::ModelStatus status=
          interface.source_at_observer_cobpoint(step,observer,source,&connectivity);
      std::cout << std::setprecision(17) << std::scientific;
      std::cout << "STATUS=" << swcme::status_code_name(status.code) << '\n';
      std::cout << "CONNECTED=" << (connectivity.connected?1:0) << '\n';
      std::cout << "SOURCE_ACTIVE=" << (source.active?1:0) << '\n';
      std::cout << "SHOCK_RADIUS_AU=" << step.common.apex.radius_m/swcme::constants::AU_M << '\n';
      if (connectivity.connected && connectivity.selected_root<connectivity.roots.size())
        std::cout << "PATH_LENGTH_AU="
                  << connectivity.roots[connectivity.selected_root].path_length_m/
                     swcme::constants::AU_M << '\n';
      else
        std::cout << "PATH_LENGTH_AU=NA\n";
      auto emit_probe=[&](const char* key,double value) {
        std::cout << key << '=';
        if (std::isfinite(value)) std::cout << value; else std::cout << "NA";
        std::cout << '\n';
      };
      emit_probe("FAST_MACH",source.fast_mach);
      emit_probe("COMPRESSION",source.compression);
      emit_probe("Q_PHASE_SPACE",source.q_phase_space);
      return status.failure()?1:0;
    }

    // Time-history mode opens its output only after the probe-mode early
    // return.  This keeps the two wire formats disjoint and also prevents a
    // probe invoked with an accidental --output option from creating an empty
    // CSV file.
    std::ofstream file;
    std::ostream* output=&std::cout;
    if (!options.output_path.empty()) {
      file.open(options.output_path);
      if (!file) throw std::runtime_error("cannot open output: "+options.output_path);
      output=&file;
    }

    // Probe output is intentionally KEY=value text, whereas histories are CSV.
    *output << std::setprecision(17) << std::scientific;
    *output << "time_s,connected,source_active,cobpoint_x_m,cobpoint_y_m,cobpoint_z_m,"
               "path_length_m,compression,theta_Bn_rad,fast_mach,Vsh_n_m_s,"
               "upstream_density_m3,upstream_B_T,q_phase_space";
    for (double E : options.energies_MeV) *output << ",relative_J_" << E << "MeV";
    *output << '\n';

    const double start_s=options.start_hours*3600.0;
    const double end_s=options.end_hours*3600.0;
    const double step_s=options.step_hours*3600.0;
    for (double time_s=start_s; time_s<=end_s+0.5*step_s; time_s+=step_s) {
      const auto step=interface.prepare(time_s);
      swcme::sep::SEPSourceState source;
      swcme3d::ConnectivityState connectivity;
      const swcme::ModelStatus status=
          interface.source_at_observer_cobpoint(step,observer,source,&connectivity);
      const bool connected=status.ok() && connectivity.connected;
      const swcme3d::ConnectivityRoot* root=nullptr;
      if (connected && connectivity.selected_root<connectivity.roots.size())
        root=&connectivity.roots[connectivity.selected_root];

      *output << time_s << ',' << (connected?1:0) << ',' << (source.active?1:0) << ',';
      if (root) {
        *output << root->position_m[0] << ',' << root->position_m[1] << ','
                << root->position_m[2] << ',' << root->path_length_m << ',';
      } else {
        *output << "NA,NA,NA,NA,";
      }
      auto emit=[&](double value) {
        if (std::isfinite(value)) *output << value;
        else *output << "NA";
      };
      emit(source.compression); *output << ',';
      emit(source.theta_Bn_rad); *output << ',';
      emit(source.fast_mach); *output << ',';
      emit(source.normal_speed_m_s); *output << ',';
      emit(source.upstream_density_m3); *output << ',';
      emit(source.upstream_B_T); *output << ',';
      emit(source.q_phase_space);
      for (double E : options.energies_MeV) {
        *output << ',';
        double shape=0.0;
        const swcme::ModelStatus spectrum_status=
            swcme::sep::relative_intensity_shape(source,E,shape);
        if (spectrum_status.ok()) emit(shape); else *output << "NA";
      }
      *output << '\n';
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "sep_reference: " << error.what() << '\n';
    return 2;
  }
}
