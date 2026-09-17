// VP03 production API driver.
//
// One public model is constructed per radial bin because the 1-D ray latitude
// affects Parker winding.  All other inputs remain untouched defaults.  The
// driver evaluates density, speed, and magnetic magnitude through the checked
// public API, then asks the production thermodynamic closure for pressure and
// sound speed.  Alfvén and perpendicular fast-mode speeds are derived only
// from those production outputs and immutable physical constants.

#include <swcme1d.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

double parse_number(const std::string& text, const char* name) {
  std::size_t used = 0;
  const double value = std::stod(text, &used);
  if (used != text.size() || !std::isfinite(value))
    throw std::runtime_error(std::string("invalid ") + name);
  return value;
}

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      std::cerr << "usage: vp03_model_driver INPUT_CSV OUTPUT_CSV\n";
      return 2;
    }
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    if (!input || !output) throw std::runtime_error("cannot open input/output CSV");
    std::string line;
    if (!std::getline(input, line) || line != "radius_au,latitude_deg")
      throw std::runtime_error("unexpected input schema");
    output << "radius_au,latitude_deg,density_cm3,speed_km_s,pressure_pa,"
              "sound_speed_km_s,alfven_speed_km_s,fast_speed_km_s\n"
           << std::setprecision(17);
    while (std::getline(input, line)) {
      if (line.empty()) continue;
      const std::size_t comma = line.find(',');
      if (comma == std::string::npos || line.find(',', comma + 1) != std::string::npos)
        throw std::runtime_error("malformed input row");
      const double radius_au = parse_number(line.substr(0, comma), "radius_au");
      const double latitude_deg = parse_number(line.substr(comma + 1), "latitude_deg");
      swcme1d::Params parameters;
      parameters.sin_theta = std::cos(latitude_deg * swcme1d::PI / 180.0);
      parameters.V0_sh_kms = parameters.V_sw_kms;
      swcme1d::Model model(parameters);
      const auto check = model.validate();
      if (!check.ok()) throw std::runtime_error(check.summary("VP03 configuration"));
      const swcme1d::StepState state = model.prepare_step(0.0);
      const double radius_m = radius_au * swcme1d::AU;
      double n = 0.0, v = 0.0, br = 0.0, bp = 0.0, bmag = 0.0, divv = 0.0;
      const auto status = model.evaluate_radii_with_B_div_checked(
          state, &radius_m, &n, &v, &br, &bp, &bmag, &divv, 1);
      if (!status.ok()) throw std::runtime_error(status.summary());
      const auto thermo = swcme::solarwind::thermodynamic_state(state.common.solar_wind, n);
      const double alfven = bmag / std::sqrt(swcme1d::MU0 * thermo.mass_density_kg_m3);
      const double fast = std::sqrt(thermo.sound_speed_m_s * thermo.sound_speed_m_s + alfven * alfven);
      output << radius_au << ',' << latitude_deg << ',' << n / 1.0e6 << ',' << v / 1000.0 << ','
             << thermo.pressure_Pa << ',' << thermo.sound_speed_m_s / 1000.0 << ','
             << alfven / 1000.0 << ',' << fast / 1000.0 << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("failed while writing output CSV");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "VP03 model driver: " << error.what() << '\n';
    return 1;
  }
}
