// VP02 production-model driver.
//
// The Python analysis owns observational filtering and plotting, but every
// value labelled SWCME is emitted by this executable through the public 1-D
// model API. Each row may have a distinct measured wind speed and latitude;
// consequently each row gets an immutable model instance whose fixed ray has
// the corresponding sin(colatitude). Sampling ahead of a stationary front
// guarantees that the returned magnetic field is the Parker background.

#include <swcme1d.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

double parse(const std::string& text, const char* field) {
  std::size_t used = 0;
  const double value = std::stod(text, &used);
  if (used != text.size() || !std::isfinite(value))
    throw std::runtime_error(std::string("invalid ") + field);
  return value;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      std::cerr << "usage: vp02_model_driver INPUT_CSV OUTPUT_CSV\n";
      return 2;
    }
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open input CSV");
    std::ofstream output(argv[2]);
    if (!output) throw std::runtime_error("cannot open output CSV");

    std::string line;
    if (!std::getline(input, line) || line != "radius_au,latitude_deg,speed_km_s")
      throw std::runtime_error("unexpected input schema");
    output << "radius_au,latitude_deg,speed_km_s,br_nt,bphi_nt,"
              "bmag_nt,parker_angle_deg\n" << std::setprecision(17);

    while (std::getline(input, line)) {
      if (line.empty()) continue;
      std::stringstream row(line);
      std::string r_text, latitude_text, speed_text;
      if (!std::getline(row, r_text, ',') ||
          !std::getline(row, latitude_text, ',') ||
          !std::getline(row, speed_text) || row.peek() != EOF)
        throw std::runtime_error("malformed input row");

      const double radius_au = parse(r_text, "radius_au");
      const double latitude_deg = parse(latitude_text, "latitude_deg");
      const double speed_km_s = parse(speed_text, "speed_km_s");
      swcme1d::Params parameters;
      parameters.V_sw_kms = speed_km_s;
      parameters.V0_sh_kms = speed_km_s;
      parameters.sin_theta = std::cos(latitude_deg * swcme1d::PI / 180.0);
      parameters.parker_radial_polarity = +1;
      swcme1d::Model model(parameters);
      const auto validation = model.validate();
      if (!validation.ok())
        throw std::runtime_error(validation.summary("VP02 configuration"));
      const swcme1d::StepState state = model.prepare_step(0.0);

      const double radius_m = radius_au * swcme1d::AU;
      double density_m3 = 0.0, velocity_m_s = 0.0;
      double br_t = 0.0, bphi_t = 0.0, bmag_t = 0.0, div_v = 0.0;
      const swcme::ModelStatus status = model.evaluate_radii_with_B_div_checked(
          state, &radius_m, &density_m3, &velocity_m_s, &br_t, &bphi_t,
          &bmag_t, &div_v, 1);
      if (!status.ok()) throw std::runtime_error(status.summary());
      const double angle_deg = std::atan2(-bphi_t, br_t) * 180.0 / swcme1d::PI;
      output << radius_au << ',' << latitude_deg << ',' << speed_km_s << ','
             << br_t * 1.0e9 << ',' << bphi_t * 1.0e9 << ','
             << bmag_t * 1.0e9 << ',' << angle_deg << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("failed while writing output CSV");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "VP02 model driver: " << error.what() << '\n';
    return 1;
  }
}
