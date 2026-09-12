// VP05 production DBM transit driver.
//
// The outer observation supplies only a target radius. Arrival time and speed
// are solved from the inner radius/speed with the production DBM and fixed
// default wind/drag parameters. Bracketed bisection avoids coupling this small
// validation adapter to any independent closed-form inversion.

#include <swcme_constants.hpp>
#include <swcme_kinematics.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string value;
  while (std::getline(stream, value, ',')) fields.push_back(value);
  return fields;
}

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      std::cerr << "usage: vp05_model_driver INPUT_CSV OUTPUT_CSV\n";
      return 2;
    }
    std::ifstream input(argv[1]);
    std::ofstream output(argv[2]);
    if (!input || !output) throw std::runtime_error("cannot open input/output CSV");
    std::string line;
    if (!std::getline(input, line) ||
        line != "pair_id,radius0_au,speed0_km_s,target_au")
      throw std::runtime_error("unexpected input schema");
    output << "pair_id,predicted_transit_s,predicted_speed_km_s\n"
           << std::setprecision(17);
    while (std::getline(input, line)) {
      if (line.empty()) continue;
      const auto fields = split(line);
      if (fields.size() != 4) throw std::runtime_error("malformed input row");
      const double radius0_au = std::stod(fields[1]);
      const double speed0_km_s = std::stod(fields[2]);
      const double target_au = std::stod(fields[3]);
      swcme::kinematics::Config config;
      config.mode = swcme::kinematics::Mode::DBM;
      config.r0_m = radius0_au * swcme::constants::AU_M;
      config.V0_m_s = speed0_km_s * 1000.0;
      config.Vsw_m_s = 400000.0;
      config.Gamma_m_inv = 1.0e-10;  // 1e-7 km^-1 converted to SI.
      if (swcme::kinematics::validate(config) != swcme::kinematics::Status::Ok)
        throw std::runtime_error("production DBM rejected pair input");
      const double target_m = target_au * swcme::constants::AU_M;
      double low = 0.0;
      double high = 4.0 * (target_m - config.r0_m) /
                    std::min(config.V0_m_s, config.Vsw_m_s);
      while (swcme::kinematics::evaluate(config, high).radius_m < target_m)
        high *= 2.0;
      for (int iteration = 0; iteration < 120; ++iteration) {
        const double middle = 0.5 * (low + high);
        const auto state = swcme::kinematics::evaluate(config, middle);
        if (state.status != swcme::kinematics::Status::Ok)
          throw std::runtime_error("production DBM failed inside arrival bracket");
        if (state.radius_m < target_m) low = middle;
        else high = middle;
      }
      const double arrival_s = 0.5 * (low + high);
      const auto state = swcme::kinematics::evaluate(config, arrival_s);
      output << fields[0] << ',' << arrival_s << ',' << state.speed_m_s / 1000.0 << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("failed while writing model CSV");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "VP05 model driver: " << error.what() << '\n';
    return 1;
  }
}
