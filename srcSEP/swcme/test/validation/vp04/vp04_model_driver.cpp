// VP04 production data-driven-kinematics driver.
//
// Fit and query files are deliberately separate: the runner can withhold
// observations without accidentally presenting them to the model. This adapter
// calls the shared production swcme::kinematics API directly. That scope is
// intentional: constructing a complete 1-D/3-D model also attempts an unrelated
// downstream Rankine-Hugoniot closure at every coronagraph height, which can
// reject a valid interpolation query before VP04 obtains the apex state it is
// designed to validate.

#include <swcme_constants.hpp>
#include <swcme_kinematics.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

std::vector<std::pair<double, double>> read_pairs(
    const char* path, const std::string& header) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(std::string("cannot open ") + path);
  std::string line;
  if (!std::getline(input, line) || line != header)
    throw std::runtime_error("unexpected CSV header");
  std::vector<std::pair<double, double>> rows;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::size_t comma = line.find(',');
    if (comma == std::string::npos ||
        line.find(',', comma + 1) != std::string::npos)
      throw std::runtime_error("malformed CSV row");
    rows.emplace_back(std::stod(line.substr(0, comma)),
                      std::stod(line.substr(comma + 1)));
  }
  if (rows.empty()) throw std::runtime_error("CSV has no data rows");
  return rows;
}

int main(int argc, char** argv) {
  try {
    if (argc != 4) {
      std::cerr << "usage: vp04_model_driver FIT_CSV QUERY_CSV OUTPUT_CSV\n";
      return 2;
    }
    const auto fit = read_pairs(argv[1], "time_s,radius_rs");
    const auto query = read_pairs(argv[2], "time_s,observed_radius_rs");

    swcme::kinematics::Config config;
    config.mode = swcme::kinematics::Mode::DataDriven;
    // The analytical validator also requires finite baseline DBM fields even
    // though DATA_DRIVEN does not consume them. Supply ordinary physical
    // values so the adapter exercises the same public validation contract.
    config.r0_m = fit.front().second * swcme::constants::SOLAR_RADIUS_M;
    config.V0_m_s = 400000.0;
    config.Vsw_m_s = 400000.0;
    config.Gamma_m_inv = 0.0;
    for (const auto& knot : fit) {
      config.data_time_s.push_back(knot.first);
      config.data_radius_m.push_back(
          knot.second * swcme::constants::SOLAR_RADIUS_M);
    }
    if (swcme::kinematics::validate(config) != swcme::kinematics::Status::Ok)
      throw std::runtime_error("production kinematics rejected the fit table");

    std::ofstream output(argv[3]);
    if (!output) throw std::runtime_error("cannot open output CSV");
    output << "time_s,swcme_radius_rs,swcme_speed_km_s\n"
           << std::setprecision(17);
    for (const auto& point : query) {
      const swcme::kinematics::State state =
          swcme::kinematics::evaluate(config, point.first);
      if (state.status != swcme::kinematics::Status::Ok)
        throw std::runtime_error(std::string("kinematics query failed: ") +
                                 swcme::kinematics::status_name(state.status));
      output << point.first << ','
             << state.radius_m / swcme::constants::SOLAR_RADIUS_M << ','
             << state.speed_m_s / 1000.0 << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("failed while writing output CSV");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "VP04 model driver: " << error.what() << '\n';
    return 1;
  }
}
