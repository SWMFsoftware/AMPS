// VP01 production-model driver.
//
// The observational analysis is written in Python, but the values labelled
// "SWCME" must come from the public C++ API rather than from a duplicated
// equation in the plotting script.  This small program reads the radial-bin
// centers, constructs an ordinary 1-D model, prepares one immutable state, and
// writes checked density evaluations in cm^-3 for independent comparison.

#include <swcme1d.hpp>

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

double parse_positive(const char* text, const char* name) {
  char* end = nullptr;
  errno = 0;
  const double value = std::strtod(text, &end);
  if (errno != 0 || end == text || *end != '\0' || !std::isfinite(value) ||
      value <= 0.0) {
    throw std::invalid_argument(std::string(name) +
                                " must be a finite positive number");
  }
  return value;
}

std::vector<double> read_bin_centers(const char* path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(std::string("cannot open ") + path);

  std::string line;
  if (!std::getline(input, line)) throw std::runtime_error("empty bin CSV");
  if (line.rfind("radius_au", 0) != 0)
    throw std::runtime_error("bin CSV must begin with radius_au");

  std::vector<double> radii_au;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::size_t comma = line.find(',');
    const std::string first = line.substr(0, comma);
    const double radius = parse_positive(first.c_str(), "radius_au");
    radii_au.push_back(radius);
  }
  if (radii_au.empty()) throw std::runtime_error("bin CSV has no data rows");
  return radii_au;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 4) {
      std::cerr << "usage: vp01_model_driver BIN_CSV OUTPUT_CSV N1AU_CM3\n";
      return 2;
    }

    const std::vector<double> radii_au = read_bin_centers(argv[1]);
    const double density_at_1au_cm3 = parse_positive(argv[3], "N1AU_CM3");

    swcme1d::Params parameters;
    parameters.n1AU_cm3 = density_at_1au_cm3;
    // Matching the nominal CME and wind speeds guarantees that the sampled
    // profile is the undisturbed background, independent of front location.
    parameters.V0_sh_kms = parameters.V_sw_kms;
    swcme1d::Model model(parameters);
    const auto validation = model.validate();
    if (!validation.ok())
      throw std::runtime_error(validation.summary("VP01 model configuration"));
    const swcme1d::StepState state = model.prepare_step(0.0);

    // Helios corefit has no alpha-particle density.  Evaluate a documented
    // 0%-to-5% composition bracket rather than silently asserting that proton
    // core density is total mass density.  The first endpoint is SWCME's exact
    // compatibility closure; the second uses the production multi-species
    // closure with n_alpha/n_proton=0.05.
    swcme1d::Params alpha_parameters = parameters;
    alpha_parameters.thermodynamic_closure =
        swcme::solarwind::ThermodynamicClosure::MultiSpecies;
    alpha_parameters.alpha_to_proton_ratio = 0.05;
    swcme1d::Model alpha_model(alpha_parameters);
    const swcme1d::StepState alpha_state = alpha_model.prepare_step(0.0);

    std::vector<double> radius_m(radii_au.size());
    std::vector<double> density_m3(radii_au.size());
    std::vector<double> speed_m_s(radii_au.size());
    for (std::size_t i = 0; i < radii_au.size(); ++i)
      radius_m[i] = radii_au[i] * swcme1d::AU;

    const swcme::ModelStatus status = model.evaluate_radii_fast_checked(
        state, radius_m.data(), density_m3.data(), speed_m_s.data(),
        radius_m.size());
    if (!status.ok())
      throw std::runtime_error(std::string("VP01 evaluation: ") + status.summary());

    std::ofstream output(argv[2]);
    if (!output) throw std::runtime_error(std::string("cannot open ") + argv[2]);
    output << "radius_au,swcme_density_cm3,"
              "swcme_mass_density_proton_only_kg_m3,"
              "swcme_mass_density_alpha5_kg_m3\n"
           << std::setprecision(17);
    for (std::size_t i = 0; i < radii_au.size(); ++i) {
      const auto proton_only = swcme::solarwind::thermodynamic_state(
          state.common.solar_wind, density_m3[i]);
      const auto alpha_five_percent = swcme::solarwind::thermodynamic_state(
          alpha_state.common.solar_wind, density_m3[i]);
      output << radii_au[i] << ',' << density_m3[i] / 1.0e6 << ','
             << proton_only.mass_density_kg_m3 << ','
             << alpha_five_percent.mass_density_kg_m3 << '\n';
    }
    output.flush();
    if (!output) throw std::runtime_error("failed while writing model CSV");
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "VP01 model driver: " << error.what() << '\n';
    return 1;
  }
}
