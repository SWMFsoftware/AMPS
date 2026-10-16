// Native MMS01 product calculator.
//
// The executable is intentionally small: it parses a prepared, authenticated table
// and calls the same FluxNumerics kernels available to Earth production code.  CDF
// interpretation and HPCA calibration metadata remain in prepare_MMS01.py, allowing
// the comparator to use the archived Level-2 flux as an independent reference.
// No trajectory, boundary mapping, energy-bin edge, or density is introduced here.

#include "../../util/FluxNumerics.h"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace FN=Earth::FluxNumerics;

// CMAD revision 4, equation 5-13/Table 5-15 uses mass number one for H+.
// This is intentionally the document's product coefficient, not the SI proton
// rest mass converted at greater precision; the latter would change the archived
// Level-2 convention this adapter is required to test.
static constexpr double kHpcaHplusMass_eV_s2_cm2=1.0453453e-12;

struct SpectrumCell {
  double energy_eV{0.0};
  std::vector<double> intensity;
  std::vector<double> weight_sr;
};

int main(int argc,char** argv) {
  if (argc!=4) {
    std::cerr << "usage: native_hpca_products INPUT DIRECTIONAL_CSV SPECTRUM_CSV\n";
    return EXIT_FAILURE;
  }
  try {
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open prepared native input");
    std::ofstream directional(argv[2]);
    std::ofstream spectrum(argv[3]);
    if (!directional || !spectrum) throw std::runtime_error("cannot create native output");
    directional << std::setprecision(17);
    spectrum << std::setprecision(17);
    directional << "distribution_index,sweep_index,start_azimuth,anode,energy_index,"
                   "energy_eV,kernel_reconstructed_differential_intensity\n";

    std::map<std::pair<int,int>,SpectrumCell> cells;
    std::string line;
    std::size_t rows=0;
    while (std::getline(input,line)) {
      if (line.empty() || line[0]=='#') continue;
      std::istringstream parser(line);
      int distributionIndex=0,sweepIndex=0,startAzimuth=0,anode=0,energyIndex=0;
      double energy_eV=0.0,psd=0.0,weight_sr=0.0;
      if (!(parser >> distributionIndex >> sweepIndex >> startAzimuth >> anode >>
            energyIndex >> energy_eV >> psd >> weight_sr))
        throw std::runtime_error("malformed prepared row at data row "+std::to_string(rows+1));
      std::string trailing;
      if (parser >> trailing)
        throw std::runtime_error("extra column in prepared row "+std::to_string(rows+1));
      const double intensity=FN::DifferentialIntensityFromPhaseSpaceDensity(
          energy_eV,kHpcaHplusMass_eV_s2_cm2,psd);
      if (!std::isfinite(intensity))
        throw std::runtime_error("invalid PSD conversion at data row "+std::to_string(rows+1));
      directional << distributionIndex << ',' << sweepIndex << ',' << startAzimuth << ','
                  << anode << ',' << energyIndex << ',' << energy_eV << ',' << intensity << '\n';
      SpectrumCell& cell=cells[std::make_pair(distributionIndex,energyIndex)];
      if (cell.intensity.empty()) cell.energy_eV=energy_eV;
      else if (cell.energy_eV!=energy_eV)
        throw std::runtime_error("energy coordinate changed inside one spectrum cell");
      cell.intensity.push_back(intensity);
      cell.weight_sr.push_back(weight_sr);
      ++rows;
    }
    if (rows==0) throw std::runtime_error("prepared native input has no data rows");

    spectrum << "distribution_index,energy_index,energy_eV,"
                "kernel_reconstructed_angular_integral\n";
    for (const auto& item:cells) {
      const double integral=FN::WeightedAngularIntegral(
          item.second.intensity,item.second.weight_sr);
      if (!std::isfinite(integral))
        throw std::runtime_error("invalid angular integral");
      spectrum << item.first.first << ',' << item.first.second << ','
               << item.second.energy_eV << ',' << integral << '\n';
    }
  }
  catch (const std::exception& error) {
    std::cerr << "MMS01 native calculation failed: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
