// Small CSPICE-independent probe used by the Stage-A3 Python tests.  It calls
// the same strict parser, compatibility check, basis, addition, and clipping
// functions as cSourceModel, while keeping failure tests fast and isolated.
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "ActivityEnhancement.h"

int main(int argc,char** argv) {
  if (argc!=8) {
    std::cerr << "usage: activity_enhancement_probe FILE BODY_FRAME LEGACY_ID x y z legacy_linear\n";
    return EXIT_FAILURE;
  }
  try {
    const ActivityEnhancement::Configuration configuration=
        ActivityEnhancement::ReadFile(argv[1]);
    ActivityEnhancement::ValidateCompatibility(configuration,argv[2],argv[3]);
    const Vec3 point(std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]));
    const double legacy=std::stod(argv[7]);
    double residual,total;
    const double applied=ActivityEnhancement::Apply(configuration,point,legacy,
        &residual,&total);
    std::cout << std::setprecision(17)
        << configuration.MaximumDegree << " " << configuration.Coefficients.size()
        << " " << residual << " " << total << " " << applied << "\n";
  }
  catch (const std::exception& exception) {
    std::cerr << "Error: " << exception.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
