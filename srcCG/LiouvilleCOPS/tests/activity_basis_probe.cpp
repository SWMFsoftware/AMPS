// Small standard-library-only driver used by test_a2_activity_basis.py.
// It keeps the cross-language comparison independent of CSPICE and the main
// executable: each input line is one Cartesian vector, each output line is the
// complete C++ basis in normative coefficient order.
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include "ActivityEnhancementBasis.h"

int main(int argc,char** argv) {
  if (argc!=3) {
    std::cerr << "usage: activity_basis_probe MIN_DEGREE MAX_DEGREE\n";
    return 2;
  }
  try {
    const int minDegree=std::stoi(argv[1]),maxDegree=std::stoi(argv[2]);
    // Validate the range even when stdin is empty.
    ActivityEnhancementBasis::CoefficientCount(minDegree,maxDegree);
    std::cout << std::setprecision(17);
    double x,y,z;
    while (std::cin >> x >> y >> z) {
      const auto values=ActivityEnhancementBasis::EvaluateCartesian(
          x,y,z,minDegree,maxDegree);
      for (std::size_t index=0;index<values.size();index++) {
        if (index) std::cout << ' ';
        std::cout << values[index];
      }
      std::cout << '\n';
    }
    if (!std::cin.eof()) throw std::invalid_argument("invalid Cartesian input");
  }
  catch (const std::exception& error) {
    std::cerr << "Error: " << error.what() << '\n';
    return 2;
  }
  return 0;
}
