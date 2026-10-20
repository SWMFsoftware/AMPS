// CSPICE-independent probe for the strict per-facet replay format.
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "ActivityEnhancement.h"
#include "Mesh.h"

int main(int argc,char** argv) {
  if (argc!=9) {
    std::cerr << "usage: probe FILE BODY_FRAME LEGACY_ID x y z area legacy\n";
    return EXIT_FAILURE;
  }
  try {
    cMesh mesh;
    cTriangle face={};
    face.Centroid=Vec3(std::stod(argv[4]),std::stod(argv[5]),std::stod(argv[6]));
    face.Area=std::stod(argv[7]);
    mesh.Faces.push_back(face);
    const auto configuration=ActivityEnhancement::ReadFacetFile(argv[1],mesh);
    ActivityEnhancement::ValidateCompatibility(configuration,argv[2],argv[3]);
    double residual,total;
    const double applied=ActivityEnhancement::Apply(
        configuration,face.Centroid,std::stod(argv[8]),&residual,&total,0);
    std::cout << std::setprecision(17) << configuration.FacetResiduals.size()
              << " " << residual << " " << total << " " << applied << "\n";
  }
  catch (const std::exception& exception) {
    std::cerr << "Error: " << exception.what() << "\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
