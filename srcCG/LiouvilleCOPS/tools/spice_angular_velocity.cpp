// Print the instantaneous angular velocity of the configured nucleus body
// frame at one or more UTC epochs.  The implementation deliberately calls the
// same cSpiceGeometry::Evaluate routine as the forward model, so a derived
// centrifugal-acceleration field cannot silently use a different frame,
// kernel set, or CK interpolation convention.
//
// Usage:
//   ./spice_angular_velocity INPUT.in UTC [UTC ...]
//
// The machine-readable records begin with "OMEGA"; cSpiceGeometry also prints
// the kernel/FOV audit that proves which SPICE inputs were loaded.

#include <cstdio>
#include <filesystem>
#include <string>

#include "Geometry.h"
#include "Input.h"
#include "Paths.h"

int main(int argc,char** argv) {
  if (argc<3) {
    std::fprintf(stderr,"usage: %s <input file> <UTC> [UTC ...]\n",argv[0]);
    return EXIT_FAILURE;
  }

  std::string inputPath=argv[1];
  if (!std::filesystem::exists(inputPath)) inputPath=CopsPaths::ProductInput(inputPath);

  cInput input;
  input.Read(inputPath.c_str());
  cSpiceGeometry geometry;
  geometry.Init(input);

  for (int index=2;index<argc;index++) {
    cObservationGeometry state;
    geometry.Evaluate(argv[index],state);
    std::printf("OMEGA %s %.17e %.17e %.17e %.17e\n",argv[index],
        state.Omega.x,state.Omega.y,state.Omega.z,Length(state.Omega));
  }
  return EXIT_SUCCESS;
}
