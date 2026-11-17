//$Id$
//The runtime input of the Orbiter application (the drag-coefficient test, DragCoefficientTest in main.cpp):
//
//   ./amps -i FILE -case NAME
//
//FILE uses the generic AMPS runtime-input grammar (src/INPUT_FILE.md, reader src/general/input/amps_input.h).  It holds
//exactly one block 'application: orbiter':
//
//   #begin application: orbiter
//     schema = 1
//     #begin upstream                          ! the free stream injected through the boundary of the domain
//       velocity_m_s = -3000, 0, 0
//       temperature_k = 293
//       #begin species: H2O
//         number_density_m-3 = 1.5e21          ! every species of the executable
//       #end species
//     #end upstream
//     #begin time_step
//       characteristic_speed_m_s = 5000        ! the time step of a block is 0.3 (cell size)/(this speed)
//     #end time_step
//     #begin run
//       iterations = 100                       ! the run ends at the first output after this number of iterations
//     #end run
//     #begin physics: gas_surface_interaction  ! owned by src/models/surface (Surface.h)
//       ...
//     #end physics
//     #begin case: NAME                        ! named cases (AMPS::Input::ResolveCase): -case selects one
//       ...
//     #end case
//   #end application
//
//When the file has cases, -case is required; an unknown case is an error.  Every key is required, an unknown key or
//block is an error, and the file is checked before the mesh is built.

#ifndef _ORBITER_RUNTIME_INPUT_H_
#define _ORBITER_RUNTIME_INPUT_H_

#include <string>

namespace Orbiter {
  namespace RuntimeInput {
    //parse FILE, select the case, check it and install it (collective: every MPI process reads the file and the
    //processes must agree on the resolved input); false and the reason on an error
    bool LoadAndApply(const std::string &path,const std::string &caseName,std::string *error);

    extern bool Loaded;
    extern std::string LoadedPath,SelectedCase;

    //run.iterations
    extern int nTotalIterations;

    //a description of the installed input (for the run log and the test result file)
    std::string Describe();
  }
}

#endif
