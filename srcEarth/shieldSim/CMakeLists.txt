#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.27)
project(B1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# See the documentation for a guide on how to enable/disable specific components
#
find_package(Geant4 REQUIRED ui_all vis_all)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Build the shieldSim application.  This executable is independent from the
# standard B1 example and demonstrates a custom shielding simulation.  The
# source file ``shieldSim.cc`` is located in the project root (not in ``src``)
# so that it will not be inadvertently included in the exampleB1 target.
#
add_executable(shieldSim shieldSim.cc)
target_link_libraries(shieldSim PRIVATE ${Geant4_LIBRARIES})

# Optionally copy the default spectrum file to the build directory so that
# users can run ``shieldSim`` without specifying ``--spectrum``.  If the
# file is absent at configure time, this step will be skipped.
if(EXISTS ${PROJECT_SOURCE_DIR}/spectrum.txt)
  configure_file(
    ${PROJECT_SOURCE_DIR}/spectrum.txt
    ${PROJECT_BINARY_DIR}/spectrum.txt
    COPYONLY
  )
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
