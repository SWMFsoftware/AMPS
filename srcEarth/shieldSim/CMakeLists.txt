cmake_minimum_required(VERSION 3.16)
project(shieldSim LANGUAGES CXX)

# Geant4 supplies compiler flags and include/library paths through Geant4_USE_FILE.
# Build example:
#   mkdir build && cd build
#   cmake .. -DGeant4_DIR=/path/to/geant4/lib/Geant4-*/
#   make -j
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Capture build/source identity in a generated header. A dirty state is
# intentionally recorded rather than rejected: development results remain
# reproducible enough to identify, while validation gates can require clean
# reviewed revisions later.
execute_process(COMMAND git -C ${PROJECT_SOURCE_DIR} rev-parse HEAD
  OUTPUT_VARIABLE SHIELDSIM_GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET RESULT_VARIABLE _git_revision_status)
if(NOT _git_revision_status EQUAL 0)
  set(SHIELDSIM_GIT_REVISION "unavailable")
endif()
execute_process(COMMAND git -C ${PROJECT_SOURCE_DIR} status --porcelain
  OUTPUT_VARIABLE _git_status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(_git_status STREQUAL "")
  set(SHIELDSIM_GIT_STATE "clean")
else()
  set(SHIELDSIM_GIT_STATE "dirty")
endif()
string(TOUPPER "${CMAKE_BUILD_TYPE}" SHIELDSIM_BUILD_TYPE_UPPER)
set(SHIELDSIM_BUILD_TYPE_FLAGS "${CMAKE_CXX_FLAGS_${SHIELDSIM_BUILD_TYPE_UPPER}}")
configure_file(${PROJECT_SOURCE_DIR}/include/BuildProvenance.hh.in
               ${PROJECT_BINARY_DIR}/generated/BuildProvenance.hh @ONLY)

add_executable(shieldSim
  shieldSim.cc
  src/AbsorberGeometry.cc
  src/AccuracyController.cc
  src/CaseDefinition.cc
  src/CLI.cc
  src/ComputedQuantities.cc
  src/ContentHash.cc
  src/DetectorConstruction.cc
  src/EventAction.cc
  src/GCRSpectrum.cc
  src/MaterialCatalog.cc
  src/MaterialLibrary.cc
  src/MaterialFactory.cc
  src/OutputUtils.cc
  src/PrimaryGeneratorAction.cc
  src/PublishedSourceAdapters.cc
  src/RunAction.cc
  src/SpectrumDataset.cc
  src/SteppingAction.cc
  src/StrictJson.cc
)

target_include_directories(shieldSim PRIVATE ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/generated)
target_link_libraries(shieldSim ${Geant4_LIBRARIES})

# Source ingestion is intentionally a standalone, Geant4-free executable. A
# successful adapter conversion is evidence about data semantics and units,
# not authorization to use that dataset as a transport boundary condition.
add_executable(shieldSimSourceTool
  sourceTool.cc
  src/ContentHash.cc
  src/PublishedSourceAdapters.cc
  src/SpectrumDataset.cc
  src/StrictJson.cc
)
target_include_directories(shieldSimSourceTool PRIVATE ${PROJECT_SOURCE_DIR}/include)

include(CTest)
if(BUILD_TESTING)
  # Analytic/observational tests write compact numerical evidence here.  A
  # dependent test renders review figures in the build tree; source fixtures
  # and the repository remain immutable during testing.
  set(SHIELDSIM_VALIDATION_EVIDENCE_DIR "${CMAKE_BINARY_DIR}/validation_evidence")
  set(SHIELDSIM_VALIDATION_FIGURE_DIR "${CMAKE_BINARY_DIR}/validation_figures")
  add_executable(stage1_core_tests
    tests/stage1_core_tests.cc
    src/AccuracyController.cc
    src/CaseDefinition.cc
    src/ContentHash.cc
    src/MaterialLibrary.cc
    src/SpectrumDataset.cc
    src/StrictJson.cc
  )
  target_include_directories(stage1_core_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
  target_compile_definitions(stage1_core_tests PRIVATE SHIELDSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME stage1_core
    COMMAND stage1_core_tests ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME stage1_case_boundaries
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage1_case_boundaries.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  # Some local Geant4 builds provide imported link targets without populating
  # Geant4_LIBRARY_DIRS. Derive the install lib directory from Geant4_DIR so the
  # validation executable can be loaded by CTest in either packaging style.
  get_filename_component(_geant4_cmake_parent "${Geant4_DIR}" DIRECTORY)
  get_filename_component(_geant4_runtime_lib "${_geant4_cmake_parent}" DIRECTORY)
  set_tests_properties(stage1_case_boundaries PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${_geant4_runtime_lib}:$ENV{LD_LIBRARY_PATH}")
  add_test(NAME stage1_transport_smoke
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage1_transport_smoke.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  set_tests_properties(stage1_transport_smoke PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${_geant4_runtime_lib}:$ENV{LD_LIBRARY_PATH}")

  add_executable(stage2_source_adapter_tests
    tests/stage2_source_adapter_tests.cc
    src/ContentHash.cc
    src/PublishedSourceAdapters.cc
    src/SpectrumDataset.cc
    src/StrictJson.cc
  )
  target_include_directories(stage2_source_adapter_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
  target_compile_definitions(stage2_source_adapter_tests PRIVATE SHIELDSIM_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
  add_test(NAME stage2_source_adapters
    COMMAND stage2_source_adapter_tests ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  add_test(NAME stage2_reference_spectra
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_reference_spectra.py)
  add_test(NAME stage2_source_tool
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage2_source_tool.py
            $<TARGET_FILE:shieldSimSourceTool> ${PROJECT_SOURCE_DIR})
  add_test(NAME validation_figures
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/generate_validation_figures.py
            ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${SHIELDSIM_VALIDATION_FIGURE_DIR} $<TARGET_FILE:shieldSimSourceTool>)
  set_tests_properties(validation_figures PROPERTIES
    DEPENDS "stage1_core;stage1_case_boundaries;stage1_transport_smoke;stage2_source_adapters;stage2_reference_spectra;stage2_source_tool")
endif()

# Copy example input files into the build tree for convenience.
file(COPY ${PROJECT_SOURCE_DIR}/examples DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})
