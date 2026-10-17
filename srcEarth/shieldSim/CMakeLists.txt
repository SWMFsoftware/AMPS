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
  src/AncestryTrackingAction.cc
  src/CaseDefinition.cc
  src/CLI.cc
  src/ComputedQuantities.cc
  src/ContentHash.cc
  src/DetectorConstruction.cc
  src/DamageResponseLibrary.cc
  src/EventAction.cc
  src/GCRSpectrum.cc
  src/GeometryTransform.cc
  src/MaterialCatalog.cc
  src/MaterialLibrary.cc
  src/MaterialFactory.cc
  src/OutputUtils.cc
  src/PrimaryGeneratorAction.cc
  src/PublishedSourceAdapters.cc
  src/RunAction.cc
  src/SpectrumDataset.cc
  src/Stage6Physics.cc
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

# Response matrices are assembled and folded outside the Geant4 executable so
# deterministic schema/linearity tests do not need to initialize transport.
add_executable(shieldSimResponseTool
  responseTool.cc
  src/ContentHash.cc
  src/ResponseMatrix.cc
  src/StrictJson.cc
)
target_include_directories(shieldSimResponseTool PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Stage-4 reference probing is deliberately its own executable. It exposes the
# Geant4 EM stopping/range/attenuation tables in NIST-compatible units without
# coupling those comparisons to stochastic slab geometry or dose scoring.
add_executable(shieldSimReferenceTool
  referenceTool.cc
  src/StrictJson.cc
)
target_include_directories(shieldSimReferenceTool PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(shieldSimReferenceTool ${Geant4_LIBRARIES})

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
    src/DamageResponseLibrary.cc
    src/GeometryTransform.cc
    src/MaterialLibrary.cc
    src/PublishedSourceAdapters.cc
    src/SpectrumDataset.cc
    src/Stage6Physics.cc
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
  # Carry the exact dataset paths exported by this Geant4 installation into
  # transport tests.  In particular FTFP_BERT_HP must fail rather than quietly
  # skipping its neutron data when G4NEUTRONHPDATA is unavailable.
  set(_geant4_test_environment "LD_LIBRARY_PATH=${_geant4_runtime_lib}:$ENV{LD_LIBRARY_PATH}")
  # Python tempfile users must stay inside the selected build/test-output tree.
  # This makes test location explicit and avoids hidden state under /tmp.
  set(SHIELDSIM_TEST_TEMP_DIR "${CMAKE_BINARY_DIR}/test_tmp")
  file(MAKE_DIRECTORY "${SHIELDSIM_TEST_TEMP_DIR}")
  list(APPEND _geant4_test_environment "TMPDIR=${SHIELDSIM_TEST_TEMP_DIR}")
  foreach(_dataset_name IN LISTS Geant4_DATASETS)
    set(_dataset_environment_name "${Geant4_DATASET_${_dataset_name}_ENVVAR}")
    set(_dataset_path "${Geant4_DATASET_${_dataset_name}_PATH}")
    list(APPEND _geant4_test_environment
      "${_dataset_environment_name}=${_dataset_path}")
  endforeach()
  set_tests_properties(stage1_case_boundaries PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}")
  add_test(NAME stage1_transport_smoke
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage1_transport_smoke.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  set_tests_properties(stage1_transport_smoke PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}")

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
  set_tests_properties(stage2_source_tool PROPERTIES
    ENVIRONMENT "TMPDIR=${SHIELDSIM_TEST_TEMP_DIR}")
  # This is the Stage-2 executable boundary: checksum-pinned published AMS-02
  # proton bins enter the normal case parser, primary generator, Geant4
  # transport, and TID/H100/10 proxy writers. It is intentionally separate
  # from adapter-only tests so ingestion cannot be mistaken for end-to-end use.
  add_test(NAME stage2_published_transport
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage2_published_transport.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  set_tests_properties(stage2_published_transport PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}")

  add_executable(stage3_core_tests
    tests/stage3_core_tests.cc
    src/AccuracyController.cc
    src/CaseDefinition.cc
    src/ContentHash.cc
    src/DamageResponseLibrary.cc
    src/GeometryTransform.cc
    src/MaterialLibrary.cc
    src/PublishedSourceAdapters.cc
    src/ResponseMatrix.cc
    src/SpectrumDataset.cc
    src/Stage6Physics.cc
    src/StrictJson.cc
  )
  target_include_directories(stage3_core_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
  add_test(NAME stage3_core
    COMMAND stage3_core_tests ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  add_test(NAME stage3_response_weights
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage3_response_weights.py
            $<TARGET_FILE:shieldSimResponseTool>
            ${PROJECT_SOURCE_DIR}/scripts/build_response_weights.py
            ${CMAKE_BINARY_DIR}/stage3_response_weight_tests
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR} ${Python3_EXECUTABLE})
  # This is the stochastic Stage-3 acceptance gate.  Its driver owns a named
  # build-tree directory so cases, Geant4 logs, matrices, folds, and evidence
  # never spill into /tmp or the source tree.
  add_test(NAME stage3_transport
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage3_transport.py
            $<TARGET_FILE:shieldSim> $<TARGET_FILE:shieldSimResponseTool>
            ${PROJECT_SOURCE_DIR}/scripts/build_response_weights.py
            ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${CMAKE_BINARY_DIR}/stage3_transport_runs full)
  set_tests_properties(stage3_transport PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}"
    TIMEOUT 600)

  add_executable(stage4_core_tests
    tests/stage4_core_tests.cc
    src/AccuracyController.cc
    src/CaseDefinition.cc
    src/ContentHash.cc
    src/DamageResponseLibrary.cc
    src/GeometryTransform.cc
    src/MaterialLibrary.cc
    src/PublishedSourceAdapters.cc
    src/SpectrumDataset.cc
    src/Stage6Physics.cc
    src/StrictJson.cc
  )
  target_include_directories(stage4_core_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
  add_test(NAME stage4_core
    COMMAND stage4_core_tests ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  add_test(NAME stage4_nist_reference
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage4_nist_reference.py
            $<TARGET_FILE:shieldSimReferenceTool> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${CMAKE_BINARY_DIR}/stage4_reference_runs)
  set_tests_properties(stage4_nist_reference PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}"
    TIMEOUT 180)
  add_test(NAME stage4_tid_transport
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage4_tid_transport.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${CMAKE_BINARY_DIR}/stage4_tid_transport_runs)
  set_tests_properties(stage4_tid_transport PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}"
    TIMEOUT 600)

  add_executable(stage5_core_tests
    tests/stage5_core_tests.cc
    src/AccuracyController.cc
    src/CaseDefinition.cc
    src/ContentHash.cc
    src/DamageResponseLibrary.cc
    src/GeometryTransform.cc
    src/MaterialLibrary.cc
    src/PublishedSourceAdapters.cc
    src/SpectrumDataset.cc
    src/Stage6Physics.cc
    src/StrictJson.cc
  )
  target_include_directories(stage5_core_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
  add_test(NAME stage5_core
    COMMAND stage5_core_tests ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  add_test(NAME stage5_transport
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage5_transport.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${CMAKE_BINARY_DIR}/stage5_transport_runs)
  set_tests_properties(stage5_transport PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}"
    TIMEOUT 600)

  add_executable(stage6_core_tests
    tests/stage6_core_tests.cc
    src/Stage6Physics.cc
    src/StrictJson.cc
  )
  target_include_directories(stage6_core_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
  add_test(NAME stage6_core
    COMMAND stage6_core_tests ${SHIELDSIM_VALIDATION_EVIDENCE_DIR})
  add_test(NAME stage6_case_boundaries
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage6_case_boundaries.py
            $<TARGET_FILE:shieldSim> ${PROJECT_SOURCE_DIR}
            ${CMAKE_BINARY_DIR}/stage6_case_boundary_runs)
  set_tests_properties(stage6_case_boundaries PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}")
  add_test(NAME stage6_reference
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage6_reference.py
            $<TARGET_FILE:shieldSimReferenceTool> ${PROJECT_SOURCE_DIR}
            ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${CMAKE_BINARY_DIR}/stage6_reference_runs)
  set_tests_properties(stage6_reference PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}"
    TIMEOUT 240)
  add_test(NAME stage6_transport
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/test_stage6_transport.py
            $<TARGET_FILE:shieldSim> $<TARGET_FILE:shieldSimReferenceTool>
            ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${CMAKE_BINARY_DIR}/stage6_transport_runs)
  set_tests_properties(stage6_transport PROPERTIES
    ENVIRONMENT "${_geant4_test_environment}"
    TIMEOUT 900)
  add_test(NAME validation_figures
    COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tests/generate_validation_figures.py
            ${PROJECT_SOURCE_DIR} ${SHIELDSIM_VALIDATION_EVIDENCE_DIR}
            ${SHIELDSIM_VALIDATION_FIGURE_DIR} $<TARGET_FILE:shieldSimSourceTool>)
  set_tests_properties(validation_figures PROPERTIES
    DEPENDS "stage1_core;stage1_case_boundaries;stage1_transport_smoke;stage2_source_adapters;stage2_reference_spectra;stage2_source_tool;stage2_published_transport;stage3_core;stage3_response_weights;stage3_transport;stage4_core;stage4_nist_reference;stage4_tid_transport;stage5_core;stage5_transport;stage6_core;stage6_case_boundaries;stage6_reference;stage6_transport")
endif()

# Copy example input files into the build tree for convenience.
file(COPY ${PROJECT_SOURCE_DIR}/examples DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})
file(COPY ${PROJECT_SOURCE_DIR}/data DESTINATION ${PROJECT_BINARY_DIR})
