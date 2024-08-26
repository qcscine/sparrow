# Load source file lists
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
import_utils_os()
include(ImportCore)
import_core()
include(ImportCereal)
import_cereal()

if(SCINE_PARALLELIZE)
  find_package(OpenMP)
endif()

if(NOT TARGET Boost::filesystem OR NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS program_options filesystem)
endif()

set_source_files_properties(
    ${SPARROW_UNITY_DISABLE} PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON
)

# Sparrow is a Scine module and is always shared, never linked against
add_library(Sparrow SHARED ${SPARROW_MODULE_FILES})

target_include_directories(Sparrow PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/Scine>
)

set_target_properties(Sparrow PROPERTIES
  PREFIX ""
  OUTPUT_NAME sparrow.module
  WINDOWS_EXPORT_ALL_SYMBOLS ON
)

if(MSVC)
  target_compile_definitions(Sparrow PUBLIC BOOST_ALL_NO_LIB)
endif()

# Add an option to repropduce the standard PM6 errors
# This is needed until a new reparametrization of PM6 is done with the right method.
option(REPRODUCE_PM6_ERRORS "Reproduce the errors of standard PM6 implementations" ON)

if(SCINE_PARALLELIZE)
  target_compile_options(Sparrow PRIVATE
    $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
  )
endif()
target_compile_definitions(Sparrow
 PRIVATE
  $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
  $<$<BOOL:${REPRODUCE_PM6_ERRORS}>:REPRODUCE_PM6_ERRORS>
)

target_link_libraries(Sparrow
  PUBLIC
    Scine::UtilsOS
    Scine::Core
  PRIVATE
    cereal
    Boost::filesystem
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

if(MSVC)
  target_compile_options(Sparrow PRIVATE /openmp- /bigobj /MP)
endif()

# Add namespaced aliases
add_library(Scine::Sparrow ALIAS Sparrow)
add_library(Scine::SparrowModule ALIAS Sparrow)

# -- Define Install
# Targets
install(TARGETS Sparrow
  EXPORT sparrowTargets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT  ScineSparrow
  EXPORT_NAME sparrowTargets
)

if(SCINE_BUILD_TESTS OR SCINE_BUILD_PYTHON_BINDINGS)
  find_package(PythonInterp REQUIRED)
endif()

# Tests
if(SCINE_BUILD_TESTS)
  include(ImportGTest)
  import_gtest()
  add_executable(Sparrow_tests ${SPARROW_TEST_FILES} $<$<NOT:$<CONFIG:Debug>>:${SPARROW_TEST_SLOW_FILES}>)
  target_link_libraries(Sparrow_tests PRIVATE
    GTest::Main
    GMock::GMock
    Scine::Sparrow
    cereal
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    ${CMAKE_DL_LIBS}
  )
  add_test(NAME Sparrow COMMAND Sparrow_tests)
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    set(TEST_SELECTION "::TestSparrowFast")
  else()
    set(TEST_SELECTION "")
  endif()

  if(MSVC)
    target_compile_options(Sparrow_tests PRIVATE /bigobj /MP /openmp-)
  endif()

  # App tests
  add_test(
    NAME Sparrow_App
    COMMAND ${PYTHON_EXECUTABLE} -B -m pytest test_sparrow.py${TEST_SELECTION} --junitxml=${CMAKE_CURRENT_BINARY_DIR}/app_report.xml
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Tests/AppTests
  )
  set_tests_properties(Sparrow_App PROPERTIES
    ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{SCINE_MODULE_PATH}
    ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}
    ENVIRONMENT PATH=${CMAKE_CURRENT_BINARY_DIR}:${PROJECT_BINARY_DIR}:$ENV{PATH}
  )
endif()

# Set the RPATH to be used when installing.
if(APPLE)
  set(CMAKE_INSTALL_RPATH "@loader_path;@loader_path/../lib")
else()
  # $ORIGIN refers to the installation directory of the binary.
  set(CMAKE_INSTALL_RPATH "\$ORIGIN;\$ORIGIN/../lib")
endif()

# Executable
add_executable(SparrowApp ${SPARROW_APP_FILES})
add_executable(Scine::SparrowApp ALIAS SparrowApp)

if(MSVC)
  target_compile_definitions(SparrowApp PUBLIC BOOST_ALL_DYN_LINK)
endif()

set_target_properties(SparrowApp PROPERTIES OUTPUT_NAME sparrow)
if(APPLE)
  set_target_properties(SparrowApp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

if(MSVC)
  target_compile_options(SparrowApp PRIVATE /bigobj /MP /openmp-)
endif()

target_link_libraries(SparrowApp PRIVATE
  Boost::program_options
  Scine::UtilsOS
  Scine::Core
  ${CMAKE_DL_LIBS}
)

add_subdirectory(Sparrow/Implementations/RealTimeSpectroscopy)

add_executable(RTSpectroscopyDriver App/RTSpectroscopyDriver.cpp)
target_link_libraries(RTSpectroscopyDriver
  PRIVATE
    Scine::UtilsOS
    Scine::Core
    Scine::RealTimeSpectroscopy
    yaml-cpp
    ${CMAKE_DL_LIBS}
)

if(MSVC)
  target_compile_options(RTSpectroscopyDriver PRIVATE /bigobj /MP /openmp-)
endif()

# Add link library specific to MinGW
if(WIN32 AND MINGW)
  target_link_libraries(SparrowApp PRIVATE ws2_32)
endif()

# Python Bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  # Figure out which targets the python bindings are going to need copied
  include(TargetLibName)
  set(_py_targets_to_copy Sparrow)
  target_lib_type(Scine::UtilsOS _utils_libtype)
  if(_utils_libtype STREQUAL "SHARED_LIBRARY")
    list(APPEND _py_targets_to_copy Scine::UtilsOS)
  endif()
  unset(_utils_libtype)

  # Generate generator expressions for each targets and figure out filenames
  # for the python setup file
  set(sparrow_PY_DEPS "")
  foreach(target ${_py_targets_to_copy})
    list(APPEND _py_target_gen_exprs "\$<TARGET_FILE:${target}>")
    target_lib_filename(${target} _target_filename)
    string(APPEND sparrow_PY_DEPS ", \"${_target_filename}\"")
  endforeach()
  message(STATUS "Targets to copy for python bindings: ${_py_targets_to_copy}")
  unset(_py_targets_to_copy)

  add_custom_command(TARGET Sparrow POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${_py_target_gen_exprs} ${CMAKE_CURRENT_BINARY_DIR}/scine_sparrow
    COMMENT "Copying required shared libraries into python package directory"
  )
  unset(_py_target_gen_exprs)

  install(CODE
    "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade --no-deps ${CMAKE_CURRENT_BINARY_DIR}
                     RESULT_VARIABLE retcode)
     if(NOT \$\{retcode\} EQUAL 0)
       message(FATAL_ERROR \"Fatal error when installing Python module using PIP.\")
     endif()"
  )

  # Copy the setup file and README.rst into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
      ${PROJECT_SOURCE_DIR}/LICENSE.txt
      ${EXTERNAL_LICENSE_FILES}
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_sparrow/
  )

  # Add setup and readme to current binary dir
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
    ${CMAKE_CURRENT_BINARY_DIR}/setup.py
    @ONLY
  )
  unset(sparrow_PY_DEPS)
  file(
    COPY ${PROJECT_SOURCE_DIR}/README.rst
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )

  # Python Binding Tests
  if(SCINE_BUILD_TESTS)
    add_test(
      NAME Sparrow_Calculator_Python_Bindings
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests/CalculatorTest.py --junitxml=${CMAKE_CURRENT_BINARY_DIR}/bindings_calculator_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_test(
      NAME Sparrow_CalculatorWithReference_Python_Bindings
      COMMAND ${PYTHON_EXECUTABLE} -B -m pytest ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests/CalculatorWithReferenceTest.py --junitxml=${CMAKE_CURRENT_BINARY_DIR}/bindings_calculator_with_reference_report.xml
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set(Sparrow_Python_Bindings
      Sparrow_Calculator_Python_Bindings
      Sparrow_CalculatorWithReference_Python_Bindings
      )
    set_tests_properties(${Sparrow_Python_Bindings} PROPERTIES
        ENVIRONMENT SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR}
        ENVIRONMENT PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}:$ENV{PYTHONPATH}
      )

    # In-tree scine utilities?
    if(scine-utils-os_BINARY_DIR)
      set(UTILS_PYTHONPATH ${scine-utils-os_BINARY_DIR}/src/Utils)
      set_tests_properties(${Sparrow_Python_Bindings} PROPERTIES
        ENVIRONMENT PYTHONPATH=${UTILS_PYTHONPATH}:$ENV{PYTHONPATH}
      )
    endif()
  endif()

  include(SphinxDocumentation)
  scine_sphinx_documentation(
    TARGET scine_sparrow
    CONFIGURATION ${CMAKE_CURRENT_SOURCE_DIR}/Python/sphinx/conf.py.in
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Python/sphinx
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/doc-py
    LINK UtilsPythonDoc
  )
endif()

add_subdirectory(Embed)

# Targets
install(TARGETS SparrowApp RUNTIME DESTINATION bin)
