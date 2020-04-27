# Load source file lists
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
import_utils_os()
include(ImportCore)
import_core()
include(ImportCereal)
import_cereal()
find_package(OpenMP)
find_package(Boost REQUIRED filesystem)
set(Boost_USE_STATIC_LIBS OFF)

# Shared Sparrow library / Sparrow module
add_library(Sparrow SHARED ${SPARROW_MODULE_FILES})
set_property(TARGET Sparrow PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(Sparrow PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>
)

if(WIN32)
  set_target_properties(Sparrow PROPERTIES PREFIX "")
  set_target_properties(Sparrow PROPERTIES IMPORT_PREFIX "")
  set_target_properties(Sparrow PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()


set_target_properties(Sparrow PROPERTIES
  OUTPUT_NAME sparrow
)

# Set the option to repropduce the standard PM6 errors
# This is needed until a new reparametrization of PM6 is done with the right method.
option(REPRODUCE_PM6_ERRORS "Reproduce the errors of standard PM6 implementations" ON)

# include it in the build / move to the appropriate folder for build
target_compile_options(Sparrow
                       PRIVATE
                       $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>
                       PUBLIC
                       $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
                       )
target_compile_definitions(Sparrow
                           PRIVATE
                           $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>
                           $<$<BOOL:${REPRODUCE_PM6_ERRORS}>:REPRODUCE_PM6_ERRORS>
                           PUBLIC
                           BOOST_ALL_DYN_LINK
                           )


target_link_libraries(Sparrow
  PUBLIC
    Scine::Core
    Scine::UtilsOS
  PRIVATE
    cereal::cereal
    Boost::filesystem
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
)

# Add namespaced aliases
add_library(Scine::Sparrow ALIAS Sparrow)
add_library(Scine::SparrowModule ALIAS Sparrow)

# -- Define Install
# Headers
install(DIRECTORY ${PROJECT_SOURCE_DIR}/${PROJECT_NAME}/
        DESTINATION include/${PROJECT_NAME}
        FILES_MATCHING PATTERN "*.h"
)

# Targets
install(TARGETS Sparrow
        EXPORT sparrowTargets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER DESTINATION include/${PROJECT_NAME})

# Symlink on install
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_INSTALL_PREFIX}/lib/${CMAKE_SHARED_LIBRARY_PREFIX}sparrow${CMAKE_SHARED_LIBRARY_SUFFIX} \
        ${CMAKE_INSTALL_PREFIX}/lib/sparrow.module${CMAKE_SHARED_LIBRARY_SUFFIX})"
)

# Symlink on build
add_custom_command(TARGET Sparrow POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E create_symlink
                   "${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_SHARED_LIBRARY_PREFIX}sparrow${CMAKE_SHARED_LIBRARY_SUFFIX}"
                   "${CMAKE_CURRENT_BINARY_DIR}/sparrow.module${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

# Parameters
install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/Sparrow/Resources/ DESTINATION resources/Parameters)
## Configure file containing the install location of the parameters
set(parameter_dir ${CMAKE_INSTALL_PREFIX}/resources/Parameters/)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Sparrow/ParametersHeader.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/Sparrow/ParametersHeader.h)
## Configure file containing the location of the parameters for the tests.
set(parameter_dir "${CMAKE_CURRENT_SOURCE_DIR}/Sparrow/Resources/")
set(target_dir "${CMAKE_CURRENT_SOURCE_DIR}/App")
configure_file(${target_dir}/parameters_location.h.in
               ${target_dir}/parameters_location.h)

# Component-level cmake package files
scine_install_component_cmake_files(
  COMPONENT  Sparrow
  EXPORT_NAME sparrowTargets
)

# Tests
if(SCINE_BUILD_TESTS)
  add_executable(Sparrow_tests ${SPARROW_TEST_FILES})
  target_link_libraries(Sparrow_tests
    PRIVATE
      gtest_main
      gmock
      Scine::Sparrow
      cereal::cereal
      $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
  )
  add_test(NAME Sparrow COMMAND Sparrow_tests)
  target_compile_options(Sparrow_tests
    PUBLIC
      $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
  )  

  # Hard-code location of the parameter sets into the test by configuring a header file for it.
  set(parameter_dir "${CMAKE_CURRENT_SOURCE_DIR}/Sparrow/Resources/")
  set(methodTests_dir "${CMAKE_CURRENT_SOURCE_DIR}/Tests/MethodsTests")
  configure_file(${methodTests_dir}/parameters_location.h.in
                 ${methodTests_dir}/parameters_location.h)
endif()


# Set the RPATH to be used when installing. $ORIGIN refers to the installation directory of the binary.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Executable
find_package(Boost COMPONENTS program_options system REQUIRED)

add_executable(SparrowApp ${SPARROW_APP_FILES})

add_executable(Scine::SparrowApp ALIAS SparrowApp)

if (NOT APPLE)
  set_target_properties(SparrowApp PROPERTIES
    OUTPUT_NAME sparrow
  )
endif()

target_compile_definitions(SparrowApp
                           PUBLIC
                           BOOST_ALL_DYN_LINK
                           )
target_link_libraries(SparrowApp
  PRIVATE
    Boost::program_options
    Scine::UtilsOS
    Scine::Core
)
target_compile_options(SparrowApp
  PUBLIC
    $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
)  


# Add link library specific for MinGW
if(WIN32)
  target_link_libraries(SparrowApp PRIVATE $<$<BOOL:MINGW>:ws2_32>)
endif()

# Python Bindings
if(SCINE_BUILD_PYTHON_BINDINGS)
  include(ImportPybind11)
  import_pybind11()
  set(PYBIND11_PYTHON_VERSION ${PYTHONVERSION})
  pybind11_add_module(scine_sparrow ${SPARROW_PYTHON_FILES})
  set_target_properties(scine_sparrow PROPERTIES SUFFIX ".so")
  set_target_properties(scine_sparrow PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/scine_sparrow)
  target_link_libraries(scine_sparrow
          PRIVATE
          Scine::Sparrow
          Scine::UtilsOS
          )
  install(CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} --upgrade ${CMAKE_CURRENT_BINARY_DIR})")

  # Python tests
  if(SCINE_BUILD_TESTS)
    add_test(
      NAME Sparrow_Python_Bindings
      COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Python/Tests/CalculationTest.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )
    set_tests_properties(Sparrow_Python_Bindings PROPERTIES
      ENVIRONMENT "SCINE_MODULE_PATH=${CMAKE_CURRENT_BINARY_DIR};PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}"
    )
  endif()
  # Copy the setup file and README.md into the binary dir so you can install
  # using pip or make a distribution with the compiled binary
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/__init__.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/scine_sparrow
  )
  file(
    COPY
      ${CMAKE_CURRENT_SOURCE_DIR}/Python/setup.py
      ${PROJECT_SOURCE_DIR}/README.md
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
  )
endif()

# Targets
install(TARGETS SparrowApp RUNTIME 
  DESTINATION bin
)

if (APPLE)
  # Rename to sparrow upon install
  install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E rename
          ${CMAKE_INSTALL_PREFIX}/bin/SparrowApp \
          ${CMAKE_INSTALL_PREFIX}/bin/sparrow)"
  )
endif()
