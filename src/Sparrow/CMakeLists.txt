# Load source file lists
include(Files.cmake)

# Find/Import dependecies
include(ImportUtilsOS)
import_utils_os()
include(ImportCereal)
import_cereal()
find_package(OpenMP)

# Shared Sparrow library / Sparrow module
add_library(${PROJECT_NAME} SHARED ${SPARROW_MODULE_FILES})
set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
target_include_directories(${PROJECT_NAME} PUBLIC
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

# include it in the build / move to the appropriate folder for build
target_compile_options(Sparrow PRIVATE $<$<BOOL:${OpenMP_CXX_FOUND}>:${OpenMP_CXX_FLAGS}>)
target_compile_definitions(Sparrow PRIVATE $<$<CONFIG:Debug>:EIGEN_INITIALIZE_MATRICES_BY_NAN>)

# Set the option to repropduce the standard PM6 errors
# This is needed until a new reparametrization of PM6 is done with the right method.
option(REPRODUCE_PM6_ERRORS "Reproduce the errors of standard PM6 implementations" ON)
target_compile_definitions(Sparrow PRIVATE $<$<BOOL:${REPRODUCE_PM6_ERRORS}>:REPRODUCE_PM6_ERRORS>)


target_link_libraries(Sparrow
  PUBLIC
    Scine::Core
    Scine::UtilsOS
  PRIVATE
    cereal::cereal
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


  # Hard-code location of the parameter sets into the test by configuring a header file for it.
  set(parameter_dir "${CMAKE_CURRENT_SOURCE_DIR}/Sparrow/Resources/")
  set(methodTests_dir "${CMAKE_CURRENT_SOURCE_DIR}/Tests/MethodsTests")
  configure_file(${methodTests_dir}/parameters_location.h.in
                 ${methodTests_dir}/parameters_location.h)
endif()

# Find/Import dependecies
include(ImportUtilsOS)
import_utils_os()

# Set the RPATH to be used when installing. $ORIGIN refers to the installation directory of the binary.
set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")

# Executable
find_package(Boost COMPONENTS program_options REQUIRED)

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

# Add link library specific for MinGW
if(WIN32)
  target_link_libraries(SparrowApp PRIVATE $<$<BOOL:MINGW>:ws2_32>)
endif()

# Targets
install(TARGETS SparrowApp RUNTIME 
  DESTINATION bin
)

if (APPLE)
  # Rename to sparrow. upon install
  install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E rename
          ${CMAKE_INSTALL_PREFIX}/bin/SparrowApp \
          ${CMAKE_INSTALL_PREFIX}/bin/sparrow)"
  )
endif()