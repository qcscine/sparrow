project(RealTimeSpectroscopyTests)

set(sources
    ARealTimeSpectroscopyTest.cpp)

add_executable(RTSpectroscopyTest ${sources})
target_link_libraries(RTSpectroscopyTest
  PRIVATE
    gtest_main
    gmock
    Scine::RealTimeSpectroscopy
    Scine::Core
    $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>
    ${CMAKE_DL_LIBS}
)

  add_test(NAME RealTimeSpectroscopy COMMAND RTSpec_tests)
  target_compile_options(RTSpectroscopyTest
    PUBLIC
      $<TARGET_PROPERTY:Scine::Core,INTERFACE_COMPILE_OPTIONS>
  )
