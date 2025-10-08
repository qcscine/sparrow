project(RealTimeSpectroscopy)

include(ImportCereal)
import_cereal()

set(headers
    SpectroscopySettings.h
    IR/IRCalculator.h
    IR/IntensitiesCalculator.h
    UvVis/UvVisCalculator.h
    Utils/Spectrum.h
    Utils/LineWidthGenerator.h
    Utils/Observers.h
    )
set(sources
    IR/IRCalculator.cpp
    IR/IntensitiesCalculator.cpp
    UvVis/UvVisCalculator.cpp
    Utils/LineWidthGenerator.cpp
    )

add_library(RealTimeSpectroscopy STATIC ${headers} ${sources})

target_link_libraries(RealTimeSpectroscopy
  PUBLIC
    Scine::Sparrow
    PRIVATE
    ${CMAKE_DL_LIBS}
  )

set_property(TARGET RealTimeSpectroscopy PROPERTY PUBLIC_HEADER ${headers})
add_library(Scine::RealTimeSpectroscopy ALIAS RealTimeSpectroscopy)

#if(SCINE_BUILD_TESTS)
#  add_subdirectory(Tests)
#endif()

install(TARGETS RealTimeSpectroscopy
          DESTINATION lib
        PUBLIC_HEADER
          DESTINATION include)
