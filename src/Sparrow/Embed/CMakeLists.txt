add_executable(Embed
  ${CMAKE_CURRENT_SOURCE_DIR}/embed.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Basisfile.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/DftbParameterSet.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/NddoParameters.cpp
)
target_include_directories(Embed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Embed
  PRIVATE Scine::Sparrow
)
