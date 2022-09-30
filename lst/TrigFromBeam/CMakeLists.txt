add_library(
  TrigFromBeam MODULE TrigFromBeam.cpp
)

target_include_directories(
  TrigFromBeam
  PUBLIC ${TriDAS_SOURCE_DIR}/DAQ/TCPU/inc
  PUBLIC ${TriDAS_SOURCE_DIR}/packages/DataFormats
)

target_link_libraries(
  TrigFromBeam
  PUBLIC Packages
)

install(
  TARGETS TrigFromBeam
  LIBRARY DESTINATION trigger_plugins
)
