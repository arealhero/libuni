cmake_minimum_required (VERSION 3.0)

set (CMAKE_CXX_STANDARD 20)

project (libuni
  VERSION 0.0.1
  LANGUAGES CXX
)

add_library (libuni INTERFACE)
target_include_directories (libuni INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
