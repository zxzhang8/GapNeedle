cmake_minimum_required(VERSION 3.21)
project(GapNeedleCpp VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(GAPNEEDLE_BUILD_GUI "Build Qt GUI" ON)
option(GAPNEEDLE_BUILD_CLI "Build CLI" ON)
option(GAPNEEDLE_BUILD_TESTS "Build tests" ON)
option(GAPNEEDLE_USE_MINIMAP2 "Enable minimap2 integration" ON)

if(GAPNEEDLE_BUILD_GUI)
  find_package(Qt6 REQUIRED COMPONENTS Core Widgets Gui Concurrent)
endif()

set(GAPNEEDLE_PUBLIC_HEADERS
  include/gapneedle/types.hpp
  include/gapneedle/aligner.hpp
  include/gapneedle/facade.hpp
  include/gapneedle/fasta_io.hpp
  include/gapneedle/paf.hpp
  include/gapneedle/mapping_service.hpp
  include/gapneedle/stitch_service.hpp
  include/gapneedle/telomere_service.hpp
)

set(GAPNEEDLE_CORE_SOURCES
  src/io/fasta_io.cpp
  src/io/paf_parser.cpp
  src/core/mapping_service.cpp
  src/core/stitch_service.cpp
  src/core/telomere_service.cpp
  src/core/facade.cpp
)

if(GAPNEEDLE_USE_MINIMAP2)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/minimap.h)
    include(cmake/minimap2.cmake)
    add_library(gapneedle_minimap2_bridge STATIC
      src/minimap2_bridge/minimap2_bridge.c
      src/minimap2_bridge/minimap2_aligner.cpp
    )
    target_include_directories(gapneedle_minimap2_bridge
      PUBLIC include
      PRIVATE third_party/minimap2
    )
    target_link_libraries(gapneedle_minimap2_bridge PUBLIC minimap2)
    target_compile_definitions(gapneedle_minimap2_bridge PUBLIC GAPNEEDLE_HAS_MINIMAP2=1)
  else()
    message(WARNING "GAPNEEDLE_USE_MINIMAP2=ON but third_party/minimap2/minimap.h is missing. Falling back to stub aligner.")
    set(GAPNEEDLE_USE_MINIMAP2 OFF)
  endif()
endif()

if(NOT GAPNEEDLE_USE_MINIMAP2)
  add_library(gapneedle_minimap2_bridge STATIC
    src/minimap2_bridge/minimap2_aligner_stub.cpp
  )
  target_include_directories(gapneedle_minimap2_bridge PUBLIC include)
  target_compile_definitions(gapneedle_minimap2_bridge PUBLIC GAPNEEDLE_HAS_MINIMAP2=0)
endif()

add_library(gapneedle_core STATIC
  ${GAPNEEDLE_PUBLIC_HEADERS}
  ${GAPNEEDLE_CORE_SOURCES}
)

target_include_directories(gapneedle_core PUBLIC include)
target_link_libraries(gapneedle_core PUBLIC gapneedle_minimap2_bridge)

if(GAPNEEDLE_BUILD_CLI)
  add_executable(gapneedle_cli src/cli/main.cpp)
  target_link_libraries(gapneedle_cli PRIVATE gapneedle_core)
endif()

if(GAPNEEDLE_BUILD_GUI)
  add_executable(gapneedle_gui
    src/gui/main.cpp
    src/gui/ui_theme.cpp
    src/gui/ui_components.cpp
    src/gui/resources/gui_resources.qrc
    src/gui/main_window.cpp
    src/gui/align_page.cpp
    src/gui/paf_viewer_page.cpp
    src/gui/manual_stitch_page.cpp
    src/gui/fasta_search_page.cpp
  )
  target_include_directories(gapneedle_gui PRIVATE src/gui)
  set_target_properties(gapneedle_gui PROPERTIES AUTOMOC ON AUTOUIC ON AUTORCC ON)
  target_link_libraries(gapneedle_gui PRIVATE gapneedle_core Qt6::Core Qt6::Gui Qt6::Widgets Qt6::Concurrent)
endif()

if(GAPNEEDLE_BUILD_TESTS)
  enable_testing()
  add_executable(gapneedle_tests tests/test_main.cpp)
  target_link_libraries(gapneedle_tests PRIVATE gapneedle_core)
  add_test(NAME gapneedle_tests COMMAND gapneedle_tests)
endif()
