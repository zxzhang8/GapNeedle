find_package(Threads REQUIRED)
find_package(ZLIB QUIET)

add_library(minimap2 STATIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/kthread.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/kalloc.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/misc.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/bseq.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/sketch.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/sdust.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/options.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/index.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/lchain.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/align.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/hit.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/seed.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/jump.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/map.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/format.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/pe.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/esterr.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/splitidx.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/ksw2_ll_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/ksw2_extz2_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/ksw2_extd2_sse.c
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2/ksw2_exts2_sse.c
)

target_include_directories(minimap2 PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/minimap2
)
target_compile_definitions(minimap2 PUBLIC
  _GNU_SOURCE
  HAVE_KALLOC
)
target_link_libraries(minimap2 PUBLIC Threads::Threads)
if(ZLIB_FOUND)
  target_link_libraries(minimap2 PUBLIC ZLIB::ZLIB)
elseif(MINGW)
  # MinGW toolchains usually provide zlib as -lz even when FindZLIB can't locate package metadata.
  target_link_libraries(minimap2 PUBLIC z)
  message(WARNING "FindZLIB failed; falling back to linking -lz for MinGW.")
else()
  message(FATAL_ERROR "zlib is required to build minimap2. Please install zlib development files.")
endif()
if(NOT WIN32)
  target_link_libraries(minimap2 PUBLIC m)
endif()
