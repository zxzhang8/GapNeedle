add_library(minimap2 STATIC
  bseq.c
  format.c
  index.c
  kalloc.c
  kthread.c
  map.c
  misc.c
  options.c
  sdust.c
  sketch.c
  splitidx.c
)

target_include_directories(minimap2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(minimap2 PUBLIC _GNU_SOURCE)
