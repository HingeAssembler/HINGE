cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=gnu++11")

add_library(DB DB.h DB.c QV.h QV.c)
add_library(LA align.h align.c)

add_library(LAInterface LAInterface.cpp LAInterface.h)
target_link_libraries(LAInterface LA DB)

add_executable(LAInterface_test LAInterface_test)
target_link_libraries(LAInterface_test LAInterface)

add_executable(LA4Falcon LA4Falcon)
target_link_libraries(LA4Falcon LA DB)

add_executable(Falcon kmer_lookup.c falcon.c DW_banded.c )

add_executable(Greedy_best_ovl Greedy_best_ovl)
target_link_libraries(Greedy_best_ovl LAInterface)

add_executable(Greedy_best_ovl2 Greedy_best_ovl_bursty_resistant)
target_link_libraries(Greedy_best_ovl2 LAInterface)


add_library(OverlapGraph OverlapGraph.cpp)

add_executable(Draft_assembly Draft_assembly)
target_link_libraries(Draft_assembly LAInterface)

add_executable(SwissCheeseScrubber Swiss_cheese_scrubber.cpp)
target_link_libraries(SwissCheeseScrubber LAInterface)

add_executable(DenovoQV DenovoQV.cpp)
target_link_libraries(DenovoQV LAInterface)

add_executable(DenovoSV DenovoSV.cpp)
target_link_libraries(DenovoSV LAInterface)