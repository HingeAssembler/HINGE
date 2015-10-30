cmake_minimum_required(VERSION 3.2)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    #set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
    message( "OS X Detected, using gcc49")
    set(CMAKE_C_COMPILER /usr/local/bin/gcc-4.9)
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-4.9)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_AUTOMOC ON)

endif()



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=gnu++11")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_library(ini ini.h ini.c INIReader.cpp INIReader.h)

add_library(DB DB.h DB.c QV.h QV.c)
add_library(LA align.h align.c)

add_library(LAInterface LAInterface.cpp LAInterface.h)
target_link_libraries(LAInterface LA DB)

add_executable(LAInterface_test LAInterface_test)
target_link_libraries(LAInterface_test LAInterface)

add_executable(LA4Falcon_test LA4Falcon)
target_link_libraries(LA4Falcon_test LA DB)

add_library(kmer_lookup kmer_lookup.c)
add_library(DW_banded DW_banded.c)
add_library(falcon falcon.c)
target_link_libraries(falcon kmer_lookup DW_banded)
add_executable(Falcon_test Falcon_test)
target_link_libraries(Falcon_test falcon)

add_executable(Greedy_best_ovl Greedy_best_ovl)
target_link_libraries(Greedy_best_ovl LAInterface)

add_executable(Greedy_best_ovl2 Greedy_best_ovl_bursty_resistant)
target_link_libraries(Greedy_best_ovl2 LAInterface ini)

add_library(OverlapGraph OverlapGraph.cpp)

add_executable(Draft_assembly Draft_assembly)
target_link_libraries(Draft_assembly LAInterface ini falcon)

add_executable(Circular_consensus Circular_consensus)
target_link_libraries(Circular_consensus LAInterface ini falcon)

add_executable(SwissCheeseScrubber Swiss_cheese_scrubber.cpp)
target_link_libraries(SwissCheeseScrubber LAInterface)

add_executable(DenovoQV DenovoQV.cpp)
target_link_libraries(DenovoQV LAInterface)

add_executable(DenovoSV DenovoSV.cpp)
target_link_libraries(DenovoSV LAInterface)

add_executable(omp_test omp_test.c)

add_executable(INIReaderTest INIReaderTest)
target_link_libraries(INIReaderTest ini)
