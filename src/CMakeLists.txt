cmake_minimum_required(VERSION 3.2)

if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set(WINDOWS TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    set(LINUX TRUE)
elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
    message( "OS X Detected, using gcc49")
    set(CMAKE_C_COMPILER /usr/local/bin/gcc-4.9)
    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-4.9)
    set(CMAKE_INCLUDE_CURRENT_DIR ON)
    set(CMAKE_AUTOMOC ON)
endif()



set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=gnu++11")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


## Libraries


add_library(ini ini.h ini.c INIReader.cpp INIReader.h)

add_library(DB DB.h DB.c QV.h QV.c)
add_library(LA align.h align.c)

add_library(LAInterface LAInterface.cpp LAInterface.h)
target_link_libraries(LAInterface LA DB)


add_library(kmer_lookup kmer_lookup.c)
add_library(DW_banded DW_banded.c)
add_library(falcon falcon.c falcon_utils.c)
target_link_libraries(falcon kmer_lookup DW_banded)




## Executables


add_executable(Greedy_best_ovl Greedy_best_ovl_bursty_resistant)
target_link_libraries(Greedy_best_ovl LAInterface ini)

add_executable(NSG NSG_v3)
target_link_libraries(NSG LAInterface ini)

add_executable(hinging hinging_v1)
target_link_libraries(hinging LAInterface ini)

add_library(OverlapGraph OverlapGraph.cpp)

add_executable(Draft_assembly Draft_assembly)
target_link_libraries(Draft_assembly LAInterface ini falcon)

add_executable(Circular_consensus Circular_consensus)
target_link_libraries(Circular_consensus LAInterface ini falcon)

add_executable(DenovoQV DenovoQV.cpp) # TBD
target_link_libraries(DenovoQV LAInterface)

add_executable(DenovoSV DenovoSV.cpp) # TBD
target_link_libraries(DenovoSV LAInterface)

add_executable(Draft_consensus Draft_consensus.cpp)
target_link_libraries(Draft_consensus LAInterface falcon ini)

add_executable(Gap_closer Gap_closer.cpp)
target_link_libraries(Gap_closer LAInterface falcon ini)

add_executable(Unitig_consensus Unitig_consensus)
target_link_libraries(Unitig_consensus LAInterface ini falcon)

## Tests
add_executable(LAInterface_test LAInterface_test)
target_link_libraries(LAInterface_test LAInterface)

add_executable(LAInterface_test_2DB LAInterface_test_2DB)
target_link_libraries(LAInterface_test_2DB LAInterface)

add_executable(LA4Falcon_test LA4Falcon)
target_link_libraries(LA4Falcon_test LA DB)

add_executable(Falcon_test Falcon_test)
target_link_libraries(Falcon_test falcon)

add_executable(omp_test omp_test.c)

add_executable(INIReaderTest INIReaderTest)
target_link_libraries(INIReaderTest ini)

add_executable(Consensus_test LAInterface_consensus_test.cpp)
target_link_libraries(Consensus_test LAInterface falcon)

add_executable(Reads_filter filter.cpp)
target_link_libraries(Reads_filter LAInterface ini)