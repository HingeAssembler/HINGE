cmake_minimum_required(VERSION 3.3)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_library(DB DB.h DB.c QV.h QV.c)
add_library(LA align.h align.c)

add_executable(DBshow DBshow)
target_link_libraries(DBshow DB)
target_link_libraries(LA DB)

add_executable(simulator simulator)
target_link_libraries(simulator DB)

add_executable(fasta2DB fasta2DB.c)
target_link_libraries(fasta2DB DB)


add_executable(LAshow LAshow.c)
target_link_libraries(LAshow LA DB)

add_library(LAInterface LAInterface.cpp LAInterface.h)
target_link_libraries(LAInterface LA DB)

add_executable(LAInterface_test LAInterface_test)
target_link_libraries(LAInterface_test LAInterface)