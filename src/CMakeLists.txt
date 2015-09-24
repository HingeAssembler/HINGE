cmake_minimum_required(VERSION 3.3)


add_library(DB DB.h DB.c QV.h QV.c)

add_executable(DBshow DBshow)
target_link_libraries(DBshow DB)

add_library(LA align.h align.c)
target_link_libraries(LA DB)

add_executable(LAshow LAshow)
target_link_libraries(LAshow LA)

add_library(LAInterface LAInterface)
target_link_libraries(LAInterface LA)

add_executable(LAInterface_test LAInterface_test)
target_link_libraries(LAInterface_test LA LAInterface)
