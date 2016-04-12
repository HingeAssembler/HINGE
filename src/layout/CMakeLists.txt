cmake_minimum_required(VERSION 3.2)

#include_directories("../ogdf/include")

set(Boost_USE_STATIC_LIBS   ON)

FIND_PACKAGE( Boost COMPONENTS graph REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(Greedy_best_ovl Greedy_best_ovl_bursty_resistant)
target_link_libraries(Greedy_best_ovl LAInterface ini)

add_executable(NSG NSG)
target_link_libraries(NSG LAInterface ini)

add_executable(hinging hinging)
target_link_libraries(hinging LAInterface ini spdlog ${Boost_LIBRARIES})
