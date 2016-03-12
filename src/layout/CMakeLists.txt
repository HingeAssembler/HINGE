cmake_minimum_required(VERSION 3.2)

add_executable(Greedy_best_ovl Greedy_best_ovl_bursty_resistant)
target_link_libraries(Greedy_best_ovl LAInterface ini)

add_executable(NSG NSG)
target_link_libraries(NSG LAInterface ini)

add_executable(hinging hinging)
target_link_libraries(hinging LAInterface ini)