cmake_minimum_required(VERSION 3.2)

add_executable(Draft_assembly Draft_assembly)
target_link_libraries(Draft_assembly LAInterface ini falcon)

add_executable(Circular_consensus Circular_consensus)
target_link_libraries(Circular_consensus LAInterface ini falcon)


add_executable(Draft_consensus Draft_consensus.cpp)
target_link_libraries(Draft_consensus LAInterface falcon ini)

add_executable(Gap_closer Gap_closer.cpp)
target_link_libraries(Gap_closer LAInterface falcon ini)

add_executable(Unitig_consensus Unitig_consensus)
target_link_libraries(Unitig_consensus LAInterface ini falcon)