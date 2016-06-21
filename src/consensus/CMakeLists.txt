cmake_minimum_required(VERSION 3.2)

add_executable(draft_assembly draft)
target_link_libraries(draft_assembly LAInterface ini falcon spdlog)

add_executable(consensus Draft_consensus.cpp)
target_link_libraries(consensus LAInterface falcon ini)

add_executable(io io_base)
target_link_libraries(io LAInterface ini spdlog)
