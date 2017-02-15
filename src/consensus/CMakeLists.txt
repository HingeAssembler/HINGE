cmake_minimum_required(VERSION 3.2)

add_executable(draft_assembly draft)
target_link_libraries(draft_assembly LAInterface ini falcon spdlog)

add_executable(consensus consensus.cpp)
target_link_libraries(consensus LAInterface falcon ini)

install(TARGETS draft_assembly consensus DESTINATION ${libexec})
