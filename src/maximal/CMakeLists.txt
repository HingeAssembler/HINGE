cmake_minimum_required(VERSION 3.2)

add_executable(get_maximal_reads maximal)
target_link_libraries(get_maximal_reads LAInterface ini spdlog)

install(TARGETS get_maximal_reads DESTINATION ${libexec})
