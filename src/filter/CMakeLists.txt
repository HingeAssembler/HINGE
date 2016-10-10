cmake_minimum_required(VERSION 3.2)

add_executable(Reads_filter filter)
target_link_libraries(Reads_filter LAInterface ini spdlog)

install(TARGETS Reads_filter DESTINATION ${libexec})
