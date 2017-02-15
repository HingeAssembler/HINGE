cmake_minimum_required(VERSION 3.2)

#if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#    set(WINDOWS TRUE)
#elseif (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#    set(LINUX TRUE)
#    message( "Linux Detected, using gcc48")
#    set(CMAKE_C_COMPILER /usr/bin/gcc-4.8)
#    set(CMAKE_CXX_COMPILER /usr/bin/g++-4.8)
#elseif (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#    set(MACOSX TRUE)
#    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/build")
#    message( "OS X Detected, using gcc49")
#    set(CMAKE_C_COMPILER /usr/local/bin/gcc-4.9)
#    set(CMAKE_CXX_COMPILER /usr/local/bin/g++-4.9)
#    set(CMAKE_INCLUDE_CURRENT_DIR ON)
#    set(CMAKE_AUTOMOC ON)
#endif()



set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=gnu++11")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

include_directories(include)

## Libraries

add_subdirectory(lib lib)

add_subdirectory(spdlog)

#add_subdirectory(ogdf)

find_program(PANDOC pandoc)

if(PANDOC)
  add_custom_command(
    OUTPUT hinge.1
    COMMAND ${PANDOC} --from markdown --to man -s ${CMAKE_CURRENT_SOURCE_DIR}/hinge.1.md -o hinge.1
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/hinge.1.md
    VERBATIM
    )

  add_custom_target(man ALL DEPENDS hinge.1)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/hinge.1 DESTINATION share/man/man1)
endif()

install(PROGRAMS hinge DESTINATION bin)

## Executables

add_subdirectory(filter filter)

add_subdirectory(maximal maximal)

add_subdirectory(consensus consensus)

add_subdirectory(layout layout)

## Tests

add_subdirectory(test)
