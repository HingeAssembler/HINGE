cmake_minimum_required(VERSION 3.2)
project(AwAssembler)

set(libexec "lib/hinge")
add_subdirectory(src bin)

install(
DIRECTORY scripts/
DESTINATION ${libexec}
USE_SOURCE_PERMISSIONS
)
