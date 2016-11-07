cmake_minimum_required(VERSION 3.2)

set(Boost_USE_STATIC_LIBS   ON)

FIND_PACKAGE( Boost COMPONENTS graph REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

add_executable(hinging hinging)
target_link_libraries(hinging LAInterface ini spdlog ${Boost_LIBRARIES})

install(TARGETS hinging DESTINATION ${libexec})
