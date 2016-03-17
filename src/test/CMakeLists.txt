cmake_minimum_required(VERSION 3.2)

enable_testing()

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(Interface_test LAInterface_test1.cpp)

target_link_libraries(Interface_test gtest gtest_main)

target_link_libraries(Interface_test LAInterface)

add_test(Interface_test_1 Interface_test)