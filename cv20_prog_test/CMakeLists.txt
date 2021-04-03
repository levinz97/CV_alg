cmake_minimum_required(VERSION 3.17)
project(cv20_prog_test C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(simple_filter simple_filter.c)
target_link_libraries(simple_filter m)