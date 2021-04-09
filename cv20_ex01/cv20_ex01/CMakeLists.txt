cmake_minimum_required(VERSION 3.17)
project(cv20_ex01 C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_FLAGS "-O3")
add_executable(gradient_map gradient_map.c)
target_link_libraries(gradient_map m)