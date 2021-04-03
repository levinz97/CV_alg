cmake_minimum_required(VERSION 3.17)
project(simple_filter)

find_package(OpenCV REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})


set(CMAKE_CXX_STANDARD 20)

add_executable(simple_filter main.cpp)

target_link_libraries(simple_filter ${OpenCV_LIBS})