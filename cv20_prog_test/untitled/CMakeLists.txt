cmake_minimum_required(VERSION 3.1)
project(untitled)
find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
set(CMAKE_CXX_STANDARD 14)

add_executable(untitled main.cpp)
target_link_libraries(untitled ${OpenCV_LIBS})