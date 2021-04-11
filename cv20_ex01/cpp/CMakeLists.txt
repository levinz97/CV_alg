cmake_minimum_required(VERSION 3.17)
project(Hough)

#set(CMAKE_CXX_STANDARD 14)
find_package(OpenCV 4.5.1 REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
        ${Eigen3_INCLUDE_DIRS}
        ${OpenCV_INCLUDE_DIRS})

add_executable(HoughTransform main.cpp)
target_link_libraries(HoughTransform ${OpenCV_LIBS})