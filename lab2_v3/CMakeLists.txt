cmake_minimum_required(VERSION 3.9)
project(lab2_v3)

set(CMAKE_CXX_STANDARD 11)

add_executable(lab2_v3 main.cpp)
target_link_libraries(lab2_v3 armadillo)