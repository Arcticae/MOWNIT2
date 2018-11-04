cmake_minimum_required(VERSION 3.9)
project(lab2)

set(CMAKE_CXX_STANDARD 14)
add_executable(zad1 jacobi_zad1_var1.cpp)
add_executable(zad2 jacobi_zad1_var2.cpp)
target_link_libraries(zad1 armadillo)
target_link_libraries(zad2 armadillo)
