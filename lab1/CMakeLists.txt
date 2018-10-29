cmake_minimum_required(VERSION 3.2)


set(CMAKE_C_STANDARD 99)



add_executable(gauss gaussian_elimination.c)
add_executable(cond cond_calculation.cpp)
target_link_libraries(gauss m)
target_link_libraries(cond armadillo)

