cmake_minimum_required(VERSION 3.19)
project(test LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_FLAGS "-Ofast -m64 -march=native -fopenmp -DLIBBSC_OPENMP_SUPPORT")
set(CMAKE_CXX_FLAGS "-Ofast -m64 -march=native -fopenmp -DLIBBSC_OPENMP_SUPPORT")

add_executable(bench benchmark.cpp)
target_link_libraries(bench PRIVATE stdc++fs)