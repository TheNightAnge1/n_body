cmake_minimum_required(VERSION 3.14)

set(CMAKE_CXX_COMPILER "icx-cl.exe")
set(CMAKE_C_COMPILER "icx-cl.exe")

add_compile_options(/QxCORE-AVX2 /Qopenmp /Qprec-sqrt)
# -Qopt-report=max -debug /QxCORE-AVX2 /Qopenmp /Qprec-sqrt

project(n_body)

add_subdirectory(ver0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)