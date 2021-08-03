cmake_minimum_required(VERSION 3.19)
project(hanm C)

set(CMAKE_C_STANDARD 99)

find_package(GSL REQUIRED)    # See below (2)
find_package(JSON-C)
#target_link_libraries(Unitsv1 GSL::gsl GSL::gslcblas)

add_executable(hanm HANM.c)