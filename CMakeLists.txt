cmake_minimum_required(VERSION 2.8.10)
project(histogram)
enable_language(Fortran)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(MPI REQUIRED)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})
link_directories(${MPI_Fortran_LIBRARIES})

file(GLOB SOURCES *.f90 *.F90 *.c)

add_executable(histogram ${SOURCES})
target_link_libraries(histogram ${MPI_Fortran_LIBRARIES})

list(REMOVE_ITEM SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.f90)
enable_testing()

add_executable(test_histogram tests/test_histogram.c ${SOURCES})
target_link_libraries(test_histogram ${MPI_Fortran_LIBRARIES})
add_test(test_histogram test_histogram)