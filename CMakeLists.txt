cmake_minimum_required(VERSION 3.13)
project(Ants C)

add_library(pso STATIC pso.c)

target_include_directories(pso
  PUBLIC
  ${CMAKE_SOURCE_DIR}
  $ENV{HOME}/miniconda3/include
)

target_link_directories(pso
  PUBLIC
  $ENV{HOME}/miniconda3/lib
)
target_link_libraries(pso
  PUBLIC
  xtb gfortran m
)

add_executable(Ants main.c)
target_include_directories(Ants
  PUBLIC
  ${CMAKE_SOURCE_DIR}
  $ENV{HOME}/miniconda3/include
)
target_link_directories(Ants
  PUBLIC
  $ENV{HOME}/miniconda3/lib
)
target_link_libraries(Ants
  PUBLIC
  pso
  xtb gfortran m
)
