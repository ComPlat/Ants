cmake_minimum_required(VERSION 3.13)
project(Ants C)

add_library(calc_energy STATIC xtb/calc_energy.c)
target_include_directories(calc_energy
  PUBLIC
  ${CMAKE_SOURCE_DIR}/xtb
  $ENV{HOME}/miniconda3/include
)
target_link_directories(calc_energy
  PUBLIC
  $ENV{HOME}/miniconda3/lib
)
target_link_libraries(calc_energy
  PUBLIC
  xtb gfortran m
)

add_library(random STATIC random/mersenne_twister.c)
target_include_directories(random
  PUBLIC
  ${CMAKE_SOURCE_DIR}/random
)

add_library(pso STATIC pso/pso.c)
target_include_directories(pso
  PUBLIC
  ${CMAKE_SOURCE_DIR}/random
  ${CMAKE_SOURCE_DIR}/pso
)
target_link_directories(pso
  PUBLIC
  random
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
  ${CMAKE_SOURCE_DIR}/random
)
target_link_libraries(Ants
  PUBLIC
  pso
  random
  calc_energy
  xtb gfortran m
)


add_executable(pso_tests pso/tests.c)
target_include_directories(pso_tests
  PUBLIC
  ${CMAKE_SOURCE_DIR}/xtb
  ${CMAKE_SOURCE_DIR}/pso
)
target_link_libraries(pso_tests
  pso
  random
)
