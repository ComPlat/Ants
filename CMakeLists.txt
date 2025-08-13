cmake_minimum_required(VERSION 3.13)
project(Ants)

set(TEST_FILES_DIR "${CMAKE_SOURCE_DIR}/test_files")

# convert zmat2xyz
add_library(convert STATIC convert/convert.c)
target_include_directories(convert
  PUBLIC
  ${CMAKE_SOURCE_DIR}/convert
)
target_link_libraries(convert
  PUBLIC
  m
)

# calc energy
add_library(calc_energy STATIC calc_energy/calc_energy.c)
target_include_directories(calc_energy
  PUBLIC
  ${CMAKE_SOURCE_DIR}/calc_energy
  ${CMAKE_SOURCE_DIR}/convert
  $ENV{HOME}/miniconda3/include
)
target_link_directories(calc_energy
  PUBLIC
  $ENV{HOME}/miniconda3/lib
)
target_link_libraries(calc_energy
  PUBLIC
  xtb gfortran m convert
)

# mersenne twister
add_library(random STATIC random/mersenne_twister.c)
target_include_directories(random
  PUBLIC
  ${CMAKE_SOURCE_DIR}/random
)

# pso
add_library(pso STATIC pso/pso.c)
target_include_directories(pso
  PUBLIC
  ${CMAKE_SOURCE_DIR}/random
  ${CMAKE_SOURCE_DIR}/pso
)
target_link_libraries(pso
  PUBLIC
  random
  m
)

# main
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
  ${CMAKE_SOURCE_DIR}/convert
)
target_link_libraries(Ants
  PUBLIC
  pso
  random
  calc_energy
  convert
  xtb gfortran m
)

# pso tests
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

# convert tests
add_executable(convert_tests convert/test_convert.c)
target_include_directories(convert_tests
  PUBLIC
  ${CMAKE_SOURCE_DIR}/convert
)
target_link_libraries(convert_tests
  convert
  m
)
add_custom_command(TARGET convert_tests POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/convert/input.zmat
    $<TARGET_FILE_DIR:convert_tests>/input.zmat
)
add_custom_command(TARGET convert_tests POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy
    ${CMAKE_SOURCE_DIR}/convert/expected.xyz
    $<TARGET_FILE_DIR:convert_tests>/expected.xyz
)

# xtb test
add_executable(xtb_test calc_energy/test_xtb_c_api.c)
target_include_directories(xtb_test
  PUBLIC
  $ENV{HOME}/miniconda3/include
)
target_link_directories(xtb_test
  PUBLIC
  $ENV{HOME}/miniconda3/lib
)
target_link_libraries(xtb_test
  PUBLIC
  xtb gfortran m convert
)
