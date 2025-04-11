cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(ODELandscaper 
        VERSION 0.3 
        LANGUAGES CXX
)

configure_file(
    "include/version.h.in"
    "version.h"
)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-O3 -pthread -fopenmp")


file(GLOB_RECURSE SOURCES src/*.cpp)
list(FILTER SOURCES EXCLUDE REGEX "src/ODESystem.cpp")

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PUBLIC
"${PROJECT_BINARY_DIR}"
"${PROJECT_SOURCE_DIR}/include"
)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)


# Test on some data
# add testing so it can be called using make test
include(CTest)
enable_testing()
message(STATUS "${PROJECT_SOURCE_DIR}")

# test NAR
add_test(
  NAME testNAR
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_NARPAR.csv" -o "${PROJECT_SOURCE_DIR}/tests/test_NAR_out.csv" -s "NAR" -t 4 -p 0.406,0.44 -w 0.05,0.05
)
# test PAR
add_test(
  NAME testPAR
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_NARPAR.csv" -o "${PROJECT_SOURCE_DIR}/tests/test_PAR_out.csv" -s "PAR" -t 4 -p 0.28,0.161 -w 0.05,0.05
)

# test FFLC1
add_test(
  NAME testFFLC1
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_FFL.csv" -o "${PROJECT_SOURCE_DIR}/tests/test_FFLC1_out.csv" -s "FFLC1" -t 4 -p 1.304,0.1,0.157 -w 0.05,0.05,0.05
)

# test FFLI1
add_test(
  NAME testFFLI1
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_FFL.csv" -o "${PROJECT_SOURCE_DIR}/tests/test_FFLI1_out.csv" -s "FFLI1" -t 4 -p 0.466,0.336,5.2 -w 0.05,0.05,0.05
)

# test FFBH
add_test(
  NAME testFFBH
  COMMAND ${PROJECT_NAME} -i "${PROJECT_SOURCE_DIR}/tests/test_FFBH.csv" -o "${PROJECT_SOURCE_DIR}/tests/test_FFBH_out.csv" -s "FFBH" -t 4 -p 1.536,0.179,0.933,0.006 -w 0.05,0.05,0.05,0.05
)
