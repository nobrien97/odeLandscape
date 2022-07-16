cmake_minimum_required(VERSION 3.1)
project(ODELandscaper)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-O3")
include_directories("include")
add_executable(ODESolver src/main.cpp)
target_include_directories(ODESolver PUBLIC include)