cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(ODELandscaper 
        VERSION 0.1 
        LANGUAGES CXX
)

configure_file(
    "include/version.h.in"
    "version.h"
)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-O3 -pthread")


add_executable(${PROJECT_NAME} src/main.cpp src/odePar.cpp src/ODESystem.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}/include"
)
message(${PROJECT_SOURCE_DIR})