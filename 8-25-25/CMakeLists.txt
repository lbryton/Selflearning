cmake_minimum_required(VERSION 3.15)
project(LinkedListProject VERSION 1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include path
include_directories(${PROJECT_SOURCE_DIR}/includes)

# Collect all source files
file(GLOB SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

# Add executable
add_executable(linkedlist ${SOURCES})