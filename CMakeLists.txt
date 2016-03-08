cmake_minimum_required(VERSION 3.3)
project(molly)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# include Eigen matrix library
find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

# make library
set(SOURCE_FILES Molecule.cpp Molecule.h)
add_library(molly ${SOURCE_FILES})
target_link_libraries(molly)

# add test folder
add_subdirectory(tests)

add_executable(mollymain main.cpp)
target_link_libraries(molly)