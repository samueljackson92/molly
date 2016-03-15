cmake_minimum_required(VERSION 3.3)
project(molly)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# include Eigen matrix library
find_package(Eigen3 REQUIRED)
include_directories(EIGEN3_INCLUDE_DIR)

# include VTK library
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# make library
set(SOURCE_FILES Molecule.cpp Molecule.h Simulation.cpp Simulation.h SystemVisualisation.cpp SystemVisualisation.h SystemParameters.cpp SystemParameters.h Property.cpp Property.h Cell.cpp Cell.h CellMatrix.cpp CellMatrix.h)
add_library(molly STATIC ${SOURCE_FILES})
target_link_libraries(molly ${VTK_LIBRARIES})

# add test folder
add_subdirectory(tests)

add_executable(mollymain main.cpp)
target_link_libraries(mollymain molly)