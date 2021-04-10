# CMake File
cmake_minimum_required(VERSION 3.0)
project(mpm_2d VERSION 0.1)

# set cpp standard
set(CMAKE_CXX_STANDARD 17)

#required on linux systems
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC ")
#message("My Message: ${CMAKE_CXX_FLAGS}")

#Use pybind11
find_package (Python3 COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

# gather cpp files
set(CPP_FILES main.cpp)

# Set include directories
include_directories(mpm_2d PUBLIC ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR} ${pybind11_INCLUDE_DIR} /home/karthik/anaconda3/include/python3.8/
)

#Build MPM Library
add_subdirectory(src)

#Build MPM Library
add_subdirectory(PythonBinding)

#Build PyBinding Library
add_executable(Hello_world Hello_world.cpp)

#Use Cmake to build python module
pybind11_add_module(MPM_2D MPM_2D_module.cpp)

#Link the module against the MPM Library
target_link_libraries(MPM_2D PRIVATE MPM_library)
