cmake_minimum_required(VERSION 3.7)
project(FastSLAM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O3")

include_directories(External/matplotlib-cpp)

find_package(Python COMPONENTS Development)
find_package(Eigen3 3.3 REQUIRED)

add_executable(FastSLAM src/main.cpp)

target_include_directories(FastSLAM PRIVATE ${Python_INCLUDE_DIRS})
target_link_libraries(FastSLAM ${Python_LIBRARIES} Eigen3::Eigen)
