cmake_minimum_required(VERSION 3.0.0)
project(percolation VERSION 0.1.0)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#include(CTest)
#enable_testing()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   message(   "Clang is used" )
   set(CMAKE_CXX_FLAGS "-O2  -Wall -Wextra ")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   message(   "gcc is used" )
   set(CMAKE_CXX_FLAGS "-O2  -march=native -Wall -Wextra -pedantic -Wshadow ")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
   message(   "intel is used" )
#elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
# using Visual Studio C++
else()
   message( FATAL_ERROR  "Not supported" )
endif()

find_package(Boost  REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

add_executable(percolation_analysis_3d main3d.cpp)
include_directories(${PRETTYPRINT})
include_directories(~/src/cxx-prettyprint)

target_link_libraries(percolation_analysis_3d Boost::program_options)



