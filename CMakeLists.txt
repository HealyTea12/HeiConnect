cmake_minimum_required(VERSION 3.20...3.30)
project(HeiConnect VERSION 0.1.0 LANGUAGES CXX)

set(CXX True)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

if("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed. Please create a separate build directory.")
endif()

option(ENABLE_VERBOSE_LOG "Enable --verbose logging flag" ON)
option(HEICONNECT_BUILD_TESTS "Build tests" ON )
option(HEICONNECT_BUILD_APP "Build application" OFF )
option(HEICONNECT_BUILD_EXPERIMENTS "Build experiments" OFF )

include(FetchContent)

find_package(GUROBI REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)
FetchContent_Declare(
  pugixml
  GIT_REPOSITORY https://github.com/zeux/pugixml
  GIT_TAG master
)
FetchContent_MakeAvailable(pugixml)

set(heiconnect-patterns-sources src/*.cpp src/set_cover/*.cpp src/min_cut/*.cpp)
file(GLOB heiconnect-patterns-sources heiconnect-sources)

add_library(heiconnect STATIC "")
target_include_directories(heiconnect
  PUBLIC 
    ${PROJECT_SOURCE_DIR}/src/include
  PRIVATE
    ${PROJECT_SOURCE_DIR}/src 
    ${PROJECT_SOURCE_DIR}/src/include/HeiConnect
    extern/viecut
)
target_sources(heiconnect 
PRIVATE 
  src/approx.cpp src/graph.cpp src/greedy.cpp src/ilp.cpp src/ring_graph.cpp src/util.cpp
  src/set_cover/set_cover.cpp  src/min_cut/simple_mincut.cpp )
target_link_libraries(heiconnect
PUBLIC
  pugixml
PRIVATE
  gurobi
  Boost::boost
  OpenMP::OpenMP_CXX
)

#add_subdirectory(src)
if(HEICONNECT_BUILD_APP)
  add_subdirectory(app)
endif()
if(HEICONNECT_BUILD_TESTS)
  add_subdirectory(tests)
endif()
# this could be completely separate
if(HEICONNECT_BUILD_EXPERIMENTS)
  add_subdirectory(experiments)
endif()