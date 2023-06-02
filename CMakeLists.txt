cmake_minimum_required(VERSION 3.2)
project(pyKmeansConcaveHull)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O0")

set(INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/build)

set(SOURCE_FILES main.cpp)

add_subdirectory(KmeansConcaveHull)

add_executable(pyKmeansConcaveHull ${SOURCE_FILES})
target_link_libraries(pyKmeansConcaveHull KmeansConcaveHull)
install(TARGETS pyKmeansConcaveHull DESTINATION ${INSTALL_BIN_DIR})