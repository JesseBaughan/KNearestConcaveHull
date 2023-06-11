cmake_minimum_required(VERSION 3.2)
project(pyKmeansConcaveHull)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -O0")
set(BINARY pyKmeansConcaveHull)

set(INSTALL_BIN_DIR ${PROJECT_SOURCE_DIR}/build)

set(SOURCE_FILES main.cpp)

add_subdirectory(KMeansConcaveHull)
add_subdirectory(test)
add_subdirectory(lib/googletest)

add_executable(${BINARY} ${SOURCE_FILES})
target_link_libraries(${BINARY} PUBLIC KMeansConcaveHull)

install(TARGETS ${BINARY} DESTINATION ${INSTALL_BIN_DIR})