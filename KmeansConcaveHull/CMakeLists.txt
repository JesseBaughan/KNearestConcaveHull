
cmake_minimum_required(VERSION 3.2)
project(KmeansConcaveHull)

set(SOURCE_FILES
    KMeansConcaveHull.hpp
    KMeansConcaveHull.cpp
)

add_library(KmeansConcaveHull SHARED STATIC ${SOURCE_FILES})
install(TARGETS KmeansConcaveHull DESTINATION ${LIB_DIR})