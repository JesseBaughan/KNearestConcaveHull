
cmake_minimum_required(VERSION 3.2)
project(pyKmeansConcaveHull)

set(SOURCE_FILES
    KMeansConcaveHull.cpp
)

add_library(KMeansConcaveHull SHARED STATIC ${SOURCE_FILES})
install(TARGETS KMeansConcaveHull DESTINATION ${LIB_DIR})