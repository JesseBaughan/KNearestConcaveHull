set(BINARY ${CMAKE_PROJECT_NAME}_tst)

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)

set(SOURCES ${TEST_SOURCES})

add_executable(${BINARY} ${TEST_SOURCES})

add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(${BINARY} PUBLIC KMeansConcaveHull gtest)

target_include_directories(${BINARY} PUBLIC ../KMeansConcaveHull)
target_include_directories(${BINARY} PUBLIC ../lib)