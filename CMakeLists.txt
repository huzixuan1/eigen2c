cmake_minimum_required(VERSION 3.10)
project(c_eigen_replacement C)

set(CMAKE_C_STANDARD 99)

include_directories(include)

file(GLOB_RECURSE SRC_FILES "src/*.c")

add_library(c_eigen STATIC ${SRC_FILES})

file(GLOB_RECURSE TEST_FILES "test/*.c")
foreach(test_file ${TEST_FILES})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})
    target_link_libraries(${test_name} c_eigen)
endforeach()
