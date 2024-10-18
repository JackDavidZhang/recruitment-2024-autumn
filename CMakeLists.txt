cmake_minimum_required(VERSION 3.22.1)
project(heptagon_24autumn_recruitment LANGUAGES CUDA CXX)

find_package(CUDA REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS on)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp")

include_directories(include)

file(GLOB SOURCE_FILES "src/*.cc" "main.cc")

cuda_add_library(CUDA_lib src/cuda.cu
        include/cuda.h)

add_executable (SmithWaterman ${SOURCE_FILES})

target_link_libraries(SmithWaterman CUDA_lib)

enable_testing()

set(dataset ${CMAKE_SOURCE_DIR}/data)

add_test(NAME "test case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/test_q.fna ${dataset}/test_t.fna ${dataset}/test_ref.txt)

add_test(NAME "tiny case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/tiny_q.fna ${dataset}/tiny_t.fna ${dataset}/tiny_ref.txt)
add_test(NAME "small case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/small_q.fna ${dataset}/small_t.fna ${dataset}/small_ref.txt)
add_test(NAME "medium case" COMMAND ${CMAKE_BINARY_DIR}/SmithWaterman ${dataset}/medium_q.fna ${dataset}/medium_t.fna ${dataset}/medium_ref.txt)
