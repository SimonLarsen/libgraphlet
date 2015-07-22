cmake_minimum_required(VERSION 3.1)
project(liborca)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(OpenCL 1.1 REQUIRED)

include_directories(
	${CMAKE_SOURCE_DIR}/compute/include
	${CMAKE_SOURCE_DIR}/include
)

add_subdirectory(src/orca)
