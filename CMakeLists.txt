cmake_minimum_required(VERSION 3.1)
project(liborca)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected. Defaulting to Release.")
	set(CMAKE_BUILD_TYPE "Release")
endif()

set(LIBGRAPHLET_WITH_OPENCL false CACHE BOOL "Compile with OpenCL")

find_package(Boost REQUIRED)

if(LIBGRAPHLET_WITH_OPENCL)
	find_package(OpenCL 1.1 REQUIRED)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")

include_directories(
	${CMAKE_SOURCE_DIR}/tclap/include
	${CMAKE_SOURCE_DIR}/libgraph/include
	${CMAKE_SOURCE_DIR}/compute/include
	${CMAKE_SOURCE_DIR}/include
	${Boost_INCLUDE_DIRS}
	${OpenCL_INCLUDE_DIRS}
)

add_executable(gdv
	${CMAKE_SOURCE_DIR}/src/GDV.cpp
)

target_link_libraries(gdv
	orca
)

add_executable(gdv_similarity
	${CMAKE_SOURCE_DIR}/src/GDVSimilarity.cpp
)

target_link_libraries(gdv_similarity
	orca
	graphlet
	pthread
)

add_executable(gdd
	${CMAKE_SOURCE_DIR}/src/GDD.cpp
)

target_link_libraries(gdd
	orca
	graphlet
)

add_executable(gdd_agreement
	${CMAKE_SOURCE_DIR}/src/GDDAgreement.cpp
)

target_link_libraries(gdd_agreement
	orca
	graphlet
)

add_subdirectory(src/orca)
add_subdirectory(src/libgraphlet)
