set(ORCA_SOURCES
	Orca.cpp
	Similarity.cpp
)

if(NOT ORCA_NO_OPENCL)
	list(APPEND ORCA_SOURCES SimilarityGPU.cpp)
endif()

add_library(orca
	${ORCA_SOURCES}
)
