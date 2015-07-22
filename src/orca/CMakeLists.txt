set(ORCA_SOURCES
	Orca.cpp
	Similarity.cpp
)

if(ORCA_WITH_OPENCL)
	list(APPEND ORCA_SOURCES SimilarityGPU.cpp)
endif()

add_library(orca
	${ORCA_SOURCES}
)
