set(LIBGRAPHLET_SOURCES
	GDD.cpp
	Similarity.cpp
)

if(LIBGRAPHLET_WITH_OPENCL)
	list(APPEND LIBGRAPHLET_SOURCES
		SimilarityGPU.cpp
	)
endif()

add_library(graphlet
	${LIBGRAPHLET_SOURCES}
)
