add_executable(gdv
	GDV.cpp
)

target_link_libraries(gdv
	orca
)

add_executable(gdv_similarity
	GDVSimilarity.cpp
)

target_link_libraries(gdv_similarity
	orca
)

add_subdirectory(orca)
