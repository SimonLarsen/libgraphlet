add_executable(gdv
	GDV.cpp
)

target_link_libraries(gdv
	orca
)

add_subdirectory(orca)
