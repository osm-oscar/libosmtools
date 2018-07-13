cmake_minimum_required(VERSION 3.9)
project(libosmtools)
find_package(CGAL REQUIRED)

set(MY_INCLUDE_DIRS
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
)

set(MY_LINK_LIBRARIES
	sserialize
	osmpbf
)

file(GLOB_RECURSE LIB_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

set(LIB_SOURCES_CPP
	src/AreaExtractor.cpp
	src/AreaExtractorFilters.cpp 
	src/types.cpp
	src/OsmTriangulationRegionStore.cpp
	src/OsmGridRegionTree.cpp
	src/CellCriteria.cpp
)

add_library(${PROJECT_NAME} STATIC
	${LIB_SOURCES_CPP}
	${LIB_HEADERS}
)
target_include_directories(${PROJECT_NAME} PUBLIC ${MY_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${MY_LINK_LIBRARIES})
