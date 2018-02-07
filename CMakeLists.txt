cmake_minimum_required(VERSION 2.8)
project(libosmtools)

set(MY_CXX_FLAGS "-std=gnu++14 -Wall -Wextra -Wstrict-overflow=3")
set(MY_C_FLAGS "-std=c99 -Wall -Wextra -Wstrict-overflow=3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MY_CXX_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MY_C_FLAGS}")

find_package(CGAL)

set(LIBOSMTOOLS_INCLUDE_DIR
	"${CMAKE_CURRENT_SOURCE_DIR}/include"
	${SSERIALIZE_INCLUDE_DIR}
	CACHE STRING "libosmtools include directories"
	FORCE
)

set(LIBOSMTOOLS_LIBRARIES
	${PROJECT_NAME}
	CACHE STRING "libosmtools libraries"
	FORCE
)

set(LIBOSMTOOLS_LINK_LIBRARIES
	${PROJECT_NAME}
	${SSERIALIZE_LINK_LIBRARIES}
	"CGAL"
	"CGAL_Core"
	"boost_thread"
	"boost_system"
	CACHE STRING "libosmtools link libraries"
	FORCE
)

mark_as_advanced(
	LIBOSMTOOLS_INCLUDE_DIR
	LIBOSMTOOLS_LIBRARIES
	LIBOSMTOOLS_LINK_LIBRARIES
)

include_directories(
	${LIBOSMTOOLS_INCLUDE_DIR}
)


file(GLOB_RECURSE INCLUDE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

set(LIB_SOURCES_CPP
	src/AreaExtractor.cpp
	src/AreaExtractorFilters.cpp 
	src/types.cpp
	src/OsmTriangulationRegionStore.cpp
)

add_library(${PROJECT_NAME} STATIC
	${LIB_SOURCES_CPP}
	${INCLUDE_HEADERS}
)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${SSERIALIZE_COMPILE_DEFINITIONS})
