cmake_minimum_required(VERSION 2.6)

project(VItA_lib)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(DOWNLOAD_DEPENDENCIES "ON if you prefere to download VTK automatically" ON)

if(DOWNLOAD_DEPENDENCIES)
	include(dependencies.cmake)
	link_directories(${VTK_LIB_DIR})
	include_directories(${VTK_INCLUDE_DIR})
#else()
#	set(VTK_DIR "/usr/local/VTK")
#	find_package(VTK REQUIRED)
#	include(${VTK_USE_FILE})
endif()

file(GLOB_RECURSE SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

file(GLOB_RECURSE HEADER_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

#	Fix this to allow shared libraries.
#option(SHARED_LIBRARY "when ON generates a .so library." OFF)
#if(SHARED_LIBRARY)
#	add_library(VItA SHARED ${SOURCE_FILES} ${HEADER_FILES})
#	set_target_properties(VItA PROPERTIES LINK_FLAGS "-fPIC")
#else()
#    add_library(VItA STATIC ${SOURCE_FILES} ${HEADER_FILES})
#	set_target_properties(VItA PROPERTIES LINK_FLAGS "-Wl,--whole-archive")
#endif()

add_library(VItA STATIC ${SOURCE_FILES} ${HEADER_FILES})

#set_target_properties(VItA PROPERTIES PUBLIC_HEADER ${HEADER_FILES})
target_link_libraries(VItA ${VTK_LIBRARIES})

if(DOWNLOAD_DEPENDENCIES)
	add_dependencies(VItA vtk)
endif()

install(TARGETS VItA DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN "build/*" EXCLUDE PATTERN "doc/*" EXCLUDE)
#install(FILES ${HEADER_FILES} DESTINATION include)
