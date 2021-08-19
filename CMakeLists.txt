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
else()
	find_path(VTK_INCLUDE_DIRS vtkAbstractArray.h
          HINTS .. ENV VITA_PATH /usr/local/include/vtk-8.1
          PATH_SUFFIXES vita_build/include/vtk-8.1)
	find_path(VTK_LIBRARY_DIRS libvtkCommonCore-8.1.so
          HINTS .. ENV VITA_PATH /usr/local/lib
          PATH_SUFFIXES vita_build/lib)
	link_directories(${VTK_LIBRARY_DIRS})
	include_directories(${VTK_INCLUDE_DIRS})
	set(VTK_LIBRARIES
			vtkCommonCore-8.1
			vtkCommonDataModel-8.1
			vtkCommonExecutionModel-8.1
			vtkFiltersFlowPaths-8.1
			vtkFiltersModeling-8.1
			vtkIOCore-8.1
			vtkIOLegacy-8.1
			vtkIOXML-8.1
			vtkIOGeometry-8.1
		)
endif()

file(GLOB_RECURSE SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

file(GLOB_RECURSE HEADER_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

#	Fix this to allow shared libraries.
option(SHARED_LIBRARY "when ON generates a .so library." OFF)
if(SHARED_LIBRARY)
	add_library(VItA SHARED ${SOURCE_FILES} ${HEADER_FILES})
	set_target_properties(VItA PROPERTIES POSITION_INDEPENDENT_CODE True)
else()
    add_library(VItA STATIC ${SOURCE_FILES} ${HEADER_FILES})
endif()

target_link_libraries(VItA ${VTK_LIBRARIES})

if(DOWNLOAD_DEPENDENCIES)
	add_dependencies(VItA vtk)
endif()

install(TARGETS VItA DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN "build/*" EXCLUDE PATTERN "doc/*" EXCLUDE)
#install(FILES ${HEADER_FILES} DESTINATION include)
