cmake_minimum_required(VERSION 2.6)

project(VItA_lib)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include(dependencies.cmake)
link_directories(${VTK_LIB_DIR})
include_directories(${VTK_INCLUDE_DIR})

file(GLOB_RECURSE SOURCE_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
)

file(GLOB_RECURSE HEADER_FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
)

# Fix this to allow shared libraries. 
option(SHARED_LIBRARY "when ON generates a .so library." OFF)
if(SHARED_LIBRARY)
  add_library(VItA SHARED ${SOURCE_FILES} ${HEADER_FILES})
  set_target_properties(VItA PROPERTIES LINK_FLAGS "-fPIC")
else()
   add_library(VItA STATIC ${SOURCE_FILES} ${HEADER_FILES})
   set_target_properties(VItA PROPERTIES LINK_FLAGS "-Wl,--whole-archive")
endif()

option(DEBUG "Enables debugging" OFF)
set(CXX_GDB "-g")
if(DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_GDB}")
endif()

option(OPENMP "Enables OpenMP parallelism" OFF)
set(CXX_OMP "-fopenmp")
if(OPENMP)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_OMP}")
endif()

option(GPROF "Enables gprof single-threaded profiling" OFF)
set(CXX_PROFILER "-pg")
if(GPROF)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_PROFILER}")
endif()

#set_target_properties(VItA PROPERTIES PUBLIC_HEADER ${HEADER_FILES})
target_link_libraries(VItA ${VTK_LIBRARIES})

install(TARGETS VItA DESTINATION lib)
install(DIRECTORY ${VItA} DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN "build/*" EXCLUDE PATTERN "doc/*" EXCLUDE)
#install(FILES ${HEADER_FILES} DESTINATION include)
