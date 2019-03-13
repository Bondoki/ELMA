cmake_minimum_required (VERSION 2.8)

# cmake .. -DLEMONADE_INCLUDE_DIR=~/LeMonADE-lib-v2.0.0/include/ -DLEMONADE_LIBRARY_DIR=~/LeMonADE-lib-v2.0.0/lib/

#SET (LEMONADE_INCLUDE_DIR "~/LeMonADE-lib-v2.0.0/include/")
#SET (LEMONADE_LIBRARY_DIR "~/LeMonADE-lib-v2.0.0/lib/")
# SET (EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

SET (PROJECTS_DIR "./")

SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECTS_DIR}/bin")

include_directories("./updater" "./analyzer" "./feature" "./utility" "./io" "./include/extern" )

if (NOT DEFINED LEMONADE_INCLUDE_DIR)
message("LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location")
endif()

MESSAGE("LEMONADE_INCLUDE_DIR is ${LEMONADE_INCLUDE_DIR}")

if (NOT DEFINED LEMONADE_LIBRARY_DIR)
message("LEMONADE_LIBRARY_DIR is not provided. If build fails, use -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE/lib/ or install to default location")
endif()

MESSAGE("LEMONADE_LIBRARY_DIR is ${LEMONADE_LIBRARY_DIR}")


include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

#define possible flags
SET (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -msse2 -mssse3 -fexpensive-optimizations -std=c++11")
SET (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -msse2 -mssse3 -fexpensive-optimizations -std=c++11")

SET (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -DDEBUG ")
SET (CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -Wall -DDEBUG ")
# SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -O2 -std=c++11")
# SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEBUG} -O2 ")



find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#define value of CMAKE_BUILD_TYPE depending on input
IF(NOT CMAKE_BUILD_TYPE)
SET (CMAKE_BUILD_TYPE "Release") #default build type is Release
ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")
SET (CMAKE_BUILD_TYPE "Release") 
ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
SET (CMAKE_BUILD_TYPE "Debug") 
ELSE(NOT CMAKE_BUILD_TYPE)
MESSAGE(FATAL_ERROR "Invalid build type ${CMAKE_BUILD_TYPE} specified.")
ENDIF(NOT CMAKE_BUILD_TYPE)

#output depending on build type
IF(CMAKE_BUILD_TYPE STREQUAL "Release")
SET (CMAKE_VERBOSE_MAKEFILE 0)
MESSAGE("Build type is ${CMAKE_BUILD_TYPE}")
MESSAGE("USING CXX COMPILER FLAGS ${CMAKE_CXX_FLAGS_RELEASE}")
MESSAGE("USING C COMPILER FLAGS ${CMAKE_C_FLAGS_RELEASE}")
ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Debug")
SET (CMAKE_VERBOSE_MAKEFILE 1)
MESSAGE("Build type is ${CMAKE_BUILD_TYPE}")
MESSAGE("USING CXX COMPILER FLAGS ${CMAKE_CXX_FLAGS_DEBUG}")
MESSAGE("USING C COMPILER FLAGS ${CMAKE_C_FLAGS_DEBUG}")
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Release")

## ###############  System Creators ############# ##
add_subdirectory(projects)

