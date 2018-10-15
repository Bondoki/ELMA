cmake_minimum_required(VERSION 3.1)

if (NOT DEFINED LEMONADE_INCLUDE_DIR)
message("LEMONADE_INCLUDE_DIR is not provided. If build fails, use -DLEMONADE_INCLUDE_DIR=/path/to/LeMonADE/headers/ or install to default location")
endif()

if (NOT DEFINED LEMONADE_LIBRARY_DIR)
message("LEMONADE_LIBRARY_DIR is not provided. If build fails, use -DLEMONADE_LIBRARY_DIR=/path/to/LeMonADE/lib/ or install to default location")
endif()

include_directories (${LEMONADE_INCLUDE_DIR})
link_directories (${LEMONADE_LIBRARY_DIR})

SET (Eigen3_DIR "/usr/include/eigen3/")

include_directories (${Eigen3_DIR})

add_executable(AnalyzerRouseMatrix mainAnalyzerRouseMatrix.cpp)

target_link_libraries(AnalyzerRouseMatrix LeMonADE )

