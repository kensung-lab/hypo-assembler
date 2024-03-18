# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)
set (CMAKE_CXX_STANDARD 11)

# Name of project and that it is C++ only.
project (overlap_finder CXX)

if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug mode")
endif()

option(optimise_for_native "Build with -march=native" OFF)

if(optimise_for_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3 -march=native")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3")
endif()

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

#find_package (ZLIB)

# include the sdsl lite installed in libs/sdsl-lite/install and htslib installed in libs/htslib-1.9/install
include_directories(SYSTEM libs/sdsl-lite/install/include
                    libs/htslib-1.9/include)
link_directories(libs/sdsl-lite/install/lib
                libs/htslib-1.9/lib)

if (NOT TARGET sdsl)
    add_subdirectory(libs/sdsl-lite EXCLUDE_FROM_ALL)    
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

# Add CXX flags
set (CMAKE_CXX_FLAGS "-pthread ${CMAKE_CXX_FLAGS} -O3")

add_executable(find_overlap overlap.cpp)
target_link_libraries(find_overlap hts z sdsl OpenMP::OpenMP_CXX)

# add_executable(find_overlap_txt overlap_txt.cpp)
# target_link_libraries(find_overlap_txt hts z OpenMP::OpenMP_CXX)

# add_executable(find_overlap_bv overlap_bitvector.cpp)
# target_link_libraries(find_overlap_bv z sdsl hts OpenMP::OpenMP_CXX)
