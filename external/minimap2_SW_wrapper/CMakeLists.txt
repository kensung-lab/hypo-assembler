cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(hypo LANGUAGES C CXX VERSION 1.0.0)

include(GNUInstallDirs)

#set( CMAKE_VERBOSE_MAKEFILE on )

if(CMAKE_BUILD_TYPE MATCHES Debug)
    message("Debug mode")
endif()

option(optimise_for_native "Build with -march=native" OFF)

if(optimise_for_native)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3 -march=native")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -O3")
endif()



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)  # Compiler-specific features are not used.

# Output directories for binaries and libraries.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Source files
set(SOURCES 
    ${PROJECT_SOURCE_DIR}/include/kalloc.h
    ${PROJECT_SOURCE_DIR}/src/kalloc.c
    ${PROJECT_SOURCE_DIR}/include/ksw2.h
    ${PROJECT_SOURCE_DIR}/src/ksw2_dispatch.c
    ${PROJECT_SOURCE_DIR}/src/ksw2_extd2_sse.c
    ${PROJECT_SOURCE_DIR}/src/ksw2_exts2_sse.c
    ${PROJECT_SOURCE_DIR}/src/ksw2_extz2_sse.c
    ${PROJECT_SOURCE_DIR}/src/ksw2_ll_sse.c
    ${PROJECT_SOURCE_DIR}/include/wrapper.hpp
    ${PROJECT_SOURCE_DIR}/src/wrapper.cpp
)

# Enable .cpp files to add headers in the include folder.
include_directories (
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

add_library(kalloc include/kalloc.h src/kalloc.c)
add_library(ksw2_extd2_sse41 src/ksw2_extd2_sse.c include/ksw2.h include/kalloc.h src/ksw2_ll_sse.c)
add_library(ksw2_extd2_sse2 src/ksw2_extd2_sse.c include/ksw2.h include/kalloc.h src/ksw2_ll_sse.c)

# Build binaries.
add_library(sw_wrapper ${SOURCES})
target_link_libraries(sw_wrapper kalloc ksw2_extd2_sse41 ksw2_extd2_sse2)
target_include_directories(sw_wrapper PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
