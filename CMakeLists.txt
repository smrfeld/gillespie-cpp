cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(gilsp)

include_directories(/usr/local/include)

AUX_SOURCE_DIRECTORY(src SRC_FILES) 

file (GLOB HEADER_FILES "include/gilsp_bits/*.hpp")
file (GLOB HEADER_FILES_FWDS "include/gilsp_bits/fwds/*.hpp")

add_library(gilsp SHARED 
	src/helpers.cpp
	src/counts.cpp
	src/counts_hist.cpp
	src/gillespie.cpp
	src/rxn.cpp
	include/gilsp_bits/helpers.hpp
	include/gilsp_bits/counts.hpp
	include/gilsp_bits/counts_hist.hpp
	include/gilsp_bits/gillespie.hpp
	include/gilsp_bits/rxn.hpp
	${HEADER_FILES_FWDS}
)

target_include_directories(gilsp PRIVATE include/gilsp_bits include/gilsp_bits/fwds)

install(TARGETS gilsp DESTINATION lib)

install(FILES 
	include/gilsp 
	DESTINATION include)

install(DIRECTORY 
	include/gilsp_bits 
	DESTINATION include)