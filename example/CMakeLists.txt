cmake_minimum_required(VERSION 3.1)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(example)

include_directories(/usr/local/include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR}/../bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR}/../bin)

find_library(GILSP_LIB gilsp HINTS /usr/local/lib/)

############################################################
# main_gillespie
############################################################

add_executable(main_gillespie src/main_gillespie.cpp)

target_link_libraries(main_gillespie PUBLIC ${GILSP_LIB})

############################################################
# main_tau_leaping
############################################################

add_executable(main_tau_leaping src/main_tau_leaping.cpp)

target_link_libraries(main_tau_leaping PUBLIC ${GILSP_LIB})

############################################################
# extra
############################################################

# If want to include install target
# install(TARGETS bmla_layer_1 RUNTIME DESTINATION bin)