cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -DDEBUG -DDEBUG_FILE -DPRINT_TO_TERMINAL -DSNIFFBOT_CONFIG -DPSO_CONFIG")
set(CMAKE_BUILD_TYPE "Debug")

# Configure Paths
set(INCLUDE_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build")

# Configure sources
aux_source_directory(${SOURCE_DIR} SRC_LIST)

# Build output configuration

# set(executables_list main_hls_model main_hls_system)
set(executables_list main_hls_model main_hls_system main_hls_pso)
foreach(exec ${executables_list})
    project(${exec} VERSION 0.1.0)
    add_executable(${PROJECT_NAME}  ${SOURCE_DIR}/${exec}.cpp)
    target_include_directories(${PROJECT_NAME} PUBLIC
        ${PROJECT_BINARY_DIR})
    target_include_directories(${PROJECT_NAME} PRIVATE
        ${INCLUDE_FOLDER})
    target_link_libraries(${PROJECT_NAME} PRIVATE 
        m)
endforeach()


# ## Configure sources PSO
# set(SOURCE_LIST hls_nonlinear_solver.cpp)

# ## Build for PSO
# set(exec main_hls_pso)
# project(${exec} VERSION 0.1.0)
# add_executable(${PROJECT_NAME} ${SOURCE_DIR}/${exec}.cpp
#     ${SOURCE_DIR}/${SOURCE_LIST}
# )
# target_include_directories(${PROJECT_NAME} PUBLIC
#     ${PROJECT_BINARY_DIR})
# target_include_directories(${PROJECT_NAME} PRIVATE
#     ${INCLUDE_FOLDER})
# target_link_libraries(${PROJECT_NAME} PRIVATE 
#     m)
