cmake_minimum_required(VERSION 3.0.0)
project(nmpc_pso VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ")

# Configure Paths
set(INCLUDE_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build")

# Configure sources
aux_source_directory(${SOURCE_DIR} SRC_LIST)

## Build output configuration
add_executable(nmpc_pso ${SOURCE_DIR}/main.cpp
    ${SRC_LIST}
    ${SOURCE_DIR}/models/inverted_pendulum.cpp
    ${SOURCE_DIR}/models/vant.cpp
    ${SOURCE_DIR}/models/sniffbot.cpp
    )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${PROJECT_BINARY_DIR})

target_include_directories(${PROJECT_NAME} PRIVATE
    ${INCLUDE_FOLDER})

target_link_libraries(${PROJECT_NAME} PRIVATE 
    m)