cmake_minimum_required(VERSION 3.0.0)
project(nmpc_pso VERSION 0.1.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DSNIFFBOT_CONFIG")

# Configure Paths
set(INCLUDE_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PROJECT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build")

# Configure sources
aux_source_directory(${SOURCE_DIR} SRC_LIST)

## Build output configuration
#    ${SOURCE_DIR}/models/inverted_pendulum.cpp
#    ${SOURCE_DIR}/models/sniffbot.cpp
add_executable(nmpc_pso ${SOURCE_DIR}/test/main.cpp
    ${SRC_LIST}
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${PROJECT_BINARY_DIR})
target_include_directories(${PROJECT_NAME} PRIVATE
    ${INCLUDE_FOLDER})
target_link_libraries(${PROJECT_NAME} PRIVATE 
    m)

#set(HLS_SOURCES ${SOURCE_DIR}/pso.cpp ${SOURCE_DIR}/inverted_pendulum.cpp ${SOURCE_DIR}/system.cpp ${SOURCE_DIR}/aux_functions.cpp)
#add_executable(nmpc_pso_hls ${SOURCE_DIR}/test/main_hls.cpp
#    ${HLS_SOURCES}
#)
#target_include_directories(nmpc_pso_hls PUBLIC
#    ${PROJECT_BINARY_DIR})
#target_include_directories(nmpc_pso_hls PRIVATE
#    ${INCLUDE_FOLDER})
#target_link_libraries(nmpc_pso_hls PRIVATE 
#    m)