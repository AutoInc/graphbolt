cmake_minimum_required(VERSION 3.16)
project(untitled)

set(CMAKE_CXX_STANDARD 14)


add_executable(SSSP ${PROJECT_SOURCE_DIR}/apps/SSSP.C)
target_compile_definitions(SSSP PUBLIC -DEDGEDATA)

add_executable(PHP ${PROJECT_SOURCE_DIR}/apps/PHP.C)
target_compile_definitions(PHP PUBLIC -DEDGEDATA)

add_executable(PANGERANK ${PROJECT_SOURCE_DIR}/apps/PageRank.C)
