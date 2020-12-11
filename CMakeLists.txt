cmake_minimum_required(VERSION 3.16)
project(untitled)

set(CMAKE_CXX_STANDARD 14)

add_definitions(-DEDGEDATA)

add_executable(PHP ${PROJECT_SOURCE_DIR}/apps/PHP.C)