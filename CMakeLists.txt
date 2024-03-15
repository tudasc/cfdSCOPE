cmake_minimum_required(VERSION 3.16)

project(
    MiniCFD
    VERSION 0.0.1
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(minicfd
    src/MiniCFD.cpp
)