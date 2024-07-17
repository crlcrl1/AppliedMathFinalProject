cmake_minimum_required(VERSION 3.20)

project(part3 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

add_executable(part3 main.cpp
        util.cpp
        util.h)