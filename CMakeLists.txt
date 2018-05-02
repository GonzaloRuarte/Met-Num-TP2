cmake_minimum_required(VERSION 3.9)
project(PageRank)

set(CMAKE_CXX_STANDARD 11)

add_executable(tp2 main.cpp util.cpp util.h
        rdtsc.h ppmloader.cpp ppmloader.h)