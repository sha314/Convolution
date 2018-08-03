cmake_minimum_required(VERSION 3.0)

project(Convolution)

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_CXX_FLAGS  "-pthread -fopenmp")

set(SOURCE_FILES
        main.cpp
        src/binomial.cpp
        src/binomial.h
        src/convolution.cpp
        src/convolution.h
        src/data_reader.cpp
        src/data_reader.h
        )

add_executable(convolution ${SOURCE_FILES})
