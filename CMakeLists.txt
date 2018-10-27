cmake_minimum_required(VERSION 3.0)

project(Convolution)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS  "-pthread -fopenmp")

set(SOURCE_FILES
        src/main.cpp
        src/convolution/binomial.cpp
        src/include/binomial.h
        src/convolution/convolution.cpp
        src/include/convolution.h
        src/data_reader.cpp
        src/include/data_reader.h
        src/data_writer.cpp
        src/include/data_writer.h
        src/tests/test1.cpp
        src/include/printer.h
        src/string_methods.cpp
        src/include/string_methods.h)

add_executable(convolution ${SOURCE_FILES})
