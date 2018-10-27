cmake_minimum_required(VERSION 3.0)

project(Convolution)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS  "-pthread -fopenmp")

set(SOURCE_FILES
        src/main.cpp
        src/convolution/binomial.cpp
        src/convolution/binomial.h
        src/convolution/convolution.cpp
        src/convolution/convolution.h
        src/data_reader.cpp
        src/data_reader.h
        src/tests/test1.cpp
        src/printer.h
        src/string_methods.cpp
        src/string_methods.h src/tests/test.h src/data_writer.cpp src/data_writer.h src/utype.cpp src/utype.h)

add_executable(convolution ${SOURCE_FILES})
