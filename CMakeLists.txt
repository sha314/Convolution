cmake_minimum_required(VERSION 3.0)

project(Convolution)

set(CMAKE_CXX_STANDARD 14)

SET(CMAKE_CXX_FLAGS  "-pthread -fopenmp") # with Open MP
#SET(CMAKE_CXX_FLAGS  "-pthread -fopenacc") # with Open ACC

set(SOURCE_FILES
        src/main.cpp
        src/convolution/binomial.cpp
        src/convolution/binomial.h
        src/convolution/convolution.cpp
        src/convolution/convolution.h
        src/io/data_reader.cpp
        src/io/data_reader.h
        src/tests/test1.cpp
        src/include/printer.h
        src/string_methods.cpp
        src/include/string_methods.h
        src/tests/test1.h
        src/io/data_writer.cpp
        src/io/data_writer.h
        src/utype.cpp
        src/include/utype.h
        src/tests/test2.cpp
        src/tests/test2.h)

add_executable(convolution ${SOURCE_FILES})
