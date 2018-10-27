cmake_minimum_required(VERSION 3.0)

project(Convolution)

set(CMAKE_CXX_STANDARD 17)

SET(CMAKE_CXX_FLAGS  "-pthread -fopenmp")

set(SOURCE_FILES
        src/main.cpp
        src/convolution/binomial.cpp
        src/include/binomial.h
        src/convolution/convolution.cpp
        src/include/convolution.h
        src/data_reader.cpp
        src/include/data_reader.h
        src/tests/test1.cpp
        src/include/printer.h
        src/string_methods.cpp
        src/include/string_methods.h
        src/include/data_writer.h
        src/include/test.h
        src/include/utype.h
        src/data_writer.cpp
        src/utype.cpp)

add_executable(convolution ${SOURCE_FILES})
