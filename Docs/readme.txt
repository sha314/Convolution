A program to perform convolution on particular columns of a data file
Then writes the original data and the convolved data to a new file with a header

# using threshold
using sufficiently small threshold (1e-9 or 1e-15)keeps the normalization constant unchanged.
convolved output seems to be unchanged also.


# For Mac OS on M1. tested on M1 Mac Pro
$ brew install llvm
$ brew install libomp

Added these two lines to CmakeLists.txt file after checking the OS. Don’t need anything else.
set(CMAKE_C_COMPILER /opt/homebrew/opt/llvm/bin/clang)
set(CMAKE_CXX_COMPILER /opt/homebrew/opt/llvm/bin/clang++)