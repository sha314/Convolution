from distutils.core import setup, Extension
import os
import sys

# so that providing no command line argument do just fine
sys.argv = ['setup.py', 'build_ext', '--inplace']

homedir = os.environ['HOME']  # get the home_directory
# print (homedir)


source_list = [
    'python/converter.cpp',
    'python/p_extension.cpp',
    'src/main.cpp',
    'src/convolution/binomial.cpp',
    'src/convolution/convolution.cpp',
    'src/data_reader.cpp',
    'src/tests/test1.cpp',
    'src/string_methods.cpp'
]

include_directory_list = [
    'python/include',
    'src/include',
    '/usr/local/include',
    homedir + '/.local/lib/python3.6/site-packages/numpy/core/include/numpy/'
]

module1 = Extension(
    'convolution',
    include_dirs = include_directory_list,
    libraries = ['pthread'],
    sources= source_list,
    language="c++",  # so that the compiler knows about the language
    extra_compile_args=["-std=c++11", "-fopenmp"],
    extra_link_args=["-std=c++11", "-fopenmp"]
)

setup(
    name = 'convolution',
    version ='1.0',
    description='Convolution Program in C++ api',
    author='M.S. Rahman',
    url='not provided',
    ext_modules=[module1]
)
