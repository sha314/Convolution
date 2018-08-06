from distutils.core import setup, Extension

module1 = Extension(
    'convolution',
    include_dirs = ['/usr/local/include'],
    libraries = ['pthread'],
    sources=['src/python/p_extension.cpp']
)

setup(
    name = 'convolution',
    version ='1.0',
    description='Python C api example',
    author='Shahnoor',
    url='not provided',
    ext_modules=[module1]
)