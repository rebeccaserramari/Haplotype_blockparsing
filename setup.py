from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext


ext_modules = [Extension("blockparser", ["blockparser.pyx", "DP_matrix.cpp"], language='c++',extra_compile_args=["-std=c++11", "-O2"], extra_link_args=["-std=c++11"])]

setup(cmdclass = {'build_ext': build_ext}, ext_modules = ext_modules)