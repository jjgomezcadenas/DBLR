from distutils.core import setup
from Cython.Build import cythonize
import os

setup(ext_modules=cythonize('*.pyx'))
os.system("rm -f *.c")

