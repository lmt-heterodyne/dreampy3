#from distutils.core import setup
from setuptools import setup, find_packages, Extension
import os, glob

__version__ = '0.1'

# punwrap2_src = glob.glob("dreampy/holography/src/punwrap/*.c")
# macros = []

# try:
#     import numpy
#     numpy_include = numpy.get_include()
#     from numpy.distutils.system_info import get_info
# except ImportError:
#     raise ImportError("""
#     The numpy package is needed to compile the holpy tools extensions but
#     is not installed. Quitting the setup.
#     """)

# holography_punwrap_extension = Extension('dreampy.holography.punwrap._punwrap2D',
#                                          punwrap2_src,
#                                          include_dirs=[numpy_include,],
#                                          libraries=['m'],
#                                          define_macros=macros)

# swig_src = glob.glob("dreampy/holography/python_cxx/*.i")
# swig_src.extend(glob.glob("dreampy/holography/python_cxx/*.cxx"))

# print swig_src
# holography_swig_extension = Extension('dreampy.holography.python_cxx._holofile',
#                                       swig_src,
#                                       swig_opts=['-c++', '-shadow'],
#                                       include_dirs = ['dreampy/holography/python_cxx'],
#                                       )

setup(
    name = "dreampy3",
    version = __version__,
    description = "DREAMPY: Data REduction and Analysis Methods in PYthon",
    author = "Gopal Narayanan",
    author_email = "gopal@astro.umass.edu",
    packages = find_packages(),
    #setup_requires=['nose>=0.11', 'sphinx'],
    #setup_requires=['nose', 'sphinx'],
    # scripts = ['bin/dreampy', 'bin/holocal', 'bin/holoxform', 'bin/holoxform2', 
    #            'bin/holophasefit', 'bin/holophasefit2', 'bin/holophasefit3',
    #            'bin/holopanelfit', 'bin/holodiff', 'bin/tle2catalog',
    #            'bin/holocontrolreport', 'bin/holoactuator', 'bin/hologetreporte',
    #            'bin/holoenableactuators', 'bin/holoavg',
    #            'bin/holoreportdiagnose', 'bin/holorestoreactuators',
    #            'bin/rsrmakecorrcal', 'bin/makepointingfile',
    #            #'bin/holoxformsp', 'bin/holophasefit2sp', 'bin/holopanelfitsp'
    #            ],
    # data_files = [('/usr/share/dreampy', ['data/ipythonrc-dreampy',
    #                                       'data/ipy_user_conf.py'])],
    # ext_modules = [holography_punwrap_extension,
    #                holography_swig_extension,
    #                ],

    )
