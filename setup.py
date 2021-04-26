from setuptools import find_packages
import pkgconfig
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
import sys

ED_libraries=pkgconfig.parse('scifor')['libraries']
ED_library_dirs=pkgconfig.parse('scifor')['library_dirs']
ED_include_dirs=pkgconfig.parse('scifor')['include_dirs']

ED_libraries.extend(pkgconfig.parse('mpi')['libraries'])
ED_library_dirs.extend(pkgconfig.parse('mpi')['library_dirs'])
ED_include_dirs.extend(pkgconfig.parse('mpi')['include_dirs'])

ED_libraries.extend(pkgconfig.parse('slave_spins')['libraries'])
ED_library_dirs.extend(pkgconfig.parse('slave_spins')['library_dirs'])
ED_include_dirs.extend(pkgconfig.parse('slave_spins')['include_dirs'])

ext1 = Extension(
    name='ss2py',
    sources=['src/SS_INPUT_VARS.f90','src/ss2py/ss2py.f90'],
    f2py_options=["--quiet"],
    libraries=ED_libraries,
    library_dirs=ED_library_dirs,
    include_dirs=ED_include_dirs,
    extra_f90_compile_args=["-O2", "-ffree-line-length-none","-cpp","-D_MPI"])


setup(
    name = "sspy",
    version = "0.0.7",
    description = "SLAVE_SPINS python API",
    author = "Adriano Amaricci",
    author_email = "amaricci@sissa.it",
    url='https://github.com/QcmPlab/SLAVE_SPINS',
    package_dir={"": "python"},
    packages=find_packages(where="python"),
    ext_modules=[ext1])


