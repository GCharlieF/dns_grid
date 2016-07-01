dns_grid
============
Code for direct numerical simulations of Navier-Stokes equations in triperiodic domain with pseudo-spectral discretization

Required Libraries
-------------------------------
FFTW

OPENMPI

P3DFFT

HDF5

export CC=icc

export F9X=ifort(mpifort)

export CXX=icpc(mpicpc?)

./configure --prefix=/usr/local/hdf5-1.8.15 --enable-fortran --enable-cxx (--enable-parallel)

add hdf5/bin to both PATH and LD_LIBRARY_PATH

It is recommended to compile all libraries with the same compiler The correct Include and
 Lib paths should be given directly in the makefile.
error in linking the libraries might arise from either wrong paths or libraries
installed with a different compiler.
Be sure that runtime LD_LIBRARY_PATH is properly set with the paths to all the
 libraries used.
If proper path to needed libraries is not found try to use the proper wrapper
compiler with -show flag, it will show the compilation options including the correct linking  /lib and /include paths.

Compilation errors and runtime errors related to undefined reference to subroutines are likely
due to errors in linking libraries and runtime linking.
