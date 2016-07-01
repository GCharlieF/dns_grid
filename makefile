.PHONY: default clean cleanall

FC = mpifort

FFTW_INC = /localhome/hi214/programming/fortran/fft_lib/include
FFTW_LIB = /localhome/hi214/programming/fortran/fft_lib/lib

HDF_LIB = /localhome/hi214/programming/fortran/fortran_libraries/hdf5-1.8.17/build_16.0.1_ompi_1.10.3/lib
HDF_INC = /localhome/hi214/programming/fortran/fortran_libraries/hdf5-1.8.17/build_16.0.1_ompi_1.10.3/include

# OMPI_INC=/opt/openmpi/openmpi_1.10.3_ifort16/include/
# OMPI_LIB=/opt/openmpi/openmpi_1.10.3_ifort16/lib

# FORTRAN_LIB = -lhdf5_fortran -lhdf5 -lz -ldl -lmpi_mpifh -limf -lifcore

# LDFLAGS = -L/opt/openmpi/openmpi_1.8.4/lib -I/opt/openmpi/openmpi_1.8.4/include/
# LDFLAGS = -L/opt/openmpi/openmpi_1.10.3_ifort16/lib -I/opt/openmpi/openmpi_1.10.3_ifort16/include/

PACKAGE = /user/hi214/programming/fortran/fortran_libraries/p3dfft-master/build/libp3dfft.a
P3DFFT_INC= /user/hi214/programming/fortran/fortran_libraries/p3dfft-master/include
P3DFFT_LIB = /user/hi214/programming/fortran/fortran_libraries/p3dfft-master/lib

INCLUDE=$(FFT_INC) $(P3DFFT_INC) #$(OMPI_INC)
LIBS=$(FFT_LIB) $(P3DFFT_LIB) #$(OMPI_LIB)

LDFLAGS =  -L$(LIBS) -I$(INCLUDE) -L$(HDF_LIB) -I$(HDF_INC) -lhdf5_fortran -lhdf5 -lz -ldl

# OP_COMP = -g -u -fpp -traceback -check bounds -convert big_endian -I$(INCLUDE) -L$(LIBS) $(LDFLAGS) $(FORTRAN_LIB) -lfftw3
 # OP_COMP =  -u -O3 -fpp -convert big_endian -I$(INCLUDE) -L$(LIBS) $(LDFLAGS) $(FORTRAN_LIB) -lfftw3 -lhdf5_fortran -lhdf5 -lz
 OP_COMP =  -u -O3 -fpp -convert big_endian -lfftw3 -lmpi_mpifh -limf -lifcore

DEP = parameters_mod.o variables_mod.o mpi_mod.o grid_forcing_mod.o stats_and_probes_mod.o hit_forcings_mod.o h5util_mod.o IO_mod.o time_advancement_mod.o dns_grid.o $(PACKAGE)


default: dns_grid

dns_grid: $(DEP)
	$(FC) $(OP_COMP) -o dns_grid $(DEP) $(LDFLAGS)

%.o: %.f90
	$(FC) $(OP_COMP) -c $< $(LDFLAGS)

clean:
	-rm  -f *.o

cleanall:
	-rm  -f dns_grid *.o *.mod
