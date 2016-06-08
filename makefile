FC = mpifort


FFTW_INC = /localhome/hi214/programming/fortran/fft_lib/include
FFTW_LIB = /localhome/hi214/programming/fortran/fft_lib/lib

FORTRAN_LIB = -lmpi_mpifh -limf -lifcore

LDFLAGS = -L/opt/openmpi/openmpi_1.8.4/lib -I/opt/openmpi/openmpi_1.8.4/include/

PACKAGE = /user/hi214/programming/fortran/fortran_libraries/p3dfft-master/build/libp3dfft.a
P3DFFT_INC= /user/hi214/programming/fortran/fortran_libraries/p3dfft-master/include
P3DFFT_LIB = /user/hi214/programming/fortran/fortran_libraries/p3dfft-master/lib

INCLUDE=$(FFT_INC) $(P3DFFT_INC)
LIBS=$(FFT_LIB) $(P3DFFT_LIB)

OP_COMP = -g -u -fpp -traceback -check bounds -convert big_endian -I$(INCLUDE) -L$(LIBS) $(LDFLAGS) $(FORTRAN_LIB) -lfftw3
 # OP_COMP =  -u -O3 -fpp -convert big_endian -I$(INCLUDE) -L$(LIBS) $(LDFLAGS) $(FORTRAN_LIB) -lfftw3

DEP = parameters_mod.o variables_mod.o mpi_mod.o fft_mod.o grid_forcing_mod.o stats_and_probes_mod.o hit_forcings_mod.o IO_mod.o time_advancement_mod.o dns_grid.o $(PACKAGE)
dns_grid: $(DEP)
	$(FC) $(OP_COMP) -o dns_grid $(DEP)



#Use this if you  want to compile all the files in the folder
#%: %.o
#	$(FC) $(OP_COMP) -o $@ $^

%.o: %.f90
	$(FC) $(OP_COMP) -c $<
clean:
	-rm  -f *.o
allclean:
	-rm  -f dns_grid *.o *.mod
