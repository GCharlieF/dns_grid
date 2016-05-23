# FC = ifort
FC = mpifort
FFT_INC = /localhome/hi214/programming/fortran/fft_lib/include
FFT_LIB = /localhome/hi214/programming/fortran/fft_lib/lib
# OP_COMP = -g -u -traceback -check bounds -convert big_endian -lfftw3 #-lfftw3_threads -traceback
 OP_COMP =  -u -O3 -convert big_endian -I$(FFT_INC) -L$(FFT_LIB)  -lfftw3  -lfftw3_mpi #-lfftw3_threads
#Loading order of modules matters!1

DEP = parameters_mod.o variables_mod.o mpi_mod.o fft_mod.o grid_forcing_mod.o stats_and_probes_mod.o hit_forcings_mod.o IO_mod.o time_advancement_mod.o dns_grid.o
dns_grid: $(DEP)
	$(FC) $(OP_COMP) -o dns_grid $(DEP) -lfftw3 -lfftw3_mpi



#Use this if you  want to compile all the files in the folder
#%: %.o
#	$(FC) $(OP_COMP) -o $@ $^

%.o: %.f90
	$(FC) $(OP_COMP) -c $< -lfftw3 -lfftw3_mpi
clean:
	-rm  -f *.o
allclean:
	-rm  -f dns_grid *.o *.mod
