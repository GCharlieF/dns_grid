FC = ifort
#OP_COMP = -g -u -check bounds -convert big_endian -lfftw3 -lfftw3_threads -traceback
OP_COMP =  -u -O3 -convert big_endian -lfftw3 -lfftw3_threads
#Loading order of modules matters!

DEP = parameters_mod.o variables_and_IO_mod.o fft_mod.o IO_mod.o time_advancement_mod.o dns_grid.o
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
