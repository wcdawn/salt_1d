# Makefile

FC = gfortran
FFLAGS = -O2 -fbacktrace -fbounds-check

.SUFFIXES : .f90 .o

default : salt_1d.exe

.f90.o :
	$(FC) $(FFLAGS) -c $*.f90

salt_1d.exe : saltprops.o salt_1d.o
	$(FC) $(FFLAGS) saltprops.o salt_1d.o -o $@

clean :
	rm -f core *.o *.mod *.exe

query :
	cvs -n update -P

update : 
	cvs update -P