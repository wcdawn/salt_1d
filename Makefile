# Makefile

FC = gfortran
FLAGS = -Wall
OBJS = saltprops.o salt_1d.o
EXEC = salt_1d.exe

$(EXEC) : $(OBJS)
	$(FC) $(FLAGS) -o $(EXEC) $(OBJS)

%.o : %.f90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.exe *.o *.mod