SHELL = /bin/bash

PROG = linretrace

INSTALL_DIR = ../bin
OBJECTS = types.o estruct.o params.o mpi_org.o \
					estruct.o root.o rootQ.o response.o main.o

include ../make_include

# delete this later
#CC = gfortran-5.1 #mpif90 #ifort #gfortran
#CC = gfortran-5.1 -I/usr/lib64/openmpi/1.4-gcc/include -m64 -pthread -I/usr/lib64/openmpi/1.4-gcc/lib -L/usr/lib64/openmpi/1.4-gcc/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -traceback -check all
#CC = gfortran-5.1 -I/usr/lib64/openmpi/1.4-gcc/include -I/usr/lib64/openmpi/1.4-gcc/lib -L/usr/lib64/openmpi/1.4-gcc/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -lquadmath
#CFLAGS = -g -O0 #-fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall -Wextra -Wconversion #-O0 #debugging
#CFLAGS = -g -O0 -traceback -check bounds -check uninit -fpe-all=3 # ifort debug
#CFLAGS = -g -fbounds-check -fbacktrace  -Wall -Wextra -Wconversion -O0 #  gfortran debugging
#CFLAGS = -O3 #-xHost -ipo  # ifort production

LIBS= ./digamma/cern_digamma.a -lquadmath $(LAPACK) $(BLAS)

.PHONY: all
all: $(INSTALL_DIR) $(INSTALL_DIR)/$(PROG)

$(INSTALL_DIR):
	if [ ! -d $(INSTALL_DIR) ]; then mkdir -p $(INSTALL_DIR) ; fi

$(INSTALL_DIR)/$(PROG): $(OBJECTS)
	$(FC) $^ -o $@ $(FFLAGS) $(LIBS)

$(PROG): $(OBJ)
	$(FC) -o $(PROG) $(FFLAGS) $(OBJ) $(LDFLAGS) $(LIBS)
	mv main linretrace

%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) $(LIBS) -c $< -o $@

estruct.o: params.o types.o
main.o: params.o mpi_org.o types.o estruct.o response.o root.o rootQ.o
mpi_org.o: types.o
response.o: mpi_org.o types.o params.o estruct.o
responseBQ.o: response.o
responseQ.o: response.o
root.o: params.o types.o
rootQ.o: params.o types.o

.PHONY: clean
clean :
	rm -f *.o *.mod

.PHONY: pristine
pristine:
	rm -f *.o *.mod
	rm -f $(INSTALL_DIR)/$(PROG)
