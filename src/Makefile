SHELL = /bin/bash

PROG = linretrace
BZ   = setupbz

INSTALL_DIR = ../bin

DIGAMMA = ./digamma/cern_digamma.a
OBJECTS = nrtype.o nrutil.o types.o estruct.o params.o mpi_org.o \
				mymath.o estruct.o root.o response.o main.o hdf5_wrapper.o input.o

BZOBJECTS = types.o estruct.o params.o setupbz.o hdf5_wrapper.o \
						lookup.o config.o

include ../make_include

# delete this later
#CC = gfortran-5.1 #mpif90 #ifort #gfortran
#CC = gfortran-5.1 -I/usr/lib64/openmpi/1.4-gcc/include -m64 -pthread -I/usr/lib64/openmpi/1.4-gcc/lib -L/usr/lib64/openmpi/1.4-gcc/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -traceback -check all
#CC = gfortran-5.1 -I/usr/lib64/openmpi/1.4-gcc/include -I/usr/lib64/openmpi/1.4-gcc/lib -L/usr/lib64/openmpi/1.4-gcc/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -lquadmath
#CFLAGS = -g -O0 #-fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall -Wextra -Wconversion #-O0 #debugging
#CFLAGS = -g -O0 -traceback -check bounds -check uninit -fpe-all=3 # ifort debug
#CFLAGS = -g -fbounds-check -fbacktrace  -Wall -Wextra -Wconversion -O0 #  gfortran debugging
#CFLAGS = -O3 #-xHost -ipo  # ifort production

LIBS= -lquadmath $(LAPACK) $(BLAS)

.PHONY: all
all: $(DIGAMMA) $(INSTALL_DIR) $(INSTALL_DIR)/$(PROG) $(INSTALL_DIR)/$(BZ)

$(DIGAMMA):
	if [ ! -f $(DIGAMMA) ]; then cd digamma; make all; fi

$(INSTALL_DIR):
	if [ ! -d $(INSTALL_DIR) ]; then mkdir -p $(INSTALL_DIR) ; fi

$(INSTALL_DIR)/$(PROG): $(OBJECTS)
	$(FC) $^ -o $@ $(FFLAGS) $(LIBS) $(HDF5) $(DIGAMMA)

$(INSTALL_DIR)/$(BZ): $(BZOBJECTS)
	$(FCDG) $^ -o $@ $(FFLAGS) $(HDF5)

%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) $(HDF5) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) $(LIBS) $(HDF5) -c $< -o $@

nrutil.o: nrtype.o
mymath.o: params.o nrtype.o nrutil.o
estruct.o: params.o types.o mymath.o
main.o: params.o mpi_org.o types.o mymath.o estruct.o response.o root.o input.o
mpi_org.o: types.o
response.o: mpi_org.o types.o params.o estruct.o
root.o: mpi_org.o params.o types.o estruct.o
input.o: params.o types.o hdf5_wrapper.o
params.o: types.o
config.o: lookup.o

setupbz.o: params.o estruct.o types.o hdf5_wrapper.o lookup.o config.o

.PHONY: clean
clean :
	rm -f *.o *.mod

.PHONY: pristine
pristine:
	rm -f *.o *.mod
	rm -f $(INSTALL_DIR)/$(PROG)
	rm -f $(INSTALL_DIR)/$(BZ)
