SHELL = /bin/bash

PROG = linretrace

INSTALL_DIR = ../bin

DIGAMMA = ./digamma/cern_digamma.a
OBJECTS = types.o params.o mpi_org.o psi_fast.o \
				main.o hdf5_wrapper.o input.o output.o fermi.o \
				auxiliary.o lookup.o config.o response.o root.o

include ../make_include

LIBS= -lquadmath

.PHONY: all
all: $(DIGAMMA) $(INSTALL_DIR) $(INSTALL_DIR)/$(PROG) $(INSTALL_DIR)/$(BZ)

$(DIGAMMA):
	if [ ! -f $(DIGAMMA) ]; then cd digamma; make all; fi

$(INSTALL_DIR):
	if [ ! -d $(INSTALL_DIR) ]; then mkdir -p $(INSTALL_DIR) ; fi

$(INSTALL_DIR)/$(PROG): $(OBJECTS)
	$(FC) $^ -o $@ $(FFLAGS) $(LIBS) $(HDF5) $(DIGAMMA)

%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) $(HDF5) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) $(LIBS) $(HDF5) -c $< -o $@

# psi_fast is a Fortran MODULE, so it must be compiled by $(FC) -- the same
# compiler as everything else -- and its .mod must land in this directory.
# It therefore does NOT go into digamma/cern_digamma.a, which is built by
# $(FCDG) and may use a different compiler; .mod files are not portable
# between compilers.  The source lives next to the routines it replaces.
psi_fast.o: digamma/psi_fast.F90
	$(FC) $(FFLAGS) $(FPPFLAGS) $(LIBS) -c $< -o $@

main.o: params.o mpi_org.o types.o input.o output.o auxiliary.o config.o response.o root.o hdf5_wrapper.o
auxiliary.o: mpi_org.o
response.o: mpi_org.o types.o params.o fermi.o psi_fast.o
root.o: mpi_org.o params.o types.o fermi.o auxiliary.o psi_fast.o
input.o: params.o types.o hdf5_wrapper.o auxiliary.o
output.o: mpi_org.o types.o params.o hdf5_wrapper.o
config.o: params.o lookup.o auxiliary.o types.o
fermi.o: params.o psi_fast.o

.PHONY: clean
clean :
	cd digamma; make clean
	rm -f *.o *.mod

.PHONY: pristine
pristine:
	cd digamma; make pristine
	rm -f *.o *.mod
	rm -f $(INSTALL_DIR)/$(PROG)
