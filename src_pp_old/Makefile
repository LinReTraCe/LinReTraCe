SHELL = /bin/bash

PROG = setupdata

INSTALL_DIR = ../bin

OBJECTS = types.o estruct.o params.o mpi_org.o \
				main.o hdf5_wrapper.o aux.o lookup.o config.o

include ../make_include

LIBS= -lquadmath $(LAPACK) $(BLAS)

.PHONY: all
all: $(DIGAMMA) $(INSTALL_DIR) $(INSTALL_DIR)/$(PROG) $(INSTALL_DIR)/$(BZ)

$(INSTALL_DIR):
	if [ ! -d $(INSTALL_DIR) ]; then mkdir -p $(INSTALL_DIR) ; fi

$(INSTALL_DIR)/$(PROG): $(OBJECTS)
	$(FC) $^ -o $@ $(FFLAGS) $(LIBS) $(HDF5) $(DIGAMMA)

$(INSTALL_DIR)/$(BZ): $(BZOBJECTS)
	$(FCDG) $^ -o $@ $(FFLAGS) $(HDF5)

%.o: %.f90
	$(FC) $(FFLAGS) $(LIBS) $(HDF5) -c $< -o $@

estruct.o: params.o types.o
main.o: params.o estruct.o types.o hdf5_wrapper.o lookup.o config.o aux.o
config.o: lookup.o

.PHONY: clean
clean :
	rm -f *.o *.mod

.PHONY: pristine
pristine:
	rm -f *.o *.mod
	rm -f $(INSTALL_DIR)/$(PROG)
