SHELL = /bin/bash

PROG = test

include ../../make_include

LIBS= -lquadmath

.PHONY: all
all: hdf5 quad

hdf5:
	$(FC) $(FFLAGS) $(HDF5) 01_test_hdf5.f90 -o $(PROG)
	./$(PROG)

quad:
	$(FC) $(FFLAGS) $(LIBS) 02_test_lquad.f90 -o $(PROG)
	./$(PROG)

.PHONY: clean
clean :
	rm -f *.hdf5
	rm -f $(PROG)
