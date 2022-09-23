SHELL = /bin/bash

PROG = tests

include ../../make_include

LIBS= -lquadmath

.PHONY: all
all:
	$(FC) $(FFLAGS) $(LIBS) $(HDF5) tests.f90 -o $(PROG)
	./$(PROG)

.PHONY: clean
clean :
	rm -f *.hdf5
	rm -f $(PROG)
