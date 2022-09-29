SHELL = /bin/bash

PROG = test

include ../../make_include

LIBS= -lquadmath

.PHONY: all
all: prec hdf5

prec:
	$(FCDG) $(FFLAGS) $(LIBS) 01_test_precision.f90 -o $(PROG)
	./$(PROG)

hdf5:
	$(FCDG) $(FFLAGS) $(HDF5) 02_test_hdf5.f90 -o $(PROG)
	./$(PROG)

.PHONY: clean
clean :
	rm -f *.hdf5
	rm -f $(PROG)
