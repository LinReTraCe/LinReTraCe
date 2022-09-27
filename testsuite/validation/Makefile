SHELL = /bin/bash

PROG = test

include ../../make_include

LIBS= -lquadmath

.PHONY: all
all: hdf5 prec

hdf5:
	$(FCDG) $(FFLAGS) $(HDF5) 01_test_hdf5.f90 -o $(PROG)
	./$(PROG)

prec:
	$(FCDG) $(FFLAGS) $(LIBS) 02_test_precision.f90 -o $(PROG)
	./$(PROG)

.PHONY: clean
clean :
	rm -f *.hdf5
	rm -f $(PROG)
