SHELL = /bin/bash

PROG = test

include ../../make_include

LIBS= -lquadmath

.PHONY: all
all: prec hdf5

prec:
	$(FCDG) $(FFLAGS) $(LIBS) -c 01_test_precision.f90 -o 01_test_precision.o
	$(FCDG) 01_test_precision.o -o $(PROG) $(FFLAGS)
	./$(PROG)

hdf5:
	$(FCDG) $(FFLAGS) $(HDF5) -c 02_test_hdf5.f90 -o 02_test_hdf5.o
	$(FCDG) 02_test_hdf5.o -o $(PROG) $(FFLAGS) $(HDF5)
	./$(PROG)

.PHONY: clean
clean :
	rm -f *.o
	rm -f *.hdf5
	rm -f $(PROG)
