include ../../make_include

# NOTE: psi_fast.F90 lives in this directory but is NOT built here.  It is a
# Fortran module, so it is compiled by the parent Makefile with $(FC) into
# src_linretrace/psi_fast.o -- module files are compiler specific and would
# not be usable from $(FCDG) builds.  See the psi_fast.o rule one level up.
#
# The archive below is the deprecated CERNLIB path.  It is still built and
# linked so that the commented-out wpsipg / wpsipghp calls in fermi.f90,
# response.F90 and root.F90 can be restored without touching the build; the
# two kmax settings only affect those routines.
kmax_double = 7
kmax_quad = 16

CPP = cpp
CPPFLAGS = -P -traditional

ifeq (${FCDG},gfortran)
	CPPFLAGS += -DGFORTRAN
endif

ARCH = ar
ARCHFLAGS = crs

# overwrite the make_include file
FFLAGS = -O3

.PHONY: all
all: cern_digamma.a

wpsipg.o: wpsipg.F77
	${CPP} ${CPPFLAGS} -D kmax=${kmax_double} $< wpsipg.f
	${FCDG} -c wpsipg.f -o $@ ${FFLAGS}
	rm -f wpsipg.f

wpsipg_hp.o: wpsipg_hp.F77
	${CPP} ${CPPFLAGS} -D kmax=${kmax_quad} $< wpsipg_hp.f
	${FCDG} -c wpsipg_hp.f -o $@ ${FFLAGS} -lquadmath
	rm -f wpsipg_hp.f

%.o: %.f
	${FCDG} -c $< -o $@ ${FFLAGS}

cern_digamma.a: wpsipg.o abend.o cpsipg.o lenocc.o mtlset.o mtlprt.o wpsipg_hp.o
	rm -f $@
	${ARCH} ${ARCHFLAGS} $@ $^

.PHONY: clean
clean :
	rm -f cern_digamma.a *.o *.mod
