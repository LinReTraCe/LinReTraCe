.SUFFIXES: .F90

#CC = ifort


CC = gfortran-5.1 #mpif90 #ifort #gfortran
CC = gfortran-5.1 -I/usr/lib64/openmpi/1.4-gcc/include -m64 -pthread -I/usr/lib64/openmpi/1.4-gcc/lib -L/usr/lib64/openmpi/1.4-gcc/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -traceback -check all

CC = mpif90

#CC = gfortran-5.1 -I/usr/lib64/openmpi/1.4-gcc/include -I/usr/lib64/openmpi/1.4-gcc/lib -L/usr/lib64/openmpi/1.4-gcc/lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -lquadmath

#CFLAGS = -g -O0 #-fbounds-check -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall -Wextra -Wconversion #-O0 #debugging

#CFLAGS = -g -O0 -traceback -check bounds -check uninit -fpe-all=3 # ifort debug
#CFLAGS = -g -fbounds-check -fbacktrace  -Wall -Wextra -Wconversion -O0 #  gfortran debugging
CFLAGS = -O3 -funroll-loops -ffast-math #-march=native -flto -fstack-arrays  # gfortran production run
#CFLAGS = -O3 #-xHost -ipo  # ifort production





CPP = cpp -P -C -traditional
# for MPI:
CPPFLAGS = -DGFORTRAN -DFLUSH -DMPI #-DGFORTRAN #-DDEBUG # GFORTRAN flag also works with ifort. gfortran does not work properly without.
# need GFORTRAN flag onn HCLM for mpif90...
LDFLAGS =

#intel
#CFLAGS = -O3 -fast -ipo -xHost #production

PROG=main
LAPACK=/home/maggio/lib/lapack/lapack-3.7.1/liblapack.a
BLAS=/home/maggio/lib/blas/BLAS-3.7.1/libblas.a

#link library
LIBS= digamma/cern_digamma.a -lquadmath $(LAPACK) $(BLAS) #-L/usr/local/lib/libquadmath.so.0.0.0

OBJ=params.o mpi_org.o estruct.o root.o rootQ.o response.o $(PROG).o


.F90.o: $(OBJ) $(PROG)
	$(CPP) $(CPPFLAGS) $< $*.f90
	$(CC) $(CFLAGS) $(FFLAGS) -c  $*.f90


all: $(OBJ) $(PROG)


$(PROG): $(OBJ)
	$(CC) -o $(PROG) $(CFLAGS) $(OBJ) $(LDFLAGS) $(LIBS)
	mv main linretrace

#	$(CC) -o $(PROG) -lquadmath $(CFLAGS) $(OBJ) $(LDFLAGS) $(LIBS)


#if more o-files needed...
locate.o: locate.f90
	$(CC) -c $(CFLAGS) $(INCLUDES) locate.f90

clean :
	rm -f *.o *.mod *.f90 *genmod* $(PROG) linretrace
