
#expansion order used in the digamma functions. 1<=kmax<=12
kmax_double=7 # <=7 . higher accuracy beyond double precision for typical cases.
kmax_quad=16 # 16 is maximum

COMPILER=ifort #gfortran
COMPILER=gfortran-5.1
CPPFLAGS="-DGFORTRAN" # use if gfortran is used
#FFLAGS="-g -O0" #-malign-double -fbounds-check -fbacktrace -Wall -Wextra -Wconversion -ffpe-trap=zero,overflow,underflow -O0 " #debugging
FFLAGS= -O3 #-malign-double -funroll-loops -ffast-math -march=native -flto -fstack-arrays #production
#FFLAGS="-O3 -xHost -ipo" # -fpmodel precise / -fpmodel except / -fpmodel strict

rm *.a

FILES=$(ls *.F77 *.F90 | grep -v hp)
for f in `echo $FILES`
do
f2="${f%.*}"
cpp -P -C -traditional $CPPFLAGS -D kmax=$kmax_double $f > $f2.f # precompiling .F?? files 
done

FILES=$(ls *.F77 *.F90 | grep hp)
for f in `echo $FILES`
do
f2="${f%.*}"
cpp -P -C -traditional $CPPFLAGS -D kmax=$kmax_quad $f > $f2.f # precompiling .F?? files 
done

FILES=$(ls *.f | grep -v hp)
for f in `echo $FILES`
do
$COMPILER $FFLAGS -c $f   #compiling all fortran files into object files
f2="${f%.*}"
xiar rv cern_digamma.a $f2.o # collect object files into library
done
FILES=$(ls *.f | grep hp)
for f in `echo $FILES`
do
$COMPILER $FFLAGS -c -lquadmath $f   #compiling all fortran files into object files
f2="${f%.*}"
xiar rv cern_digamma.a $f2.o # collect object files into library
done






