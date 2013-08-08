cd `dirname $0`


rm -r build
mkdir build

cp src/* build/

cd build

gfortran -O3 -c nrtype.f90
gfortran -O3 -c nr.f90
gfortran -O3 -c large_array.f90
gfortran -O3 -o grninfer.out nr.o nrtype.o large_array.o regnw.f90 dsvdcmp.for dpythag.for SumAbsMin.f90 draw.f90 lpr.f

cp grninfer.out ../

rm grninfer.out
