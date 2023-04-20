#! /bin/bash
BUILDDIR=temp
DLSODE=/home/david/libf/libodepack.a
cd $BUILDDIR
gfortran -fPIC -c ../cheby.f90 ../vcheby.f90 ../f3breg.f90 ../tensorops.f90 ../odes.f90 
ar -r libreg.a cheby.o \
			    vcheby.o \
			    f3breg.o \
			    tensorops.o \
			    odes.o 
ar -rcT liball.a $DLSODE ./libreg.a
mv liball.a ../
mv libreg.a ../
cd ..
rm $BUILDDIR/*
