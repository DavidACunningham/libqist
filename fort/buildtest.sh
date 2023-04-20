#! /bin/bash

NAME="f3breg"


gfortran -c $NAME.f90 tensorops.f90
gfortran *.o

./a.out
