#! /bin/bash

NAME="f3breg"

pushd ~/wrk/scipy-env && source ./bin/activate && popd

f2py $NAME.f90 -m $NAME -h $NAME.pyf --overwrite-signature
f2py -c $NAME.f90 $NAME.pyf tensorops.f90 odes.f90

mv $NAME.pyf ../
mv $NAME.cpython* ../
