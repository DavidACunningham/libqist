#! /bin/bash
HOME=/home/david
UTILLIBS=$HOME/libf
SRCDIR=src
LIBDIR=./lib
TEMPDIR=temp
DOP853=$UTILLIBS/libdop853.a
FRK=$UTILLIBS/libfrk_light.a
SPICE=$UTILLIBS/spicelib.a
GENOBJECTS="globals.o \
  		tensorops.o \
	    makemodel.o\
		genqist.o"
RUNOBJECTS="globals.o \
			tensorops.o \
			qist.o\
			q_inter.o"
mkdir $TEMPDIR
cd $SRCDIR
cp -v $RUNOBJECTS ../$TEMPDIR
cd ../$TEMPDIR
ar -x $FRK
ar -r libqist.a *.o
rm -v $RUNOBJECTS
cd ../$SRCDIR
pwd
cp -v $GENOBJECTS ../$TEMPDIR
cd  ../$TEMPDIR
pwd
ar -x $SPICE
ar -r libgenqist.a *.o
cd ..
mv -v $TEMPDIR/libgenqist.a $LIBDIR
mv -v $TEMPDIR/libqist.a $LIBDIR
rm -r $TEMPDIR
