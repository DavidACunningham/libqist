HOME := /home/david/
FC := gfortran
PFC := f2py
UTILLIBS := $(HOME)libf/
SRCDIR := ./
DOP853 := $(UTILLIBS)libdop853.a
FRK := $(UTILLIBS)libfrk_light.a
EXT_SUFFIX := $(shell python3-config --extension-suffix)

OBJECTS :=  globals.o \
			tensorops.o \
			packmod.o \
			qist.o\
			q_inter.o\
		    $(DOP853) \
		    $(FRK)


PFFLAGS := -lgomp 
FFLAGS := -O3 -g -mcmodel=large \
	      -Wl,--no-relax\
		  -fPIC\
		  -ffast-math \
		  -march=native \
		  -ffree-line-length-none\
		  -fcheck=bounds\
		  -fopenmp
			 
pq: pyqist$(EXT_SUFFIX)
pqftest: pqftest.f90 $(OBJECTS) pyqist.o
	$(FC) $(FFLAGS) $^ -o $@

%$(EXT_SUFFIX) : $(SRCDIR)%.f90 %.pyf q_inter.o $(OBJECTS)
	$(PFC) --f90flags='$(FFLAGS)' $(PFFLAGS) -c $^ 

%.pyf : %.f90
	$(PFC) $^ -m $* -h $@ --overwrite-signature

%.o : $(SRCDIR)%.f90 
	$(FC) $(FFLAGS) -c -o $@ $^ 

clean:
	@echo " Cleaning";
	-rm --verbose *.o *.mod *.pyf *$(EXT_SUFFIX)

%.o: %.mod
