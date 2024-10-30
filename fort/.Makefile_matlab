HOME := /home/david/
FC := gfortran
PFC := f2py
UTILLIBS := $(HOME)libf/
SRCDIR := ./
LIBDIR := ../lib/
DOP853 := $(UTILLIBS)libdop853.a
FRK := $(UTILLIBS)libfrk_light.a
SPICE := $(UTILLIBS)spicelib.a
EXT_SUFFIX := $(shell python3-config --extension-suffix)

OBJECTS :=  globals.o \
			tensorops.o \
			qist.o\
			q_inter.o\
		    makemodel.o\
			genqist.o\
			frkmin.o \
			frkmin_q.o \
			tinysh.o \
			denselight.o \
			cheby.o \
			$(SPICE)

PFFLAGS := -lgomp 
FFLAGS := -O3 -g -mcmodel=large \
	      -Wl,--no-relax\
		  -fPIC\
		  -ffast-math \
		  -march=native \
		  -ffree-line-length-none\
		  -fcheck=bounds\
		  -fopenmp
MFLAGS := -O3 -mcmodel=large \
	      -Wl,--no-relax\
		  -fPIC\
		  -ffast-math \
		  -march=native \
		  -ffree-line-length-none\
		  -fcheck=bounds\
		  -fopenmp\
		  -fno-underscoring\
		  -shared
			 
mat: mqist.a
pq: pyqist$(EXT_SUFFIX)
pqftest: pqftest.f90 $(OBJECTS) pyqist.o
	$(FC) $(FFLAGS) $^ -o $@

%.a : $(SRCDIR)%.f90 $(OBJECTS)
	$(FC) $(MFLAGS) $^ -o $@

%$(EXT_SUFFIX) : $(SRCDIR)%.f90 %.pyf q_inter.o $(OBJECTS)
	$(PFC) --f90flags='$(FFLAGS)' $(PFFLAGS) -c $^ 
	-cp -v $@ $(LIBDIR)

%.pyf : %.f90
	$(PFC) $^ -m $* -h $@ --overwrite-signature

%.o : $(SRCDIR)%.f90 
	$(FC) $(FFLAGS) -c -o $@ $^ 

clean:
	@echo " Cleaning";
	-rm --verbose *.o *.mod *.pyf *$(EXT_SUFFIX)

%.o: %.mod
