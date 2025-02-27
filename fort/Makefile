.SUFFIXES:
SPICE := /home/david/libf/spicelib.a
FC := gfortran
PFC := f2py
SRCDIR := ./src/
OBJDIR := ./build/
TESTOBJDIR := ./test/build/
TESTANCDIR := ./test/anc/
LIBDIR := ./lib/
MODDIR := ./mod/
EXEDIR := ./exe/
TESTDIR := ./test/
ANCDIR := ./anc/
EXT_SUFFIX := $(shell python3-config --extension-suffix)
FOREXT := f90

define source-to-extension
  $(strip \
	  $(foreach ext,$(FOREXT),\
	  	$(subst .$(ext),.$2,$(filter %.$(ext),$1))))
endef
define changefolder
  $(strip \
	  $(foreach sourcefile,$1,\
	  	$(subst $3,$2,$1)))
endef

SOURCES := $(SRCDIR)globals.f90 \
		   $(SRCDIR)tensorops.f90 \
	       $(SRCDIR)frkmin.f90 \
		   $(SRCDIR)denselight.f90
OBJECTS :=$(call source-to-extension,$(SOURCES),o)
ANCHORS :=$(call source-to-extension,$(SOURCES),anc)
OBJECTS := $(call changefolder,$(OBJECTS),$(OBJDIR),$(SRCDIR))
ANCHORS := $(call changefolder,$(ANCHORS),$(ANCDIR),$(SRCDIR))

GQ_SOURCES := $(SOURCES) \
		    $(SRCDIR)cheby.f90 \
	        $(SRCDIR)quat.f90 \
	        $(SRCDIR)subspice.f90 \
	        $(SRCDIR)tinysh.f90 \
	        $(SRCDIR)makemodel.f90 \
	        $(SRCDIR)frkmin_q.f90 \
	        $(SRCDIR)genqist.f90
GQ_OBJECTS :=$(call source-to-extension,$(GQ_SOURCES),o)
GQ_ANCHORS :=$(call source-to-extension,$(GQ_SOURCES),anc)
GQ_OBJECTS := $(call changefolder,$(GQ_OBJECTS),$(OBJDIR),$(SRCDIR))
GQ_ANCHORS := $(call changefolder,$(GQ_ANCHORS),$(ANCDIR),$(SRCDIR))
Q_SOURCES := $(SOURCES) \
			 $(SRCDIR)qist.f90
Q_OBJECTS :=$(call source-to-extension,$(Q_SOURCES),o)
Q_ANCHORS :=$(call source-to-extension,$(Q_SOURCES),anc)
Q_OBJECTS := $(call changefolder,$(Q_OBJECTS),$(OBJDIR),$(SRCDIR))
Q_ANCHORS := $(call changefolder,$(Q_ANCHORS),$(ANCDIR),$(SRCDIR))

I_SOURCES := $(Q_SOURCES) \
			 $(SRCDIR)q_inter.f90
I_OBJECTS :=$(call source-to-extension,$(I_SOURCES),o)
I_ANCHORS :=$(call source-to-extension,$(I_SOURCES),anc)
I_OBJECTS := $(call changefolder,$(I_OBJECTS),$(OBJDIR),$(SRCDIR))
I_ANCHORS := $(call changefolder,$(I_ANCHORS),$(ANCDIR),$(SRCDIR))

T_SOURCES := $(GQ_SOURCES) \
			 $(TESTDIR)test_globals.f90\
			 $(TESTDIR)test_util.f90\
			 $(SRCDIR)qist.f90\
			 $(TESTDIR)findiffmod.f90\
			 $(TESTDIR)test_sh.f90\
			 $(TESTDIR)test_cheby.f90\
			 $(TESTDIR)test_quat.f90\
			 $(TESTDIR)test_tensorops.f90\
			 $(TESTDIR)test_frkmin.f90\
			 $(TESTDIR)test_frkmin_q.f90\
			 $(TESTDIR)test_genqist.f90\
			 $(TESTDIR)test_makemodel.f90
T_OBJECTS :=$(call source-to-extension,$(T_SOURCES),o)
T_ANCHORS :=$(call source-to-extension,$(T_SOURCES),anc)
T_OBJECTS := $(call changefolder,$(T_OBJECTS),$(OBJDIR),$(SRCDIR))
T_ANCHORS := $(call changefolder,$(T_ANCHORS),$(ANCDIR),$(SRCDIR))
T_OBJECTS := $(call changefolder,$(T_OBJECTS),$(TESTOBJDIR),$(TESTDIR))
T_ANCHORS := $(call changefolder,$(T_ANCHORS),$(TESTANCDIR),$(TESTDIR))

ifneq ($(MODDIR),)
  $(shell test -d $(MODDIR) || mkdir -p $(MODDIR))
	FFLAGS+= -J$(MODDIR)
endif

PFFLAGS := -lgomp \
		   --build-dir $(MODDIR)
PFLAGS := -O3 -mcmodel=large \
	      -Wl,--no-relax\
		  -fPIC\
		  -ffast-math \
		  -march=native \
		  -ffree-line-length-none\
		  -fcheck=bounds
		  
FFLAGS := -O3  -mcmodel=large \
	      -Wl,--no-relax\
		  -fPIC\
		  -ffast-math \
		  -march=native \
		  -ffree-line-length-none\
		  -fcheck=bounds\
		  -J$(MODDIR)
MFLAGS := -O3 -mcmodel=large \
	      -Wl,--no-relax\
		  -fPIC\
		  -ffast-math \
		  -march=native \
		  -ffree-line-length-none\
		  -fcheck=bounds\
		  -fno-underscoring\
		  -shared\
		  -J$(MODDIR)

			 
MAKEMOD.F08 = $(FC) $(FFLAGS) -fsyntax-only -c 
COMPILE.F08 = $(FC) $(FFLAGS)  -c 

mat: mqist.a
test: $(TESTDIR)unit_test_run 
genqist: make_qist_existing_kernel
genkernel: make_new_kernel
gencheby: make_cheby_kernel
pq: pyqist$(EXT_SUFFIX)
pqtest: pyqist.a
qist: qist.a
modfiles: $(T_ANCHORS) $(I_ANCHORS)

.PHONY: clean
clean:
	@echo " Cleaning";
	-rm --verbose $(SRCDIR)*.o *.pyf *.anc *$(EXT_SUFFIX)
	-test -d $(MODDIR) && rm --verbose $(MODDIR)*
	-test -d $(ANCDIR) && rm --verbose $(ANCDIR)*
	-test -d $(OBJDIR) && rm --verbose $(OBJDIR)*
	-test -d $(TESTOBJDIR) && rm --verbose $(TESTOBJDIR)*
	-test -d $(TESTANCDIR) && rm --verbose $(TESTANCDIR)*
	-test -d $(TESTDIR) && rm --verbose $(TESTDIR)*.subspice
	-test -d $(TESTDIR) && rm --verbose $(TESTDIR)*.odesolution
	-test -d $(TESTDIR) && rm --verbose $(TESTDIR)*.qist
	-test -d $(TESTDIR) && rm --verbose $(TESTDIR)*.rot
	-test -d $(TESTDIR) && rm --verbose $(TESTDIR)unit_test_run

make_resample: $(GQ_OBJECTS)
	$(FC) $(FFLAGS)  $(EXEDIR)$@.f90 $^ $(SPICE) -o $@
	-mv -v $@ $(EXEDIR)
make_rot: $(GQ_OBJECTS)
	$(FC) $(FFLAGS)  $(EXEDIR)$@.f90 $^ $(SPICE) -o $@
	-mv -v $@ $(EXEDIR)
make_cheby_kernel: $(GQ_OBJECTS)
	$(FC) $(FFLAGS)  $(EXEDIR)$@.f90 $^ $(SPICE) -o $@
	-mv -v $@ $(EXEDIR)
make_new_kernel: $(GQ_OBJECTS)
	$(FC) $(FFLAGS)  $(EXEDIR)$@.f90 $^ $(SPICE) -o $@
	-mv -v $@ $(EXEDIR)
make_qist_existing_kernel: $(GQ_OBJECTS)
	$(FC) $(FFLAGS)  $(EXEDIR)$@.f90 $^ $(SPICE) -o $@
	-mv -v $@ $(EXEDIR)
$(TESTDIR)unit_test_run: $(T_OBJECTS)
	$(FC) $(FFLAGS)  $@.f90 $^ $(SPICE) -o $@
qist.a: $(Q_OBJECTS)
	$(FC) $(FFLAGS) -shared $^ -o $@
	-cp -v $@ $(LIBDIR)
%.a : $(SRCDIR)%.f90 $(I_OBJECTS)
	$(FC) $(MFLAGS) $(SPICE) $^ -o $@
%$(EXT_SUFFIX) : $(SRCDIR)%.f90 $(I_OBJECTS)
	$(PFC) -m $* --f90flags='$(PFLAGS)' $(PFFLAGS) -c $^ 
	-cp -v $@ $(LIBDIR)
$(ANCDIR)%.anc: $(SRCDIR)%.f90
	$(MAKEMOD.F08) $<
	@touch $@
$(OBJDIR)%.o : $(ANCDIR)%.anc 
	$(FC) $(FFLAGS) -c -o $(OBJDIR)$*.o $(SRCDIR)$*.f90
	@touch $@
$(TESTANCDIR)%.anc: $(TESTDIR)%.f90
	$(MAKEMOD.F08) $<
	@touch $@
$(TESTOBJDIR)%.o : $(TESTANCDIR)%.anc 
	$(FC) $(FFLAGS) -c -o $(TESTOBJDIR)$*.o $(TESTDIR)$*.f90
	@touch $@
