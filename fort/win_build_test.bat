@ECHO OFF 

SET SPICE="C:\Users\tester\Downloads\toolkit\toolkit\lib\spicelib.lib"

SET PRODUCT=unit_test_run

SET HERE="%cd%"
SET SRCDIR=%HERE%\src\
SET OBJDIR=%HERE%\build\
SET TESTOBJDIR=%HERE%\test\build\
SET LIBDIR=%HERE%\lib\
SET MODDIR=%HERE%\mod\
SET EXEDIR=%HERE%\exe\
SET TESTDIR=%HERE%\test\

SET SOURCES=^
             %SRCDIR%globals.f90^
             %SRCDIR%tensorops.f90^
             %SRCDIR%frkmin_q.f90^
             %SRCDIR%frkmin.f90^
             %SRCDIR%denselight.f90^
             %SRCDIR%qist.f90^
             %SRCDIR%cheby.f90^
             %SRCDIR%quat.f90^
             %SRCDIR%subspice.f90^
             %SRCDIR%tinysh.f90^
             %SRCDIR%makemodel.f90^
             %SRCDIR%genqist.f90^
             %TESTDIR%test_globals.f90^
             %TESTDIR%test_util.f90^
             %TESTDIR%findiffmod.f90^
             %TESTDIR%test_sh.f90^
             %TESTDIR%test_cheby.f90^
             %TESTDIR%test_quat.f90^
             %TESTDIR%test_tensorops.f90^
             %TESTDIR%test_frkmin.f90^
             %TESTDIR%test_frkmin_q.f90^
             %TESTDIR%test_genqist.f90^
             %TESTDIR%test_makemodel.f90

SET OBJECTS=^
             %OBJDIR%globals.obj^
             %OBJDIR%tensorops.obj^
             %OBJDIR%frkmin_q.obj^
             %OBJDIR%frkmin.obj^
             %OBJDIR%denselight.obj^
             %OBJDIR%qist.obj^
             %OBJDIR%cheby.obj^
             %OBJDIR%quat.obj^
             %OBJDIR%subspice.obj^
             %OBJDIR%tinysh.obj^
             %OBJDIR%makemodel.obj^
             %OBJDIR%genqist.obj^
             %OBJDIR%test_globals.obj^
             %OBJDIR%test_util.obj^
             %OBJDIR%findiffmod.obj^
             %OBJDIR%test_sh.obj^
             %OBJDIR%test_cheby.obj^
             %OBJDIR%test_quat.obj^
             %OBJDIR%test_tensorops.obj^
             %OBJDIR%test_frkmin.obj^
             %OBJDIR%test_frkmin_q.obj^
             %OBJDIR%test_genqist.obj^
             %OBJDIR%test_makemodel.obj

SET FC=ifx

SET FFLAGS_MOD= /syntax-only^
                /module:%MODDIR%^
                /c

SET FFLAGS=^
           /check:bounds^
           /module:%MODDIR%^
           /heap-arrays:10240^
           /Od^
           /debug:full^
           /traceback

%FC% %FFLAGS_MOD% %SOURCES%
%FC% %FFLAGS% /c %SOURCES%
MOVE *.obj %OBJDIR%
DEL %TESTDIR%%PRODUCT%.exe
%FC% %FFLAGS% %TESTDIR%%PRODUCT%.f90 %SPICE% %OBJECTS%
MOVE %PRODUCT%.exe %TESTDIR%%PRODUCT%.exe
DEL %OBJECTS%
ECHO unit_test_run.exe generated.
