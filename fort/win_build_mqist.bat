@ECHO OFF
REM
REM   win_build_mqist.bat
REM
REM      Builds matlab QIST relative motion runtime for the Intel Fortran Environment 64bit
REM      makenv-tag: PC-WINDOWS-64BIT-IFORT
REM
REM   Version
REM
REM      win_build_qist.bat Script Version 1.0.0 5-Apr-2025
REM   

REM
REM   Product name
REM
SET PRODUCT=mqist

REM
REM   Directory info
REM 
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
             %SRCDIR%frkmin.f90^
             %SRCDIR%denselight.f90^
             %SRCDIR%qist.f90^
             %SRCDIR%q_inter.f90^
             %SRCDIR%mqist.f90
SET OBJECTS=^
             %OBJDIR%globals.obj^
             %OBJDIR%tensorops.obj^
             %OBJDIR%frkmin.obj^
             %OBJDIR%denselight.obj^
             %OBJDIR%qist.obj^
             %OBJDIR%q_inter.obj^
             %OBJDIR%mqist.obj

REM
REM   Set the FC variable to contain the name of the compiler that will
REM   be used to build the toolkit.
REM      Intel Fortran          - ifx
REM
SET FC=ifx

REM
REM   Set the FFLAGS environment variable.
REM
REM      /nodebug         - Production level code, no debug info.
REM      /check:bounds    - Enable array bounds checking
REM
SET FFLAGS_MOD= /syntax-only^
                /module:%MODDIR%^
                /c

SET FFLAGS=  /MP^
	     /heap-arrays0^
	     /warn:interfaces^
	     /O2^
	     /module:%MODDIR%^
             /nodebug^
	     /check:bounds^
             /c
REM
REM   Set the LD environment variable to define the linker we will use
REM   to link.
REM
SET LD=link

REM
REM   Set the LDOPS environment variable to define the options we will
REM   pass to the linker.
REM
SET LDOPS= /nologo /verbose /out:%LIBDIR%%PRODUCT%.dll

%FC% %FFLAGS_MOD% %SOURCES%
%FC% %FFLAGS% %SOURCES%
move *.obj %OBJDIR%
%LD% %LDOPS% %OBJECTS%

ECHO %PRODUCT%.dll generated.
