# QIST
> Version 0.1 (beta)
>
> Quadratic Interpolated State Transition: A library for fast propagation of
> second-order relative motion around arbitrary SPICE kernels.
> Authors: David Cunningham and Ryan P. Russell.
> Copyright 2025.

Current and future space mission planning can require relative motion in regimes 
that are highly non-Keplerian and potentially non-linear.
QIST is an attempt to build a generalized relative motion framework for a broad class
of gravitational dynamics. 

QIST allows:
* Numerical integration of STMs and STTs around arbitrary SPICE SPK trajectories from a beginning simulation time ``t0`` to an end simulation time ``tf`` with accelerations including:
	* Ephemeris gravitational dynamics from other SPICE bodies
	* Nonsingular spherical harmonics gravitation of the central body (Pines formulation)
* Fast propagation of trajectories between any two points ``ta`` and ``tb`` such that ``t0 <= ta <= tb <= tf`` relative to the reference trajectory with quadratic accuracy
* Propagation of relative covariance between ta and tb through direct access to the STM and STT.


## Installation

Currently, QIST must be built from source.

**Terminology for this section:**

* The root directory of the QIST source: ``$LIBQIST``
* Shell prompt: ``user@system$``

Currently, QIST compilation is documented and tested for Linux systems with the gfortran compiler installed. Installation for OS X should be nearly identical with a configured gfortran compiler, with the exception of some file extension changes (i.e. shared library extension) Installation instructions for Windows computers will be supplied later.

### System/Software Requirements
* Mandatory:
    * GNU make
    * gfortran
    * SPICE
* Optional:
    * Python 3 (versions from 3.8 to 3.11 tested) including packages:
        * numpy [https://pypi.org/project/numpy/](https://pypi.org/project/numpy/)
        * matplotlib [https://pypi.org/project/matplotlib/](https://pypi.org/project/matplotlib/)
        * spiceypy [https://pypi.org/project/spiceypy/](https://pypi.org/project/spiceypy/)
    * MATLAB (versions from R2017 to R2024 tested)

NOTE: f2py for Python 3.12 switched to the meson build backend, which has introduced changes to the build process. These changes have not yet been incorporated into the QIST build system. If this makes sense to you, you probably have a good shot at updating things yourself–just have a look at the Makefile.

### OS X & Linux Installation

#### Pre-install setup (Mandatory for all users)

1. Download the source and unzip it into the ``$LIBQIST`` directory. 
2. In the ``$LIBQIST/fort`` directory, edit line 2 of the file Makefile to point to your local SPICE build's shared library file, usually called spicelib.a. This is the only external dependency of QIST. If you don't have SPICE, you can download and build it from the JPL NAIF website at [https://naif.jpl.nasa.gov/naif/toolkit.html](https://naif.jpl.nasa.gov/naif/toolkit.html). QIST has been tested with the Fortran version of SPICE.
3. (OPTIONAL–only for Python wrapper users) In the ``$LIBQIST/fort`` directory, edit line 3 of the file Makefile to point to the desired python 3 virtual environment activate script.
4. OS-X users: alpha testing shows that the removal of the following two gfortran compiler options in the makefile may be necessary:
    * If you’re unable to generate ``\*.o`` files, remove ``\-mcmodel=large``
    * If you get a linker error, remove ``\-Wl,`` ``–no-relax``

   

#### Temporary Build Notes

This software is still largely in prototype phase, particularly with respect to the build system and filepath references. I (David C) am in the process of adapting the build files to remove references to absolute filepaths on the user system. As of now, the makefile and some SPICE metakernels contain some references to the path /home/david/wrk/nstgro/. You can find these using, e.g., grep. They need to be replaced with the equivalent path value on your system. Fortran source code files should be free of these path references. Explicitly, the first value in the ``PATH_VALUES`` variable should be where you keep your SPICE ephemerides, with optional additional ``PATH_VALUES``.

#### Recommended SPICE Setup

Some familiarity with SPICE is useful for using QIST. SPICE documentation is available at [https://naif.jpl.nasa.gov/pub/naif/toolkit\_docs/C/index.html](https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/index.html).  
The following NAIF kernels are needed to run the QIST test suite:

* de440.bsp ([https://naif.jpl.nasa.gov/pub/naif/generic\_kernels/spk/planets/de440.bsp](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp))  
* naif0012.tls  ([https://naif.jpl.nasa.gov/pub/naif/generic\_kernels/lsk/naif0012.tls](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls) or, if using MS Windows, [https://naif.jpl.nasa.gov/pub/naif/generic\_kernels/lsk/naif0012.tls.pc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls.pc))  
* moon\_pa\_de440\_200625.bpc ([https://naif.jpl.nasa.gov/pub/naif/generic\_kernels/pck/moon\_pa\_de440\_200625.bpc](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/moon_pa_de440_200625.bpc))  
* Moon\_de440\_220930.tf ([https://naif.jpl.nasa.gov/pub/naif/generic\_kernels/fk/satellites/moon\_de440\_220930.tf](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/fk/satellites/moon_de440_220930.tf))

The following NAIF kernel is only needed if building the included Deimos example

* mar097.bsp ([https://naif.jpl.nasa.gov/pub/naif/generic\_kernels/spk/satellites/mar097.bsp](https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/mar097.bsp)) 

If you are using SPICE for the first time along with QIST, the shell script ``$LIBQIST/kernels/get_naif_kernels.sh`` can be invoked to download all necessary kernels using wget on machines where this utility is available (e.g. \*nix, MacOS, git bash, etc.). The file ``$LIBQIST/kernels/kernel_links`` is a text file containing newline-separated links of all kernels necessary for running QIST tests and examples.

#### Building and running unit tests

Since QIST is in active development, it is recommended for all new users to build and run the unit test suite to confirm QIST functionality.

To build unit tests, run 
```sh
user@system:$LIBQIST/fort$ make test_build
```
NOTE: Depending on system performance/number of cores available, building can take up to 30 seconds on consumer hardware.  
Before running the unit tests:

* Make sure the user has a compiled copy of Fortran SPICE (all versions for compatible systems can be downloaded and compiled or installed from [https://naif.jpl.nasa.gov/naif/toolkit\_FORTRAN.html](https://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html))  
*  the user must edit the files ``$LIBQIST/kernels/mk_test.tf`` and ``$LIBQIST/kernels/mk_test_withorbit.tf`` to reflect the SPICE setup on the user’s machine.  
* *NOTE:* The ``.bsp`` file ``test_orbit.bsp`` is generated as part of the unit test process, and should not be available on the user system before running unit tests.

To run unit tests, change to the ``$LIBQIST/fort/test`` directory and run:
```sh
user@system:$LIBQIST/fort/test$ unit_test_run
```
NOTE: Since running the unit tests involves checking the performance of analytic derivatives and numerical integration, running the unit test suite can take up to 3 minutes on consumer hardware.  
Any output with the word ``FAIL`` means the tests have failed.

#### Building the relative motion runtime library

To propagate relative motion models, you need to build a version of the QIST  
relative motion propagation runtime.  
NOTE: Depending on system performance/number of cores available, building can take up to 30 seconds on consumer hardware.  
If you want to build the Python interface using f2py, run:  
```sh
user@system:$LIBQIST/fort$ make python_wrapper  
```
The ``j`` argument is optional but recommended, enabling multi core compilation.  
If you want to build the Matlab interface, run:  
```sh
user@system:$LIBQIST/fort$ make matlab_wrapper  
```
If you want to build the native Fortran shared library, run:  
```sh
user@system:$LIBQIST/fort$ make native_qist
```

#### Building the model generation executables

If you'll be generating your own QIST models, you need to build the model  
generation facilities. To do this, run:  
```sh
user@system:$LIBQIST/fort$ make genqist
```

## Usage 

see the ``$LIBQIST/doc`` folder for documentation, including a User's Guide.

Currently, the user's guide can also be accessed in a live version on Google Docs [here](https://docs.google.com/document/d/1CWJEghDJ7qMbSs2IA_GUeV96SGAyMXTjxKXaUPOHEKM/edit?usp=sharing)

After building, see the examples subdirectory for minimal usage examples.

### Overview
To run a QIST model, you need:
* The appropriate (i.e. Python, Fortran or Matlab) QIST binary built from the instructions above
* The model, usually stored as a ``.qist`` file
* A calling program
* Optional but highly recommended, the QIST configuration parameters should be stored in a configuration namelist.

To create a QIST model, you need:
* The appropriate qist generation executables built from the instructions above
* SPICE
* A set of SPICE kernels containing all the data needed for your trajectories
* a configuration namelist. An example of this configuration namelist is in the ``$LIBQIST/examples`` folder. It is recommended to generate these namelists automatically, e.g. with Python, to minimize the introduction of errors. An example namelist generation script has been included in the ``$LIBQIST/examples`` folder.

## Release History
* 0.1
    * Initial beta release

## Support and Contact Info

Please contact David Cunningham for questions about usage, bug reports, and distribution information.

David Cunningham – david.cunningham@utexas.edu

## LICENSING
Copyright 2025. This work is openly licensed via Creative Commons BY-NC 4.0 [https://creativecommons.org/licenses/by-nc/4.0/](https://creativecommons.org/licenses/by-nc/4.0/).

For a commercial license for $1000USD, please follow the instructions at: [https://discoveries.utexas.edu/EULA/v1.0.](https://discoveries.utexas.edu/EULA/v1.0.).
