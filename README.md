# QIST
> Version 0.1 (beta)
>
> Quadratic Interpolated State Transition: A library for fast propagation of
> second-order relative motion around arbitrary SPICE kernels.

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


## Building from source

OS X & Linux:

Download the source and unzip it into a directory I'll refer to as ``$LIBQIST``.
First, in the ``$LIBQIST/fort`` directory, edit the second line of the file ``Makefile``
to point to your local SPICE build's shared library file, usually called ``spicelib.a``.
This is the only external dependency of QIST.
If you don't have SPICE, you can download and build it from the JPL NAIF website at 
(<https://naif.jpl.nasa.gov/naif/toolkit.html>).
QIST has been tested with the Fortran version of SPICE.

Next, build the module files:
```sh
user@system:$LIBQIST/fort$ make modfiles
```
### Building the relative motion runtime
To propagate relative motion models, you need to build a version of the QIST
relative motion propagation runtime.

If you want to build the Python interface using f2py, run:
```sh
user@system:$LIBQIST/fort$ make -j pq
```
The 'j' argument is optional but recommended, enabling multi core compilation.

If you want to build the Matlab interface, run:
```sh
user@system:$LIBQIST/fort$ make -j mqist
```
If you want to build the native Fortran shared library, run:
```sh
user@system:$LIBQIST/fort$ make -j qist
```

### Model generation setup

If you'll be generating your own QIST models, you need to build the model
generation facilities. To do this, run:
```sh
user@system:$LIBQIST/fort$ make -j make_resample
user@system:$LIBQIST/fort$ make -j make_rot
user@system:$LIBQIST/fort$ make -j genqist
user@system:$LIBQIST/fort$ make -j make_qist_existing_kernel
user@system:$LIBQIST/fort$ make -j make_new_kernel
```

### Test build

To build unit tests, run:
```sh
user@system:$LIBQIST/fort$ make -j test
```

To build unit tests, run:
```sh
user@system:$LIBQIST/fort$ make -j test
```
Before running the unit tests:
 * the user must edit the files ``$LIBQIST/kernels/mk_test.tf`` and ``$LIBQIST/kernels/mk_test_withorbit.tf`` to reflect the SPICE setup on the user’s machine.
 * NOTE: The ``.bsp`` file ``test_orbit.bsp`` is generated as part of the unit test process, and should not be available on the user system before running unit tests.
To run unit tests, change to the ``$LIBQIST/fort/test`` directory and run:
```sh
user@system:$LIBQIST/fort/test$ unit_test_run
```

Any output with the word “FAIL” means the tests have failed.


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
