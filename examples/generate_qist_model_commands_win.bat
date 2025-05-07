@echo off
REM These two commands use the example namelist that is in this directory to 
REM     generate a new QIST model for an example reference orbit around Deimos.

REM The arguments for both programs are are:
REM     1. A string parameter containing the filename of the namelist file
REM   2. A logical parameter to select whether to resample the SPICE 
REM	      ephemeris (must be True if not done before)
REM	  3. A logical parameter to select whether to make a new interpolated 
REM         rotation matrix (must be True if not done before)

REM The following command uses the initial conditions and dynamics parameters in 
REM     the namelist file to numerically integrate the reference orbit and write
REM     the resulting trajectory to a SPICE kernel.

..\generate_kernel\x64\Release\generate_kernel.exe "./curve_deimos_config_namelist_2026Nov26120000002026Nov2617450000.nml" True True
..\generate_qist_model\x64\Release\generate_qist_model.exe "./curve_deimos_config_namelist_2026Nov26120000002026Nov2617450000.nml" True True
pause