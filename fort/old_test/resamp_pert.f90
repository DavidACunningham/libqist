! Resamples SPICE to generate ephemeris of celestial bodies
! along with reference trajectory
program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use genqist, only: make_spice_subset
    call make_spice_subset("/home/david/wrk/nstgro/qist/libqist/fort/test/buildconfig.nml")
end program main
