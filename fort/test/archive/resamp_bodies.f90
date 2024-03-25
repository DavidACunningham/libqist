! Resamples SPICE to generate ephemeris of celestial bodies ONLY 
! For use in generating reference trajectories
program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cheby, only: spice_subset
    implicit none
    type(spice_subset)            :: subspice
    real(dp)                      :: a, b 
    integer                       :: bodlist(9), deg
    character(len=:), allocatable :: kernel, filename
    a = 0._dp
    b = 2._dp * 24._dp * 3600._dp
    deg = 2000
    kernel = "/home/david/wrk/nstgro/qist/kernels/mk_noorb.tf"
    filename = "./resampled_celestial_bodies.subspice"
    bodlist = [1, 2, 4, 5, 6, 7, 8, 301, 10] ! All planet system barycenters, moon, sun
    call subspice%init( &
                       kernel, & ! spice kernel
                       399, &    ! central body
                       bodlist, & ! list of bodies to resample
                       a, & ! epoch start
                       b, & ! epoch end
                       deg & ! degree of fit
                      )
    open(file=filename,unit=73,access="stream",status="replace")
    call subspice%write(73)
    close(73)
end program main
