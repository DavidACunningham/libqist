! Resamples SPICE to generate ephemeris of celestial bodies
! along with reference trajectory
program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cheby, only: spice_subset
    implicit none
    type(spice_subset)            :: subspice
    real(dp)                      :: a, b 
    integer                       :: bodlist(4), deg
    character(len=:), allocatable :: kernel, filename

    a = 0._dp
    b = 2._dp * 24._dp * 3600._dp
    deg = 2000
    kernel = "/home/david/wrk/nstgro/qist/kernels/mk.tf"
    filename = "./perturbed_reference.subspice"
    bodlist = [10, 301, 5, -998] ! sun, moon, jupiter barycenter, reference
    call subspice%init( &
                       kernel, & ! Spice kernel
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
