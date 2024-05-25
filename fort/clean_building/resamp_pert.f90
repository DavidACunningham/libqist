! Resamples SPICE to generate ephemeris of celestial bodies
! along with reference trajectory
program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use genqist, only: make_spice_subset
    character(len=1000) :: arg
    call get_command_argument(1,arg)
    call make_spice_subset(trim(adjustl(arg)))
end program main
