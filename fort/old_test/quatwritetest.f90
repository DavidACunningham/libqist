program main
    use, intrinsic :: iso_fortran_env, only: dp =>real64
    use genqist, only: make_rotation
    use quat, only: rothist
    implicit none
    integer, parameter       :: nnodes=1000
    type(rothist)            :: rotb
    real(dp) t0
    integer num
    call make_rotation("/home/david/wrk/nstgro/qist/libqist/fort/data/moon_PA_config.nml")
    open(file="../data/moon_PA_one_orbit.rot", status="old", access="stream", newunit=num)
        call rotb%read(num)
    close(num)
    t0 = 769269009.185_dp
    print *, rotb%call(t0+100._dp)
    print *, rotb%calldot(t0+100._dp)
    print *, rotb%callddot(t0+100._dp)
end program main
