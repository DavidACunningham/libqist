program main
    use qist, only: itraj
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    type(itraj)         :: it
    real(dp), parameter :: t0=0._dp, tf=2._dp*24._dp*3600._dp, tof = tf
    real(dp)            :: itraj_stm(8,8), &
                         & itraj_stt(8,8,8)
    integer i,j

    call it%init(t0, 1._dp, &
                 "/home/david/wrk/nstgro/qist/libqist/fort/test/",&
                 "qist_sol_packed.lightsol")
    itraj_stm = it%stm(1._dp)
    itraj_stt = it%stt(1._dp)
print *, "ITRAJ STM"
    do i = 1,8
        print *, real(itraj_stm(i,:),4)
    end do
    print *, "ITRAJ STT"
    do i = 1,8
        print *, "PAGE", i
    do j = 1,8
        print *, real(itraj_stt(i,j,:),4)
    end do
    end do

end program main
