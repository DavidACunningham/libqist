program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use qist, only: itraj
    implicit none
    character(len=1000) :: arg
    type(Itraj)         :: it
    real(dp)            :: stm(8,8), stt(8,8,8), time, timeb, xa(8), xb(8)
    integer i, j
    call get_command_argument(1,arg)
    print *, trim(adjustl(arg))
    call it%init(trim(adjustl(arg)))
    xa = 0._dp
    xa(:3) = [10., 10., 10.]
    time = 769269009.185_dp + 5._dp*24._dp*3600._dp
    timeb = time + 3600._dp*24._dp*1._dp
    call it%stts_ab(time, timeb, stm, stt)
    ! stm = it%stm(time)
    ! stt = it%stt(time)
    print *, "STM"
    do i=1,8
        print *, real(stm(i,:),4)
    enddo
    print *, "STT"
    do i = 1,8
    print *, "PAGE ", i
    do j = 1,8
            print *, real(stt(i,j,:),4)
    end do
    end do
    xb = it%prop(time, timeb, xa, 2)
    print *, "XB"
    print *, real(xb, 4)
end program main
