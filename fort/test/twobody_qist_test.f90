program main
    use cheby, only: spice_subset
    use genqist, only: gqist
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    use frkmain, only: Odesolution
    implicit none
    type(gqist)        :: qist
    type(odesolution)  :: sol
    real(dp)           :: cdum(2,2), sdum(2,2)

    call qist%init(0._qp, 2._qp*24._qp*3600._qp, &
                 & "./twobody_resample.subspice", &
                 & -999, &
                 &  399, &
                 &  [10, 301, 5], &
                 &  398600.5_dp, &
                 &  6378.137_dp, &
                 & [1.32712440041279419e11_dp, &
                 &  4902.800118_dp, &
                 &  126712764.1_dp], &
                 &  .False., &
                 &  cdum, sdum, .true.)
    sol = qist%integrate(3600._qp, 5*3600._qp)
    

    open(unit=75, file="test_tbod_sol.odesolution", access="stream", status="replace")
    call sol%write(75,dp)
    close(75)

    print *, sol%call(3600._qp*2)


    contains
        subroutine print_to_file(fname, var)
            integer io,j
            real(dp), intent(in) :: var(:)
            character(len=*) :: fname
            open(newunit=io, file=trim(adjustl(fname))//".txt")
            do j = 1,size(var)
            write(io,*) var(j)
            end do
            close(io)
        end subroutine print_to_file
end program main
