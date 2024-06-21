module test_cheby
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    use cheby, only: vectorcheb, &
                     spice_subset, &
                     chcall, &
                     chnodes, &
                     chfit, &
                     chderiv, &
                     pi
    implicit none
    contains
    subroutine test_scalar_cheby(testpass)
        logical, intent(inout) :: testpass
        character(len=100) :: msg
        integer, parameter  :: deg=80, tlen=1000
        real(dp), parameter :: tol=1.e-15_dp, dtol=1.e-14_dp
        real(dp), parameter :: a=0._dp, b=2._dp*pi
        real(dp)            :: nodes(deg), fnodes(deg), testpoints(tlen), &
                               truth(tlen), dtruth(tlen), &
                               cheb(tlen), dcheb(tlen), &
                               coeffs(deg), dcoeffs(deg-1), &
                               diffnorm, ddiffnorm
        integer i

        nodes = chnodes(deg, a, b)
        testpoints = [(a + i*b/tlen, i=1,tlen)]
        fnodes = fitfun(nodes)
        truth = fitfun(testpoints)
        dtruth = dfitfun(testpoints)
        coeffs = chfit(deg,fnodes)
        dcoeffs = chderiv(coeffs,a,b)
        cheb = chcall(a,b,coeffs,testpoints)
        dcheb = chcall(a,b,dcoeffs,testpoints)
        diffnorm = norm2(truth-cheb)/tlen
        ddiffnorm = norm2(dtruth-dcheb)/tlen
        testpass = .true.
        msg = ""
        if (diffnorm.ge.tol) then
            testpass=.false.
            msg = trim(msg)//" Scalar Chebyshev interpolant test fail."
        endif
        if (ddiffnorm.ge.dtol) then
            testpass=.false.
            msg = trim(msg)//new_line("a")//" Scalar Chebyshev interpolant derivative test fail."
        endif
        if (testpass) msg = "Scalar Chebyshev interpolation tests passed."
        write (*,*) trim(msg)
        contains
            elemental function fitfun(x) result(res)
                real(dp), intent(in) :: x
                real(dp)             :: res
                res = sin(x)*cos(3._dp*x) + sin(5._dp*x)*cos(7._dp*x)
            end function
            elemental function dfitfun(x) result(res)
                real(dp), intent(in) :: x
                real(dp)             :: res
                res = cos(x)*cos(3._dp*x) &
                    - 3._dp*sin(x)*sin(3._dp*x) &
                    + 5._dp*cos(5._dp*x)*cos(7._dp*x) & 
                    - 7._dp*sin(5._dp*x)*sin(7._dp*x)
            end function
    end subroutine test_scalar_cheby

    subroutine test_spice_subset(metakernel_filepath,testpass)
        logical, intent(inout) :: testpass

        character(len=*), intent(in)  :: metakernel_filepath
        character(len=500)  :: msg
        integer, parameter  :: tlen=200, indvar=1
        type(spice_subset)  :: spiceb, spice
        integer, parameter  :: deg=20
        real(dp), parameter :: a=772254069.1843683_dp, & ! 2024 Jun 21 15:00
                             & b=774846069.1835366_dp, & ! 2024 Jul 21 15:00
                             & tol = 1.e-12_dp, &
                             & dtol = 1.e-10_dp
        real(dp)            :: lt_dum, &
                               spkgeo_out(6)
        real(dp)            :: testpoints(tlen), &
                               truth(tlen,6,3), dtruth(tlen,6,3), &
                               cheb(tlen,6,3), dcheb(tlen,6,3), &
                               pnorm, vnorm
        integer             :: i, j, bodlist(6), stat, num

        bodlist = [10, 2, 301, 4, 5, 6]
        msg = ""
        call spiceb%init(trim(adjustl(metakernel_filepath)), 399, bodlist, a, b, deg)

        write (*,*) "Writing to temporary subspice file. . ."
        open(file="./test_resamp.subspice",newunit=num,access="stream",status="replace",iostat=stat)
        if (stat==0) then
            write (*,*) "Done"
        else
            testpass = .false.
            msg = trim(msg)//new_line("a")//" Failed to write subspice file."
        end if 
        call spiceb%write(num)
        close(num)
        open(file="./test_resamp.subspice",newunit=num,access="stream",status="old",iostat=stat)
        call spice%read(num)
        close(num)
        write (*,*) "Reading from temporary subspice file. . ."
        if (stat==0) then
            write (*,*) "Done"
        else
            testpass = .false.
            msg = trim(msg)//new_line("a")//" Failed to write subspice file."
        end if 

        testpoints = [(a + i*(b-a)/tlen, i=1,tlen)]
        do i=1,size(testpoints)
            do j=1,size(bodlist)
                call spkgeo(bodlist(j), testpoints(i), "J2000", 399, &
                                  & spkgeo_out, lt_dum)
                truth(i,j,1) = spkgeo_out(1)
                truth(i,j,2) = spkgeo_out(2)
                truth(i,j,3) = spkgeo_out(3)
                dtruth(i,j,1) = spkgeo_out(4)
                dtruth(i,j,2) = spkgeo_out(5)
                dtruth(i,j,3) = spkgeo_out(6)
            end do
        end do
        do j=1,size(bodlist)
            cheb(:,j,:) = spice%call(testpoints,bodlist(j),"p")
            dcheb(:,j,:) = spice%call(testpoints,bodlist(j),"v")
        end do
        pnorm = norm2(truth-cheb)/(norm2(truth)*tlen)
        vnorm = norm2(dtruth-dcheb)/(norm2(dtruth)*tlen)
        if (pnorm.ge.tol) then
            testpass=.false.
            msg = trim(msg)//" SPICE resample position test fail."
        endif
        if (vnorm.ge.dtol) then
            testpass=.false.
            msg = trim(msg)//new_line("a")//" SPICE resample velocity test fail."
        endif
        write(*,*) "Removing temporary subspice file. . ."
        open(file="./test_resamp.subspice",newunit=num,access="stream",status="old", iostat=stat)
        if (stat == 0) then
            close(num, status='delete')
            write (*,*) "Done."
        endif
        if (testpass) msg = "Spice resample tests passed."
        write (*,*) trim(msg)
    end subroutine test_spice_subset

end module test_cheby
