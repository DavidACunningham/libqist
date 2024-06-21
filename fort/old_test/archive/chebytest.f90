program main
    use cheby
    implicit none
    integer, parameter  :: deg=30, tlen=1000
    real(dp), parameter :: a=0._dp, b=2._dp*pi
    real(dp)            :: nodes(deg), fnodes(deg), testpoints(tlen), &
                           truth(tlen), dtruth(tlen), ddtruth(tlen), &
                           cheb(tlen), dcheb(tlen), ddcheb(tlen), &
                           coeffs(deg), dcoeffs(deg-1),ddcoeffs(deg-2)
    integer i

    nodes = chnodes(deg, a, b)
    testpoints = [(a + i*b/tlen, i=1,tlen)]
    fnodes = fitfun(nodes)
    truth = fitfun(testpoints)
    dtruth = dfitfun(testpoints)
    ddtruth = ddfitfun(testpoints)
    coeffs = chfit(deg,fnodes)
    dcoeffs = chderiv(coeffs,a,b)
    ddcoeffs = chderiv(dcoeffs,a,b)
    cheb = chcall(a,b,coeffs,testpoints)
    dcheb = chcall(a,b,dcoeffs,testpoints)
    ddcheb = chcall(a,b,ddcoeffs,testpoints)
    call print_to_file("truth", truth)
    call print_to_file("dtruth", dtruth)
    call print_to_file("ddtruth", ddtruth)
    call print_to_file("coeffs", coeffs)
    call print_to_file("dcoeffs", dcoeffs)
    call print_to_file("cheb", cheb)
    call print_to_file("dcheb", dcheb)
    call print_to_file("ddcheb", ddcheb)
    call print_to_file("nodes", nodes)
    call print_to_file("testpoints", testpoints)
    contains
        elemental function fitfun(x) result(res)
            real(dp), intent(in) :: x
            real(dp)             :: res
            ! res = x**3 + 2.*x**2 - 4.*x - 2.
            res = sin(x)
        end function
        elemental function dfitfun(x) result(res)
            real(dp), intent(in) :: x
            real(dp)             :: res
            res = cos(x)
        end function
        elemental function ddfitfun(x) result(res)
            real(dp), intent(in) :: x
            real(dp)             :: res
            res = -sin(x)
        end function
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
