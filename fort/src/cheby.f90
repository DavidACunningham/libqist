module cheby
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)
    interface chcall
        pure function chcall_s(a,b,coeffs,x)
            import
            real(dp), intent(in) :: a, b, x, coeffs(:)
            real(dp) :: chcall_s
        end function chcall_s
        pure function chcall_v(a,b,coeffs,x)
            import
            real(dp), intent(in) :: a, b, x(:), coeffs(:)
            real(dp) :: chcall_v(size(x))
        end function chcall_v
    end interface chcall
    interface chnodes
        pure function chnodes(n,a,b) result(res)
            import
            integer, intent(in) :: n
            real(dp), intent(in) :: a, b
            real(dp)             :: res(n)
        end function
    end interface chnodes
    interface chcall_mult
        pure function chcall_mult(a,b,coeffs,x) result(res)
            import
            real(dp), intent(in) :: a, b, x, coeffs(:)
            real(dp) :: res(size(coeffs))
        end function chcall_mult
    end interface chcall_mult
    interface chfit
        pure function chfit_s(n,fi) result(res)
            import
            integer, intent(in) :: n
            real(dp), intent(in) :: fi(n)
            real(dp) :: res(n)
        end function chfit_s
    end interface chfit
    interface chderiv
        pure function chderiv_s(coeffs,a,b) result(res)
            import
            real(dp), intent(in) :: coeffs(:), a, b
            real(dp) :: res(size(coeffs)-1)
        end function chderiv_s
    end interface chderiv

    type vectorcheb
        integer               :: ndim, ndeg
        real(dp), allocatable :: coeffs(:,:) !(ndeg, ndim)
        real(dp)              :: a,b
        contains
            procedure fit
            generic, public :: call => call_s, call_v
            procedure, private :: call_s
            procedure, private :: call_v
            procedure deriv
    end type vectorcheb
    contains
        subroutine fit(me, fi, a, b)
            class(vectorcheb), intent(inout) :: me
            real(dp),          intent(in)    :: fi(:,:), &
                                                ! (n x ndim)
                                              & a, b
            integer i
            me%ndim = size(fi,2)
            me%ndeg = size(fi,1)
            me%a    = a
            me%b    = b
            if (allocated(me%coeffs)) deallocate(me%coeffs)
            allocate(me%coeffs(me%ndeg,me%ndim))
            do i=1,me%ndim
                me%coeffs(:,i) = chfit(me%ndeg, fi(:,i))
            end do
        end subroutine fit
        function call_s(me, x) result(res)
            class(vectorcheb), intent(in) :: me
            real(dp),          intent(in)    :: x
            real(dp)                         :: res(me%ndim)
            integer                          :: i
            do i=1,me%ndim
                res(i) = chcall(me%a, me%b, me%coeffs(:,i), x)
            enddo
        end function call_s
        function call_v(me, x) result(res)
            class(vectorcheb), intent(in) :: me
            real(dp),          intent(in)    :: x(:)
            real(dp)                         :: res(size(x),me%ndim)
            integer                          :: i
            do i=1,me%ndim
                res(:,i) = chcall(me%a, me%b, me%coeffs(:,i), x)
            enddo
        end function call_v
        function deriv(me) result(res)
            class(vectorcheb), intent(in) :: me
            type(vectorcheb)             :: res
            integer                          :: i
            res%ndim = me%ndim
            res%ndeg = me%ndeg-1
            res%a    = me%a
            res%b    = me%b
            if (allocated(res%coeffs)) deallocate(res%coeffs)
            allocate(res%coeffs(res%ndeg,res%ndim))
            do i=1,res%ndim
                res%coeffs(:,i) = chderiv(me%coeffs(:,i), me%a, me%b)
            enddo
        end function deriv
end module cheby
pure function chnodes(n,a,b) result(res)
    use cheby, only: pi, dp
    integer, intent(in) :: n
    integer :: i
    real(dp), intent(in) :: a, b
    real(dp) :: bma, bpa
    real(dp) :: res(n)
    ! n: Order of approximation
    ! a: beginning of range
    ! b: end of range
    ! nodearray: n-dimensional array of chebyshev nodes

    ! return the list of chebyshev nodes at which to evaluate 
    ! a function for chebyshev approximation
    bma=0.5_dp*(b-a)
    bpa = 0.5_dp*(b+a)
    res = [(cos(pi*(i-0.5_dp)/n)*bma+bpa, i=n,1,-1)]
end function chnodes

pure function chcall_v(a,b,coeffs,x)
    use cheby, only: chnodes, pi, dp
    implicit none
    integer :: n,j
    real(dp), intent(in) :: a, b, x(:), coeffs(:)
    real(dp) :: u(size(x)), Tjm1(size(x)), Tj(size(x)), Tjp1(size(x))
    real(dp) :: chcall_v(size(x))
    ! a: beginning of range
    ! b: end of range
    ! val: value(s) at which to evaluate function (independent variable)
    ! coeffs: the list of chebyshev coefficients for the interpolant

    ! return the value(s) of a chebyshev polynomial (defined by coeffs)
    ! at a value val

    n=size(coeffs)
    u = (2.0_dp*x-a-b)/(b-a)
    Tjm1 = 1.0_dp
    Tj = u
    chcall_v = Tjm1*coeffs(1)
    do j=2,n
        Tjp1 = 2.0_dp*u*Tj - Tjm1
        chcall_v = chcall_v + Tj*coeffs(j)
        Tjm1 = Tj
        Tj = Tjp1
    end do
end function chcall_v

pure function chcall_s(a,b,coeffs,x)
    use cheby, only: chnodes, pi, dp
    implicit none
    integer :: n,j
    real(dp), intent(in) :: a, b, x, coeffs(:)
    real(dp) :: u, Tjm1, Tj, Tjp1
    real(dp) :: chcall_s
    ! a: beginning of range
    ! b: end of range
    ! val: value(s) at which to evaluate function (independent variable)
    ! coeffs: the list of chebyshev coefficients for the interpolant

    ! return the value(s) of a chebyshev polynomial (defined by coeffs)
    ! at a value val

    n=size(coeffs)
    u = (2.0_dp*x-a-b)/(b-a)
    Tjm1 = 1.0_dp
    Tj = u
    chcall_s = Tjm1*coeffs(1)
    do j=2,n
        Tjp1 = 2.0_dp*u*Tj - Tjm1
        chcall_s = chcall_s + Tj*coeffs(j)
        Tjm1 = Tj
        Tj = Tjp1
    end do
end function chcall_s

pure function chcall_mult(a,b,coeffs,x) result(res)
    use cheby, only: chnodes, pi, dp
    implicit none
    integer :: n,j
    real(dp), intent(in) :: a, b, x, coeffs(:)
    real(dp) :: u, Tjm1, Tj, Tjp1
    real(dp) :: res(size(coeffs))
    ! a: beginning of range
    ! b: end of range
    ! val: value(s) at which to evaluate function (independent variable)
    ! coeffs: the list of chebyshev coefficients for the interpolant

    ! return the value of a chebyshev polynomial (defined by coeffs)
    ! at a value val

    n=size(coeffs)
    u = (2.0_dp*x-a-b)/(b-a)
    Tjm1 = 1.0_dp
    Tj = u
    res(1) = Tjm1*coeffs(1)
    do j=2,n
        Tjp1 = 2.0_dp*u*Tj - Tjm1
        res(j) = res(j-1) + Tj*coeffs(j)
        Tjm1 = Tj
        Tj = Tjp1
    end do
end function chcall_mult

pure function chfit_s(n,fi) result(res)
    use cheby, only: chnodes, pi, dp
    implicit none
    integer, intent(in) :: n
    integer :: j,i
    real(dp), intent(in) :: fi(n)
    real(dp) :: res(n), u(n), polyarray(n,n)
    ! n: Order of approximation
    ! a: beginning of range
    ! b: end of range
    ! fi: n-dimensional array of function evaluated at chebyshev nodes 
    !     (see above)

    ! return the list of chebyshev coefficients approximating the function
    ! encoded by the values of the xi
    u = [(cos((i-0.5_dp)*pi/n), i=n,1,-1)]
    if (n>0) polyarray(:,1) = 1.0_dp
    if (n>1) polyarray(:,2) = u
    if (n>2) then
        do j=3,n 
            polyarray(:,j) = [(2.0_dp*u(i)*polyarray(i,j-1)-polyarray(i,j-2), i=1,n)]
        end do
    end if
    if (n>1) polyarray(:,2:) = 2.0_dp*polyarray(:,2:)
    res = matmul(fi,polyarray)/n
end function chfit_s

pure function chderiv_s(coeffs,a, b) result(res)
    ! Source: "Modern Computing Methods", Goodwin, E. T., ed., 1961, p. 79
    use cheby, only: dp
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    real(dp), intent(in) :: coeffs(:), a, b
    real(dp)             :: res(size(coeffs)-1)
    real(dp)             :: primecoeffs(size(coeffs)-1)
    real(qp)             :: quadp(size(coeffs)+1)
    integer              :: deg, i
    deg = size(coeffs)
    quadp = 0._qp
    do i=deg,2,-1
        quadp(i-1) = quadp(i+1) + 2._qp*(i-1)*real(coeffs(i),qp)
    end do
    res = real(quadp(:deg-1)*(2._qp/(b-a)),dp)
    res(1) = res(1)/2.
    ! primecoeffs = 0._dp
    ! do i=deg,2,-1
    !     primecoeffs(i-1) = primecoeffs(i+1) + 2._dp*(i-1)*coeffs(i)
    ! end do
    ! res = primecoeffs(:deg-1)*(2._dp/(b-a))
    ! res(1) = res(1)/2.
end function chderiv_s
