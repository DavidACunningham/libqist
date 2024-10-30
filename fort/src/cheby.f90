! Title: cheby.f90 
! Description:
!     types, functions, and subroutines implementing Chebyshev
!     interpolation of functions of one variable.
!     
!
! References:
!     Goodwin, E.T., ed., Modern Computing Methods, Her Majesty's Stationery
!         Office, 1961, p. 79ff
!     Press, W.H., Teukolsky, S.A., Flannery, B.P., Vetterling, W.T.: Numerical
!         Recipes in Fortran 77: The Art of Scientific Computing, vol. 1. 
!         Cambridge University Press, Cambridge, 1992
! 
! author: David Cunningham
! Last edited: See git log
module cheby
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)
    ! See functions below for documentation of these interfaces
    interface chcall
        function chcall_s(a,b,coeffs,x)
            import
            real(dp), intent(in) :: a, b, x, coeffs(:)
            real(dp) :: chcall_s
        end function chcall_s
        function chcall_v(a,b,coeffs,x)
            import
            real(dp), intent(in) :: a, b, x(:), coeffs(:)
            real(dp) :: chcall_v(size(x))
        end function chcall_v
    end interface chcall
    interface chnodes
        function chnodes(n,a,b) result(res)
           import
           integer, intent(in) :: n
           real(dp), intent(in) :: a, b
           real(dp)             :: res(n)
        end function
    end interface chnodes
    interface chfit
        function chfit_s(n,fi) result(res)
           import
           integer, intent(in) :: n
           real(dp), intent(in) :: fi(n)
           real(dp) :: res(n)
        end function chfit_s
    end interface chfit
    interface chderiv
        function chderiv_s(coeffs,a,b) result(res)
           import
           real(dp), intent(in) :: coeffs(:), a, b
           real(dp) :: res(size(coeffs)-1)
        end function chderiv_s
    end interface chderiv
    type vectorcheb
        ! vectorcheb: type for computing, manipulating, and storing
        !    Chebyshev interpolants of functions g: R -> R^n
        integer               :: ndim, ndeg
        real(dp), allocatable :: coeffs(:,:) !Shape will be (ndeg, ndim)
        real(dp)              :: a,b ! Beginning and end of interpolant domain
        contains
            procedure fit
            generic, public :: call => call_s, call_v
            procedure, private :: call_s
            procedure, private :: call_v
            procedure deriv
            procedure, public  :: write => vectorchebwrite
            procedure, public  :: read => vectorchebread
    end type vectorcheb
    contains
        subroutine vectorchebwrite(me,unit_num)
            ! vectorchebwrite:
            ! Write the calling instance of the type to disk
            ! NOTE: as a type bound procedure, the signature of this function
            ! will be:
            !   CALL me%write(unit_num)
            ! INPUTS:
            ! NAME     TYPE        DESCRIPTION
            ! me       vectorcheb  Calling instance of the type
            ! unitnum  integer     unit number/file handle of open
            !                      file in which type will be written.
            !                      file must be opened with the ``ACCESS''
            !                      flag set to ``stream'' and the ``STATUS''
            !                      flag set to ``new'' or ``replace''
            class(vectorcheb), intent(inout) :: me
            integer,             intent(in) :: unit_num
            write(unit_num) me%ndim
            write(unit_num) me%ndeg
            write(unit_num) me%a
            write(unit_num) me%b
            write(unit_num) me%coeffs
        end subroutine vectorchebwrite
        subroutine vectorchebread(me,unit_num)
            ! vectorchebread:
            ! Read a stored vectorcheb type instance from disk into
            ! the calling instance
            ! NOTE: as a type bound procedure, the signature of this function
            ! will be:
            !   CALL me%read(unit_num)
            ! NOTE: This function modifies the state of the calling instance
            ! INPUTS:
            ! NAME     TYPE        DESCRIPTION
            ! me       vectorcheb  Calling instance of the type
            ! unitnum  integer     unit number/file handle of open
            !                      file in which type will be written.
            !                      file must be opened with the ``ACCESS''
            !                      flag set to ``stream'' and the ``STATUS''
            !                      flag set to ``old''
            class(vectorcheb), intent(inout) :: me
            integer,             intent(in) :: unit_num
            read(unit_num) me%ndim
            read(unit_num) me%ndeg
            read(unit_num) me%a
            read(unit_num) me%b
            allocate(me%coeffs(me%ndeg,me%ndim))
            read(unit_num) me%coeffs
        end subroutine vectorchebread
        subroutine fit(me, fi, a, b)
            ! fit:
            ! compute and store a set of Chebyshev coefficients in the calling
            ! instance
            ! NOTE: as a type bound procedure, the signature of this function
            ! will be:
            !   CALL me%fit(fi, a, b)
            ! NOTE: This function modifies the state of the calling instance
            ! INPUTS:
            ! NAME     TYPE           DESCRIPTION
            ! me       vectorcheb     Calling instance of the type
            ! unitnum  integer        unit number/file handle of open
            !                         file in which type will be written.
            !                         file must be opened with the ``ACCESS''
            !                         flag set to ``stream'' and the ``STATUS''
            !                         flag set to ``old''
            ! fi       real (n,ndim)  Values of function to interpolate at 
            !                         Chebyshev nodes of order n 
            !                         (i.e. evaluated at the output of 
            !                         chnodes(n,a,b))
            ! a        real           beginning of interpolant domain
            ! b        real           end of interpolant domain
            class(vectorcheb), intent(inout) :: me
            real(dp),          intent(in)    :: fi(:,:), a, b
            integer i
            me%ndeg  = size(fi,1)
            me%ndim  = size(fi,2)
            me%a     = a
            me%b     = b
            if (allocated(me%coeffs)) deallocate(me%coeffs)
            allocate(me%coeffs(me%ndeg, me%ndim))
            do i=1, me%ndim
                me%coeffs(:,i) = chfit(me%ndeg,fi(:,i))
            end do
        end subroutine fit
        function call_s(me, x) result(res)
            ! call_s:
            ! Return the value of a Chebyshev polynomial at a value specified by x.
            ! This function returns a _single_ value of an interpolant
            ! NOTE: as a type bound procedure, the signature of this function
            ! will be:
            !   res = me%call(x)
            ! INPUTS:
            ! NAME     TYPE        DESCRIPTION
            ! me       vectorcheb  Calling instance of the type
            ! x        real        value at which to evaluate 
            !                      interpolant
            ! OUTPUTS:
            ! NAME     TYPE   DESCRIPTION
            ! res      real   value of interpolant at x
            class(vectorcheb), intent(in) :: me
            real(dp),          intent(in)    :: x
            real(dp)                         :: res(me%ndim)
            integer                          :: i
            do i=1,me%ndim
                res(i) = chcall(me%a, me%b, me%coeffs(:,i), x)
            enddo
        end function call_s
        function call_v(me, x) result(res)
            ! call_v:
            ! Return the values of a Chebyshev polynomial at values specified by x.
            ! This function returns _multiple_ value of an interpolant
            ! NOTE: as a type bound procedure, the signature of this function
            ! will be:
            !   res = me%call(x)
            ! INPUTS:
            ! NAME     TYPE            DESCRIPTION
            ! me       vectorcheb      Calling instance of the type
            ! x        real (:)        values at which to evaluate 
            !                          interpolant
            ! OUTPUTS:
            ! NAME     TYPE            DESCRIPTION
            ! res      real (size(x))  values of interpolant at elements of x
            class(vectorcheb), intent(in) :: me
            real(dp),          intent(in)    :: x(:)
            real(dp)                         :: res(size(x),me%ndim)
            integer                          :: i
            do i=1,me%ndim
                res(:,i) = chcall(me%a, me%b, me%coeffs(:,i), x)
            enddo
        end function call_v
        function deriv(me) result(res)
            ! deriv:
            ! Return the list of Chebyshev coefficients that interpolate a function
            ! evaluated at the Chebyshev nodes of a given degree
            ! NOTE: as a type bound procedure, the signature of this function
            ! will be:
            !   res = me%deriv()
            ! INPUTS:
            ! chderiv:
            ! INPUTS:
            ! NAME     TYPE        DESCRIPTION
            ! me       vectorcheb  Calling instance of the type
            ! OUTPUTS:
            ! NAME     TYPE        DESCRIPTION
            ! res      vectorcheb  a new instance of vectorcheb
            !                      storing the elementwise derivative of
            !                      ``me''
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
function chnodes(n,a,b) result(res)
    ! chnodes:
    ! return the list of Chebyshev nodes at which to evaluate
    ! a function for Chebyshev approximation
    ! INPUTS:
    ! NAME     TYPE     DESCRIPTION
    ! n        integer  Degree of interpolating polynomial
    ! a        real     beginning of interpolant domain
    ! b        real     end of interpolant domain
    ! OUTPUTS:
    ! NAME     TYPE     DESCRIPTION
    ! res      real (n) n-dimensional array of Chebyshev nodes
    use cheby, only: pi, dp
    integer, intent(in) :: n
    integer :: i
    real(dp), intent(in) :: a, b
    real(dp) :: bma, bpa
    real(dp) :: res(n)
    bma=0.5_dp*(b-a)
    bpa = 0.5_dp*(b+a)
    res = [(cos(pi*(i-0.5_dp)/n)*bma+bpa, i=n,1,-1)]
end function chnodes

function chcall_v(a,b,coeffs,x)
    ! chcall_v:
    ! Return the values of a Chebyshev polynomial at values specified by x. 
    ! This function returns _multiple_ values of an interpolant
    ! INPUTS:
    ! NAME      TYPE           DESCRIPTION
    ! a         real           beginning of interpolant domain
    ! b         real           end of interpolant domain
    ! coeffs    real (:)       list of Chebyshev coefficients for the 
    !                          interpolant
    ! x         real (:)       list of values at which to evaluate
    !                          interpolant
    ! OUTPUTS:
    ! NAME      TYPE           DESCRIPTION
    ! chcall_v  real (size(x)) values of interpolant at each element of x      
    use cheby, only: chnodes, pi, dp
    implicit none
    integer :: n,j
    real(dp), intent(in) :: a, b, x(:), coeffs(:)
    real(dp) :: u(size(x)), Tjm1(size(x)), Tj(size(x)), Tjp1(size(x))
    real(dp) :: chcall_v(size(x))
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

function chcall_s(a,b,coeffs,x)
    ! chcall_s:
    ! Return the value of a Chebyshev polynomial at a value specified by x. 
    ! This function returns a _single_ value of an interpolant
    ! INPUTS:
    ! NAME      TYPE   DESCRIPTION
    ! a         real   beginning of interpolant domain
    ! b         real   end of interpolant domain
    ! coeffs    real   list of Chebyshev coefficients for the interpolant
    ! x         real   value at which to evaluate interpolant
    ! OUTPUTS:
    ! NAME      TYPE   DESCRIPTION
    ! chcall_s  real   value of interpolant at x
    use cheby, only: chnodes, pi, dp
    implicit none
    integer :: n,j
    real(dp), intent(in) :: a, b, x, coeffs(:)
    real(dp) :: u, Tjm1, Tj, Tjp1
    real(dp) :: chcall_s
    j = 0
    n=size(coeffs)
    u = (2.0_dp*x-a-b)/(b-a)
    Tjm1 = 1.0_dp
    Tj = u
    chcall_s = Tjm1*coeffs(1)
    if (size(coeffs).lt.j) then
        print *, j
        print * , size(coeffs)
        error stop
    end if
    do j=2,n
        Tjp1 = 2.0_dp*u*Tj - Tjm1
        chcall_s = chcall_s + Tj*coeffs(j)
        Tjm1 = Tj
        Tj = Tjp1
    end do
end function chcall_s

function chfit_s(n,fi) result(res)
    ! chfit_s:
    ! return the list of Chebyshev coefficients that interpolate a function
    ! evaluated at the Chebyshev nodes of a given degree
    ! INPUTS:
    ! NAME     TYPE     DESCRIPTION
    ! n        integer  Degree of interpolating polynomial
    ! fi       real (n) Values of function to interpolate at Chebyshev nodes
    !                   of order n (i.e. evaluated at the output of 
    !                   chnodes(n,a,b))
    ! OUTPUTS:
    ! NAME     TYPE     DESCRIPTION
    ! res      real (n) list of Chebyshev coefficients approximating the
    !                   function used to build fi.
    use cheby, only: chnodes, pi, dp
    implicit none
    integer, intent(in) :: n
    integer :: j,i
    real(dp), intent(in) :: fi(n)
    real(dp) :: res(n), u(n), polyarray(n,n)
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

function chderiv_s(coeffs,a, b) result(res)
    ! chderiv_s:
    ! return the list of Chebyshev coefficients that interpolate a function
    ! evaluated at the Chebyshev nodes of a given degree
    ! INPUTS:
    ! NAME     TYPE     DESCRIPTION
    ! a        real     beginning of interpolant domain
    ! b        real     end of interpolant domain
    ! coeffs   real (:) list of Chebyshev coefficients of interpolant
    !                   
    ! OUTPUTS:
    ! NAME     TYPE                   DESCRIPTION
    ! res      real (size(coeffs)-1)  list of Chebyshev coefficients 
    !                                 approximating the derivative of the
    !                                 interpolant specified by the input
    !                                 coefficients
    use cheby, only: dp
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    real(dp), intent(in) :: coeffs(:), a, b
    real(dp)             :: res(size(coeffs)-1)
    real(qp)             :: quadp(size(coeffs)+1)
    integer              :: deg, i
    deg = size(coeffs)
    quadp = 0._qp
    do i=deg,2,-1
        quadp(i-1) = quadp(i+1) + 2._qp*(i-1)*real(coeffs(i),qp)
    end do
    res = real(quadp(:deg-1)*(2._qp/(b-a)),dp)
    res(1) = res(1)/2.
end function chderiv_s
