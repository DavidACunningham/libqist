module cheby
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)
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
    interface chcall_mult
         function chcall_mult(a,b,coeffs,x) result(res)
            import
            real(dp), intent(in) :: a, b, x, coeffs(:)
            real(dp) :: res(size(coeffs))
        end function chcall_mult
    end interface chcall_mult
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
    type spice_subset
        integer               :: nbods, ndeg, central_body
        integer, allocatable  :: bodlist(:)
        real(dp)              :: a, b
        real(dp), allocatable :: pcoeffs(:,:), & !(ndeg,3*nbods)
                               & vcoeffs(:,:), & !(ndeg,3*nbods)
                               & acoeffs(:,:)    !(ndeg,3*nbods)
        contains
        procedure          :: init => fitspice
        procedure          :: write => spicewrite
        procedure          :: read => spiceread
        procedure, private :: call_spice_one
        procedure, private :: call_spice_sev
        generic, public    :: call => call_spice_one, call_spice_sev
    end type spice_subset
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
            procedure, public  :: write => vectorchebwrite
            procedure, public  :: read => vectorchebread
    end type vectorcheb
    contains
        subroutine vectorchebwrite(me,unit_num)
            class(vectorcheb), intent(inout) :: me
            integer,             intent(in) :: unit_num
            write(unit_num) me%ndim
            write(unit_num) me%ndeg
            write(unit_num) me%a
            write(unit_num) me%b
            write(unit_num) me%coeffs
        end subroutine vectorchebwrite
        subroutine vectorchebread(me,unit_num)
            class(vectorcheb), intent(inout) :: me
            integer,             intent(in) :: unit_num
            read(unit_num) me%ndim
            read(unit_num) me%ndeg
            read(unit_num) me%a
            read(unit_num) me%b
            allocate(me%coeffs(me%ndeg,me%ndim))
            read(unit_num) me%coeffs
        end subroutine vectorchebread
        subroutine spicewrite(me,unit_num)
            class(spice_subset), intent(in) :: me
            integer,             intent(in) :: unit_num
            write(unit_num) me%nbods
            write(unit_num) me%ndeg
            write(unit_num) me%central_body
            write(unit_num) me%bodlist
            write(unit_num) me%a
            write(unit_num) me%b
            write(unit_num) me%pcoeffs
            write(unit_num) me%vcoeffs
            write(unit_num) me%acoeffs
        end subroutine spicewrite
        subroutine spiceread(me,unit_num)
            class(spice_subset), intent(inout) :: me
            integer,             intent(in) :: unit_num
            read(unit_num) me%nbods
            read(unit_num) me%ndeg
            read(unit_num) me%central_body
            allocate(me%bodlist(me%nbods))
            read(unit_num) me%bodlist
            read(unit_num) me%a
            read(unit_num) me%b
            allocate(me%pcoeffs(me%ndeg,3*me%nbods), &
                   & me%vcoeffs(me%ndeg,3*me%nbods), &
                   & me%acoeffs(me%ndeg,3*me%nbods))
            read(unit_num) me%pcoeffs
            read(unit_num) me%vcoeffs
            read(unit_num) me%acoeffs
            me%acoeffs(me%ndeg-1:,:) = 0._dp
        end subroutine spiceread
        subroutine fitspice(me, kernelfile, central_body, bodlist, a, b, ndeg)
            class(spice_subset), intent(inout) :: me
            character(len=*),  intent(in)    :: kernelfile
            integer,             intent(in)    :: central_body, &
                                                & bodlist(:), &
                                                & ndeg
            real(dp),            intent(in)    :: a,b
            real(dp)                           :: nodes(ndeg), &
                                                & pfi(ndeg,3*(size(bodlist)+1)), &
                                                & vfi(ndeg,3*(size(bodlist)+1)), &
                                                & spkgeo_out(6), lt_dum
            integer i,j
            ! load kernel
            call FURNSH(trim(adjustl(kernelfile)))
            if (allocated(me%pcoeffs)) deallocate(me%pcoeffs)
            if (allocated(me%vcoeffs)) deallocate(me%vcoeffs)
            if (allocated(me%acoeffs)) deallocate(me%acoeffs)
            if (allocated(me%bodlist)) deallocate(me%bodlist)
            ! assign type properties
            me%nbods = size(bodlist) + 1
            me%ndeg  = ndeg
            me%a     = a
            me%b     = b
            me%central_body = central_body
            ! allocate bodlist, p/v/acoeffs
            allocate(me%pcoeffs(ndeg,3*me%nbods), &
                     me%vcoeffs(ndeg,3*me%nbods), &
                     me%acoeffs(ndeg,3*me%nbods), &
                     me%bodlist(me%nbods))
            ! assign bodies
            me%bodlist(1) = central_body
            me%bodlist(2:) = bodlist
            ! get chebyshev nodes
            nodes = chnodes(me%ndeg, me%a, me%b)
            !  for each node, do:
            do j = 1,me%ndeg
                ! SPKGEO call at time node(i) for pos/vel of central_body 
                ! with observer solar_system_barycenter
                call spkgeo(me%bodlist(1), nodes(j), "J2000", 0, &
                                  & spkgeo_out, lt_dum)
                pfi(j,:3) = spkgeo_out(:3)
                vfi(j,:3) = spkgeo_out(4:)
            end do
            ! for each body in bodlist, do:
            do i=2,me%nbods
            !  for each node, do:
                do j=1,me%ndeg
                    ! SPKGEO call at time node(i) for pos/vel of body(i) with
                    ! observer central_body
                    call spkgeo(me%bodlist(i), nodes(j), "J2000", me%bodlist(1), &
                                      & spkgeo_out, lt_dum)
                    pfi(j,3*(i-1)+1:3*i) = spkgeo_out(:3)
                    vfi(j,3*(i-1)+1:3*i) = spkgeo_out(4:)
                end do
            end do
            ! Fit all coeffs
            do i=1, 3*me%nbods
                me%pcoeffs(:,i) = chfit(me%ndeg,pfi(:,i))
                me%vcoeffs(:,i) = chfit(me%ndeg,vfi(:,i))
                me%acoeffs(:,i) = chderiv(me%vcoeffs(:,i),me%a,me%b)
            end do
        end subroutine fitspice
        function call_spice_one(me,x, bod_id,mode) result(res)
            class(spice_subset), intent(in) :: me
            real(dp),            intent(in) :: x
            integer,             intent(in) :: bod_id
            character(len=1),    intent(in) :: mode
            real(dp)                        :: res(3)
            integer                         :: i
            i = findloc(me%bodlist,bod_id,dim=1) - 1
            select case(mode)
            case ("p")
                res(1) = chcall(me%a,me%b,me%pcoeffs(:,3*i+1),x)
                res(2) = chcall(me%a,me%b,me%pcoeffs(:,3*i+2),x)
                res(3) = chcall(me%a,me%b,me%pcoeffs(:,3*i+3),x)
            case ("v")
                res(1) = chcall(me%a,me%b,me%vcoeffs(:,3*i+1),x)
                res(2) = chcall(me%a,me%b,me%vcoeffs(:,3*i+2),x)
                res(3) = chcall(me%a,me%b,me%vcoeffs(:,3*i+3),x)
            case ("a")
                res(1) = chcall(me%a,me%b,me%acoeffs(:,3*i+1),x)
                res(2) = chcall(me%a,me%b,me%acoeffs(:,3*i+2),x)
                res(3) = chcall(me%a,me%b,me%acoeffs(:,3*i+3),x)
            end select
        end function
        function call_spice_sev(me,x,bod_id,mode) result(res)
            class(spice_subset), intent(in) :: me
            real(dp),            intent(in) :: x(:)
            integer,             intent(in) :: bod_id
            character(len=1),    intent(in) :: mode
            real(dp)                        :: res(size(x),3)
            integer                         :: i
            i = findloc(me%bodlist,bod_id,dim=1) - 1
            select case(mode)
            case ("p")
                res(:,1) = chcall(me%a,me%b,me%pcoeffs(:,3*i+1),x)
                res(:,2) = chcall(me%a,me%b,me%pcoeffs(:,3*i+2),x)
                res(:,3) = chcall(me%a,me%b,me%pcoeffs(:,3*i+3),x)
            case ("v")
                res(:,1) = chcall(me%a,me%b,me%vcoeffs(:,3*i+1),x)
                res(:,2) = chcall(me%a,me%b,me%vcoeffs(:,3*i+2),x)
                res(:,3) = chcall(me%a,me%b,me%vcoeffs(:,3*i+3),x)
            case ("a")
                res(:,1) = chcall(me%a,me%b,me%acoeffs(:,3*i+1),x)
                res(:,2) = chcall(me%a,me%b,me%acoeffs(:,3*i+2),x)
                res(:,3) = chcall(me%a,me%b,me%acoeffs(:,3*i+3),x)
            end select
        end function call_spice_sev
        subroutine fit(me, fi, a, b)
            class(vectorcheb), intent(inout) :: me
            real(dp),          intent(in)    :: fi(:,:), &
                                                ! (n x ndim)
                                              & a, b
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
 function chnodes(n,a,b) result(res)
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

 function chcall_v(a,b,coeffs,x)
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

function chcall_s(a,b,coeffs,x)
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

 function chcall_mult(a,b,coeffs,x) result(res)
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

 function chfit_s(n,fi) result(res)
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

 function chderiv_s(coeffs,a, b) result(res)
    ! Source: "Modern Computing Methods", Goodwin, E. T., ed., 1961, p. 79
    ! This is done in quad based on advice from the source above
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
