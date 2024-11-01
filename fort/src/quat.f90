! Title: cheby.f90 
! Description:
!     implements two types: 
!           ``quaternion'': a quaternion type with algebraic multiplications defined
!           ``rothist'' : A type for storing continuous interpolated rotations
!                         using Chebyshev interpolants on the quaternion
!                         elements and appropriate renormalizations to compute 
!                         the derivatives of the quaternions and make sure they
!                         stay in T*SO(3).
!     Only double precision is currently implemented.
!
! References:
!   None
! 
! author: David Cunningham
! Last edited: See git log
module quat
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use cheby, only: vectorcheb, chnodes
    implicit none
    integer, parameter :: wp = dp
    real(wp), parameter :: pi = 4._wp*atan(1._wp)
    

    type :: quaternion
        ! A type implementing quaternion algebra
        ! And operations to convert the underlying unit quaternion
        ! to and from DCM and axis angle rotations
        ! The only attribute is the four-vector q which 
        ! stores the elements of the quaternion.
        real(dp) :: q(4)

        contains
        procedure norm
        procedure fromdcm
        procedure asdcm
        procedure axang
        procedure rotate_vec
        procedure :: renorm => quat_renorm
    end type quaternion
    type  rothist
        ! A type implementing the interpolation of rotations
        ! and the computation of time derivatives of those 
        ! rotations
        !
        ! Attributes:
        !  els: vectorcheb
        !   An array of four Chebyshev interpolants, interpolating
        !   the time history of the elements of the quaternion.
        !  elsdot: vectorcheb
        !   An array of four Chebyshev interpolants, interpolating
        !   the time history of the time derivative of the elements
        !   of the quaternion, computed in R^4, not constrained to be in T*SO(3)
        !  elsddot: vectorcheb
        !   An array of four Chebyshev interpolants, interpolating
        !   the time history of the second time derivative of the elements
        !   of the quaternion, computed in R^4, not constrained to be in T*SO(3)
        !  qstat: quaternion
        !   the static part of the rotation, i.e. we break the rotation down into two
        !   parts: J2000 -> BF(t_0) -> BF(t_f)
        !   so qstat only captures the first rotation. This is done since most celestial
        !   bodies are simple spinners.
        !  degree: integer
        !   the degree of the polynomial interpolation of els
        !  t0: real
        !   initial time
        !  tf: real
        !   final time
        !  fit: function handle
        !   see fit_func signature. The function that takes in times and returns a quaternion
        !   for the rotation BF(t_0) -> BF(t_f)
        type(vectorcheb)             :: els, elsdot, elsddot
        type(quaternion)             :: qstat
        integer                      :: degree
        real(dp)                     :: t0, tf
        procedure(fit_func), pointer :: fit => null()

        contains
            procedure :: read => rotread
            procedure :: write => rotwrite
            generic, public :: init => init_rot_q, init_rot_dcm
            procedure :: callq => getrotquat
            procedure :: callqdot => getrotquatdot
            procedure :: callqddot => getrotquatddot
            procedure :: call => getdcm
            procedure :: calldot => getdcmdot
            procedure :: callddot => getdcmddot
            procedure, private :: renorm
            procedure, private :: init_rot_q
            procedure, private :: init_rot_dcm
    end type rothist

    abstract interface 
        function fit_func(me, ta, tb) result(res)
           ! Function wrapping the interpolant target
           ! Inputs:
           ! me: The calling instance
           ! ta: real- The initial time for the rotation
           ! tb: real- The final time for the rotation
           ! Outputs:
           ! res: real (4) the rotation quaternion from
           !      ta to tb
            import :: wp, rothist
            implicit none
            class(rothist), intent(inout) :: me
            real(wp), intent(in) :: ta, tb
            real(wp)             :: res(4)
        end function fit_func
    end interface
    interface operator ( * )
        ! Binary operator overloading to enable
        ! quaternion multiplications
        ! by scalars and other quaternions.
        ! Quaternion composition is not a separate
        ! operator, so multiplication of one object
        ! of quaternion type by another object of quaternion
        ! type automatically performs quaternion composition.
        module procedure qmul
        module procedure smul1
        module procedure smul2
        module procedure imul1
        module procedure imul2
    end interface
    interface operator ( + )
        ! Binary operator overloading
        ! to enable seamless quaternion addition
        module procedure qadd
    end interface
    interface operator ( / )
        ! Binary operator overloading
        ! to enable quaternion division by a scalar
        module procedure sdiv
    end interface
    interface operator ( - )
        ! Unary/binary operator overloading
        ! to enable subtracting and computing
        ! additive inverses of quaternion objects
        module procedure qsub
        module procedure qainv
    end interface
    interface operator ( .T. )
        ! Unary operator to enable transpose
        ! of quaternions--more properly termed adjoint but
        ! .A. didn't seem as intuitive.
        module procedure qconj
    end interface


    contains
        subroutine rotwrite(me,unit_num)
            ! Write the calling rotation interpolation history
            ! to disk at unit unit_num
            class(rothist), intent(inout) :: me
            integer,        intent(in)    :: unit_num
            write(unit_num) me%t0
            write(unit_num) me%tf
            write(unit_num) me%degree
            write(unit_num) me%qstat%q
            call me%els%write(unit_num)
        end subroutine rotwrite
        subroutine rotread(me,unit_num)
            ! Read a rotation interpolation history 
            ! from disk at unit unit_num into the calling
            ! instance. Requires ``STREAM'' access.
            class(rothist), intent(inout) :: me
            integer,        intent(in)    :: unit_num
            read(unit_num) me%t0
            read(unit_num) me%tf
            read(unit_num) me%degree
            read(unit_num) me%qstat%q
            call me%els%read(unit_num)
            me%elsdot = me%els%deriv()
            me%elsddot = me%elsdot%deriv()
        end subroutine rotread
        subroutine init_rot_q(me, fit, t0, tf, order, qstat)
            ! Initializer for rotation history from a quaternion
            ! INPUTS:
            ! NAME:         TYPE:             DESCRIPTION:
            ! me            rothist           calling istance (STATE WILL CHANGE)
            ! fit           function pointer  function that generates rotation 
            !                                 quaternions
            ! t0            real              initial time for interpolant
            ! tf            real              final time for interpolant
            ! order         integer           degree of interpolating polynomial
            ! qstat         real (4)          4-d array containing unit quaternion
            !                                 of static part of rotation
            class(rothist),     intent(inout) :: me
            procedure(fit_func)               :: fit
            real(wp),           intent(in)    :: t0,tf, qstat(4)
            integer,            intent(in)    :: order
            real(wp)                          :: nodes(order), &
                                               & vals(order,4), &
                                               & first(4)
            integer i
            me%fit => fit
            me%degree = order
            me%t0 = t0
            me%tf = tf
            nodes = chnodes(order, t0, tf)
            first = me%fit(t0, t0+(nodes(1)-t0)/10._wp)
            vals(1,:) = fitwrap(t0, nodes(1), first)
            do i = 2, order
                vals(i,:) = fitwrap(t0,nodes(i), vals(i-1,:))
            end do
            call me%els%fit(vals, t0,tf)
            me%qstat%q = qstat
            me%elsdot = me%els%deriv()
            me%elsddot = me%elsdot%deriv()
                contains
                    function fitwrap(a, b, prev) result(res)
                        real(wp), intent(in) :: a, b, prev(4)
                        real(wp)             :: res(4), cur(4), dot
                        cur = me%fit(a, b)
                        dot = dot_product(prev,cur)
                        if (dot.lt.0._wp) then
                            res = -cur
                        else
                            res = cur
                        end if
                    end function fitwrap
        end subroutine
        subroutine init_rot_dcm(me, fit, t0, tf, order, dcmstat)
            ! Initializer for rotation history from a direction cosine matrix
            ! INPUTS:
            ! NAME:         TYPE:             DESCRIPTION:
            ! me            rothist           calling istance (STATE WILL CHANGE)
            ! fit           function pointer  function that generates rotation 
            !                                 quaternions
            ! t0            real              initial time for interpolant
            ! tf            real              final time for interpolant
            ! order         integer           degree of interpolating polynomial
            ! qstat         real (3,3)        array containing DCM for
            !                                 static part of rotation
            class(rothist),     intent(inout) :: me
            procedure(fit_func)               :: fit
            real(wp),           intent(in)    :: t0,tf, dcmstat(3,3)
            integer,            intent(in)    :: order
            real(wp)                          :: nodes(order), &
                                               & vals(order,4), &
                                               & first(4)
            integer i
            me%fit => fit
            me%degree = order
            me%t0 = t0
            me%tf = tf
            nodes = chnodes(order, t0, tf)
            first = me%fit(t0, t0+(nodes(1)-t0)/10._wp)
            vals(1,:) = fitwrap(nodes(1), first)
            do i = 2, order
                vals(i,:) = fitwrap(nodes(i), vals(i-1,:))
            end do
            call me%els%fit(vals, t0,tf)
            call me%qstat%fromdcm(dcmstat)
            me%elsdot = me%els%deriv()
            me%elsddot = me%elsdot%deriv()
                contains
                    function fitwrap(t, prev) result(res)
                        real(wp), intent(in) :: t, prev(4)
                        real(wp)             :: res(4), cur(4), dot
                        cur = me%fit(t0, t)
                        dot = dot_product(prev,cur)
                        if (dot.lt.0._wp) then
                            res = -cur
                        else
                            res = cur
                        end if
                    end function fitwrap
        end subroutine
        function getdcm(me, t) result(res)
            ! Return the DCM representation
            ! of the rotation at time t
            class(rothist), intent(inout) :: me
            type(quaternion)              :: qtot, qdyn
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(3,3)
            real(wp)                      :: q(4)
            q = me%els%call(t)
            call me%renorm(q)
            qdyn%q = q
            qtot = qdyn*me%qstat
            res = transpose(qtot%asdcm())
        end function
        function getrotquat(me, t) result(res)
            ! Return the quaternion representation
            ! of the rotation at time t
            class(rothist), intent(inout) :: me
            type(quaternion)              :: qtot, qdyn
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(4)
            real(wp)                      :: q(4)
            q = me%els%call(t)
            call me%renorm(q)
            qdyn%q = q
            qtot = qdyn*me%qstat
            res = qtot%q
        end function
        function getdcmdot(me, t) result(res)
            ! Return the matrix representation
            ! of the time derivative of the rotation at time t
            class(rothist), intent(inout) :: me
            real(wp)                      :: q(4), qdot(4)
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(3,3)
            q = me%callq(t)
            qdot = me%callqdot(t)
            associate (q0 => q(1), q1 => q(2), q2 => q(3), q3 => q(4), &
                     & q0dot =>qdot(1), q1dot => qdot(2), &
                     & q2dot => qdot(3), q3dot => qdot(4))
            res(1,1) = 2*q0*q0dot + 2*q1*q1dot - 2*q2*q2dot - 2*q3*q3dot
            res(1,2) = 2*q0*q3dot + 2*q0dot*q3 + 2*q1*q2dot + 2*q1dot*q2
            res(1,3) = -2*q0*q2dot - 2*q0dot*q2 + 2*q1*q3dot + 2*q1dot*q3
            res(2,1) = -2*q0*q3dot - 2*q0dot*q3 + 2*q1*q2dot + 2*q1dot*q2
            res(2,2) = 2*q0*q0dot - 2*q1*q1dot + 2*q2*q2dot - 2*q3*q3dot
            res(2,3) = 2*q0*q1dot + 2*q0dot*q1 + 2*q2*q3dot + 2*q2dot*q3
            res(3,1) = 2*q0*q2dot + 2*q0dot*q2 + 2*q1*q3dot + 2*q1dot*q3
            res(3,2) = -2*q0*q1dot - 2*q0dot*q1 + 2*q2*q3dot + 2*q2dot*q3
            res(3,3) = 2*q0*q0dot - 2*q1*q1dot - 2*q2*q2dot + 2*q3*q3dot
            end associate
        end function
        function getdcmddot(me, t) result(res)
            ! Return the matrix representation
            ! of the second time derivative of the rotation at time t
            class(rothist), intent(inout) :: me
            real(wp)                      :: q(4), qdot(4), qddot(4)
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(3,3)
            q = me%callq(t)
            qdot = me%callqdot(t)
            qddot = me%callqddot(t)
            associate (q0 => q(1), q1 => q(2), q2 => q(3), q3 => q(4), &
                     & q0dot =>qdot(1), q1dot => qdot(2), &
                     & q2dot => qdot(3), q3dot => qdot(4), &
                     & q0ddot =>qddot(1), q1ddot => qddot(2), &
                     & q2ddot => qddot(3), q3ddot => qddot(4))
                    res(1,1) = 2*q0*q0ddot + 2*q0dot**2 + 2*q1*q1ddot &
                           & + 2*q1dot**2 - 2*q2*q2ddot - 2*q2dot**2 &
                           & - 2*q3*q3ddot - 2*q3dot**2
                    res(1,2) = 2*q0*q3ddot + 2*q0ddot*q3 + 4*q0dot*q3dot &
                        & + 2*q1*q2ddot + 2*q1ddot*q2 + 4*q1dot*q2dot
                    res(1,3) = -2*q0*q2ddot - 2*q0ddot*q2 - 4*q0dot*q2dot &
                        & + 2*q1*q3ddot + 2*q1ddot*q3 + 4*q1dot*q3dot
                    res(2,1) = -2*q0*q3ddot - 2*q0ddot*q3 - 4*q0dot*q3dot &
                        & + 2*q1*q2ddot + 2*q1ddot*q2 + 4*q1dot*q2dot
                    res(2,2) = 2*q0*q0ddot + 2*q0dot**2 - 2*q1*q1ddot &
                        & - 2*q1dot**2 + 2*q2*q2ddot + 2*q2dot**2 &
                        & - 2*q3*q3ddot - 2*q3dot**2
                    res(2,3) = 2*q0*q1ddot + 2*q0ddot*q1 + 4*q0dot*q1dot &
                        & + 2*q2*q3ddot + 2*q2ddot*q3 + 4*q2dot*q3dot
                    res(3,1) = 2*q0*q2ddot + 2*q0ddot*q2 + 4*q0dot*q2dot &
                        & + 2*q1*q3ddot + 2*q1ddot*q3 + 4*q1dot*q3dot
                    res(3,2) = -2*q0*q1ddot - 2*q0ddot*q1 - 4*q0dot*q1dot &
                        & + 2*q2*q3ddot + 2*q2ddot*q3 + 4*q2dot*q3dot
                    res(3,3) = 2*q0*q0ddot + 2*q0dot**2 - 2*q1*q1ddot &
                        & - 2*q1dot**2 - 2*q2*q2ddot - 2*q2dot**2 &
                        & + 2*q3*q3ddot + 2*q3dot**2
            end associate
        end function
        function getrotquatdot(me, t) result(res)
            ! Return the four-vector representation
            ! of the time derivative of the rotation at time t
            class(rothist), intent(inout) :: me
            type(quaternion)              :: q, qdot, dummy
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(4), elsdot(4), &
                                             qnorm, qnorm2, qhatdot(4), &
                                             onebyqnorm, onebyqnorm3, &
                                             onebyqnorm2
            elsdot = me%elsdot%call(t)
            dummy%q = elsdot
            qdot = dummy*me%qstat
            q%q = me%callq(t)
            qnorm2 = sum(q%q*q%q)
            qnorm = sqrt(qnorm2)
            onebyqnorm = 1._wp/qnorm
            onebyqnorm2 = 1._wp/qnorm2
            onebyqnorm3 = 1._wp/qnorm**3
            qhatdot = onebyqnorm*qdot%q - sum(q%q*qdot%q)*onebyqnorm3*q%q
            res = qhatdot
        end function
        function getrotquatddot(me, t) result(res)
            ! Return the four-vector representation
            ! of the second time derivative of the rotation at time t
            class(rothist), intent(inout) :: me
            type(quaternion)              :: q, qdot, qddot, dummy
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(4), elsdot(4), elsddot(4),&
                                             qnorm, qnorm2, &
                                             onebyqnorm, onebyqnorm3, &
                                             onebyqnorm2, qddothat(4)
            elsdot = me%elsdot%call(t)
            elsddot = me%elsddot%call(t)
            dummy%q = elsdot
            qdot = dummy*me%qstat
            dummy%q = elsddot
            qddot = dummy*me%qstat
            q%q = me%callq(t)
            qnorm2 = sum(q%q*q%q)
            qnorm = sqrt(qnorm2)
            onebyqnorm = 1._wp/qnorm
            onebyqnorm2 = 1._wp/qnorm2
            onebyqnorm3 = 1._wp/qnorm**3
            qddothat = &
                    (3._wp*sum(q%q*qdot%q)**2*onebyqnorm2 &
                     - sum(q%q*qddot%q) &
                     - sum(qdot%q**2))*onebyqnorm3*q%q &
                  - 2._wp*sum(q%q*qdot%q)*onebyqnorm3*qdot%q &
                  + onebyqnorm*qddot%q
            res = qddothat
        end function
        pure subroutine quat_renorm(me)
            ! Renormalize the quaternion. Modifies
            ! the state of the calling instance.
            class(quaternion), intent(inout) :: me
            real(dp)                   :: q(4)
            real(dp)                   :: q_raw(4), norm
            q_raw = me%q
            norm = sqrt(sum(q_raw**2))
            q = q_raw/norm
            me%q = q
        end subroutine
        subroutine renorm(me, q)
            ! Renormalize the quaternion stored as elements in q. 
            class(rothist), intent(in) :: me
            real(dp),    intent(inout) :: q(4)
            real(dp)                   :: q_raw(4), norm
            q_raw = q
            norm = sqrt(sum(q_raw**2))
            q = q_raw/norm
        end subroutine
        pure function norm(q) result(res)
            ! Return the norm of the calling instance
            class(quaternion), intent(in) :: q
            real(wp) :: res
            res = sqrt(sum(q%q*q%q))
        end function norm
        elemental function smul1(s,q) result(res)
            type(quaternion), intent(in) :: q
            real(wp), intent(in)         :: s
            type(quaternion)             :: res
                res = quaternion(s*q%q)
        end function smul1
        elemental function smul2(q,s) result(res)
            type(quaternion), intent(in) :: q
            real(wp), intent(in)         :: s
            type(quaternion)             :: res
                res = quaternion(s*q%q)
        end function smul2
        elemental function imul1(s,q) result(res)
            type(quaternion), intent(in) :: q
            integer, intent(in)          :: s
            type(quaternion)             :: res
                res = real(s,wp)*q
        end function imul1
        elemental function imul2(q,s) result(res)
            type(quaternion), intent(in) :: q
            integer, intent(in)          :: s
            type(quaternion)             :: res
                res = real(s,wp)*q
        end function imul2
        elemental function qmul(q1,q2) result(res)
            type(quaternion), intent(in) :: q1,q2
            real(wp) :: prod(4), mat(4,4)
            type(quaternion) :: res
            integer i
            associate (q => q1%q, p => q2%q)
                ! This is not a traditional matrix multiply
                ! it is hard-coded to attempt to minimize
                ! striding and maximize parallelizability
                mat(:,1) = [q(1), -q(2), -q(3), -q(4)]
                mat(:,2) = [q(2),  q(1), -q(4),  q(3)]
                mat(:,3) = [q(3),  q(4),  q(1), -q(2)]
                mat(:,4) = [q(4), -q(3),  q(2),  q(1)]
                forall (i=1:4)
                    prod(i) = mat(1,i)*p(1) &
                            + mat(2,i)*p(2) &
                            + mat(3,i)*p(3) &
                            + mat(4,i)*p(4)
                end forall
            end associate
            ! prod = prod/norm2(prod)
            res = quaternion(prod)
        end function qmul
        elemental function qadd(q1,q2) result(res)
            type(quaternion), intent(in) :: q1,q2
            type(quaternion) :: res
            associate (a => q1%q, b => q2%q)
                res = quaternion([ a(1) + b(1), &
                                   a(2) + b(2), &
                                   a(3) + b(3), &
                                   a(4) + b(4)])
            end associate
        end function qadd
        elemental function sdiv(q,s) result(res)
            type(quaternion), intent(in) :: q
            real(wp), intent(in) :: s
            type(quaternion) :: res
            res = (1._wp/s)*q
        end function sdiv
        elemental function qsub(q1,q2) result(res)
            type(quaternion), intent(in) :: q1,q2
            type(quaternion) :: res
            associate (q => q1%q, p => q2%q)
                res = quaternion(q-p)
            end associate
        end function qsub
        elemental function qainv(q1) result(res)
            type(quaternion), intent(in) :: q1
            type(quaternion) :: res
                res = quaternion(-q1%q)
        end function qainv
        elemental function qconj(q1) result(res)
            type(quaternion), intent(in) :: q1
            type(quaternion) :: res
            associate (q => q1%q)
                res = quaternion([q(1), -q(2),-q(3), -q(4)])
            end associate
        end function qconj
        pure function asdcm(q) result (res)
            ! Return the DCM representation of the calling instance
            class(quaternion), intent(in) :: q
            real(wp) :: sq(0:3), q12, q30, q13, q20, q23, q10
            real(wp) :: res(3,3)
            integer i 
            sq = [(q%q(i)**2, i=1,4)]
            associate (q0 => q%q(1), q1 => q%q(2), q2 => q%q(3), q3 => q%q(4))
                q12 = q1*q2; q30 = q3*q0; q13 = q1*q3
                q20 = q2*q0; q23 = q2*q3; q10 = q1*q0
            end associate
            res(:,1) = [sq(1) - sq(2) -sq(3) + sq(0), &
                      & 2*(q12 + q30), &
                      & 2*(q13 - q20)]
            res(:,2) = [2*(q12 - q30), &
                      & -sq(1) + sq(2) - sq(3) + sq(0), &
                      & 2*(q23 + q10)]
            res(:,3) = [2*(q13 + q20), &
                      & 2*(q23 - q10), &
                      & -sq(1) - sq(2) + sq(3) + sq(0)]
        end function asdcm
        pure subroutine fromdcm(q,dcm) 
            ! Store a quaternion representation of dcm
            ! in the calling instance
            class(quaternion), intent(inout) :: q
            real(wp), intent(in) :: dcm(3,3)
            integer i
            q%q(1) = 0.5_wp*sqrt(1+sum([(dcm(i,i), i=1,3)]))
            q%q(2) = (dcm(3,2)-dcm(2,3))/(4*q%q(1))
            q%q(3) = (dcm(1,3)-dcm(3,1))/(4*q%q(1))
            q%q(4) = (dcm(2,1)-dcm(1,2))/(4*q%q(1))
            call q%renorm()
        end subroutine fromdcm
        pure subroutine axang(q,ax, ang)
            ! Store the quaternion representation
            ! of (ax,ang) in the calling instance
            class(quaternion), intent(inout) :: q
            real(wp), intent(in) :: ax(3), ang
            q%q(1) = cos(ang/2)
            q%q(2) = ax(1)*sin(ang/2)
            q%q(3) = ax(2)*sin(ang/2)
            q%q(4) = ax(3)*sin(ang/2)
            q%q = q%q/q%norm()
        end subroutine axang
        pure function ax(q) result (res)
            ! Return the rotation axis of the calling instance
            class(quaternion), intent(in) :: q
            real(wp) :: res(3)
            res = q%q(2:4)/sqrt(sum(q%q(2:4)**2))
        end function ax
        pure function ang(q) result (res)
            ! Return the rotation angle of the calling instance
            class(quaternion), intent(in) :: q
            real(wp) :: res
            res = 2._wp*atan2(sqrt(sum(q%q(2:4)**2)),q%q(1))
        end function ang
        pure function rotate_vec(q,v) result(res)
            ! Perform a vector rotation through the quaternion
            ! Stored in the calling instance
            class(quaternion), intent(in) :: q
            real(wp),          intent(in) :: v(3)
            real(wp)                      :: res(3), norm
            type(quaternion)              :: vquat
            norm = sqrt(sum(v**2))
            vquat = quaternion([0._wp, v])
            vquat = .T.q*vquat*q
            res = norm*vquat%q(2:)
        end function rotate_vec
end module quat
