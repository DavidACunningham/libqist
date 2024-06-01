module quat
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use cheby, only: vectorcheb, chnodes
    implicit none
    integer, parameter :: wp = dp
    real(wp), parameter :: pi = 4._wp*atan(1._wp)
    

    type :: quaternion
        real(dp) :: q(4)

        contains
        procedure norm
        procedure fromdcm
        procedure asdcm
        procedure axang
        procedure rotate_vec
    end type quaternion
    type  rothist
        type(vectorcheb)             :: els, elsdot, elsddot
        type(quaternion)             :: qstat
        integer                      :: degree
        real(dp)                     :: t0, tf
        procedure(fit_func), pointer :: fit => null()

        contains
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
            import :: wp, rothist
            implicit none
            class(rothist), intent(inout) :: me
            real(wp), intent(in) :: ta, tb
            real(wp)             :: res(4)
        end function fit_func
    end interface
    interface operator (*)
        module procedure qmul
        module procedure smul1
        module procedure smul2
        module procedure imul1
        module procedure imul2
    end interface
    interface operator (+)
        module procedure qadd
    end interface
    interface operator (/)
        module procedure sdiv
    end interface
    interface operator (-)
        module procedure qsub
        module procedure qainv
    end interface
    interface operator (.T.)
        module procedure qconj
    end interface


    contains
        subroutine init_rot_q(me, fit, t0, tf, order, qstat)
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
            first = me%fit(t0, (t0+nodes(1))/10._wp)
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
            first = me%fit(t0, (t0+nodes(1))/10._wp)
            vals(1,:) = fitwrap(t0, nodes(1), first)
            do i = 2, order
                vals(i,:) = fitwrap(t0,nodes(i), vals(i-1,:))
            end do
            call me%els%fit(vals, t0,tf)
            call me%qstat%fromdcm(dcmstat)
            me%elsdot = me%els%deriv()
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
        function getdcm(me, t) result(res)
            class(rothist), intent(inout) :: me
            type(quaternion)              :: qtot, qdyn
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(3,3)
            real(wp)                      :: q(4)
            q = me%els%call(t)
            call me%renorm(q)
            qdyn%q = q
            qtot = qdyn*me%qstat
            res = qtot%asdcm()
        end function
        function getrotquat(me, t) result(res)
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
            class(rothist), intent(inout) :: me
            type(quaternion)              :: qdot
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(4)
            real(wp)                      :: q(4)
            q = me%elsdot%call(t)
            call me%renorm(q)
            qdot%q = q
            res = qdot%q
        end function
        function getrotquatddot(me, t) result(res)
            class(rothist), intent(inout) :: me
            type(quaternion)              :: qddot
            real(wp),       intent(in)    :: t
            real(wp)                      :: res(4)
            real(wp)                      :: q(4)
            q = me%elsddot%call(t)
            call me%renorm(q)
            qddot%q = q
            res = qddot%q
        end function
        subroutine renorm(me, q)
            class(rothist), intent(in) :: me
            real(dp),    intent(inout) :: q(4)
            real(dp)                   :: q_raw(4), norm
            q_raw = q
            norm = sqrt(sum(q_raw**2))
            q = q_raw/norm
        end subroutine
        pure function norm(q) result(res)
            class(quaternion), intent(in) :: q
            real(wp) :: res
            res = norm2(q%q)
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
            prod = prod/norm2(prod)
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
                      & 2*(q12 - q30), &
                      & 2*(q13 + q20)]
            res(:,2) = [2*(q12 + q30), &
                      & -sq(1) + sq(2) - sq(3) + sq(0), &
                      & 2*(q23 - q10)]
            res(:,3) = [2*(q13 - q20), &
                      & 2*(q23 + q10), &
                      & -sq(1) - sq(2) + sq(3) + sq(0)]
        end function asdcm
        pure subroutine fromdcm(q,dcm) 
            class(quaternion), intent(inout) :: q
            real(wp), intent(in) :: dcm(3,3)
            integer i
            q%q(1) = 0.5_wp*sqrt(1+sum([(dcm(i,i), i=1,3)]))
            q%q(2) = (dcm(3,2)-dcm(2,3))/(4*q%q(1))
            q%q(3) = (dcm(1,3)-dcm(3,1))/(4*q%q(1))
            q%q(4) = (dcm(2,1)-dcm(1,2))/(4*q%q(1))
        end subroutine fromdcm
        pure subroutine axang(q,ax, ang)
            class(quaternion), intent(inout) :: q
            real(wp), intent(in) :: ax(3), ang
            q%q(1) = cos(ang/2)
            q%q(2) = ax(1)*sin(ang/2)
            q%q(3) = ax(2)*sin(ang/2)
            q%q(4) = ax(3)*sin(ang/2)
            q%q = q%q/norm2(q%q)
        end subroutine axang
        pure function ax(q) result (res)
            class(quaternion), intent(in) :: q
            real(wp) :: res(3)
            res = q%q(2:4)/norm2(q%q(2:4))
        end function ax
        pure function ang(q) result (res)
            class(quaternion), intent(in) :: q
            real(wp) :: res
            res = 2._wp*atan2(norm2(q%q(2:4)),q%q(1))
        end function ang
        pure function rotate_vec(q,v) result(res)
            class(quaternion), intent(in) :: q
            real(wp),          intent(in) :: v(3)
            real(wp)                      :: res(3), norm
            type(quaternion)              :: vquat
            norm = norm2(v)
            vquat = quaternion([0._wp, v])
            vquat = .T.q*vquat*q
            res = norm*vquat%q(2:)
        end function rotate_vec
end module quat
