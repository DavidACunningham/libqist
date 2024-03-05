module quat
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
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
