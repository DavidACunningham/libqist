module test_quat
    use quat, only: quaternion, rothist
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use findiffmod
    implicit none
    contains
    subroutine rot_hist_test(testpass)
        logical, intent(inout) :: testpass
        integer, parameter :: nnodes=80
        real(dp), parameter :: pi = asin(1._dp)*2._dp
        real(dp), dimension(3,3) :: rotmat_comp
        real(dp)                 :: t0, tf
        real(dp)                 :: qdot(4), qdot_fd(4), &
                                    dcmdot(3,3), dcmddot(3,3), dcmdot_fd(3,3), dcmddot_fd(3,3), &
                                    qdum(4,1), midpoint_time
        type(rothist)            :: rot
        integer num, i, j, k
        !! ALL THAT'S NEEDED TO SETUP ROTS
        t0 = 0.12_dp
        tf = 4._dp*pi
        midpoint_time = (tf-t0)/2._dp + t0
        rotmat_comp(1,:) = [1._dp,         0._dp,         0._dp]
        rotmat_comp(2,:) = [0._dp, cos(pi/4._dp), sin(pi/4._dp)]
        rotmat_comp(3,:) = [0._dp,-sin(pi/4._dp), cos(pi/4._dp)]
        call rot%init(fitfun, t0, tf, nnodes, rotmat_comp)
        !! ALL THAT'S NEEDED TO SETUP ROTS

        dcmdot    = rot%calldot(midpoint_time)
        dcmddot = rot%callddot(midpoint_time)
        qdot = rot%callqdot(midpoint_time)
        qdum = real(findiffmat_td(fd_quat, real(midpoint_time,qp), (tf-t0)*1.e-7_qp, 9,real(qdum,qp) ),dp)
        qdot_fd = qdum(:,1)
        dcmdot =  rot%calldot(real(midpoint_time,dp))
        dcmddot = rot%callddot(real(midpoint_time,dp))
        dcmdot_fd = real(findiffmat_td(fd_rot, real(midpoint_time,qp), (tf-t0)*1.e-7_qp, 9, real(dcmdot,qp)),dp)
        dcmddot_fd = real(findiffmat_td(fd_rotdot, real(midpoint_time,qp), (tf-t0)*1.e-7_qp, 9, real(dcmdot,qp)),dp)

        print *, "TIME"
        print *, midpoint_time
        print *, "ANALYTIC QDOT"
            print *, real(qdot,4)
        print *, "FINITE DIFF QDOT"
            print *, real(qdot_fd,4)
        print *, "QDOT DIFFERENCE"
            print *, real(qdot_fd-qdot,4)
        print *, "ANALYTIC DCMDOT"
        do i=1,3
            print *, real(dcmdot(i,:),4)
        end do
        print *, "FINITE DIFF DCMDOT"
        do i=1,3
            print *, real(dcmdot_fd(i,:),4)
        end do
        print *, "DCMDOT DIFFERENCE"
        do i=1,3
            print *, real(dcmdot_fd(i,:)-dcmdot(i,:),4)
        end do
        print *, "ANALYTIC DCMDOTDOT"
        do i=1,3
            print *, real(dcmddot(i,:),4)
        end do
        print *, "FINITE DIFF DCMDOTDOT"
        do i=1,3
            print *, real(dcmddot_fd(i,:),4)
        end do
        print *, "DCMDOTDOT DIFFERENCE"
        do i=1,3
            print *, real(dcmddot_fd(i,:)-dcmddot(i,:),4)
        end do
        contains
        function fitfun(me, ta,tb) result(res)
            class(rothist), intent(inout) :: me
            real(dp), intent(in)          :: ta, tb
            real(dp)                      :: res(4), mat(3,3), x
            type(quaternion)              :: qclass
            x = tb
            mat(1,:) = [ cos(x*pi), sin(x*pi), 0._dp]
            mat(2,:) = [-sin(x*pi), cos(x*pi), 0._dp]
            mat(3,:) = [     0._dp,     0._dp, 1._dp]
            call qclass%fromdcm(mat)
            res = qclass%q
        end function
        function fd_quat(x) result(res)
            real(qp), intent(in) :: x
            real(qp), allocatable :: res(:,:)
            allocate(res(4,1))
            res(:,1) = real(rot%callq(real(x,dp)),qp)
        end function fd_quat
        function fd_rot(x) result(res)
            real(qp), intent(in) :: x
            real(qp), allocatable :: res(:,:)
            allocate(res(3,3))
            res = real(rot%call(real(x,dp)),qp)
        end function fd_rot
        function fd_rotdot(x) result(res)
            real(qp), intent(in) :: x
            real(qp), allocatable :: res(:,:)
            allocate(res(3,3))
            res = real(rot%calldot(real(x,dp)),qp)
        end function fd_rotdot
        end subroutine rot_hist_test

end module test_quat
