program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use findiffmod
    use quat, only: rothist, quaternion
    use cheby, only: spice_subset
    implicit none
    integer, parameter :: nnodes=50
    real(dp), dimension(3,3) :: rotmat_comp
    real(dp)                 :: t0, tf
    real(dp)                 :: qdot(4), qddot(4), qdot_fd(4), qddot_fd(4), &
                                q_hist(4,nnodes), qdot_hist(4,nnodes), qdot_fd_hist(4,nnodes), &
                                dcmdot(3,3), dcmddot(3,3), dcmdot_fd(3,3), dcmddot_fd(3,3), &
                                fdsteps(20), qdum(4,1), test_times(nnodes), midpoint_time
    type(rothist)            :: rot
    integer num, i, j, k
    !! ALL THAT'S NEEDED TO SETUP ROTS
    fdsteps = 0._qp
    fdsteps = [(10._dp**real(i,dp), i=3,-16, -1)]
    t0 = 769269009.185_dp
    tf = 769845939.185_dp
    tf = 10._dp*(tf-t0) + t0
    midpoint_time = (tf-t0)/2 + t0
    test_times = [(t0+ i*(tf-t0)/(nnodes-1), i=0,nnodes-1)]
    call furnsh("/home/david/wrk/nstgro/qist/kernels/mk_gw.tf")
    call pxform('MOON_PA','J2000',t0,rotmat_comp)
    call rot%init(fitfun, t0, tf, nnodes, rotmat_comp)
    !! ALL THAT'S NEEDED TO SETUP ROTS

    do i=1,nnodes
        q_hist(:,i) = rot%callq(test_times(i))
        qdot_hist(:,i) = rot%callqdot(test_times(i))
        qdum = real(findiffmat_td(fd_quat, &
                             real(test_times(i),qp), &
                             real((tf-t0)*10.e-7_dp,qp), &
                             9, &
                             real(qdum,qp) &
                             ),dp)
        qdot_fd_hist(:,i) = qdum(:,1)
    end do

    call print_to_file("q0",real(q_hist(1,:),dp))
    call print_to_file("q1",real(q_hist(2,:),dp))
    call print_to_file("q2",real(q_hist(3,:),dp))
    call print_to_file("q3",real(q_hist(4,:),dp))
    call print_to_file("qdot0",real(qdot_hist(1,:),dp))
    call print_to_file("qdot1",real(qdot_hist(2,:),dp))
    call print_to_file("qdot2",real(qdot_hist(3,:),dp))
    call print_to_file("qdot3",real(qdot_hist(4,:),dp))
    call print_to_file("qdotfd0",real(qdot_fd_hist(1,:),dp))
    call print_to_file("qdotfd1",real(qdot_fd_hist(2,:),dp))
    call print_to_file("qdotfd2",real(qdot_fd_hist(3,:),dp))
    call print_to_file("qdotfd3",real(qdot_fd_hist(4,:),dp))
    call print_to_file("testtimes",real(test_times,dp))

    qdot = rot%callqdot(real(midpoint_time,dp))
    qddot = rot%callqddot(real(midpoint_time,dp))
    qdum = real(findiffmat_td(fd_quat, real(midpoint_time,qp), 1._qp, 9,real(qdum,qp) ),dp)
    qdot_fd = qdum(:,1)
    dcmdot =  rot%calldot(real(midpoint_time,dp))
    dcmddot = rot%callddot(real(midpoint_time,dp))
    dcmdot_fd = real(findiffmat_td(fd_rot, real(midpoint_time,qp), 1._qp, 9, real(dcmdot,qp)),dp)
    dcmddot_fd = real(findiffmat_td(fd_rotdot, real(midpoint_time,qp), 1._qp, 9, real(dcmdot,qp)),dp)
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
        real(dp)                      :: res(4), mat(3,3)
        type(quaternion)              :: qclass
        call pxfrm2('MOON_PA','MOON_PA',ta,tb,mat)
        call qclass%fromdcm(mat)
        res = qclass%q
    end function
    function fd_sin(x) result(res)
        real(qp), intent(in) :: x
        real(qp), allocatable :: res(:,:)
        allocate(res(4,1))
        res(:,1) = sin(x)
    end function fd_sin
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
