program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use findiffmod
    use quat, only: rothist, quaternion
    use cheby, only: vectorcheb, chnodes
    implicit none
    integer, parameter :: nnodes=500
    real(qp)                 :: rtol, atol, t0, tf
    real(qp)                 :: q_hist(4,nnodes), qdot_hist(4,nnodes), qdot_fd_hist(4,nnodes), &
                                fdsteps(20), qdum(4,1), test_times(nnodes), midpoint_time, &
                                fi(4,nnodes), dfi(nnodes), chebnodes(nnodes)
    type(vectorcheb)            :: rot, rotdot
    integer i
    !! ALL THAT'S NEEDED TO SETUP ROTS
    fdsteps = 0._qp
    fdsteps = [(10._qp**real(i,qp), i=3,-16, -1)]
    rtol = 1.e-13
    atol = 1.e-20
    ! t0 = 769269009.185_qp
    ! tf = 769845939.185_qp
    ! tf = 20*(tf-t0) + t0
    t0 = 0._qp
    tf = 500._qp*2*3.1415926_qp
    midpoint_time = (tf-t0)/2 + t0
    test_times = [(t0+ i*(tf-t0)/(nnodes-1), i=0,nnodes-1)]
    chebnodes = chnodes(nnodes,real(t0,dp),real(tf,dp))
    do i=1,4
        fi(i,:) = sin(chebnodes/500.)
    end do
    call rot%fit(real(transpose(fi),dp), real(t0,dp), real(tf,dp))
    rotdot = rot%deriv()
    !! ALL THAT'S NEEDED TO SETUP ROTS


    do i=1,nnodes
        q_hist(:,i) = real(rot%call(real(test_times(i),dp)),qp)
        qdot_hist(:,i) = real(rotdot%call(real(test_times(i),dp)),qp)
        qdum = findiffmat_td(fd_quat, test_times(i), 1.e-5_qp, 9,qdum )
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

    ! qdot = real(rot%callqdot(real(midpoint_time,dp)),qp)
    ! qddot = real(rot%callqddot(real(midpoint_time,dp)),qp)
    ! qdum = findiffmat_td(fd_quat, midpoint_time, 1._qp, 9,qdum )
    ! qdot_fd = qdum(:,1)
    ! dcmdot = real(rot%calldot(real(midpoint_time,dp)),qp)
    ! dcmddot = real(rot%callddot(real(midpoint_time,dp)),qp)
    ! dcmdot_fd = findiffmat_td(fd_rot, midpoint_time, 1._qp, 9, dcmdot)
    ! print *, "TIME"
    ! print *, midpoint_time
    ! print *, "ANALYTIC QDOT"
    !     print *, real(qdot,4)
    ! print *, "FINITE DIFF QDOT"
    !     print *, real(qdot_fd,4)
    ! print *, "QDOT DIFFERENCE"
    !     print *, real(qdot_fd-qdot,4)
    ! print *, "ANALYTIC DCMDOT"
    ! do i=1,3
    !     print *, real(dcmdot(i,:),4)
    ! end do
    ! print *, "FINITE DIFF DCMDOT"
    ! do i=1,3
    !     print *, real(dcmdot_fd(i,:),4)
    ! end do
    ! print *, "DCMDOT DIFFERENCE"
    ! do i=1,3
    !     print *, real(dcmdot_fd(i,:)-dcmdot(i,:),4)
    ! end do
    contains
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
        res(:,1) = real(rot%call(real(x,dp)),qp)
    end function fd_quat
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
