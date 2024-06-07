program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel, only: dynamicsModel, accelonly_sh
    use findiffmod
    use quat, only: rothist, quaternion
    use cheby, only: spice_subset
    implicit none
    integer, parameter :: nnodes=1000
    real(dp), dimension(3,3) :: rotmat_comp
    real(qp)                 :: rtol, atol, t0, tf
    real(qp)                 :: y(8), acc(8), jac(8,8), hes(8,8,8), C(3,3), S(3,3), &
                                jac_fd(8,8), hes_fd(8,8,8), thisjacfd(8,8), &
                                qdot(4), qddot(4), qdot_fd(4), qddot_fd(4), &
                                dcmdot(3,3), dcmddot(3,3), dcmdot_fd(3,3), dcomddot_fd(3,3), &
                                fdsteps(20), qdum(4,1)
    type(rothist)            :: rot
    type(dynamicsModel)      :: dyn
    type(spice_subset)       :: subspice
    integer num, i, j, k
    !! ALL THAT'S NEEDED TO SETUP ROTS
    open(file="/home/david/wrk/nstgro/qist/libqist/fort/data/20240524_gw_resample.subspice", &
         newunit=num, status="old", access="stream")
    call subspice%read(num)
    close(num)
    fdsteps = 0._qp
    fdsteps = [(10._qp**real(i,qp), i=3,-16, -1)]
    rtol = 1.e-13
    atol = 1.e-20
    t0 = 769269009.185_qp
    tf = 769845939.185_qp
    call furnsh("/home/david/wrk/nstgro/qist/kernels/mk_gw.tf")
    call pxform('MOON_PA','J2000',real(t0,dp),rotmat_comp)
    call rot%init(fitfun, real(t0,dp), real(tf,dp), nnodes, rotmat_comp)
    !! ALL THAT'S NEEDED TO SETUP ROTS
    C(1,:) = [1._qp,    0._qp, 0._qp]
    C(2,:) = [0._qp,    0._qp, 0._qp]
    ! C(3,:) = [0._qp,    0._qp, 0._qp]
    C(3,:) = [1.e-3_qp, 0._qp, 1.e-6_qp]
    S = 0._qp

    call dyn%init(subspice, &
                  -60000, &
                  301, &
                  [399, &
                   10, &
                   5], &
                  4902.800066_qp, &
                  1740._qp, &
                  [398600.5_qp, &
                   132712440041.9394_qp, &
                   126712764.8_qp], &
                  .true., &
                  .true., &
                  rot, &
                  C, &
                  S &
                 )

    y = [dyn%trajstate((tf-t0)/2 + t0), (tf-t0)/2 + t0, 1._qp]
    call dyn%allderivs_sh((tf-t0)/2 + t0, y, acc, jac, hes)
    jac_fd = 1.e9_qp
    do i=1,20
        thisjacfd = findiff(fd_acc,y, fdsteps(i), 9)
        do j = 1,8
            do k = 1,8
                if (abs(thisjacfd(j,k)-jac(j,k))<abs(jac_fd(j,k)-jac(j,k))) then
                    jac_fd(j,k) = thisjacfd(j,k)
                end if
            end do
        end do
    end do

    qdot = real(rot%callqdot(real((tf-t0)/2 + t0,dp)),qp)
    qddot = real(rot%callqddot(real((tf-t0)/2 + t0,dp)),qp)
    qdum = findiffmat_td(fd_quat, (tf-t0)/2 + t0, 1._qp, 9,qdum )
    qdot_fd = qdum(:,1)
    dcmdot = real(rot%calldot(real((tf-t0)/2 + t0,dp)),qp)
    dcmddot = real(rot%callddot(real((tf-t0)/2 + t0,dp)),qp)
    dcmdot_fd = findiffmat_td(fd_rot, (tf-t0)/2 + t0, 1._qp, 9, dcmdot)
    print *, "TIME"
    print *, (tf-t0)/2 + t0
    print *, "GW STATE"
    print *, real(dyn%trajstate((tf-t0)/2 + t0),4)
    print *, "INERTIAL ACCELERATION"
    print *, real(acc,4)
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
    print *, "INERTIAL JACOBIAN"
    do i=1,8
        print *, real(jac(i,:),4)
    end do
    print *, "FINITE DIFF JACOBIAN"
    do i=1,8
        print *, real(jac_fd(i,:),4)
    end do
    print *, "JACOBIAN DIFFERENCE"
    do i=1,8
        print *, real(jac(i,:)-jac_fd(i,:),4)
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
    function fd_acc(x) result(res)
        real(qp), intent(in) :: x(:)
        real(qp)             ::  res(size(x))
        call accelonly_sh(dyn,(tf-t0)/2 + t0, x, res)
    end function fd_acc
end program main
