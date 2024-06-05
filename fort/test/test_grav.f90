program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel, only: dynamicsModel
    use findiffmod
    use quat, only: rothist, quaternion
    use cheby, only: spice_subset
    implicit none
    real(dp), parameter :: epoch=10000._dp, epoch2=336000._dp
    integer, parameter :: nnodes=300
    real(dp), dimension(3,3) :: rotmat_comp
    real(qp)                 :: rtol, atol, t0, tf
    real(qp)                 :: y(8), acc(8), jac(8,8), hes(8,8,8), C(3,3), S(3,3), &
                                jac_fd(8,8), hes_fd(8,8,8), thisjacfd(8,8), &
                                fdsteps(17)
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
    fdsteps = [(10._qp**real(i,qp), i=3,-13, -1)]
    rtol = 1.e-13
    atol = 1.e-20
    t0 = 769269009.185_qp
    tf = 769845939.185_qp
    call furnsh("/home/david/wrk/nstgro/qist/kernels/mk_gw.tf")
    call pxform('MOON_PA','J2000',epoch,rotmat_comp)
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
    do i=1,17
        thisjacfd = findiff(fd_acc,y, fdsteps(i), 9)
        do j = 1,8
            do k = 1,8
                if (abs(thisjacfd(j,k)-jac(j,k))<abs(jac_fd(j,k)-jac(j,k))) then
                    jac_fd(j,k) = thisjacfd(j,k)
                end if
            end do
        end do
    end do
    print *, "TIME"
    print *, (tf-t0)/2 + t0
    print *, "GW STATE"
    print *, real(dyn%trajstate((tf-t0)/2 + t0),4)
    print *, "INERTIAL ACCELERATION"
    print *, real(acc,4)
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
    function fd_acc(x) result(res)
        real(qp), intent(in) :: x(:)
        real(qp)             ::  res(size(x)), &
                                 jacdum(8,8), &
                                 hesdum(8,8,8)
        call dyn%allderivs_sh((tf-t0)/2 + t0, x, res, jacdum, hesdum)

    end function fd_acc
end program main
