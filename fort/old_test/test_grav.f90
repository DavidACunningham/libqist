program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel, only: dynamicsModel, accelonly_sh
    use genqist, only: load_gravity_model
    use findiffmod
    use quat, only: rothist, quaternion
    use cheby, only: spice_subset
    implicit none
    integer, parameter :: nnodes=1000
    character(len=1000)      :: namefile
    real(qp)                 :: t0, tf
    real(qp)                 :: y(8), acc(8), jac(8,8), hes(8,8,8), &
                                jac_fd(8,8), hes_fd(8,8,8), &
                                dcmdot(3,3), dcmddot(3,3), dcmdot_fd(3,3), &
                                fdsteps(8), rpert, vpert, tpert, &
                                mu, ref_radius
    real(qp), allocatable    :: C(:,:), S(:,:)
    type(rothist)            :: rot
    type(dynamicsModel)      :: dyn
    type(spice_subset)       :: subspice
    integer num, i, j, stat
    
    namefile= "/home/david/wrk/nstgro/qist/libqist/fort/data/moon8by8.nml"
    call load_gravity_model(namefile,ref_radius, mu, C, S, rot)
    acc = 0._qp
    jac = 0._qp
    hes = 0._qp
    rpert = 1._qp ! km
    vpert = .01_qp ! km/s
    tpert =  1000._qp ! s
    fdsteps(:3) = rpert
    fdsteps(4:6) = vpert
    fdsteps(7:8) = tpert
    !! ALL THAT'S NEEDED TO SETUP ROTS
    open(file="/home/david/wrk/nstgro/qist/libqist/fort/data/20240524_gw_resample.subspice", &
         newunit=num, status="old", access="stream", iostat=stat)
    call subspice%read(num)
    close(num)
    t0 = 769269009.185_qp
    tf = 769845939.185_qp
    !! ALL THAT'S NEEDED TO SETUP ROTS

    call dyn%init(subspice, &
                  -60000, &
                  301, &
                  [399, &
                   10, &
                   5], &
                  mu, &
                  ref_radius, &
                  [398600.5_qp, &
                   132712440041.9394_qp, &
                   126712764.8_qp], &
                  .true., &
                  .true., &
                  rot, &
                  C, &
                  S &
                 )

    y = [dyn%trajstate((tf-t0)/2 + t0), (tf-t0)/2 + t0, tf-t0]
    call dyn%allderivs_sh((tf-t0)/2 + t0, y, acc, jac, hes)
    jac_fd = 0._qp
    jac_fd = findiff_multiscale(fd_acc,y, fdsteps, 9)
    hes_fd = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
    dcmddot = real(rot%callddot(real((tf-t0)/2 + t0,dp)),qp)
    dcmdot_fd = findiffmat_td(fd_rot, (tf-t0)/2 + t0, 1._qp, 9, dcmdot)
    print *, "TIME"
    print *, (tf-t0)/2 + t0
    print *, "GW STATE"
    print *, real(dyn%trajstate((tf-t0)/2 + t0),4)
    print *, "FINITE DIFF STEPS"
    print *, real(fdsteps,4)
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
    print *, "INERTIAL HESSIAN"
    do i=1,8
        print *, "PAGE ", i
        do j=1,8
        print *, real(hes(i,j,:),4)
        end do
    end do
    print *, "FINITE DIFF HESSIAN"
    do i=1,8
        print *, "PAGE ", i
        do j=1,8
        print *, real(hes_fd(i,j,:),4)
        end do
    end do
    print *, "HESSIAN DIFFERENCE"
    do i=1,8
        print *, "PAGE ", i
        do j=1,8
            print *, real(hes_fd(i,j,:)-hes(i,j,:),4)
        end do
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
        real(qp)             :: res(size(x)), &
                              & jac(size(x),size(x)), &
                              & hes(size(x),size(x),size(x))
        res = 0._qp
        call dyn%allderivs_sh(x(7), x, res,jac,hes)
    end function fd_acc
    function fd_jac(x) result(res)
        real(qp), intent(in) :: x(:)
        real(qp)             :: res(size(x),size(x)), &
                              & acc(size(x)), &
                              & hes(size(x),size(x),size(x))
        res = 0._qp
        call dyn%allderivs_sh(x(7), x, acc,res,hes)
    end function fd_jac
end program main
