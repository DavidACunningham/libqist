program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel, only: dynamicsModel
    use findiffmod
    use cheby, only: spice_subset
    implicit none
    integer, parameter :: nnodes=100
    real(qp)                 :: t0, tf
    real(qp)                 :: y(8), acc(8), jac(8,8), hes(8,8,8), &
                                jac_fd(8,8), hes_fd(8,8,8), &
                                pos_tb(3), vel_tb(3), acc_tb(3), &
                                vel_tb_fd(3), acc_tb_fd(3), &
                                fdsteps(8), rpert, vpert, tpert, &
                                mu, ref_radius, calltime, pdum(3,1), &
                                testtimes(nnodes), hes_fd_hist(8,8,8,nnodes), &
                                hes_hist(8,8,8,nnodes), hes_diff_hist(8,8,8,nnodes)
    type(dynamicsModel)      :: dyn
    type(spice_subset)       :: subspice
    integer num, i, j, stat
    
    acc = 0._qp
    jac = 0._qp
    hes = 0._qp
    rpert = .01_qp ! km
    vpert = .001_qp ! km/s
    tpert = 1000._qp ! s
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
    calltime =  (tf-t0)/2.001_qp + t0
    testtimes = [(t0 + i*(tf-t0)/nnodes, i=1,nnodes)]
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
                  .false., &
                  .true.  &
                 )

    do i = 1,nnodes
        pos_tb = real(dyn%bod_db%call(real(testtimes(i),dp), &
                                           & dyn%bodylist(1),'p'), qp)
        vel_tb = real(dyn%bod_db%call(real(testtimes(i),dp), &
                                           & dyn%bodylist(1),'v'), qp)
        acc_tb = real(dyn%bod_db%call(real(testtimes(i),dp), &
                                           & dyn%bodylist(1),'a'), qp)
        y = [dyn%trajstate(testtimes(i)), testtimes(i), tf-t0]
        call dyn%allderivs_thirdbody(dyn%nbody_mus(1), y, &
                                     & pos_tb, vel_tb, acc_tb, &
                                     & acc, jac, hes)
        hes_hist(:,:,:,i) = hes
        hes_fd_hist(:,:,:,i) = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
    end do
    hes_diff_hist = hes_hist - hes_fd_hist
    call print_to_file("times_hes", real(testtimes,dp))
    call print_to_file("hes477", real(hes_hist(4,7,7,:),dp))
    call print_to_file("hes477_fd", real(hes_fd_hist(4,7,7,:),dp))
    call print_to_file("hes477_diff", real(hes_diff_hist(4,7,7,:),dp))
    pos_tb = real(dyn%bod_db%call(real(calltime,dp), &
                                       & dyn%bodylist(1),'p'), qp)
    vel_tb = real(dyn%bod_db%call(real(calltime,dp), &
                                       & dyn%bodylist(1),'v'), qp)
    acc_tb = real(dyn%bod_db%call(real(calltime,dp), &
                                       & dyn%bodylist(1),'a'), qp)
    y = [dyn%trajstate(calltime), calltime, tf-t0]

    call dyn%allderivs_thirdbody(dyn%nbody_mus(1), y, &
                                 & pos_tb, vel_tb, acc_tb, &
                                 & acc, jac, hes)
    jac_fd = 0._qp
    jac_fd = findiff_multiscale(fd_acc,y, fdsteps, 9)
    hes_fd = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
    pdum = findiffmat_td(fd_pos, calltime, tpert, 9, pdum)
    vel_tb_fd = pdum(:,1)
    pdum = findiffmat_td(fd_vel, calltime, tpert, 9, pdum)
    acc_tb_fd = pdum(:,1)

    
    print *, "TIME"
    print *, calltime
    print *, "GW STATE"
    print *, real(dyn%trajstate(calltime),4)
    print *, "FINITE DIFF STEPS"
    print *, real(fdsteps,4)
    print *, "INTERPOLATED VELOCITY"
    print *, real(vel_tb,4)
    print *, "FINITE DIFF VELOCITY"
    print *, real(vel_tb_fd,4)
    print *, "VELOCITY DIFFERENCE"
    print *, real(vel_tb_fd - vel_tb,4)
    print *, "INTERPOLATED ACCEL"
    print *, real(acc_tb,4)
    print *, "FINITE DIFF ACCEL"
    print *, real(acc_tb_fd,4)
    print *, "ACCEL DIFFERENCE"
    print *, real(acc_tb_fd - acc_tb,4)
    print *, "ANALYTIC ACCELERATION"
    print *, real(acc,4)
    print *, "ANALYTIC JACOBIAN"
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
    print *, "ANALYTIC HESSIAN"
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
    print *, "ANLAYTIC HESSIAN 4,7,1:3; 4,1:3,7, AND DIFF"
    print *, hes(4,7,1:3)
    print *, hes(4,1:3,7)
    print *, hes(4,7,1:3) - hes(4,1:3,7)
    print *, "FINITE DIFF HESSIAN 4,7,1:3; 4,1:3,7, AND DIFF"
    print *, hes_fd(4,7,1:3)
    print *, hes_fd(4,1:3,7)
    print *, hes_fd(4,7,1:3) - hes_fd(4,1:3,7)
    contains
    function fd_pos(x) result(res)
        real(qp), intent(in)  :: x
        real(qp), allocatable :: res(:,:)
        allocate(res, mold=pdum)
        res(:,1) = real(dyn%bod_db%call(real(x,dp), &
                   & dyn%bodylist(1),'p'), qp)
    end function
    function fd_vel(x) result(res)
        real(qp), intent(in)  :: x
        real(qp), allocatable :: res(:,:)
        allocate(res, mold=pdum)
        res(:,1) = real(dyn%bod_db%call(real(x,dp), &
                   & dyn%bodylist(1),'v'), qp)
    end function
    function fd_acc(x) result(res)
        real(qp), intent(in) :: x(:)
        real(qp)             :: res(size(x)), &
                              & jac(size(x),size(x)), &
                              & hes(size(x),size(x),size(x)), &
                              & thispos(3), thisvel(3), thisacc(3)
        res = 0._qp
        thispos = real(dyn%bod_db%call(real(x(7),dp), &
                                           & dyn%bodylist(1),'p'), qp)
        thisvel = real(dyn%bod_db%call(real(x(7),dp), &
                                           & dyn%bodylist(1),'v'), qp)
        thisacc = real(dyn%bod_db%call(real(x(7),dp), &
                                           & dyn%bodylist(1),'a'), qp)
        call dyn%allderivs_thirdbody(dyn%nbody_mus(1), x, thispos, thisvel, thisacc, res,jac,hes)
    end function fd_acc
    function fd_jac(x) result(res)
        real(qp), intent(in) :: x(:)
        real(qp)             :: res(size(x),size(x)), &
                              & acc(size(x)), &
                              & hes(size(x),size(x),size(x)), &
                              & thispos(3), thisvel(3), thisacc(3)
        res = 0._qp
        thispos = real(dyn%bod_db%call(real(x(7),dp), &
                                           & dyn%bodylist(1),'p'), qp)
        thisvel = real(dyn%bod_db%call(real(x(7),dp), &
                                           & dyn%bodylist(1),'v'), qp)
        thisacc = real(dyn%bod_db%call(real(x(7),dp), &
                                           & dyn%bodylist(1),'a'), qp)
        call dyn%allderivs_thirdbody(dyn%nbody_mus(1), x, thispos, thisvel, thisacc, acc,res,hes)
    end function fd_jac
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
