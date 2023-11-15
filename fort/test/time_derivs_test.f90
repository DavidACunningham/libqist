program timederivtest
    use, intrinsic :: iso_fortran_env, only: wp=>real64, dp=>real64, qp=>real128
    use cheby, only: spice_subset
    use makemodel, only: dynamicsmodel
    use findiffmod
    implicit none
    type(dynamicsModel)  :: dyn
    type(spice_subset)   :: subspice
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, test_time=0.879302_qp*tf, tof = 0.943939*tf
    integer, parameter   :: traj_id = -999, & 
                            central_body = 399, &
                            bodylist(3)= [10,301,5]
    logical, parameter   :: shgrav = .false.
    real(dp), parameter  :: central_body_ref_radius=6381.137_dp, &
                            central_body_mu=398600.5_dp, &
                            mu_list(3)= [1.32712440041279419e11_dp, &
                                      &  4902.800118_dp, &
                                      &  126712764.1_dp]
    real(dp), parameter  :: Cbar(2,2) = 0._dp, &
                            Sbar(2,2) = 0._dp

    real(qp)             :: acc_analytic(8), &
                          & jac_analytic(8,8), &
                          & hes_analytic(8,8,8), &
                          & jac_findiff(8,8), &
                          & hes_findiff(8,8,8), &
                          & test_state(8)
    integer i

    open(file=trim(adjustl("/home/david/wrk/nstgro/qist/libqist/fort/test/twobody_resample.subspice")),unit=73, &
       & access="stream", status="old")
    call subspice%read(73)
    close(73)
    call dyn%init(subspice, traj_id, central_body, bodylist, &
                & central_body_mu, central_body_ref_radius,  &
                & mu_list, shgrav, Cbar, Sbar,.true.)
    dyn%tof = tof
    call dyn%get_derivs(test_time, acc_analytic, jac_analytic, hes_analytic)

    test_state = [dyn%trajstate(test_time), test_time, tof]
    hes_analytic = heswrap(test_state)
    jac_findiff = findiff(kepwrap, test_state, 300._qp , 9) + &
                & findiff(tbwrap,  test_state, 300._qp , 9)
    hes_findiff = findiffhes(kepjacwrap, test_state, 300._qp , 9) + &
                & findiffhes(tbjacwrap,  test_state, 300._qp , 9)

    print *, "STATE"
    print *, real(test_state,4)
    print *, "ACCEL"
    print *, real(acc_analytic,4)
    print *, "FD ACCEL"
    print *, real(kepwrap(test_state) + tbwrap(test_state),4)
    print *, "ANALYTIC JACOBIAN"
    do i = 1,8
        print *, real(jac_analytic(i,:),4)
        print *,  ""
    end do
    print *, "FINITE DIFF JACOBIAN"
    do i = 1,8
        print *, real(jac_findiff(i,:),4)
        print *,  ""
    end do
    print *, "JACOBIAN ERROR"
    do i = 1,8
        print *, real(jac_findiff(i,:) - jac_analytic(i,:),4)
        print *,  ""
    end do
    print *, "JAC ERROR FROB"
    print *, real(maxval(abs(jac_findiff - jac_analytic)),4)
    print *, "HES ERROR FROB"
    print *, real(maxval(abs(hes_findiff - hes_analytic)),4)
    contains
        function tbwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            res = dyn%fd_acc_nbody(x)
            res(:3) = 0._qp
        end function tbwrap
        function kepwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            res = dyn%fd_acc_kepler(x)
        end function kepwrap
        function tbjacwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x))
            res = dyn%fd_jac_nbody(x)
            res(:3,:) = 0._qp
        end function tbjacwrap
        function kepjacwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x))
            res = dyn%fd_jac_kepler(x)
        end function kepjacwrap
        function heswrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(8,8,8)
            res = dyn%fd_hes_kepler(x)
            res = res + dyn%fd_hes_nbody(x)
        end function heswrap
        
! function findiffhes(j, xstar,eps,order) result(res)
! function findiff(f, xstar,eps,order) result(res)
end program timederivtest
