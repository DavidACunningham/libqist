program timederivtest
    use, intrinsic :: iso_fortran_env, only: wp=>real64, dp=>real64, qp=>real128
    use cheby, only: spice_subset
    use makemodel, only: dynamicsmodel
    use findiffmod
    implicit none
    type(dynamicsModel)  :: dyn
    type(spice_subset)   :: subspice
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, test_time=0.879302_qp*tf
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
    call dyn%get_derivs(test_time, acc_analytic, jac_analytic, hes_analytic)

    test_state = [dyn%trajstate(test_time), test_time, 1.5_qp*tf]
    jac_findiff = findiff(kepwrap, test_state, 1._qp , 9) + &
                & findiff(tbwrap,  test_state, 1._qp , 9)

    do i = 1,8
        print *, jac_findiff(i,:) - jac_analytic(i,:)
        print *, ""
    end do
    contains
        function tbwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            res = dyn%fd_acc_nbody(x)
        end function tbwrap
        function kepwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            res = dyn%fd_acc_kepler(x)
        end function kepwrap
        
! function findiffhes(j, xstar,eps,order) result(res)
! function findiff(f, xstar,eps,order) result(res)
end program timederivtest
