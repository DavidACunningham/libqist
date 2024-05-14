program main
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use findiffmod
    use genqist, only: gqist
    use tensorops, only: mattens, quad
    use makemodel, only: dynamicsmodel
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(gqist)         :: qist
    type(spice_subset)   :: subspice
    type(dynamicsModel) :: dyn
    type(odesolution)   :: base_sol, qistsol
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, tof = tf,&
                            rtol = 1.e-17_qp, atol = 1.e-20_qp
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
    real(qp)             :: init_state(8), eye(8,8), fd_stm(8,8),&
                          & analytic_stm(8,8), analytic_stt(8,8,8), & 
                          & divisor(8,8), init_stt(8**3), fd_stt(8,8,8), &
                          & stt_divisor(8,8,8)
    integer i, j

    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    open(file=trim(adjustl("/home/david/wrk/nstgro/qist/libqist/fort/test/twobody_resample.subspice")),unit=73, &
       & access="stream", status="old")
    call subspice%read(73)
    close(73)
    call dyn%init(subspice, traj_id, central_body, bodylist, &
                & central_body_mu, central_body_ref_radius,  &
                & mu_list, shgrav, Cbar, Sbar,.true.)
    dyn%tof = tof

    init_state = [dyn%trajstate(t0), t0, tof]
    print *, "Integrating base case"
    base_sol = solve_ivp(fd_eoms,&
                  & [0._qp, 1._qp], &
                  & [init_state, &
                     reshape(eye,[8**2]), &
                     init_stt], &
                  & method="DOP853",&
                  & dense_output=.false.,&
                  & rtol=rtol, &
                  & atol=atol, istep=0.5_qp)
    print *, "Done."

    analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
    analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8])

    print *, "Integrating finite diff first order"
    fd_stm = findiff(fd_integrate, init_state, 1.e-4_qp, 9)
    print *, "Done"

    print *, "Integrating finite diff second order"
    fd_stt = findiffhes(fd_integrate_jac, init_state, 1.e-4_qp, 9)
    print *, "Done"

    divisor = 1._qp
    where (abs(analytic_stm).ge.1.e-14_qp)
        divisor = analytic_stm
    end where
    stt_divisor = 1._qp
    where (abs(analytic_stt).ge.1.e-10_qp)
        stt_divisor = analytic_stt
    end where
    print *, "FD THEN ANALYTIC STM THEN NORMALIZED ERROR"
    do i = 1,8
        print *, real(fd_stm(i,:),4)
        print *, real(analytic_stm(i,:),4)
        print *, real((analytic_stm(i,:) - fd_stm(i,:))/divisor(i,:),4)
        print *,  ""
    end do
    print *, "FD THEN ANALYTIC STT THEN NORMALIZED ERROR"
    do i = 1,8
    print *, "PAGE", i
    do j = 1,8
        print *, real(fd_stt(i,j,:),4)
        print *, real(analytic_stt(i,j,:),4)
        print *, real((analytic_stt(i,j,:) - fd_stt(i,j,:))/stt_divisor(i,j,:),4)
        print *,  ""
    end do
    end do
    ! IF FD and analytic STM agree, integrate reference trajectory with dynamics under test, analytic STM, finite diff STT

    ! Integrate reference trajectory with dynamics under test and analytic STM and analytic STT

    ! Compare FD and analytic STM and analytic STT

    ! Save reference trajectory to SPICE kernel

    ! Integrate QIST trajectory using SPICE kernel

    ! Compare solutions


    !!! OLD TEST, BASIC FUNCTIONALITY ONLY
    ! call qist%init(0._qp, 2._qp*24._qp*3600._qp, &
    !              & "./twobody_resample.subspice", &
    !              & -999, &
    !              &  399, &
    !              &  [10, 301, 5], &
    !              &  398600.5_dp, &
    !              &  6378.137_dp, &
    !              & [1.32712440041279419e11_dp, &
    !              &  4902.800118_dp, &
    !              &  126712764.1_dp], &
    !              &  .False., &
    !              &  cdum, sdum, .true.)
    ! sol = qist%integrate(3600._qp, 5*3600._qp)
    

    ! open(unit=75, file="test_tbod_sol.odesolution", access="stream", status="replace")
    ! call sol%write(75,dp)
    ! close(75)

    ! print *, sol%call(0.33_qp*2)


    contains
        function fd_integrate(x) result(res)
            type(odesolution) :: fd_sol
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            fd_sol = solve_ivp(fd_eoms,&
                          & [0._qp, 1._qp], &
                          & x,&
                          & dense_output=.false.,&
                          & rtol=rtol, &
                          & atol=atol, istep=0.5_qp)
            res = fd_sol%ys(:,size(fd_sol%ts))
        end function fd_integrate
        function fd_integrate_jac(x) result(res)
            type(odesolution) :: fd_sol
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x)), &
                                  & finalstate(8+8**2)
            fd_sol = solve_ivp(fd_eoms,&
                          & [0._qp, 1._qp], &
                          & [x, reshape(eye,[8**2])],&
                          & dense_output=.false.,&
                          & rtol=rtol, &
                          & atol=atol, istep=0.5_qp)
            finalstate = fd_sol%ys(:,size(fd_sol%ts))
            res = reshape(finalstate(9:), [8,8])
        end function fd_integrate_jac
        function fdwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            res = dyn%fd_acc_nbody(x)
            res(:3) = 0._qp
            res = res + dyn%fd_acc_kepler(x)
        end function fdwrap
        function fdjacwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x))
            res = dyn%fd_jac_nbody(x)
            res(:3,:) = 0._qp
            res = res + dyn%fd_jac_kepler(x)
        end function fdjacwrap

        function heswrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(8,8,8)
            res = dyn%fd_hes_kepler(x)
            res = res + dyn%fd_hes_nbody(x)
        end function heswrap
        function fd_eoms(me, x, y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp),          intent(in)    :: x, y(:)
            real(qp)                         :: res(size(y))
            real(qp)                         :: stm(8,8), stmdot(8,8), &
                                                jac(8,8), hes(8,8,8), &
                                                stt(8,8,8), sttdot(8,8,8)
            if (size(y)==8) then

                res = fdwrap(y)
            else if (size(y)==8 + 8**2) then
                res(:8) = fdwrap(y)
                stm = reshape(y(9:),[8,8])
                stmdot = matmul(fdjacwrap(y(:8)),stm)
                res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
            else
                res(:8) = fdwrap(y)
                stm = reshape(y(9:8+8**2),[8,8])
                stt = reshape(y(9+8**2:), [8,8,8])
                jac = fdjacwrap(y(:8))
                hes = heswrap(y(:8))
                stmdot = matmul(jac,stm)
                sttdot = mattens(jac,stt,8) + quad(stm,hes,8)
                res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
                res(9+8**2:) = reshape(sttdot, [8**3])
            endif
        end function fd_eoms
        
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
