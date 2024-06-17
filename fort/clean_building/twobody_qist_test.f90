program main
    use genqist, only: make_qist_model, gqist
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use qist, only: itraj
    use findiffmod
    use tensorops, only: mattens, quad
    use makemodel, only: dynamicsmodel
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(gqist)         :: gq
    type(itraj)         :: it
    type(odesolution)   :: base_sol !, qistsol
    character(len=1000) :: arg
    real(qp)            :: t0, tf, tof, &
                           rtol, atol
    integer             :: stat, num
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: init_state(8), eye(8,8), &
                          & analytic_stm(8,8), analytic_stt(8,8,8), & 
                          & analytic_final_state(8), &
                          & divisor(8,8), init_stt(8**3), &
                          & stt_divisor(8,8,8)
    real(dp)             :: itraj_stm(8,8), itraj_stt(8,8,8)
    integer i, j

    call get_command_argument(1,arg)
    ! call make_qist_model(trim(adjustl(arg)))
    call gq%init(trim(adjustl(arg)))
    call it%init("/home/david/wrk/nstgro/qist/libqist/fort/data/gw_itraj_params.nml")
    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    tof = gq%tf-gq%t0
    ! To integrate in real time, set tof to 1.
    gq%dynmod%tof = 1._qp
    init_state = [gq%dynmod%trajstate(gq%t0), gq%t0, tof]
    gq%dynmod%tgt_on_rails = .false.
    gq%dynmod%state = init_state
    rtol = 1.e-10_qp
    atol = 1.e-20_qp
    print *, "Integrating base case"
    base_sol = solve_ivp(fd_eoms,&
                       & [0._qp, 1._qp], &
                       & [init_state, &
                          reshape(eye,[8**2]), &
                          init_stt], &
                       & dense_output=.false.,&
                       & rtol=rtol, &
                       & atol=atol, &
                       & istep=12._qp*3600._qp &
                      & )
    print *, "Done."

    analytic_final_state = base_sol%ys(:8,size(base_sol%ts))
    analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
    analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8])

    call it%stts_ab(it%t0,it%tf,itraj_stm, itraj_stt)
    divisor = 1._qp
    where (abs(analytic_stm).ge.1.e-14_qp)
        divisor = analytic_stm
    end where
    stt_divisor = 1._qp
    where (abs(analytic_stt).ge.1.e-10_qp)
        stt_divisor = analytic_stt
    end where
    print *, "ITRAJ THEN ANALYTIC STM THEN NORMALIZED ERROR"
    do i = 1,8
        print *, real(itraj_stm(i,:),4)
        print *, real(analytic_stm(i,:),4)
        print *, real((analytic_stm(i,:) - real(itraj_stm(i,:),qp))/divisor(i,:),4)
        print *,  ""
    end do
    print *, "ITRAJ THEN ANALYTIC STT THEN NORMALIZED ERROR"
    do i = 1,8
    print *, "PAGE", i
    do j = 1,8
        print *, real(itraj_stt(i,j,:),4)
        print *, real(analytic_stt(i,j,:),4)
        print *, real((analytic_stt(i,j,:) - real(itraj_stt(i,j,:),qp))/stt_divisor(i,j,:),4)
        print *,  ""
    end do
    end do
    contains
        function fd_wrap_acc(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            real(qp)             :: jac(8,8), hes(8,8,8), &
                                    acc(8)

            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = x
            call gq%dynmod%get_derivs(x(7), acc, jac, hes)
            res = acc
        end function
        function fd_wrap_jac(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x))
            real(qp)             :: jac(8,8), hes(8,8,8), &
                                    acc(8)

            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = x
            call gq%dynmod%get_derivs(x(7), acc, jac, hes)
            res = jac
        end function
        function fd_integrate(x) result(res)
            type(odesolution) :: fd_sol
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x))
            fd_sol = solve_ivp(fd_eoms,&
                             & [0._qp, 1._qp], &
                             & x(:8), &
                             & dense_output=.false., &
                             & rtol=rtol, &
                             & atol=atol, &
                             & istep=0.5_qp &
                            & )
            res = fd_sol%ys(:,size(fd_sol%ts))
        end function fd_integrate
        function fd_integrate_jac(x) result(res)
            type(odesolution) :: fd_sol
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x)), &
                                  & finalstate(8+8**2)
            fd_sol = solve_ivp(fd_eoms,&
                             & [0._qp, 1._qp], &
                             & [x(:8), reshape(eye,[8**2])],&
                             & dense_output=.false.,&
                             & rtol=rtol, &
                             & atol=atol, &
                             & istep=0.5_qp &
                            & )
            finalstate = fd_sol%ys(:,size(fd_sol%ts))
            res = reshape(finalstate(9:), [8,8])
        end function fd_integrate_jac
        function fd_eoms(me, x, y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp),          intent(in)    :: x, y(:)
            real(qp)                         :: res(size(y))
            real(qp)                         :: stm(8,8), stmdot(8,8), &
                                                jac(8,8), hes(8,8,8), &
                                                stt(8,8,8), sttdot(8,8,8), &
                                                acc(8)

            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = y
            call gq%dynmod%get_derivs(y(7), acc, jac, hes)
            if (size(y)==8) then
                res = acc
            else if (size(y)==8 + 8**2) then
                res(:8) = acc
                stm = reshape(y(9:),[8,8])
                stmdot = matmul(jac,stm)
                res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
            else
                res(:8) = acc
                stm = reshape(y(9:8+8**2),[8,8])
                stt = reshape(y(9+8**2:), [8,8,8])
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
