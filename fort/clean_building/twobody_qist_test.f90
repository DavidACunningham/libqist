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
    real(qp)            :: tof, &
                           rtol, atol
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: init_state(8), eye(8,8), &
                          & analytic_stm(8,8), analytic_stt(8,8,8), & 
                          & analytic_final_state(8), &
                          & divisor(8,8), init_stt(8**3), &
                          & stt_divisor(8,8,8), test_t0, test_tf
    real(dp)             :: itraj_stm0f(8,8), itraj_stt0f(8,8,8), itraj_stm_f(8,8), itraj_stt_f(8,8,8), &
                          & itraj_stm_0(8,8), itraj_stt_0(8,8,8)
    integer i, j

    call get_command_argument(1,arg)
    ! call make_qist_model(trim(adjustl(arg)))
    call gq%init(trim(adjustl(arg)))
    call it%init(trim(adjustl(arg)))
    ! call it%stts_ab(it%t0+24._dp*3600._dp,it%t0+48._dp*3600._dp,itraj_stmab, itraj_sttab)
    call it%stts_ab(it%t0,it%tf,itraj_stm0f, itraj_stt0f)
    itraj_stm_0 = it%stm(it%t0)
    itraj_stt_0 = it%stt(it%t0)
    itraj_stm_f = it%stm(it%tf)
    itraj_stt_f = it%stt(it%tf)
    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    test_t0 = gq%t0
    test_tf = gq%tf
    tof = test_tf-test_t0
    ! To integrate in real time, set tof to 1.
    ! tof = 1._qp
    gq%dynmod%tof = tof
    init_state = [gq%dynmod%trajstate(gq%t0), test_t0, tof]
    gq%dynmod%tgt_on_rails = .false.
    gq%dynmod%state = init_state
    print *, "Integrating base case"
    base_sol = solve_ivp(fd_eoms,&
                       & [0._qp, 1._qp], &
                       & [init_state, &
                          reshape(eye,[8**2]), &
                          init_stt], &
                       & dense_output=.false.,&
                       & rtol=gq%rtol, &
                       & atol=gq%atol, &
                       & istep=0.5_qp &
                      & )
    print *, "Done."

    analytic_final_state = base_sol%ys(:8,size(base_sol%ts))
    analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
    analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8])

    divisor = 1._qp
    where (abs(analytic_stm).ge.1.e-14_qp)
        divisor = analytic_stm
    end where
    stt_divisor = 1._qp
    where (abs(analytic_stt).ge.1.e-14_qp)
        stt_divisor = analytic_stt
    end where
    print *, "ITRAJ STM0"
    do i = 1,8
        print *, real(itraj_stm_0(i,:),4)
    end do
    print *, "ITRAJ THEN ANALYTIC STM THEN NORMALIZED ERROR"
    do i = 1,8
        print *, real(itraj_stm_f(i,:),4)
        print *, real(analytic_stm(i,:),4)
        print *, real((analytic_stm(i,:) - real(itraj_stm_f(i,:),qp))/divisor(i,:),4)
        print *,  ""
    end do
    print *, "ITRAJ THEN ANALYTIC STT THEN NORMALIZED ERROR"
    do i = 1,8
    print *, "PAGE", i
    do j = 1,8
        print *, real(itraj_stt_f(i,j,:),4)
        print *, real(analytic_stt(i,j,:),4)
        print *, real((analytic_stt(i,j,:) - real(itraj_stt_f(i,j,:),qp))/stt_divisor(i,j,:),4)
        print *,  ""
    end do
    end do
    print *, "CHAINED ITRAJ THEN ANALYTIC STM THEN NORMALIZED ERROR"
    do i = 1,8
        print *, real(itraj_stm0f(i,:),4)
        print *, real(analytic_stm(i,:),4)
        print *, real((analytic_stm(i,:) - real(itraj_stm0f(i,:),qp))/divisor(i,:),4)
        print *,  ""
    end do
    print *, "CHAINED ITRAJ THEN ANALYTIC STT THEN NORMALIZED ERROR"
    do i = 1,8
    print *, "PAGE", i
    do j = 1,8
        print *, real(itraj_stt0f(i,j,:),4)
        print *, real(analytic_stt(i,j,:),4)
        print *, real((analytic_stt(i,j,:) - real(itraj_stt0f(i,j,:),qp))/stt_divisor(i,j,:),4)
        print *,  ""
    end do
    end do
    contains
        function fd_eoms(me, x, y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp),          intent(in)    :: x, y(:)
            real(qp)                         :: res(size(y))
            real(qp)                         :: stm(8,8), stmdot(8,8), &
                                                jac(8,8), hes(8,8,8), &
                                                stt(8,8,8), sttdot(8,8,8), &
                                                acc(8)

            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = y(:8)
            call gq%dynmod%get_derivs(y(7), acc, jac, hes)
            res(:8) = acc
            stm = reshape(y(9:8+8**2),[8,8])
            stt = reshape(y(9+8**2:), [8,8,8])
            stmdot = matmul(jac,stm)
            sttdot = mattens(jac,stt,8) + quad(stm,hes,8)
            res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
            res(9+8**2:) = reshape(sttdot, [8**3])
        end function fd_eoms
end program main
