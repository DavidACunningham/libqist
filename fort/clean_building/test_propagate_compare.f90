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
    real(qp)            :: tof
    real(qp)             :: init_state_ref(6), init_state_chaser(6), &
                          & test_t0, test_tf, dx0(6)
    real(dp),allocatable  :: itraj_state_hist(:,:)
    real(qp), allocatable :: rel_state_hist(:,:)
                          
    integer i

    call get_command_argument(1,arg)
    call gq%init(trim(adjustl(arg)))
    call it%init(trim(adjustl(arg)))
    dx0 = 5.e2_qp*[-1._qp, -1._qp, 1._qp, 0._qp, 0._qp, 0._qp]
    test_t0 = (gq%tf - gq%t0)/2._qp + gq%t0
    test_tf = test_t0 + 2._qp*24._qp*3600._qp
    tof = test_tf-test_t0
    ! To integrate in real time, set tof to 1.
    tof = 1._qp
    gq%dynmod%tof = tof
    init_state_ref = gq%dynmod%trajstate(test_t0)
    init_state_chaser = init_state_ref + dx0
    gq%dynmod%tgt_on_rails = .false.
    print *, "Integrating truth"
    base_sol = solve_ivp(ephem_eoms_two_sc,&
                       & [test_t0, test_tf], &
                       & [init_state_ref, &
                       &  init_state_chaser], &
                       & dense_output=.false.,&
                       & rtol=gq%rtol, &
                       & atol=gq%atol, &
                       & istep=0.5_qp &
                      & )
    allocate(rel_state_hist(6,size(base_sol%ts)), &
             itraj_state_hist(8,size(base_sol%ts)))
    print *, "Done."
    rel_state_hist = base_sol%ys(7:12,:) - base_sol%ys(:6,:)
    do i =1, size(base_sol%ts)
        itraj_state_hist(:,i) = it%prop(real(test_t0,dp), &
                                        real(base_sol%ts(i),dp), &
                                        [real(dx0,dp), 0._dp, 0._dp], &
                                        2)
    end do

    print *, real((itraj_state_hist(:6,size(base_sol%ts)) - &
                  rel_state_hist(:,size(base_sol%ts)))/rel_state_hist(:,size(base_sol%ts)),4)

    contains
        function ephem_eoms_two_sc(me, x, y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp),          intent(in)    :: x, y(:)
            real(qp)                         :: res(size(y))
            real(qp)                         :: jac(8,8), hes(8,8,8), &
                                                acc(8)
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = [y(:6), x, 1._qp]
            call gq%dynmod%get_derivs(x, acc, jac, hes)
            res(:6) = acc(:6)
            gq%dynmod%state = [y(7:12), x, 1._qp]
            call gq%dynmod%get_derivs(x, acc, jac, hes)
            res(7:12) = acc(:6)
        end function ephem_eoms_two_sc
end program main
