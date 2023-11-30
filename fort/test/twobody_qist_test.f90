program main
    use frkmain, only: solve_ivp, Odesolution, RungeKutta
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
    type(odesolution)   :: base_sol, qist_sol
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, tof = tf,&
                            rtol = 1.e-17_qp, atol = 1.e-20_qp
    integer, parameter   :: traj_id = -998, & 
                            central_body = 399, &
                            bodylist(3)= [10,301,5]
    logical, parameter   :: shgrav = .false.
<<<<<<< HEAD
    real(qp), parameter  :: central_body_ref_radius=6381.137_qp, &
                            central_body_mu=398600.5_qp, &
                            mu_list(3)= [1.32712440041279419e11_qp, &
                                      &  4902.800118_qp, &
                                      &  126712764.1_qp]
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: init_state(8), eye(8,8), qist_stm(8,8),&
                          & base_stm(8,8), base_stt(8,8,8), & 
                          & divisor(8,8), init_stt(8**3), qist_stt(8,8,8), &
                          & stt_divisor(8,8,8), qsolbuf(1+8**2+8**3), bsolbuf(8+8**2+8**3), &
                          & testjac_q(8,8), testhes_q(8,8,8), testjac_b(8,8), testhes_b(8,8,8), testacc(8)
    integer i, j, spkhand, nsteps
    logical, parameter :: debug = .true.
=======

    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    open(file=trim(adjustl("/home/david/wrk/nstgro/qist/libqist/fort/test/perturbed_reference.subspice")),unit=73, &
       & access="stream", status="old")
    call subspice%read(73)
    close(73)
    call dyn%init(subspice, traj_id, central_body, bodylist, &
                & central_body_mu, central_body_ref_radius,  &
                & mu_list, shgrav, Cbar, Sbar,.false.)
    call qist%init(t0, tf, &
                 & "./perturbed_reference.subspice", &
                 &  traj_id, &
                 &  399, &
                 &  bodylist, &
                 &  central_body_mu, &
                 &  central_body_ref_radius, &
                 &   mu_list, &
                 &  shgrav, &
                 &  Cbar, Sbar, .true.)
    dyn%tof = tof
    init_state = [dyn%trajstate(t0), t0, tof]

    qist%dynmod%tof = tof
    call qist%dynmod%get_derivs(t0,testacc, testjac_q, testhes_q)
    dyn%state = init_state
    call dyn%get_derivs(t0, testacc, testjac_b,testhes_b)
    print *, "JAC"
    do i=1,8
        print *, real(testjac_q(i,:) - testjac_b(i,:),4)
        end do
    print *, ""
    print *, "HES"
    do i=1,8
    do j=1,8
    print *, real(testhes_q(i,j,:) - testhes_b(i,j,:),4)
    print *, ""
    end do 
    end do
    
    

    if (.not.debug) then
    print *, init_state
    print *, "Integrating base case"
    base_sol = solve_ivp(base_eoms,&
                  & [0._qp, 1._qp], &
                  & [init_state, &
                     reshape(eye,[8**2]), &
                     init_stt], &
                  & method="DOP853",&
                  & dense_output=.true.,&
                  & rtol=rtol, &
                  & atol=atol, istep=0.5_qp)
    nsteps = size(base_sol%ts)
    print *, "Done."
    endif

    print *, "Integrating QIST"
    qist%rtol = 1.e-17_qp
    qist%atol = 1.e-20_qp
    qist_sol = qist%integrate(t0, tf)
    

    open(unit=75, file="qist_sol_perturbed.odesolution", access="stream", status="replace")
    call qist_sol%write(75,qp)
    close(75)

    print *, "DONE"
    qsolbuf = qist_sol%call(0.33_qp*2)
    bsolbuf = base_sol%call(0.33_qp*2)
    qist_stm = reshape(qsolbuf(1+1:1+8**2), [8,8])
    base_stm = reshape(bsolbuf(9:8+8**2), [8,8])
    qist_stt = reshape(qsolbuf(1+1+8**2:), [8,8,8])
    base_stt = reshape(bsolbuf(9+8**2:), [8,8,8])
    print *, "BASE STM"
    do i = 1,8
        print *, real(base_stm(i,:),4)
    end do
    print *, "BASE STT"
    do i = 1,8
        print *, "PAGE", i
    do j = 1,8
        print *, real(base_stt(i,j,:),4)
    end do
    end do
    print *, "QIST STM"
    do i = 1,8
        print *, real(qist_stm(i,:),4)
    end do
    print *, "QIST STT"
    do i = 1,8
        print *, "PAGE", i
    do j = 1,8
        print *, real(qist_stt(i,j,:),4)
    end do
    end do
    print *, "STM Diff"
    do i = 1,8
        print *, real(qist_stm(i,:) - base_stm(i,:),4)
    end do
    print *, "STT Diff"
    do i = 1,8
        print *, "PAGE", i
    do j = 1,8
        print *, real(qist_stt(i,j,:) - base_stt(i,j,:),4)
    end do
    end do
    contains
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
        function base_eoms(me, x, y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp),          intent(in)    :: x, y(:)
            real(qp)                         :: res(size(y))
            res = dyn%eoms(x,y)
        end function base_eoms
        
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
