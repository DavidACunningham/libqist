program main
    use frkmain, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use denselight, only: lightsol
    use genqist, only: gqist
    use qist, only: itraj
    use tensorops, only: mattens, quad
    use makemodel, only: dynamicsmodel
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(gqist)         :: qist_i
    type(itraj)         :: it
    type(spice_subset)  :: subspice
    type(dynamicsModel) :: dyn
    type(odesolution)   :: base_sol, qist_sol
    type(lightsol)      :: packedsol
    character(len=12)   :: arg
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, tof = tf,&
                            rtol = 1.e-14_qp, atol = 1.e-20_qp
    integer, parameter   :: traj_id = -998, & 
                            central_body = 399, &
                            bodylist(3)= [10,301,5], &
                            ntest = 1000
    logical, parameter   :: shgrav = .false.
    integer              :: offset
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
                          & stt_divisor(8,8,8), &
                          & bsolbuf(8+8**2+8**3), &
                          & testjac_q(8,8), testhes_q(8,8,8), &
                          & testjac_b(8,8), testhes_b(8,8,8), &
                          & stmdiff(8,8), sttdiff(8,8,8), &
                          & testacc(8), &
                          & trand, &
                          & write_state(ntest,6), &
                          & times(ntest), &
                          & itraj_stm(8,8), &
                          & itraj_stt(8,8,8), &
                          & itraj_diff1(8,8), &
                          & itraj_diff2(8,8,8)
                        
    real(qp), allocatable :: qsolbuf(:)
    integer i, j, nsteps
    logical              :: run_base, run_qist, rails

    ! call random_number(trand)
    ! trand = trand*tof
    trand = 80343.0243863275644373278711349111887_qp
    call get_command_argument(1,arg)
    read(arg,*) run_base
    call get_command_argument(2,arg)
    read(arg,*) run_qist
    call get_command_argument(3,arg)
    read(arg,*) rails
    if (rails) then
        offset = 1
    else
        offset = 8
    endif
    allocate(qsolbuf(offset+8**2+8**3))

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
    call qist_i%init(t0, tf, &
                 & "./perturbed_reference.subspice", &
                 &  traj_id, &
                 &  399, &
                 &  bodylist, &
                 &  central_body_mu, &
                 &  central_body_ref_radius, &
                 &  mu_list, &
                 &  shgrav, &
                 &  Cbar, Sbar, rails,.true.)
    dyn%tof = tof

    qist_i%dynmod%tof = tof
    testjac_b = 0._qp
    testhes_b = 0._qp
    testjac_q = 0._qp
    testhes_q = 0._qp
    qist_i%dynmod%state = [qist_i%dynmod%trajstate(trand), trand, tof]
    dyn%state = [dyn%trajstate(trand), trand, tof]
    call qist_i%dynmod%get_derivs(trand,testacc, testjac_q, testhes_q)
    ! testacc =0._qp
    call dyn%get_derivs(trand, testacc, testjac_b,testhes_b)

    init_state = [ &
              4.898587196589413e-13_qp, &
              8000.0_qp, &
              0.0_qp, &
              -2.956795170934981_qp, &
              1.8105148709099377e-16_qp, &
              8.123727966096327_qp, &
                  t0, tof]
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
    dyn%state = init_state
    qist_i%dynmod%state = init_state
    

    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    
    print *, "Initial State"
    print *, real(init_state,4)

    if (run_base) then
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

    if (run_qist) then
        print *, "Integrating QIST"
        qist_i%rtol = 1.e-14_qp
        qist_i%atol = 1.e-20_qp
        qist_sol = qist_i%integrate(t0, tf)
        print *, "DONE"
        

        print *, "Writing out QIST solution"
        open(unit=75, file="qist_sol_perturbed.odesolution", access="stream", status="replace")
        call qist_sol%write(75,dp)
        close(75)
        print *, "DONE"
        print *, "Packing solution and writing out"
        call qist_i%pack("qist_sol_perturbed.odesolution",packedsol)
        print *, "DONE"
        open(unit=59, file="qist_sol_packed.lightsol", access='stream', status='replace')
        call packedsol%write(59)
        print *, "DONE"
        close(59)

        print *, "Initializing itraj"
        call it%init(real(t0,dp), 1._dp, &
                     "/home/david/wrk/nstgro/qist/libqist/fort/test/",&
                     "qist_sol_packed.lightsol")
        print *, "DONE"
        itraj_stm = real(it%stm(1._dp),qp)
        itraj_stt = real(it%stt(1._dp),qp)

        qsolbuf = qist_sol%call(1._qp)
        print *, "UNPACKED DENSE SOLUTION"
        do i=1,size(qsolbuf)
            print *, qsolbuf(i)
        end do
        print *, "PACKED DENSE SOLUTION"
        print *, it%call(1._dp)
    end if

    if (run_base) then
        bsolbuf = base_sol%call(1._qp)
        base_stm = reshape(bsolbuf(9:8+8**2), [8,8])
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
        times = [((tf-t0)/(ntest-1)*i + t0, i=0,ntest-1)]
        do i = 1, ntest
           bsolbuf = base_sol%call(times(i)/tf)
           write_state(i,:) = bsolbuf(:6)
        end do
        call print_to_file("base_sol_x",real(write_state(:,1),dp))
        call print_to_file("base_sol_y",real(write_state(:,2),dp))
        call print_to_file("base_sol_z",real(write_state(:,3),dp))
    end if

    if (run_qist) then
        qist_stm = reshape(qsolbuf(offset+1:offset+8**2), [8,8])
        qist_stt = reshape(qsolbuf(offset+1+8**2:), [8,8,8])
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
        print *, "ITRAJ STM"
        do i = 1,8
            print *, real(itraj_stm(i,:),4)
        end do
        print *, "ITRAJ STT"
        do i = 1,8
            print *, "PAGE", i
        do j = 1,8
            print *, real(itraj_stt(i,j,:),4)
        end do
        end do
    end if
    if (run_base.and.run_qist)then
        stmdiff = qist_stm - base_stm
        sttdiff = qist_stt - base_stt
        itraj_diff1 = itraj_stm - qist_stm
        itraj_diff2 = itraj_stt - qist_stt
        where (base_stm.ne.0._qp) stmdiff=stmdiff/base_stm
        where (base_stt.ne.0._qp) sttdiff=sttdiff/base_stt
        where (qist_stm.ne.0._qp) itraj_diff1=itraj_diff1/qist_stm
        where (qist_stt.ne.0._qp) itraj_diff2=itraj_diff2/qist_stt
        print *, "STM Diff"
        do i = 1,8
            print *, real(stmdiff(i,:),4)
        end do
        print *, "STT Diff"
        do i = 1,8
            print *, "PAGE", i
        do j = 1,8
            print *, real(sttdiff(i,j,:),4)
        end do
        end do
        print *, "ITRAJ STM Diff"
        do i = 1,8
            print *, real(itraj_diff1(i,:),4)
        end do
        print *, "ITRAJ STT Diff"
        do i = 1,8
            print *, "PAGE", i
        do j = 1,8
            print *, real(itraj_diff2(i,j,:),4)
        end do
        end do
    end if
    contains
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
