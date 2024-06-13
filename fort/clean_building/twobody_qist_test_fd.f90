program main
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use findiffmod
    use genqist, only: gqist
    use tensorops, only: mattens, quad
    use makemodel, only: dynamicsmodel
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(spice_subset)  :: subspice
    type(gqist)         :: gq
    type(odesolution)   :: base_sol !, qistsol
    character(len=1000) :: qist_config_file, arg, metakernel_filepath
    real(qp)            :: t0, tf, tof, &
                           rtol, atol, x0(6), fdrstep, fdvstep, fdtstep, &
                           epsvec(8)
    integer             :: reference_trajectory_id, & 
                           central_body_id, &
                           body_list(30), &
                           fdord
    logical             :: shgrav
    integer             :: stat, num, n_bodies
    real(qp)            :: central_body_ref_radius, &
                           central_body_mu, & 
                           mu_list(30)
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: init_state(8), eye(8,8), fd_stm(8,8), &
                          & analytic_stm(8,8), analytic_stt(8,8,8), & 
                          & analytic_final_state(8), &
                          & divisor(8,8), init_stt(8**3), fd_stt(8,8,8), &
                          & stt_divisor(8,8,8), jac(8,8), fd_jac(8,8), &
                          & hes(8,8,8), fd_hes(8,8,8)

    real(dp), allocatable :: spice_state(:,:)
    real(dp)              :: lt_dum
    integer i, j
    namelist /FD_CONFIG/   metakernel_filepath, &
                           qist_config_file, &
                           x0, &
                           t0, &
                           tf, &
                           rtol, &
                           atol, &
                           fdrstep, &
                           fdvstep, &
                           fdtstep, &
                           fdord

    call get_command_argument(1,arg)
    ! Read in namelist
    inquire(file=trim(adjustl(arg)), iostat=stat)
    if (stat .ne. 0) then 
        print *, "ERROR: Bad FD config namelist filename"
        stop
    end if
    open(file=trim(adjustl(arg)), status="old", &
         iostat=stat,newunit=num)
    read(unit=num, nml=FD_CONFIG, iostat=stat)
    if (stat .ne. 0) then 
        print *, "ERROR: bad FD config namelist format"
        print *, trim(adjustl(arg))
        print *, stat
        stop
    end if
    close(num)
    epsvec(1:3) = fdrstep
    epsvec(4:6) = fdvstep
    epsvec(7:8) = fdtstep
    call gq%init(qist_config_file )
    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    tof = 1._qp
    tof = tf-t0
    ! To integrate in real time, set tof to 1.
    gq%dynmod%tof = 1._qp
    init_state = [x0, t0, tof]
    ! gq%dynmod%tgt_on_rails = .false.
    gq%dynmod%shgrav = .false.
    jac = fdjacwrap(init_state)
    hes = heswrap(init_state)
    fd_jac = findiff_multiscale(fdwrap, init_state, epsvec, fdord)
    fd_hes = findiffhes_multiscale(fdjacwrap, init_state, epsvec, fdord)
    divisor = 1._qp
    where (abs(jac).ge.1.e-14_qp)
        divisor = analytic_stm
    end where
    stt_divisor = 1._qp
    where (abs(hes).ge.1.e-10_qp)
        stt_divisor = analytic_stt
    end where
    print *, "FD THEN ANALYTIC JACOBIAN THEN NORMALIZED ERROR"
    do i = 1,8
        print *, real(fd_jac(i,:),4)
        print *, real(jac(i,:),4)
        print *, real((jac(i,:) - fd_jac(i,:))/divisor(i,:),4)
        print *,  ""
    end do
    print *, "FD THEN ANALYTIC HESSIAN THEN NORMALIZED ERROR"
    do i = 1,8
    print *, "PAGE", i
    do j = 1,8
        print *, real(fd_hes(i,j,:),4)
        print *, real(hes(i,j,:),4)
        print *, real((hes(i,j,:) - fd_hes(i,j,:))/stt_divisor(i,j,:),4)
        print *,  ""
    end do
    end do




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

    call furnsh(trim(adjustl(metakernel_filepath)))
    allocate(spice_state(6,size(base_sol%ts)))
    do i=1,size(base_sol%ts)
        call spkgeo(gq%dynmod%traj_id, &
                    real(base_sol%ts(i)*tof+ t0,dp), &
                    "J2000", &
                    gq%dynmod%central_body, &
                    spice_state(:,i), &
                    lt_dum &
                   )

    end do
    analytic_final_state = base_sol%ys(:8,size(base_sol%ts))
    analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
    analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8])


    print *, "Integrating finite diff first order with step ", real(epsvec,4)
    fd_stm = findiff_multiscale(fd_integrate, init_state, epsvec, fdord)
    print *, "Done."

    print *, "Integrating finite diff second order with step ", real(epsvec,4)
    fd_stt = findiffhes_multiscale(fd_integrate_jac, init_state, epsvec, fdord)
    print *, "Done."

    divisor = 1._qp
    where (abs(analytic_stm).ge.1.e-14_qp)
        divisor = analytic_stm
    end where
    stt_divisor = 1._qp
    where (abs(analytic_stt).ge.1.e-10_qp)
        stt_divisor = analytic_stt
    end where
    print *, "FINAL TIME"
    print *, real(tf,8)
    print *, "FINAL STATE"
    print *, real(analytic_final_state,8)
    print *, "FINAL STATE ERROR"
    print *, real(analytic_final_state(:6) - spice_state(:,size(base_sol%ts)),8)
    print *, "FINAL POSITION ERROR"
    print *, real(norm2(analytic_final_state(:3) - spice_state(:3,size(base_sol%ts))),8)
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
    contains
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
        function fdwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x)), &
                                  & acc(8), jacdum(8,8), hesdum(8,8,8)
            res = gq%dynmod%fd_acc_nbody(x)
            res(:3) = 0._qp
            if (.not.gq%dynmod%shgrav) then
                res = res + gq%dynmod%fd_acc_kepler(x)
            else
                call gq%dynmod%allderivs_sh(x(7), x(:8), acc, jacdum, hesdum)
                res = res +  acc
            endif
        end function fdwrap
        function fdjacwrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(size(x),size(x)), &
                                  & accdum(8), jac(8,8), hesdum(8,8,8)
            res = gq%dynmod%fd_jac_nbody(x)
            res(:3,:) = 0._qp
            if (.not.gq%dynmod%shgrav) then
                res = res + gq%dynmod%fd_jac_kepler(x)
            else
                call gq%dynmod%allderivs_sh(x(7), x(:8), accdum, jac, hesdum)
                res = res + jac
            endif
        end function fdjacwrap

        function heswrap(x) result(res)
            real(qp), intent(in) :: x(:)
            real(qp)             :: res(8,8,8), &
                                  & accdum(8), jacdum(8,8), hes(8,8,8)
            res = gq%dynmod%fd_hes_nbody(x)
            res(:3,:,:) = 0._qp
            if (.not.gq%dynmod%shgrav) then
                res = gq%dynmod%fd_hes_kepler(x)
            else
                call gq%dynmod%allderivs_sh(x(7), x(:8), accdum, jacdum, hes)
                res = res + hes
            endif
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
