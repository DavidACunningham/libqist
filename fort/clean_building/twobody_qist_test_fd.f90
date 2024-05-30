program main
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use findiffmod
    use genqist, only: gqist
    use tensorops, only: mattens, quad
    use makemodel, only: dynamicsmodel
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    ! type(gqist)         :: qist
    type(spice_subset)   :: subspice
    type(dynamicsModel) :: dyn
    type(odesolution)   :: base_sol !, qistsol
    character(len=1000) :: resample_filepath, metakernel_filepath, arg
    real(qp)            :: t0, tf, tof, &
                            rtol, atol, x0(6), fdstep
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
                          & stt_divisor(8,8,8)

    real(dp), allocatable :: spice_state(:,:)
    real(dp)              :: lt_dum
    integer i, j
    namelist /FD_CONFIG/   resample_filepath, &
                           metakernel_filepath, &
                           x0, &
                           t0, &
                           tf, &
                           rtol, &
                           atol, &
                           reference_trajectory_id, &
                           central_body_id, &
                           central_body_ref_radius, &
                           central_body_mu, &
                           body_list, &
                           shgrav, &
                           mu_list, &
                           fdstep, &
                           fdord

    call get_command_argument(1,arg)
    mu_list = 0._qp
    body_list = 0._qp
    rtol = 1.e-10_qp
    metakernel_filepath = ""
    atol = 1.e-12_qp
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
    ! Set number of bodies
    n_bodies = findloc(body_list,0,dim=1)-1
    tof = 1._qp
    eye = 0._qp
    init_stt = 0._qp
    do i=1,8
        eye(i,i) = 1._qp
    end do
    open(file=trim(adjustl(resample_filepath)), status="old", &
         iostat=stat,newunit=num,access="stream")
    if (stat .ne. 0) then 
        print *, "ERROR: bad spice resample config file"
        stop
    end if
    call subspice%read(num)
    close(num)
    call dyn%init(subspice, &
                & reference_trajectory_id, &
                & central_body_id, &
                & body_list(:n_bodies), &
                & central_body_mu, &
                & central_body_ref_radius,  &
                & mu_list(:n_bodies), &
                & shgrav, &
                & Cbar, &
                & Sbar, &
                & .true. &
               & )
    ! To integrate in real time, set tof to 1.
    dyn%tof = 1._qp

    init_state = [x0, t0, tof]
    print *, "Integrating base case"
    base_sol = solve_ivp(fd_eoms,&
                       & [t0, tf], &
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
        call spkgeo(reference_trajectory_id, real(base_sol%ts(i),dp), "J2000", central_body_id, &
                     & spice_state(:,i), lt_dum)

    end do
    analytic_final_state = base_sol%ys(:8,size(base_sol%ts))
    analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
    analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8])


    print *, "Integrating finite diff first order with step ", fdstep
    fd_stm = findiff(fd_integrate, init_state, fdstep, fdord)
    print *, "Done"

    print *, "Integrating finite diff second order with step ", fdstep
    fd_stt = findiffhes(fd_integrate_jac, init_state, fdstep, fdord)
    print *, "Done"

    divisor = 1._qp
    where (abs(analytic_stm).ge.1.e-14_qp)
        divisor = analytic_stm
    end where
    stt_divisor = 1._qp
    where (abs(analytic_stt).ge.1.e-10_qp)
        stt_divisor = analytic_stt
    end where
    open(file="times.strm", access="stream", status="replace",newunit=num, iostat=stat)
    write(num) real(base_sol%ts,dp)
    close(num)
    open(file="base_sol_ys.strm", access="stream",status="replace", newunit=num, iostat=stat)
    write(num) real(base_sol%ys(:6,:),dp)
    close(num)
    open(file="spice_state.strm", access="stream", status="replace",newunit=num, iostat=stat)
    write(num) spice_state
    close(num)
    print *, "FINAL TIME"
    print *, real(tf,8)
    print *, "FINAL STATE"
    print *, real(analytic_final_state,8)
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
                             & [t0, tf], &
                             & [x(:6), t0, 1._qp], &
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
                             & [t0, tf], &
                             & [x(:6), t0, 1._qp, reshape(eye,[8**2])],&
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
