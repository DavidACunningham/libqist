program main
    use frkmain, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
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
    real(qp), parameter  :: central_body_ref_radius=6381.137_qp, &
                            central_body_mu=398600.5_qp, &
                            mu_list(3)= [1.32712440041279419e11_qp, &
                                      &  4902.800118_qp, &
                                      &  126712764.1_qp]
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: init_state(8), eye(8,8), qist_stm(8,8),&
                          & analytic_stm(8,8), analytic_stt(8,8,8), & 
                          & divisor(8,8), init_stt(8**3), qist_stt(8,8,8), &
                          & stt_divisor(8,8,8)
    integer i, j, ntest
    real(qp), allocatable :: states(:,:), time(:)
    real(dp), allocatable :: dtime(:)

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
    ! print *, "Integrating base case"
    ! base_sol = solve_ivp(fd_eoms,&
    !               & [0._qp, 1._qp], &
    !               & [init_state, &
    !                  reshape(eye,[8**2]), &
    !                  init_stt], &
    !               & method="DOP853",&
    !               & dense_output=.false.,&
    !               & rtol=rtol, &
    !               & atol=atol, istep=0.5_qp)
    ! print *, "Done."
    print *,  "Initializing QIST"
    call qist%init(t0, tf, &
                 & "./perturbed_reference.subspice", &
                 &  traj_id, &
                 &  399, &
                 &  bodylist, &
                 &  central_body_mu, &
                 &  central_body_ref_radius, &
                 &  mu_list, &
                 &  shgrav, &
                 &  Cbar, Sbar, .true., .true.)
    qist%dynmod%tof = tof
    ntest = size(base_sol%ys,2)
    ! time = base_sol%ys(7,:)
    call read_from_file("twobody_integrated_t",dtime)
    time = real(dtime,qp)
    ntest = size(dtime)
    allocate(states(6,ntest))
    do i = 1, ntest
        states(:3,i) = real(qist%dynmod%bod_db%call(real(time(i),dp),qist%dynmod%traj_id,'p'),qp)
        states(4:,i) = real(qist%dynmod%bod_db%call(real(time(i),dp),qist%dynmod%traj_id,'v'),qp)
    end do
    call print_to_file("testx",real(states(1,:),8))
    call print_to_file("testy",real(states(2,:),8))
    call print_to_file("testz",real(states(3,:),8))
    ! call print_to_file("twobody_integrated_x", real(base_sol%ys(1,:),8))
    ! call print_to_file("twobody_integrated_y", real(base_sol%ys(2,:),8))
    ! call print_to_file("twobody_integrated_z", real(base_sol%ys(3,:),8))
    ! call print_to_file("twobody_integrated_xdot", real(base_sol%ys(4,:),8))
    ! call print_to_file("twobody_integrated_ydot", real(base_sol%ys(5,:),8))
    ! call print_to_file("twobody_integrated_zdot", real(base_sol%ys(6,:),8))
    ! call print_to_file("twobody_integrated_t", real(base_sol%ys(7,:),8))
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
        subroutine read_from_file(fname, var)
            integer io,j
            real(dp), intent(inout),allocatable :: var(:)
            real(dp)        :: dummy
            character(len=*) :: fname
            integer stat, nlines
            open(newunit=io, file=trim(adjustl(fname))//".txt")
            nlines = 0
            stat = 0
            do while (stat.eq.0)
                read(io,*,iostat=stat)
                nlines = nlines + 1
            end do
            rewind(io)
            if (allocated(var)) deallocate(var)
            allocate(var(nlines-1))
            do j = 1,nlines-1
                read(io,*) var(j)
                ! write(*,*) var(j)
            end do
            close(io)
        end subroutine read_from_file
end program main
