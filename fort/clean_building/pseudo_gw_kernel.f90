program main
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use findiffmod
    use genqist, only: gqist
    use tensorops, only: mattens, quad
    use makemodel, only: dynamicsmodel
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(gqist)         :: gq
    type(odesolution)   :: base_sol !, qistsol
    character(len=1000) :: qist_config_file, arg, metakernel_filepath, &
                           output_kernel_filename
    real(qp)            :: t0, tf, tof, &
                           rtol, atol
    integer             :: stat, num, nnodes
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: analytic_final_state(6), init_state(6)

    real(dp), allocatable :: spice_state(:,:)
    real(qp), allocatable :: kernel_times(:), &
                           & kernelstates(:,:)
    real(dp)              :: lt_dum, x0(6)
    real(dp), allocatable :: x(:,:), kernel_times_double(:)
    integer i, traj_id
    namelist /KERNEL_CONFIG/   metakernel_filepath, &
                               qist_config_file, &
                               output_kernel_filename, &
                               traj_id, &
                               x0, &
                               t0, &
                               tf, &
                               rtol, &
                               atol, &
                               nnodes

    call get_command_argument(1,arg)
    ! Read in namelist
    inquire(file=trim(adjustl(arg)), iostat=stat)
    if (stat .ne. 0) then 
        print *, "ERROR: Bad kernel config namelist filename"
        stop
    end if
    open(file=trim(adjustl(arg)), status="old", &
         iostat=stat,newunit=num)
    read(unit=num, nml=KERNEL_CONFIG, iostat=stat)
    if (stat .ne. 0) then 
        print *, "ERROR: bad kernel config namelist format"
        print *, trim(adjustl(arg))
        print *, stat
        stop
    end if
    close(num)
    call gq%init(qist_config_file )
    call furnsh(trim(adjustl(metakernel_filepath)))
    call spkgeo(gq%dynmod%traj_id, &
                real(t0,dp), &
                "J2000", &
                gq%dynmod%central_body, &
                x0, &
                lt_dum &
               )
    tof = 1._qp
    gq%dynmod%tof = 1._qp
    init_state = real(x0,qp)
    gq%dynmod%tgt_on_rails = .false.
    gq%dynmod%state = [init_state, t0, 1._qp]


    print *, "Integrating kernel trajectory"
    base_sol = solve_ivp(stateonly_eoms,&
                       & [t0, tf], &
                       & init_state, &
                       & dense_output=.true.,&
                       & rtol=rtol, &
                       & atol=atol, &
                       & istep=24._qp*3600._qp &
                      & )
    print *, "Done."
    allocate(spice_state(6,size(base_sol%ts)))
    do i=1,size(base_sol%ts)
        call spkgeo(gq%dynmod%traj_id, &
                    real(base_sol%ts(i),dp), &
                    "J2000", &
                    gq%dynmod%central_body, &
                    spice_state(:,i), &
                    lt_dum &
                   )
    end do
    analytic_final_state = base_sol%ys(:,size(base_sol%ts))

    print *, "FINAL TIME"
    print *, real(tf,8)
    print *, "FINAL STATE"
    print *, real(analytic_final_state,8)
    print *, "FINAL STATE ERROR"
    print *, real(analytic_final_state(:6) - spice_state(:,size(base_sol%ts)),8)
    print *, "FINAL POSITION ERROR"
    print *, real(norm2(analytic_final_state(:3) - spice_state(:3,size(base_sol%ts))),8)
    print *, "WRITING KERNEL. . ."

    allocate(kernel_times(nnodes), &
           & kernelstates(6,nnodes), &
             x(6,nnodes), &
             kernel_times_double(nnodes))

    kernel_times = [(t0 + i*(tf-t0)/(nnodes-1), i=0,nnodes-1)]
    do i = 1, nnodes
        kernelstates(:,i) = base_sol%call(kernel_times(i))
    end do
    x = real(kernelstates,dp)
    kernel_times_double = real(kernel_times,dp)

    call spkopn(trim(adjustl(output_kernel_filename)), "pseudo_gw", 100, num)
    call spkw13( &
                num, &
                traj_id, &
                301, &
                'J2000', &
                kernel_times_double(1), &
                kernel_times_double(nnodes), &
                "pseudo-gateway", &
                5, &
                nnodes, &
                x, &
                kernel_times_double &
                )
    call spkcls(num)
    contains
        function stateonly_eoms(me, x, y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp),          intent(in)    :: x, y(:)
            real(qp)                         :: res(size(y))
            real(qp)                         :: jac(8,8), hes(8,8,8), &
                                                acc(8)
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%tof = 1._qp
            gq%dynmod%state = [y, x, 1._qp]
            call gq%dynmod%get_derivs(x, acc, jac, hes)
            res = acc(:6)
        end function stateonly_eoms
        
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
