module genqist
    use globals
    use denselight, only: lightsol, qistpack
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use quat, only: rothist, quaternion
    implicit none

    type gqist
        type(dynamicsmodel)   :: dynmod
        real(qp), allocatable :: initstate(:)
        real(qp)              :: rtol, atol, t0, tf
        character(len=1000)   :: filepath_var
        contains
        generic, public :: init => gq_namelist_init, var_init
        procedure integrate
        procedure :: pack => packsol
        procedure, private :: gq_namelist_init
        procedure, private :: var_init
    end type gqist

    contains

    subroutine make_spice_subset(namefile)
        character(len=*), intent(in)  :: namefile
        character(len=1000)           :: resample_filepath, metakernel_filepath
        real(dp)                      :: t0,tf
        integer                       :: central_body, deg, n_bodies
        integer                       :: body_list(30)
        type(spice_subset)            :: subspice
        integer                       :: stat, num
        logical dasein
        namelist /RESAMPLE_CONFIG/ metakernel_filepath, &
                                   resample_filepath, &
                                   t0, &
                                   tf, &
                                   central_body, &
                                   body_list, &
                                   deg

        ! Defaults go here
        body_list = 0
        ! Read in namelist
        inquire(file=namefile, iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: Bad spice resample namelist filename"
            stop
        end if
        open(file=namefile, status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=RESAMPLE_CONFIG, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: bad build config namelist format"
            stop
        end if
        close(num)
        ! Set number of bodies
        n_bodies = findloc(body_list,0,dim=1)-1
        ! Print status
        print *, "Resampling SPICE bodies ", body_list(:n_bodies)
        print *, "Relative to body ", central_body
        print *, "from J2000 + ", t0
        print *, "to J2000 + ", tf
        print *, "degree of fit: ", deg
        ! Pass to fitter
        call subspice%init( &
                           trim(adjustl(metakernel_filepath)), & ! Spice kernel
                           central_body, &    ! central body
                           body_list(:n_bodies), & ! list of bodies to resample
                           t0, & ! epoch start
                           tf, & ! epoch end
                           deg & ! degree of fit
                          )
        print *, "writing resampled SPICE to ", adjustl(trim(resample_filepath))
        open(file=adjustl(trim(resample_filepath)),newunit=num,access="stream",&
             status="replace",iostat=stat)
        call subspice%write(num)
        if (stat .ne. 0) then 
            print *, "ERROR: resample file not written"
            stop
        end if
        close(num)
    end subroutine make_spice_subset

    subroutine make_rotation(namefile)
        character(len=*), intent(in) :: namefile
        character(len=16)            :: inertial_frame_string, &
                                      & rotating_frame_string
        character(len=1)             :: yn
        integer                      :: fitdeg
        character(len=1000)          :: rotfile_write, kernelfile_read
        real(dp), dimension(3,3)     :: rotmat_comp
        real(qp)                     :: t0, tf
        type(rothist)                :: rot
        integer num, stat
        logical dasein
        namelist /ROTCONFIG/ kernelfile_read, &
                             rotfile_write, &
                             fitdeg, &
                             t0, &
                             tf, &
                             inertial_frame_string, &
                             rotating_frame_string
        ! Read in namelist file
        inquire(file=trim(adjustl(namefile)), iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: "
            print *, "Bad rotation history configuration namelist filename"
            stop
        end if
        open(file=trim(adjustl(namefile)), status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=ROTCONFIG, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: "
            print *, "Bad rotation history configuration namelist format"
            stop
        end if
        close(num)
        
        call furnsh(trim(adjustl(kernelfile_read)))
        call pxform(inertial_frame_string, &
                  & rotating_frame_string, &
                  & real(t0,dp), &
                  & rotmat_comp &
                   )
        inquire(file=trim(adjustl(rotfile_write)), exist=dasein)
        if (dasein) then 
            print *, "WARNING: "
            print *, trim(adjustl(rotfile_write)), " already exists"
            print *, "Continuing will overwrite. Continue?"
            yn = "X"
            do while (yn.ne."Y".and.yn.ne."N")
                write(*,'(A)',advance='no') "Y/N: "
                read(*,*) yn
                print *, ""
            end do
            if (yn.eq."N") then
                print *, "Not continuing."
                stop
            end if
        end if
        call rot%init(fitfun, real(t0,dp), real(tf,dp), fitdeg, rotmat_comp)
        open(file=trim(adjustl(rotfile_write)), &
           & status="replace", &
           & access="stream", &
           & newunit=num &
            )
            call rot%write(num)
        close(num)
        contains
        function fitfun(me, ta,tb) result(res)
            class(rothist), intent(inout) :: me
            real(dp), intent(in)          :: ta, tb
            real(dp)                      :: res(4), mat(3,3)
            type(quaternion)              :: qclass
            call pxfrm2(trim(adjustl(rotating_frame_string)), &
                      & trim(adjustl(rotating_frame_string)), &
                      & ta, &
                      & tb, &
                      & mat &
                       )
            call qclass%fromdcm(mat)
            res = qclass%q
        end function
    end subroutine make_rotation

    subroutine generate_kernel(namefile)
        type(gqist)                  :: gq
        type(odesolution)            :: base_sol !, qistsol
        character(len=*), intent(in) :: namefile
        character(len=1000)          :: qist_config_file, &
                                        metakernel_filepath, &
                                        output_kernel_filename
        real(qp)                     :: t0, tf, tof, &
                                        rtol, atol
        integer                      :: stat, num, nnodes
        real(qp), parameter          :: Cbar(2,2) = 0._qp, &
                                        Sbar(2,2) = 0._qp
        real(qp)                     :: init_state(6)
        real(qp), allocatable        :: kernel_times(:), &
                                      & kernelstates(:,:)
        real(dp)                     :: x0(6)
        real(dp), allocatable        :: x(:,:), &
                                        kernel_times_double(:)
        logical dasein
        character(len=1) :: yn
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
        ! Read in namelist
        inquire(file=trim(adjustl(namefile)), iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: Bad kernel config namelist filename"
            stop
        end if
        open(file=trim(adjustl(namefile)), status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=KERNEL_CONFIG, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: bad kernel config namelist format"
            print *, trim(adjustl(namefile))
            print *, stat
            stop
        end if
        close(num)
        inquire(file=trim(adjustl(output_kernel_filename)), exist=dasein)
        if (dasein) then 
            print *, "WARNING: "
            print *, trim(adjustl(output_kernel_filename)), " already exists"
            print *, "Continuing will overwrite. Continue?"
            yn = "X"
            do while (yn.ne."Y".and.yn.ne."N")
                write(*,'(A)',advance='no') "Y/N: "
                read(*,*) yn
                print *, ""
            end do
            if (yn.eq."N") then
                print *, "Not continuing."
                stop
            else
                open(newunit=num,file=trim(adjustl(output_kernel_filename)))
                close(num,status="delete")
            end if
        end if
        call gq%init(qist_config_file)
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
                           & istep=(tf-t0)/2._qp &
                          & )
        print *, "Done."
        print *, "Writing kernel. . ."

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
        call spkopn(trim(adjustl(output_kernel_filename)), "SPK_file", 100, num)
        call spkw13( &
                    num, &
                    traj_id, &
                    gq%dynmod%central_body, &
                    'J2000', &
                    kernel_times_double(1), &
                    kernel_times_double(nnodes), &
                    "trajectory", &
                    5, &
                    nnodes, &
                    x, &
                    kernel_times_double &
                    )
        call spkcls(num)
        print *, "Done."
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
    end subroutine generate_kernel

    subroutine load_gravity_model(gravity_file, &
                                & ref_radius, &
                                & mu, &
                                & C, &
                                & S, &
                                & rot &
                               & )
        character(len=*),      intent(in)  :: gravity_file
        real(qp),              intent(out) :: ref_radius, &
                                            & mu
        real(qp), allocatable, intent(out) :: C(:,:), S(:,:)
        real(qp), allocatable              :: Cbar(:,:), &
                                            & Sbar(:,:)
        type(rothist),         intent(out) :: rot
        integer                            :: degord, stat, num
        character(len=1000)                :: rotfile
        logical dasein
        namelist /SPHERHARM/ rotfile, &
                             degord, &
                             ref_radius, &
                             mu, &
                             Cbar, &
                             Sbar 
        ! Read in namelist file
        allocate(Cbar(0:100,0:100),Sbar(0:100, 0:100))
        Cbar = 0._qp
        Sbar = 0._qp
        inquire(file=trim(adjustl(gravity_file)), iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: Bad gravity file namelist filename"
            stop
        end if
        open(file=gravity_file, status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=SPHERHARM, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: bad gravity file namelist format"
            stop
        end if
        close(num)
        allocate(C(degord+1,degord+1), &
                 S(degord+1,degord+1))
        C = Cbar(0:degord,0:degord)
        S = Sbar(0:degord,0:degord)
        inquire(file=trim(adjustl(rotfile)), iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: Bad body frame filename"
            stop
        end if
        open(file=rotfile, status="old", &
             iostat=stat, access="stream", newunit=num)
            call rot%read(num)
        close(num)
    end subroutine load_gravity_model

    subroutine gq_namelist_init(me,namefile)
        character(len=*), intent(in)    :: namefile
        class(gqist),     intent(inout) :: me
        character(len=1000)             :: resample_filepath, &
                                         & qist_filepath, &
                                         & gravity_filepath
        real(qp)                        :: t0, tf, tof, &
                                           rtol, atol
        integer                         :: reference_trajectory_id, & 
                                         & central_body_id, &
                                         & body_list(30)
        logical                         :: shgrav
        integer                         :: stat, num, n_bodies
        real(qp)                        :: central_body_ref_radius, &
                                         & central_body_mu, & 
                                         & mu_list(30)
        real(qp), allocatable           :: C(:,:), S(:,:)
        logical                         :: rails, dasein
        type(rothist)                   :: rot
        namelist /QIST_CONFIG/ resample_filepath, &
                               qist_filepath, &
                               gravity_filepath, &
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
                               mu_list
        ! Defaults go here
        rails = .true.
        mu_list = 0._qp
        body_list = 0._qp
        rtol = 1.e-10_qp
        atol = 1.e-12_qp
        gravity_filepath = "NONE"
        ! Read in namelist
        inquire(file=namefile, iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: Bad QIST config namelist filename"
            stop
        end if
        open(file=namefile, status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=QIST_CONFIG, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: bad QIST config namelist format"
            stop
        end if
        close(num)
        if (trim(adjustl(gravity_filepath)).ne."NONE") then
            call load_gravity_model(gravity_filepath, &
                                    & central_body_ref_radius, &
                                    & central_body_mu, &
                                    & C, &
                                    & S, &
                                    & rot &
                                   & )
            shgrav = .true.
        else
            central_body_ref_radius = 0._qp
            allocate(C(2,2), S(2,2))
            C = 0._qp
            S = 0._qp
            shgrav = .false.
        endif 
        ! Set number of bodies
        n_bodies = findloc(body_list,0,dim=1)-1
        tof = tf - t0
        ! Print status
        print *, "Loading QIST model for reference trajectory ", reference_trajectory_id
        print *, "Using ephemeris data from ", trim(adjustl(resample_filepath))
        print *, "Relative to body ", central_body_id
        print *, "With perturbations from bodies ", body_list(:n_bodies)
        print *, "With spherical harmonics gravity? ", shgrav
        print *, "From J2000 + ", t0
        print *, "To J2000 + ", tf
        print *, ""
        ! Initialize QIST
        call me%init(t0, tf, &
                     &  trim(adjustl(resample_filepath)), &
                     &  reference_trajectory_id, &
                     &  central_body_id, &
                     &  body_list(:n_bodies), &
                     &  central_body_mu, &
                     &  central_body_ref_radius, &
                     &  mu_list(:n_bodies), &
                     &  rails, &
                     &  shgrav, &
                     &  rot, &
                     &  C, S &
                     )
        me%filepath_var = qist_filepath
        me%dynmod%tof = tof
        me%rtol = rtol
        me%atol = atol
    end subroutine gq_namelist_init
    subroutine make_qist_model(namefile)
        character(len=*), intent(in)  :: namefile
        type(gqist)         :: qist_i
        type(odesolution)   :: qist_sol
        type(lightsol)      :: packedsol
        integer             :: stat, num, i
        call qist_i%init(namefile)
        print *, "Integrating QIST. . ."
        ! Integrate model
        qist_sol = qist_i%integrate(qist_i%t0, qist_i%tf)
        print *, "DONE"
        print *, "Writing QIST solution. . ."
        open(newunit=num, file="temp_qist_sol.odesolution", iostat=stat, &
             access="stream", status="replace")
        call qist_sol%write(num,dp)
        close(num)
        print *, "DONE"
        print *, "Packing solution and writing out. . ."
        call qist_i%pack("temp_qist_sol.odesolution",packedsol)
        print *, "DONE"
        print *, "Writing QIST model to file ", trim(adjustl(qist_i%filepath_var)), ". . ."
        open(newunit=num, file=trim(adjustl(qist_i%filepath_var)), iostat=stat, &
             access='stream', status='replace')
        call packedsol%write(num)
        print *, "DONE"
        close(num)
        print *, "Cleaning up temp file. . ."
        open(newunit=num, file="temp_qist_sol.odesolution", iostat=stat, &
             status="old")
        close(num, status="delete")
        print *, "DONE"
    end subroutine make_qist_model

    subroutine var_init(me,t0, tf, subspicefile, traj_id, central_body, bodylist, &
                  & central_body_mu, central_body_ref_radius, mu_list, &
                  & shgrav, rails, rot, Cbar, Sbar)
        ! init_dm: method to initialize dynamicsModel object
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! t0             real           lowest possible time in seconds past 
        !                               J2000 at which model integration may 
        !                               start
        ! tf             real           highest possible time in seconds past 
        !                               J2000 at which model integration may 
        !                               end
        ! subspicefile   character      The absolute directory path for the 
        !                               spice_subset object
        ! traj_id        integer        the integer ID of the reference 
        !                               trajectory in the SPICE kernel
        ! central_body   integer        the integer ID of the selected central
        !                               body in the SPICE kernel
        ! bodylist       integer (:)    list of integer IDs of the gravitating
        !                               bodies in the SPICE kernel, other than
        !                               the central body
        ! central_body_mu . . .
        !                real          gravitational parameter of central body
        ! central_body_ref_radius . . . 
        !                real          reference radius for SH data
        ! mu_list        real          gravitational parameters of bodies
        !                               sorted in the same order ad bodylist
        ! shgrav         logical        TRUE if central body will be modeled
        !                               by spherical harmonics
        ! rails          logical        is the reference state on rails
        ! check          logical        whether to compare integrated and rails state
        !                               during integration
        ! rot            rothist        the central body-fixed rotation matrix 
        !                               interpolation. Must be R_IU: Take
        !                               vectors FROM the body-fixed (U) frame
        !                               TO the inertial frame.
        ! Cbar           real (:,:)     4 pi (Kaula) normalized cosine Stokes 
        !                               for central body
        ! Sbar           real (:,:)     4 pi (Kaula) normalized sine Stokes 
        !                               for central body
        ! OUTPUTS:
        ! NONE
        class(gqist),      intent(inout)        :: me
        type(spice_subset)                      :: subspice
        type(rothist),     intent(in),              optional :: rot
        real(qp),          intent(in)                        :: t0, tf
        integer,           intent(in)                        :: traj_id, & 
                                                                central_body, &
                                                                bodylist(:)
        logical,           intent(in)                        :: shgrav, rails
        real(qp),          intent(in)                        :: central_body_ref_radius, &
                                                                central_body_mu, &
                                                                mu_list(:)
        real(qp),          intent(in), allocatable, optional :: Cbar(:,:), &
                                                                Sbar(:,:)
        character(len=*),  intent(in)                        :: subspicefile

        open(file=trim(adjustl(subspicefile)),unit=73, &
           & access="stream", status="old")
        call subspice%read(73)
        close(73)
        if (shgrav.and.present(rot).and.present(Cbar).and.present(Sbar)) then
        call me%dynmod%init(subspice, traj_id, central_body, bodylist, &
                          & central_body_mu, central_body_ref_radius,  &
                          & mu_list, shgrav, rails, rot, Cbar, Sbar)
        else
        call me%dynmod%init(subspice, traj_id, central_body, bodylist, &
                          & central_body_mu, central_body_ref_radius,  &
                          & mu_list, shgrav, rails)
        endif
        me%t0 = t0
        me%tf = tf
        me%rtol = 1.e-10
        me%atol = 1.e-12
    end subroutine
    subroutine model_accuracy_check(me_qist, t0, tf, &
                                    bodylist, mulist, &
                                    shgrav, Cbars, Sbars,&
                                    spicepoints, testpoints, testtimes)
        integer,      parameter     :: npoints=500
        class(gqist), intent(inout) :: me_qist
        integer,      intent(in)    :: bodylist(:)
        real(dp),     intent(in)    :: t0, tf
        real(dp),     intent(in)    :: mulist(:), &
                                       Cbars(:,:), &
                                       Sbars(:,:)
        logical,      intent(in)    :: shgrav
        real(dp),     intent(out)   :: testpoints(6,npoints), &
                                       spicepoints(6,npoints), &
                                       testtimes(npoints)
        type(ODESolution)           :: testsol
        integer i

        me_qist%dynmod%tgt_on_rails = .false.
        ! Change bodylist and SHGrav Parameters
        call me_qist%dynmod%new_bodies(bodylist, mulist)
        call me_qist%dynmod%new_gravstatus(shgrav, Cbars,Sbars)
        me_qist%dynmod%tof = 1._qp
        me_qist%dynmod%tgt_on_rails = .false.
        testsol = solve_ivp(test_eoms,&
                      & [real(t0,qp), real(tf,qp)], &
                        me_qist%dynmod%trajstate(real(t0,qp)), &
                      & dense_output=.true.,&
                      & rtol=me_qist%rtol, &
                      & atol=me_qist%atol, istep=0.5_qp)
        testtimes = [(t0 + (tf-t0)/npoints * i, i=0,npoints-1)]
        do i=1,npoints
            testpoints(:,i) = real(testsol%call( &
                                     real(testtimes(i),qp)), dp)
            spicepoints(:,i) = real(me_qist%dynmod%trajstate( &
                                     real(testtimes(i),qp)), dp)
        end do
        contains
            function test_eoms(me, x, y) result(res)
                class(RungeKutta), intent(inout) :: me
                real(qp),          intent(in)    :: x, y(:)
                real(qp)                         :: acc(8), jac(8,8), hes(8,8,8)
                real(qp)                         :: res(size(y))
                me_qist%dynmod%state = [y, x, 1._qp]
                call me_qist%dynmod%get_derivs(x,acc,jac,hes)
                res = acc(:6)
            end function test_eoms
        ! Compare integrated and spice solutions
    end subroutine model_accuracy_check
    function integrate(thisqist, t0, tf) result(res)
        ! integrate: integrate a QIST model
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! t0             real           integration start time in seconds past 
        !                               J2000 
        ! tf             real           integration stop time in seconds past 
        !                               J2000
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            ODESolution    the dense solution, ready for packing
        class(gqist), intent(inout) :: thisqist
        type(ODESolution)           :: res
        real(qp), intent(in)        :: t0, tf
        real(qp)                    :: eye(8,8), hess_init(8**3)
        real(qp), allocatable       :: init_state(:)
        integer i
        eye = 0._qp
        do i=1,8
            eye(i,i) = 1._qp
        end do
        hess_init = 0._qp
        thisqist%dynmod%tof = tf-t0
        if (thisqist%dynmod%tgt_on_rails) then
            allocate(init_state(1))
            init_state = [t0]
        else
            allocate(init_state(8))
            init_state(:6) = thisqist%dynmod%trajstate(t0)
            init_state(7:) = [t0, tf - t0]
        endif
        ! thisqist%dynmod%state = [thisqist%dynmod%trajstate(t0), t0, tf-t0]
        res = solve_ivp(myint_eoms,&
                      & [0._qp, 1._qp], &
                      & [ init_state, &
                         reshape(eye,[8**2]), &
                         hess_init], &
                      & dense_output=.true.,&
                      & rtol=thisqist%rtol, &
                      & atol=thisqist%atol, istep=0.5_qp)
        contains
            function myint_eoms(me, x, y) result(res)
                class(RungeKutta), intent(inout) :: me
                real(qp),          intent(in)    :: x, y(:)
                real(qp)                         :: res(size(y))
                if (thisqist%dynmod%tgt_on_rails) then
                    res = thisqist%dynmod%eoms_rails(x,y)
                else
                    res = thisqist%dynmod%eoms(x,y)
                endif
            end function myint_eoms
        end function integrate
    subroutine packsol(me,solfile,light) 
        class(gqist),   intent(inout)  :: me
        character(len=*), intent(in)   :: solfile
        type(lightSol),    intent(out) :: light
        call light%convert_from_file_and_pack(solfile,qistpack)
    end subroutine packsol
end module genqist
