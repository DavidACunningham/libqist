module genqist
    use globals
    use denselight, only: lightsol, qistpack
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    use quat, only: rothist, quaternion
    implicit none
    type configdata
        character(len=1000)   :: config_filename="", &
                                 metakernel_filename_with_trajectory="", &
                                 metakernel_filename_no_trajectory="", &
                                 resample_filename_with_trajectory="", &
                                 resample_filename_no_trajectory="", &
                                 rotation_filename="", &
                                 output_kernel_filename="", &
                                 qist_filename="", &
                                 kvtau_filename=""
        character(len=16)     :: inertial_frame_string="", &
                                 rotating_frame_string=""
        real(qp)              :: t0_resamp = 0._qp, &
                                 tf_resamp = 0._qp, &
                                 t0_qist = 0._qp, &
                                 tf_qist = 0._qp, &
                                 mu_list(30), &
                                 rtol_kernel=0._qp, &
                                 atol_kernel=0._qp, &
                                 x0(6)=0._qp, &
                                 rtol_qist=0._qp, &
                                 atol_qist=0._qp, &
                                 central_body_ref_radius = 0._qp, &
                                 central_body_mu = 0._qp
        real(qp), allocatable :: Cbar(:,:), &
                                 Sbar(:,:)
        integer               :: central_body_id=0, &
                                 body_list(30)=0, &
                                 n_bodies=0, &
                                 resamp_fit_deg=0, &
                                 rot_fit_deg=0, &
                                 traj_id=0, &
                                 spherical_harmonics_degree=0, &
                                 nnodes_kernel=0
        logical               :: regularize=.false.
        contains
            procedure :: init =>namelist_init
    end type configdata
    type gqist
        type(dynamicsmodel)   :: dynmod
        real(qp), allocatable :: initstate(:)
        real(qp)              :: rtol, atol, t0, tf
        character(len=1000)   :: filepath_var, &
                               & kvtau_filename
        contains
        generic, public :: init => gq_namelist_init, var_init
        procedure integrate
        procedure time_regularized_integrate
        procedure :: pack => packsol
        procedure, private :: gq_namelist_init
        procedure, private :: var_init
    end type gqist
    contains
    subroutine namelist_init(me, namefile)
        class(configdata), intent(inout) :: me
        character(len=*), intent(in)  :: namefile
        logical :: dasein
        integer :: num, stat
        !!! NAMELIST VARIABLES (with defaults)
        character(len=1000)   :: metakernel_filename_with_trajectory, &
                                 metakernel_filename_no_trajectory, &
                                 resample_filename_with_trajectory, &
                                 resample_filename_no_trajectory, &
                                 rotation_filename, &
                                 output_kernel_filename, &
                                 qist_filename, &
                                 kvtau_filename
        character(len=16)     :: inertial_frame_string, &
                                 rotating_frame_string
        real(qp)              :: t0_resamp, &
                                 tf_resamp, &
                                 t0_qist, &
                                 tf_qist, &
                                 mu_list(30), &
                                 rtol_kernel, &
                                 atol_kernel, &
                                 x0(6), &
                                 rtol_qist, &
                                 atol_qist, &
                                 central_body_ref_radius, &
                                 central_body_mu
        real(qp), allocatable :: Cbar(:,:), &
                                 Sbar(:,:)
        integer               :: central_body_id, &
                                 body_list(30), &
                                 n_bodies, &
                                 resamp_fit_deg, &
                                 rot_fit_deg, &
                                 traj_id, &
                                 spherical_harmonics_degree, &
                                 nnodes_kernel
        logical               :: regularize
        !!! MASTER CONFIGURATION NAMELIST
        namelist /QIST_CONFIGURATION/ metakernel_filename_with_trajectory, &
                                      metakernel_filename_no_trajectory, &
                                      resample_filename_with_trajectory, &
                                      resample_filename_no_trajectory, &
                                      rotation_filename, &
                                      qist_filename, &
                                      kvtau_filename, &
                                      output_kernel_filename, &
                                      inertial_frame_string, &
                                      rotating_frame_string, &
                                      t0_resamp, &
                                      tf_resamp, &
                                      t0_qist, &
                                      tf_qist, &
                                      mu_list, &
                                      rtol_kernel, &
                                      atol_kernel, &
                                      x0, &
                                      rtol_qist, &
                                      atol_qist, &
                                      central_body_ref_radius, &
                                      central_body_mu, &
                                      Cbar, &
                                      Sbar, &
                                      central_body_id, &
                                      body_list, &
                                      n_bodies, &
                                      resamp_fit_deg, &
                                      rot_fit_deg, &
                                      traj_id, &
                                      spherical_harmonics_degree, &
                                      nnodes_kernel, &
                                      regularize

        ! set defaults
        allocate(Cbar(0:100,0:100), Sbar(0:100,0:100))

        ! character defaults
        metakernel_filename_with_trajectory = ""
        metakernel_filename_no_trajectory   = ""
        resample_filename_with_trajectory   = ""
        resample_filename_no_trajectory     = ""
        rotation_filename                   = ""
        output_kernel_filename              = ""
        qist_filename                       = ""
        kvtau_filename                      = ""
        inertial_frame_string               = ""
        rotating_frame_string               = ""

        ! real defaults
        t0_resamp                           = 0._qp
        tf_resamp                           = 0._qp
        t0_qist                             = 0._qp
        tf_qist                             = 0._qp
        mu_list                             = 0._qp
        rtol_kernel                         = 0._qp
        atol_kernel                         = 0._qp
        x0                                  = 0._qp
        rtol_qist                           = 0._qp
        atol_qist                           = 0._qp
        central_body_ref_radius             = 0._qp
        central_body_mu                     = 0._qp
        Cbar                                = 0._qp
        Sbar                                = 0._qp

        ! integer defaults
        central_body_id                     = 0
        body_list                           = 0
        n_bodies                            = 0
        resamp_fit_deg                      = 0
        rot_fit_deg                         = 0
        traj_id                             = 0
        spherical_harmonics_degree          = 0
        nnodes_kernel                       = 0

        !logical defaults
        regularize                          = .false.
        ! Read in namelist
        inquire(file=namefile, iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: Bad QIST confiugration filename"
            stop
        end if
        open(file=namefile, status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=QIST_CONFIGURATION, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: bad QIST configuration namelist format"
            stop
        end if
        close(num)

        ! STORE IN DATA TYPE
        allocate(me%cbar(0:100,0:100), me%sbar(0:100,0:100))
        ! integers
        me%central_body_id                       = central_body_id
        me%body_list                             = body_list
        me%n_bodies                              = n_bodies
        me%resamp_fit_deg                        = resamp_fit_deg
        me%rot_fit_deg                           = rot_fit_deg
        me%traj_id                               = traj_id
        me%spherical_harmonics_degree            = spherical_harmonics_degree
        me%nnodes_kernel                         = nnodes_kernel

        !logicals
        me%regularize                            = regularize

        ! characters
        me%config_filename                       = namefile
        me%metakernel_filename_with_trajectory   = metakernel_filename_with_trajectory
        me%metakernel_filename_no_trajectory     = metakernel_filename_no_trajectory
        me%resample_filename_with_trajectory     = resample_filename_with_trajectory
        me%resample_filename_no_trajectory       = resample_filename_no_trajectory
        me%output_kernel_filename                = output_kernel_filename
        me%rotation_filename                     = rotation_filename
        me%inertial_frame_string                 = inertial_frame_string
        me%rotating_frame_string                 = rotating_frame_string
        me%qist_filename                         = qist_filename
        me%kvtau_filename                        = kvtau_filename

        ! reals
        me%t0_resamp                             = t0_resamp
        me%tf_resamp                             = tf_resamp
        me%t0_qist                               = t0_qist
        me%tf_qist                               = tf_qist
        me%mu_list                               = mu_list
        me%x0                                    = x0
        me%rtol_kernel                           = rtol_kernel
        me%atol_kernel                           = atol_kernel
        me%rtol_qist                             = rtol_qist
        me%atol_qist                             = atol_qist
        me%central_body_ref_radius               = central_body_ref_radius
        me%central_body_mu                       = central_body_mu
        me%Cbar                                  = Cbar
        me%Sbar                                  = Sbar
    end subroutine namelist_init
    subroutine make_spice_subset(namefile, traj_exist)
        type(configdata)              :: cd
        character(len=*), intent(in)  :: namefile
        logical, intent(in), optional :: traj_exist
        logical                       :: te
        character(len=1000)           :: resample_filepath, metakernel_filepath
        integer                       :: n_bodies
        type(spice_subset)            :: subspice
        integer                       :: stat, num
        call cd%init(namefile)
        if (present(traj_exist)) then
            te = traj_exist
        else
            te = .false.
        endif
        if (te) then
            metakernel_filepath = cd%metakernel_filename_with_trajectory
            resample_filepath = cd%resample_filename_with_trajectory
        else
            metakernel_filepath = cd%metakernel_filename_no_trajectory
            resample_filepath = cd%resample_filename_no_trajectory
        end if
        ! Set number of bodies
        n_bodies = findloc(cd%body_list,0,dim=1)-1
        if (te) then
            n_bodies = n_bodies + 1
            cd%body_list(n_bodies) = cd%traj_id
        end if
        ! Print status
        print *, "Resampling SPICE bodies ", cd%body_list(:n_bodies)
        print *, "Relative to body ", cd%central_body_id
        print *, "from J2000 + ", real(cd%t0_resamp,dp)
        print *, "to J2000 + ", real(cd%tf_resamp,dp)
        print *, "degree of fit: ", cd%resamp_fit_deg
        ! Pass to fitter
        call subspice%init( &
                           trim(adjustl(metakernel_filepath)), & ! Spice kernel
                           cd%central_body_id, &    ! central body
                           cd%body_list(:n_bodies), & ! list of bodies to resample
                           real(cd%t0_resamp,dp), & ! epoch start
                           real(cd%tf_resamp,dp), & ! epoch end
                           cd%resamp_fit_deg & ! degree of fit
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
        type(configdata)             :: cd
        character(len=*), intent(in) :: namefile
        character(len=1)             :: yn
        character(len=1000)          :: kernelfile_read
        real(dp), dimension(3,3)     :: rotmat_comp
        type(rothist)                :: rot
        integer num
        logical dasein
        call cd%init(namefile)
        if (cd%metakernel_filename_no_trajectory=="") then
            kernelfile_read = cd%metakernel_filename_with_trajectory
        else
            kernelfile_read = cd%metakernel_filename_no_trajectory
        end if
        call furnsh(trim(adjustl(kernelfile_read)))
        call pxform(cd%inertial_frame_string, &
                  & cd%rotating_frame_string, &
                  & real(cd%t0_resamp,dp), &
                  & rotmat_comp &
                   )
        inquire(file=trim(adjustl(cd%rotation_filename)), exist=dasein)
        if (dasein) then 
            print *, "WARNING: "
            print *, trim(adjustl(cd%rotation_filename)), " already exists"
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
        call rot%init( &
                    & fitfun, &
                    & real(cd%t0_resamp,dp), &
                    & real(cd%tf_resamp,dp), &
                    & cd%rot_fit_deg, &
                    & rotmat_comp &
                    &)
        open(file=trim(adjustl(cd%rotation_filename)), &
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
            call pxfrm2(trim(adjustl(cd%rotating_frame_string)), &
                      & trim(adjustl(cd%rotating_frame_string)), &
                      & ta, &
                      & tb, &
                      & mat &
                       )
            call qclass%fromdcm(mat)
            res = qclass%q
        end function
    end subroutine make_rotation
    subroutine generate_kernel(namefile)
        type(configdata)             :: cd
        type(gqist)                  :: gq
        type(odesolution)            :: base_sol !, qistsol
        character(len=*), intent(in) :: namefile
        character(len=1000)          :: qist_config_file, &
                                        output_kernel_filename
        real(qp)                     :: t0, tf, tof, &
                                        rtol, atol
        integer                      :: num, nnodes
        real(qp)                     :: init_state(6)
        real(qp), allocatable        :: kernel_times(:), &
                                      & kernelstates(:,:)
        real(dp)                     :: x0(6)
        real(dp), allocatable        :: x(:,:), &
                                        kernel_times_double(:)
        logical dasein
        character(len=1) :: yn
        integer i, traj_id
        call cd%init(namefile)
        qist_config_file = cd%config_filename
        output_kernel_filename = cd%output_kernel_filename
        traj_id = cd%traj_id
        x0 = real(cd%x0,dp)
        t0 = cd%t0_resamp
        tf = cd%tf_resamp
        rtol = cd%rtol_kernel
        atol = cd%atol_kernel
        nnodes = cd%nnodes_kernel
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
        ! Make sure nonphysical integration stuff is turned off
        gq%dynmod%regularize = .false.
        tof = 1._qp
        gq%dynmod%tof = 1._qp
        gq%dynmod%tgt_on_rails = .false.
        init_state = real(x0,qp)
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
    subroutine gq_namelist_init(me,namefile,traj_exist)
        type(configdata)                :: cd
        character(len=*), intent(in)    :: namefile
        class(gqist),     intent(inout) :: me
        logical, intent(in), optional   :: traj_exist
        logical                         :: te
        character(len=1000)             :: resample_filepath, &
                                         & rotfile
        real(qp)                        :: tof
        logical                         :: shgrav
        integer                         :: stat, num, n_bodies, degord
        real(qp), allocatable           :: C(:,:), S(:,:)
        logical                         :: rails, dasein
        type(rothist)                   :: rot

        call cd%init(namefile)
        if (present(traj_exist)) then
            te = traj_exist
        else
            te = .false.
        endif
        if(te) then
            resample_filepath = cd%resample_filename_with_trajectory
        else
            resample_filepath = cd%resample_filename_no_trajectory
        endif
        degord = cd%spherical_harmonics_degree
        if (degord>0) then
            shgrav = .true.
            allocate(C(degord+1,degord+1), &
                     S(degord+1,degord+1))
            C = cd%Cbar(0:degord,0:degord)
            S = cd%Sbar(0:degord,0:degord)
            rotfile = cd%rotation_filename
        else
            shgrav = .false.
            allocate(C(2,2),S(2,2))
            C = 0._qp
            S = 0._qp
        end if
        inquire(file=trim(adjustl(rotfile)), iostat=stat, exist=dasein)
        if (stat .ne. 0 .or. .not.dasein) then 
            print *, "ERROR: Bad body frame filename"
            stop
        end if
        open(file=rotfile, status="old", &
             iostat=stat, access="stream", newunit=num)
            call rot%read(num)
        close(num)
        shgrav = .true.
        rails = .true.
        ! Set number of bodies
        n_bodies = findloc(cd%body_list,0,dim=1)-1
        tof = cd%tf_qist - cd%t0_qist
        ! Print status
        print *, "Loading dynamics model for reference trajectory ", cd%traj_id
        print *, "Using ephemeris data from ", trim(adjustl(resample_filepath))
        print *, "Relative to body ", cd%central_body_id
        print *, "With perturbations from bodies ", cd%body_list(:n_bodies)
        print *, "With spherical harmonics gravity? ", shgrav
        print *, "From J2000 + ", cd%t0_qist
        print *, "To J2000 + ", cd%tf_qist
        print *, ""
        ! Initialize QIST
        call me%init(cd%t0_qist, cd%tf_qist, &
                     &  trim(adjustl(resample_filepath)), &
                     &  cd%traj_id, &
                     &  cd%central_body_id, &
                     &  cd%body_list(:n_bodies), &
                     &  cd%central_body_mu, &
                     &  cd%central_body_ref_radius, &
                     &  cd%mu_list(:n_bodies), &
                     &  rails, &
                     &  shgrav, &
                     &  rot, &
                     &  C, S, &
                     &  cd%regularize, &
                     &  cd%kvtau_filename &
                     )
        me%filepath_var = cd%qist_filename
        me%dynmod%tof = tof
        me%rtol = cd%rtol_qist
        me%atol = cd%atol_qist
    end subroutine gq_namelist_init
    subroutine make_qist_model(namefile)
        character(len=*), intent(in)  :: namefile
        type(gqist)                   :: qist_i
        type(odesolution)             :: qist_sol, kvtau
        type(lightsol)                :: packedsol
        integer                       :: stat, num
        real(qp)                      :: kfinal
        call qist_i%init(namefile,.true.)
        print *, "Integrating QIST. . ."
        ! Integrate model
        if (qist_i%dynmod%regularize) then
            call qist_i%time_regularized_integrate(kvtau, kfinal)
            qist_sol = qist_i%integrate(qist_i%t0, qist_i%tf,kfinal)
            if (qist_i%kvtau_filename=="") then
                print *, "Warning: time regularization filename not set"
            end if
        else
            qist_sol = qist_i%integrate(qist_i%t0, qist_i%tf)
        end if
        print *, "DONE"
        if (qist_i%dynmod%regularize) then
            print *, "Writing k vs tau to file ", &
                     trim(adjustl(qist_i%kvtau_filename)), &
                     ". . ."
            open(newunit=num, file=qist_i%kvtau_filename, iostat=stat, &
                 access="stream", status="replace")
                call kvtau%write(num,dp)
            close(num)
            print *, "DONE"
        end if
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
                  & shgrav, rails, rot, Cbar, Sbar,reg,kvtau_filename)
        ! init_dm: method to initialize dynamicsModel object
        ! INPUTS:
        ! NAME            TYPE           DESCRIPTION
        ! t0              real           lowest possible time in seconds past 
        !                                J2000 at which model integration may 
        !                                start
        ! tf              real           highest possible time in seconds past 
        !                                J2000 at which model integration may 
        !                                end
        ! subspicefile    character      The absolute directory path for the 
        !                                spice_subset object
        ! traj_id         integer        the integer ID of the reference 
        !                                trajectory in the SPICE kernel
        ! central_body    integer        the integer ID of the selected central
        !                                body in the SPICE kernel
        ! bodylist        integer (:)    list of integer IDs of the gravitating
        !                                bodies in the SPICE kernel, other than
        !                                the central body
        ! central_body_mu . . .
        !                  real          gravitational parameter of central body
        ! central_body_ref_radius . . . 
        !                  real          reference radius for SH data
        ! mu_list          real          gravitational parameters of bodies
        !                                sorted in the same order ad bodylist
        ! shgrav          logical        TRUE if central body will be modeled
        !                                by spherical harmonics
        ! rails           logical        is the reference state on rails
        ! check           logical        whether to compare integrated and rails state
        !                                during integration
        ! rot             rothist        the central body-fixed rotation matrix 
        !                                interpolation. Must be R_IU: Take
        !                                vectors FROM the body-fixed (U) frame
        !                                TO the inertial frame.
        ! Cbar            real (:,:)     4 pi (Kaula) normalized cosine Stokes 
        !                                for central body
        ! Sbar            real (:,:)     4 pi (Kaula) normalized sine Stokes 
        !                                for central body
        ! reg            logical         OPTIONAL: whether to use regularization
        ! kvtau_filename character       OPTIONAL: The absolute directory path for the 
        !                                k vs tau fit
        ! OUTPUTS:
        ! NONE
        class(gqist),      intent(inout)        :: me
        type(spice_subset)                      :: subspice
        type(rothist),     intent(in),              optional :: rot
        logical,           intent(in),              optional :: reg
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
        character(len=*),  intent(in),              optional :: kvtau_filename

        open(file=trim(adjustl(subspicefile)),unit=73, &
           & access="stream", status="old")
        call subspice%read(73)
        close(73)
        if (shgrav.and.present(rot).and.present(Cbar).and.present(Sbar)) then
        call me%dynmod%init(subspice, traj_id, central_body, bodylist, &
                          & central_body_mu, central_body_ref_radius,  &
                          & mu_list, shgrav, rails, rot, Cbar, Sbar,reg)
        else
        call me%dynmod%init(subspice, traj_id, central_body, bodylist, &
                          & central_body_mu, central_body_ref_radius,  &
                          & mu_list, shgrav, rails,reg=reg)
        endif
        me%t0 = t0
        me%tf = tf
        me%rtol = 1.e-10
        me%atol = 1.e-12
        me%kvtau_filename = trim(adjustl(kvtau_filename))
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
    function integrate(thisqist, t0, tf,kfinal) result(res)
        ! integrate: integrate a QIST model
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! t0             real           integration start time in seconds past 
        !                               J2000 
        ! tf             real           integration stop time in seconds past 
        !                               J2000
        ! kfinal         real           OPTIONAL: the final regularized independent variable
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            ODESolution    the dense solution, ready for packing
        class(gqist), intent(inout)    :: thisqist
        type(ODESolution)              :: res
        real(qp), intent(in)           :: t0, tf
        real(qp), intent(in), optional :: kfinal
        real(qp)                       :: eye(8,8), hess_init(8**3), stopvar
        real(qp), allocatable          :: init_state(:)
        integer i
        eye = 0._qp
        do i=1,8
            eye(i,i) = 1._qp
        end do
        hess_init = 0._qp
        thisqist%dynmod%tof = tf-t0
        if (present(kfinal)) then
            stopvar = kfinal
            thisqist%dynmod%regularize=.true.
        else
            stopvar = 1._qp
            thisqist%dynmod%regularize=.false.
        end if
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
                      & [0._qp, stopvar], &
                      & [ init_state, &
                         reshape(eye,[8**2]), &
                         hess_init], &
                      & dense_output=.true.,&
                      & rtol=thisqist%rtol, &
                      & atol=thisqist%atol, &
                      & istep=stopvar*0.5_qp)
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
    subroutine time_regularized_integrate(thisqist,kvtau,kfinal)
        class(gqist),      intent(inout)    :: thisqist
        type(odesolution), intent(out)      :: kvtau
        real(qp),          intent(out)      :: kfinal

        ! dtau = c*r**alpha d k
        kvtau = solve_ivp(time_eom,&
                      & [0._qp, 1._qp], &
                      & [0._qp] , &
                      & dense_output=.true.,&
                      & rtol=1.e-16_qp, &
                      & atol=1.e-20_qp, istep=0.5_qp)
        kfinal = kvtau%ys(1,size(kvtau%ts))
        contains
            function time_eom(me, x, y) result(res)
                class(RungeKutta), intent(inout) :: me
                real(qp),          intent(in)    :: x, y(:)
                real(qp)                         :: res(size(y))
                real(qp)                         :: state(6), r,tau, t, kprime, tof
                tof = (thisqist%tf-thisqist%t0)
                t = x*tof + thisqist%t0 ! physical time
                tau = x
                state = thisqist%dynmod%trajstate(t)
                r = sqrt(sum(state(:3)**2))
                kprime = r**(-1.5_qp)
                res(1) = kprime
            end function time_eom
    end subroutine time_regularized_integrate
end module genqist
