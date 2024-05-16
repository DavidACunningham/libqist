module genqist
    use globals
    use denselight, only: lightsol, qistpack
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    implicit none

    type gqist
        type(dynamicsmodel)   :: dynmod
        real(qp), allocatable :: initstate(:)
        real(qp)              :: rtol, atol, t0, tf
        contains
        procedure init
        procedure integrate
        procedure :: pack => packsol
        procedure :: bounds
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
        inquire(file=namefile, iostat=stat)
        if (stat .ne. 0) then 
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
                           adjustl(trim(metakernel_filepath)), & ! Spice kernel
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

    subroutine make_qist_model(namefile)
        character(len=*), intent(in)  :: namefile
        type(gqist)         :: qist_i
        type(odesolution)   :: qist_sol
        type(lightsol)      :: packedsol
        character(len=1000) :: resample_filepath, qist_filepath
        real(qp)            :: t0, tf, tof, &
                                rtol, atol
        integer             :: reference_trajectory_id, & 
                               central_body_id, &
                               body_list(30)
        logical             :: shgrav
        integer             :: stat, num, n_bodies
        real(qp)            :: central_body_ref_radius, &
                               central_body_mu, & 
                               mu_list(30)
        real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                                Sbar(2,2) = 0._qp
        logical              :: rails
        namelist /QIST_CONFIG/ resample_filepath, &
                               qist_filepath, &
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
        ! Read in namelist
        inquire(file=namefile, iostat=stat)
        if (stat .ne. 0) then 
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
        ! Set number of bodies
        n_bodies = findloc(body_list,0,dim=1)-1
        tof = tf
        ! Print status
        print *, "Generating QIST model for reference trajectory ", reference_trajectory_id
        print *, "Using ephemeris data from ", trim(adjustl(resample_filepath))
        print *, "Relative to body ", central_body_id
        print *, "With perturbations from bodies ", body_list(:n_bodies)
        print *, "From J2000 + ", t0
        print *, "To J2000 + ", tf
        print *, ""
        ! Initialize QIST
        call qist_i%init(t0, tf, &
                     &  trim(adjustl(resample_filepath)), &
                     &  reference_trajectory_id, &
                     &  central_body_id, &
                     &  body_list(:n_bodies), &
                     &  central_body_mu, &
                     &  central_body_ref_radius, &
                     &  mu_list(:n_bodies), &
                     &  shgrav, &
                     &  Cbar, Sbar,rails)
        qist_i%dynmod%tof = tof
        print *, "Integrating QIST. . ."
        qist_i%rtol = rtol
        qist_i%atol = atol
        ! Integrate model
        qist_sol = qist_i%integrate(t0, tf)
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
        print *, "Writing QIST model to file ", trim(adjustl(qist_filepath)), ". . ."
        open(newunit=num, file=trim(adjustl(qist_filepath)), iostat=stat, &
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

    subroutine init(me,t0, tf, subspicefile, traj_id, central_body, bodylist, &
                  & central_body_mu, central_body_ref_radius, mu_list, &
                  & shgrav, Cbar, Sbar,rails)
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
        ! Cbar           real (:,:)     4 pi (Kaula) normalized cosine Stokes 
        !                               for central body
        ! Sbar           real (:,:)     4 pi (Kaula) normalized sine Stokes 
        !                               for central body
        ! rails          logical        is the reference state on rails
        ! check          logical        whether to compare integrated and rails state
        !                               during integration
        ! OUTPUTS:
        ! NONE
        class(gqist), intent(inout) :: me
        type(spice_subset)          :: subspice
        real(qp),     intent(in)    :: t0, tf
        integer,              intent(in)    :: traj_id, & 
                                               central_body, &
                                               bodylist(:)
        logical,              intent(in)    :: shgrav, rails
        real(qp),             intent(in)    :: central_body_ref_radius, &
                                               central_body_mu, &
                                               mu_list(:)
        real(qp),             intent(in)    :: Cbar(:,:), &
                                               Sbar(:,:)
        character(len=*),   intent(in)    :: subspicefile

        open(file=trim(adjustl(subspicefile)),unit=73, &
           & access="stream", status="old")
        call subspice%read(73)
        close(73)
        call me%dynmod%init(subspice, traj_id, central_body, bodylist, &
                          & central_body_mu, central_body_ref_radius,  &
                          & mu_list, shgrav, Cbar, Sbar,rails)
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
        logical                     :: boundscheck
        real(dp),     intent(out)   :: testpoints(6,npoints), &
                                       spicepoints(6,npoints), &
                                       testtimes(npoints)
        type(ODESolution)           :: testsol
        integer i

        me_qist%dynmod%tgt_on_rails = .false.
        ! Change bodylist and SHGrav Parameters
        call me_qist%dynmod%new_bodies(bodylist, mulist)
        call me_qist%dynmod%new_gravstatus(shgrav, Cbars,Sbars)
        boundscheck = me_qist%bounds(real(t0,qp), real(tf,qp))
        me_qist%dynmod%tof = tf-t0
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
                real(qp)                         :: acc(6), jac(6,6), hes(6,6,6)
                real(qp)                         :: res(size(y))
                call me_qist%dynmod%get_derivs(x,acc,jac,hes)
                res = acc
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
        logical                     :: boundscheck
        real(qp)                    :: eye(8,8), hess_init(8**3)
        real(qp), allocatable       :: init_state(:)
        integer i
        boundscheck = thisqist%bounds(t0,tf)
        if (boundscheck) then
            error stop
        end if
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
        thisqist%dynmod%state = [thisqist%dynmod%trajstate(t0), t0, tf-t0]
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
    function bounds(me, t0, tf) result(res)
        class(gqist), intent(in) :: me
        real(qp),     intent(in) :: t0, tf
        logical                  :: res
        res = .false.
        if (t0<me%t0) then
            write(*,42,advance='no') "integration start set to ", t0
            write(*,42) " but must be greater than ", me%t0
            res = .true.
        end if
        if (tf>me%tf) then
            write(*,42,advance='no') "integration end set to ", tf
            write(*,42) " but must be less than ", me%tf
            res = .true.
        end if
        42 format (a,f5.3)
    end function bounds
end module genqist
