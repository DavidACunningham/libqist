module genqist
    use globals
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel
    use frkmain, only: solve_ivp, Odesolution, RungeKutta
    use cheby, only: spice_subset
    implicit none

    type gqist
        type(dynamicsmodel)   :: dynmod
        real(qp), allocatable :: initstate(:)
        real(qp)              :: rtol, atol, t0, tf
        logical               :: check
        contains
        procedure init
        procedure integrate
    end type gqist

    contains

    subroutine make_spice_subset(ta,tb, central_body, bod_list, metakernel_filename, output_filename)
        real(dp), intent(in) :: ta,tb
        integer,  intent(in) :: central_body, bod_list(:)
        character(len=*)     :: output_filename, metakernel_filename

    end subroutine make_spice_subset

    subroutine init(me,t0, tf, subspicefile, traj_id, central_body, bodylist, &
                  & central_body_mu, central_body_ref_radius, mu_list, &
                  & shgrav, Cbar, Sbar,rails,check)
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
        logical,              intent(in)    :: shgrav, rails, check
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
        me%rtol = 1.e-14
        me%atol = 1.e-14
        me%check = check
    end subroutine
    function integrate(meqist, t0, tf) result(res)
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
        class(gqist), intent(inout) :: meqist
        type(ODESolution)           :: res
        real(qp), intent(in)        :: t0, tf
        logical                     :: boundscheck
        real(qp)                    :: eye(8,8), hess_init(8**3)
        real(qp), allocatable       :: init_state(:)
        integer i
        boundscheck = .false.
        if (t0<meqist%t0) then
            write(*,42,advance='no') "integration start set to ", t0
            write(*,42) " but must be greater than ", meqist%t0
            boundscheck = .true.
        end if
        if (tf>meqist%tf) then
            write(*,42,advance='no') "integration end set to ", tf
            write(*,42) " but must be less than ", meqist%tf
            boundscheck = .true.
        end if
        42 format (a,f5.3)
        if (boundscheck) then
            error stop
        end if
        eye = 0._qp
        do i=1,8
            eye(i,i) = 1._qp
        end do
        hess_init = 0._qp
        meqist%dynmod%tof = tf-t0
        if (meqist%dynmod%tgt_on_rails) then
            allocate(init_state(1))
            init_state = [t0]
        else
            allocate(init_state(8))
            init_state(:6) = meqist%dynmod%trajstate(t0)
            init_state(7:) = [t0, tf - t0]
        endif
        print *, init_state
        res = solve_ivp(myint_eoms,&
                      & [0._qp, 1._qp], &
                      & [ init_state, &
                         reshape(eye,[8**2]), &
                         hess_init], &
                      & dense_output=.true.,&
                      & rtol=meqist%rtol, &
                      & atol=meqist%atol, istep=0.5_qp)
        contains
            function myint_eoms(me, x, y) result(res)
                class(RungeKutta), intent(inout) :: me
                real(qp),          intent(in)    :: x, y(:)
                real(qp)                         :: res(size(y))
                if (meqist%dynmod%tgt_on_rails) then
                    res = meqist%dynmod%eoms_rails(x,y)
                else
                    res = meqist%dynmod%eoms(x,y)
                endif
            end function myint_eoms
        end function integrate
end module genqist
