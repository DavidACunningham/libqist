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
                  & shgrav, Cbar, Sbar)
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
        ! OUTPUTS:
        ! NONE
        class(gqist), intent(inout) :: me
        type(spice_subset)          :: subspice
        real(qp),     intent(in)    :: t0, tf
        integer,              intent(in)    :: traj_id, & 
                                               central_body, &
                                               bodylist(:)
        logical,              intent(in)    :: shgrav
        real(dp),             intent(in)    :: central_body_ref_radius, &
                                               central_body_mu, &
                                               mu_list(:)
        real(dp),             intent(in)    :: Cbar(:,:), &
                                               Sbar(:,:)
        character(len=*),   intent(in)    :: subspicefile

        open(file=trim(adjustl(subspicefile)),unit=73, &
           & access="stream", status="old")
        call subspice%read(73)
        close(73)
        call me%dynmod%init(subspice, traj_id, central_body, bodylist, &
                          & central_body_mu, central_body_ref_radius,  &
                          & mu_list, shgrav, Cbar, Sbar)
        me%t0 = t0
        me%tf = tf
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
        boundscheck = .false.
        if (t0<=meqist%t0) then
            write(*,42,advance='no') "integration start set to ", t0
            write(*,42) " but must be greater than ", meqist%t0
            boundscheck = .true.
        end if
        if (t0<=meqist%t0) then
            write(*,42,advance='no') "integration end set to ", tf
            write(*,42) " but must be less than ", meqist%tf
            boundscheck = .true.
        end if
        42 format (a,f5.3)
        if (boundscheck) then
            error stop
        end if
        res = solve_ivp(myint_eoms,&
                      & [t0, tf], &
                      & meqist%dynmod%trajstate(t0),&
                      & method="DOP853",&
                      & dense_output=.true.,&
                      & rtol=meqist%rtol*1.E-2_qp, &
                      & atol=meqist%atol, istep=1._qp)
        contains
            function myint_eoms(me, x, y) result(res)
                class(RungeKutta), intent(inout) :: me
                real(qp),          intent(in)    :: x, y(:)
                real(qp)                         :: res(size(y))
                res = meqist%dynmod%eoms(x,y)
            end function myint_eoms
        end function integrate
end module genqist
