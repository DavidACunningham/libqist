module makemodel
    use globals
    use frkmain, only: solve_ivp, Odesolution
    use astkindmodule
    use shdatamodule
    use pinesmodule
    implicit none

    type :: dynamicsModel
        integer               :: num_bodies, &
                               & central_body, &
                               & traj_id
        integer,  allocatable :: bodylist(:) ! num_bodies
        logical               :: shgrav
        real(dp)              :: central_body_ref_radius, &
                               & central_body_mu, &
                               & central_body_rad(3)
        real(wp), allocatable :: nbody_mus(:), &      ! (num_bodies)
                               & nbody_radii(:,:)     ! (3,num_bodies)
        type(SHarmData)       :: shdat
        type(PinesData)       :: pdat
        contains
            procedure :: get_derivs
            procedure :: acc_kepler
            procedure :: jac_kepler
            procedure :: hes_kepler
            procedure :: acc_nbody
            procedure :: jac_nbody
            procedure :: hes_nbody
            procedure :: allderivs_sh
    end type 

    contains

    subroutine init_dm(me, traj_id, central_body, bodylist, shgrav, Cbar, Sbar)
        class(dynamicsModel), intent(inout) :: me
        integer,              intent(in)    :: traj_id, & 
                                               central_body, &
                                               bodylist(:)
        logical,              intent(in)    :: shgrav
        real(wp),             intent(in)    :: Cbar(:,:), &
                                               Sbar(:,:)
        real(dp)                            :: mudum(1)
        real(dp)                            :: radii(3)
        integer                             :: i,N
        me%num_bodies = size(bodylist)
        me%shgrav = shgrav
        me%traj_id = traj_id
        allocate(me%bodylist(me%num_bodies),  &
                 me%nbody_mus(me%num_bodies), &
                 me%nbody_radii(3,me%num_bodies) &
                )
        me%bodylist = bodylist
        call bodvcd(central_body, 'RADII', 3, N, radii)
        call bodvcd(central_body, 'GM', 1, N, mudum)
        me%central_body_ref_radius = radii(1)
        me%central_body_mu = mudum(1)
        do i=1,me%num_bodies
            call bodvcd(bodylist(i), 'GM', 1, N, mudum)
            me%nbody_mus(i) = real(mudum(1), wp)
        end do

        if (me%shgrav) then
            call shdat_from_table( me%central_body_ref_radius, & 
                                 & me%central_body_mu, &
                                 & real(Cbar,dp), real(Sbar,dp), &
                                 & size(Cbar,1),size(Cbar,2), me%shdat &
                                 )
            call pinesinit(size(Cbar,1),me%shdat%Cml,me%shdat%Sml,me%pdat)
        endif
    end subroutine init_dm
    subroutine get_derivs(me, time, acc, jac, hes)
        class(dynamicsModel), intent(inout) :: me
        real(wp),             intent(in)    :: time
        real(wp),             intent(out)   :: acc(3), jac(3,3), hes(3,3,3)
        real(dp)                            :: traj_state(6), lt_dum
        integer i
        ! at time t (in s past j2000)
            ! For central body (i=1) NOTE: SPICE IS NOT THREAD SAFE, DON'T DO THIS IN ||
                ! get r from traj to body i (SPKGPS)
        call spkgeo( me%traj_id, time, "J2000", me%central_body, &
                   & traj_state, lt_dum)
            ! for bodies i = 2. . . N NOTE: SPICE IS NOT THREAD SAFE, DON'T DO THIS IN ||
        do i=1,me%num_bodies
            ! get r from traj to body i (SPKGPS)
            call spkgps( me%traj_id, time, "J2000", me%central_body, &
                               & me%nbody_radii(:,i), &     ! (3,num_bodies)
                               & traj_state, lt_dum)
            ! get r from central body to body i (SPKGPS)

        end do

            ! For central body (i=1)
                ! Choose point mass or SH model
                ! Compute accel, jac, hess for i=1
            ! for bodies i = 2. . . N NOTE: CAN BE DONE IN ||
                ! Compute point mass accel, jac, hess for i=1

            ! Sum results to get final acc, jac, hes
    end subroutine get_derivs

    function acc_kepler(me, mu, r) result (res)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3)
        real(wp)                         :: res(3)

    end function

    function jac_kepler(me, mu, r) result (res)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3)
        real(wp)                         :: res(3,3)
    end function
    
    function hes_kepler(me, mu, r) result (res)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3)
        real(wp)                         :: res(3,3,3)
    end function

    function acc_nbody(me, mu, r, rbods) result (res)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3), &
                                            rbods(3)
        real(wp)                         :: res(3)

    end function

    function jac_nbody(me, mu, r, rbods) result (res)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3), &
                                            rbods(3)
        real(wp)                         :: res(3,3)
    end function
    
    function hes_nbody(me, mu, r, rbods) result (res)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3), &
                                            rbods(3)
        real(wp)                         :: res(3,3,3)
    end function

    subroutine allderivs_sh(me, time, r, acc, jac, hes)
        class(dynamicsModel), intent(inout)  :: me
        real(wp),             intent(in)     :: time, &
                                              & r(3)
        real(dp)                             :: pot_d, &
                                              & acc_d(3), &
                                              & jac_d(3,3), &
                                              & hes_d(3,3,3)
        real(wp),             intent(out)    :: acc(3), &
                                              & jac(3,3), &
                                              & hes(3,3,3)
        real(dp)                             :: xform_mat(3,3), &
                                              & r_bf(3)

        ! TODO: Call transform to get transformation matrix 
        ! to/from body fixed frame in xform_mat
        r_bf = mmult(xform_mat,r_bf)

        call shpines( me%shdat%Rbody, me%shdat%GM, me%pdat, &
                    & me%shdat%maxdeg, me%shdat%maxorder, &
                    & real(r_bf,dp), pot_d, acc_d, jac_d, hes_d &
                    & )
        acc = real(acc_d,wp)
        jac = real(jac_d,wp)
        hes = real(hes_d,wp)
    end subroutine
end module makemodel
