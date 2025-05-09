! Title: makemodel.f90 
! Description:
!     A module for managing gravitational dynamics based on ephemeris data.
!     The main purpose is the computation of equations of motion of spacecraft
!     trajectories and state transition tensors. Additional machinery is
!     required in order to accomodate trajectory information stored in 
!     Chebyshev interpolants.
!
! References:
!   None
! 
! author: David Cunningham
! Last edited: See git log
module makemodel
    use, intrinsic :: iso_fortran_env, only: wp=>real64, dp=>real64, qp=>real128
    use tinysh, only: PinesData, pinesinit, shpines
    use tensorops, only: mattens, quad, vectens3, eyemat, vectensquad
    use subspice, only: spice_subset
    use quat, only: rothist
    use globals, only: mmult
    implicit none
    integer :: m=8

    type :: dynamicsModel
        ! dynamicsModel:
        !   derived type containing dynamics model information
        !   and the ability to compute the derivative of the state
        !   along with the Jacobian and Hessian of that state
        !   for ephemeris gravity as well as spherical harmonics
        !   gravity for the central body.
        !
        ! Attributes:
        !    bod_db: spice_subset
        !       The database containing reinterpolated ephemerides
        !       of the celestial bodies and spacecraft required
        !       for computing the trajectory derivatives.
        !    rot: rothist
        !       The interpolated rotation and rotational derivatives
        !       needed for computing spherical harmonics accelerations,
        !       Jacobians, and Hessians
        !   num_bodies: integer
        !       The number of bodies (other than the central body)
        !       gravitationally acting on the reference trajectory
        !   central_body: integer
        !       NAIF integer ID code of the central body for integration
        !   traj_id: integer
        !       NAIF integer ID code of the reference trajectory
        !   degord: integer
        !       The degree and order of the spherical harmonics model
        !       used for computing the gravity of the central body
        !   bodylist: integer (num_bodies)
        !       NAIF integer ID codes of the gravitating bodies
        !       other than the central body
        !   shgrav: logical
        !       Is the central body gravitation using spherical harmonics
        !       (TRUE) or keplerian (FALSE)
        !   tgt_on_rails: logical
        !       Is the reference trajectory precomputed from a SPICE SPK
        !       (TRUE) or being integrated along with the STM/STT (FALSE)
        !   regularize: logical
        !       Is the independent variable sundman regularized relative
        !       to the central body (TRUE) or not (FALSE)
        !   central_body_ref_radius: real
        !       The reference radius of the central body in km
        !   state: real (8)
        !       The augmented dynamical state of the reference trajectory:
        !       state = [x, y, z, xdot, ydot, zdot, t, TOF]
        !   central_body_mu: real
        !       The gravitational parameter of the central body in km3/s2
        !   tof: real
        !       The time of flight of the trajectory in seconds
        !   nbody_mus: real (num_bodies)
        !       The gravitational parameter of the perturbing bodies in km3/s2
        !   pdat: PinesData
        !       The derived type for storing the spherical harmonics computation
        !       data
        type(spice_subset)    :: bod_db
        type(rothist)         :: rot
        integer               :: num_bodies, &
                               & central_body, &
                               & traj_id, &
                               & degord
        integer,  allocatable :: bodylist(:) ! num_bodies
        logical               :: shgrav, &
                               & tgt_on_rails, &
                               & regularize
        real(qp)              :: central_body_ref_radius, &
                               & state(8), &
                               & central_body_mu, &
                               & tof
        real(qp), allocatable :: nbody_mus(:)     ! (num_bodies)
        type(PinesData)       :: pdat
        contains
            procedure:: init => init_dm
            procedure :: eoms
            procedure :: eoms_rails
            procedure :: get_derivs
            procedure :: allderivs_kepler
            procedure :: allderivs_thirdbody
            procedure :: allderivs_sh
            procedure :: trajstate
            procedure :: new_bodies
            procedure :: new_gravstatus
    end type dynamicsModel
    contains
    subroutine init_dm(me, subspice, traj_id, central_body, bodylist, &
                     & central_body_mu, central_body_ref_radius, mu_list, &
                     & shgrav, rails,rot, Cbar, Sbar,reg)
        ! init_dm: method to initialize dynamicsModel object
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! subspice       spice_subset   spice subset object containing 
        !                               necessary objects for the model--
        !                               the traj_id, central_body, and bodylist
        !                               must all be members of 
        !                               spice_subset%bodlist
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
        ! tgt_on_rails   logical        TRUE if target trajectory should not
        !                               be integrated
        ! rot            rothist        OPTIONAL: the central body-fixed 
        !                               rotation matrix interpolation. 
        !                               Must be R_IU: Take
        !                               vectors FROM the body-fixed (U) frame
        !                               TO the inertial frame.
        ! Cbar           real (:,:)     OPTIONAL: 4 pi (Kaula) normalized 
        !                               cosine Stokes coefficients for 
        !                               central body
        ! Sbar           real (:,:)     OPTIONAL: 4 pi (Kaula) normalized 
        !                               sine Stokes coefficients for
        !                               central body
        ! reg           logical         OPTIONAL: whether to use regularization
        ! OUTPUTS:
        ! NONE
        class(dynamicsModel), intent(inout)           :: me
        type(spice_subset),   intent(in)              :: subspice
        type(rothist),        intent(in), optional    :: rot
        integer,              intent(in)              :: traj_id, & 
                                                         central_body, &
                                                         bodylist(:)
        real(qp),             intent(in)              :: central_body_ref_radius, &
                                                         central_body_mu, &
                                                         mu_list(:)
        logical,              intent(in)              :: shgrav, rails
        logical,              intent(in), optional    :: reg
        real(qp),             intent(in), optional    :: Cbar(:,:), &
                                                         Sbar(:,:)
        real(dp),             allocatable             :: Cbaralloc(:,:), &
                                                         Sbaralloc(:,:)
        me%num_bodies = size(bodylist)
        me%shgrav = shgrav
        me%traj_id = traj_id
        me%bod_db = subspice
        me%tgt_on_rails = rails
        allocate(me%bodylist(me%num_bodies),  &
                 me%nbody_mus(me%num_bodies))
        me%bodylist = bodylist
        me%central_body = central_body
        me%central_body_ref_radius = real(central_body_ref_radius,qp)
        me%central_body_mu = real(central_body_mu,qp)
        me%nbody_mus = real(mu_list,qp)

        if (me%shgrav) then
            if (present(rot)) then
                me%rot = rot
            else
                print *, "Error: rotation history required. None provided"
            endif
            if (present(Cbar)) then
                me%degord = size(Cbar,1)
            else
                print *, "Error: Stokes coefficients (Cbar) required."
                print *, "None provided."
                error stop
            end if
            if (present(Sbar)) then
                allocate(Cbaralloc(me%degord,me%degord))
                allocate(Sbaralloc(me%degord,me%degord))
                Cbaralloc=transpose(real(Cbar,dp))
                Sbaralloc=transpose(real(Sbar,dp))
                call pinesinit(me%degord,Cbaralloc, &
                             & Sbaralloc,me%pdat)
            else
                print *, "Error: Stokes coefficients (Sbar) required."
                print *, "None provided."
                error stop
            end if
        endif
        if (present(reg)) then
            me%regularize = reg
        else
            me%regularize = .false.
        endif
    end subroutine init_dm
    subroutine new_bodies(me, bodylist, mulist)
        ! new_bodies: method to update list of bodies in use
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! me             dynamicsModel  calling instance
        ! bodylist       integer (:)    list of integer IDs of the gravitating
        !                               bodies in the SPICE kernel, other than
        !                               the central body
        ! mu_list        real           gravitational parameters of bodies
        !                               sorted in the same order as bodylist
        !                               in km3/s2
        ! OUTPUTS:
        ! NONE
        class(dynamicsModel), intent(inout) :: me
        integer,              intent(in)    :: bodylist(:)
        real(dp),             intent(in)    :: mulist(:)
        if (allocated(me%bodylist)) deallocate(me%bodylist)
        if (allocated(me%nbody_mus)) deallocate(me%nbody_mus)
        me%num_bodies = size(bodylist)
        allocate(me%bodylist(me%num_bodies),  &
                 me%nbody_mus(me%num_bodies))
        me%bodylist = bodylist
        me%nbody_mus = real(mulist,qp)
    end subroutine new_bodies
    subroutine new_gravstatus(me, shgrav, Cbar, Sbar,rot)
        ! new_gravstatus: method to update central body gravity model
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! me             dynamicsModel  calling instance
        ! shgrav         logical        whether to use spherical harmonics
        ! Cbar           real (double)  nxm array of normalized cosine Stokes 
        !                               coefficients
        ! Sbar           real (double)  nxm array of normalized sine Stokes 
        !                               coefficients
        ! rot            rothist        transformation from rotating to inertial
        !                               reference frames
        ! OUTPUTS:
        ! NONE
        class(dynamicsModel), intent(inout)        :: me
        type(rothist),        intent(in), optional :: rot
        logical,              intent(in)           :: shgrav
        real(dp),             intent(in)           :: Cbar(:,:), &
                                                    & Sbar(:,:)
        real(dp), allocatable                      :: Cbaralloc(:,:), &
                                                    & Sbaralloc(:,:)
        me%shgrav = shgrav
        if (me%shgrav) then
            me%degord = size(Cbar,1)
            allocate(Cbaralloc,Sbaralloc,mold=Cbar)
            Cbaralloc = transpose(Cbar)
            Sbaralloc = transpose(Sbar)
            call pinesinit(me%degord,Cbaralloc, &
                         & Sbaralloc,me%pdat)
            if (present(rot)) then
                me%rot = rot
            end if
        endif
    end subroutine new_gravstatus
    subroutine get_derivs(me, time, acc, jac, hes)
        ! get_derivs: method to generate dynamics and dynamics partials around 
        !             reference trajectory
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! time           real           time in sec past J2000 at which
        !                               partials are evaluated
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! acc            float (6)      time derivative of 6-d state,
        !                               i.e. the dynamics for the system:
        !                               (xdot, ydot, zdot, xddot, yddot, zddot)
        ! jac            float (6,6)    jacobian of dynamics: 
        !                               jac(i,j) = df_i/dx_j
        ! hes            float (6,6,6)  hessian of dynamics: 
        !                               hes(i,j,k) = d^2f_i/(dx_j*dx_k)
        class(dynamicsModel), intent(inout) :: me
        real(qp),             intent(in)    :: time
        real(qp),             intent(out)   :: acc(m), jac(m,m), hes(m,m,m)
        real(qp)                            :: acc_2b(m), &
                                             & jac_2b(m,m), &
                                             & hes_2b(m,m,m), &
                                             & acc_nb(m), &
                                             & jac_nb(m,m), &
                                             & hes_nb(m,m,m), &
                                             & acc_thisbody(m), &
                                             & jac_thisbody(m,m), &
                                             & hes_thisbody(m,m,m), &
                                             & nbody_radii(3,me%num_bodies), &
                                             & nbody_vels(3,me%num_bodies), &
                                             & nbody_accs(3,me%num_bodies), &
                                             & y(m), pos(3), r
        integer i
        ! at time t (in s past j2000)
            ! get r from central body to traj
        if (me%tgt_on_rails) then
            y(:6) = me%trajstate(time)
            y(7) = time
            y(8) = me%tof
        else
            y = me%state
        endif
        ! for bodies i = 2. .
        ! if num_bodies==0 loop will not run (Fortran standard behavior)
        do i=1,me%num_bodies 
            ! get r from central body to body i
            nbody_radii(:,i) = real(me%bod_db%call(real(time,dp), &
                                               & me%bodylist(i),'p'), qp)
            nbody_vels(:,i) = real(me%bod_db%call(real(time,dp), &
                                               & me%bodylist(i),'v'), qp)
            nbody_accs(:,i) = real(me%bod_db%call(real(time,dp), &
                                               & me%bodylist(i),'a'), qp)
        end do
        ! For central body
            ! Choose point mass or SH model
        if (me%shgrav) then
            call me%allderivs_sh(time, y, acc_2b, jac_2b, hes_2b)
        else
            call me%allderivs_kepler(me%central_body_mu, y, &
                                   & acc_2b, &
                                   & jac_2b, &
                                   & hes_2b &
                                    )
        endif
        ! for bodies i = 2. . . N NOTE: CAN BE DONE IN ||
        ! Compute point mass accel, jac, hess for i=1
        acc_nb = 0._qp
        jac_nb = 0._qp
        hes_nb = 0._qp
        do i=1,me%num_bodies
            call me%allderivs_thirdbody(me%nbody_mus(i), y, &
                                      & nbody_radii(:,i), &
                                      & nbody_vels(:,i), &
                                      & nbody_accs(:,i), &
                                      & acc_thisbody, &
                                      & jac_thisbody, &
                                      & hes_thisbody &
                                       )
            acc_nb = acc_nb + acc_thisbody
            jac_nb = jac_nb + jac_thisbody
            hes_nb = hes_nb + hes_thisbody
        end do
        acc = acc_2b
        acc(4:6) = acc_2b(4:6) + acc_nb(4:6)
        jac = jac_2b
        jac(4:6,:) = jac(4:6,:) + jac_nb(4:6,:)
        hes = hes_2b
        hes(4:6,:,:) = hes(4:6,:,:) + hes_nb(4:6,:,:)
        !! EXTRA REGULARIZATION
        if (me%regularize) then
            pos = y(:3)
            r = sqrt(sum(pos**2))
            acc = acc*r**(1.5_qp)
            jac = jac*r**(1.5_qp)
            hes = hes*r**(1.5_qp)
        end if
        contains
            function sme(y) result(res)
                real(qp), intent(in) :: y(:)
                real(qp)             :: res
                real(qp)             :: v, r
                r = norm2(y(:3))
                v = norm2(y(4:6))
                res = v**2/2._qp - me%central_body_mu/r
            end function sme
    end subroutine get_derivs
    subroutine allderivs_kepler(me, mu, y, acc, jac, hes)
        ! acc_kepler: method to compute Keplerian derivatives
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! y              real (:)       Augmented dynamical state vector
        !                               field point in first 3 elements,
        !                               could be state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8)       Augmented acceleration vector
        class(dynamicsModel), intent(in)  :: me
        real(qp),             intent(in)  :: mu
        real(qp),             intent(in)  :: y(:)
        real(qp),             intent(out) :: acc(size(y)), &
                                          & jac(size(y), size(y)), &
                                          & hes(size(y),size(y),size(y))
        real(qp)                          :: pos(3), rdotr, norm, onebyr, &
                                           & onebyr2, onebyr3, onebyr5, &
                                           & onebyr7, &
                                           F(3), J(3,3), H(3,3,3), &
                                           eye(3,3), rr(3,3), rrr(3,3,3), &
                                           ridjk(3,3,3), rjdki(3,3,3), &
                                           rkdij(3,3,3)
        integer i
        eye = 0._qp
        do i=1,3
            eye(i,i) = 1._qp
        end do
        pos = y(:3)
        rdotr = sum(pos**2)
        norm = sqrt(rdotr)
        onebyr = 1._qp/norm
        onebyr2 = onebyr*onebyr
        onebyr3 = onebyr2*onebyr
        onebyr5 = onebyr3*onebyr2
        onebyr7 = onebyr5*onebyr2
        rr = outer_2vec(pos,pos)
        rrr = outer_3vec(pos,pos,pos)
        ridjk = outer_kron(pos,1)
        rjdki = outer_kron(pos,2)
        rkdij = outer_kron(pos,3)
        F = -mu*pos*onebyr3
        J = mu*(3._qp*onebyr5*rr - onebyr3*eye)
        H = 3._qp*mu*(-5._qp*onebyr7*rrr + onebyr5*(ridjk + rjdki + rkdij))
        acc = 0._qp
        jac = 0._qp
        hes = 0._qp
        !! Assign accel
        acc(:3) = y(4:6)
        acc(4:6) = F
        acc(7) = 1._qp
        acc(8) = 0._qp
        !! Regularization
        acc = acc*y(8)
        !! Jacobian assignment
        jac(1:3,4:6) = y(8) * eye
        jac(1:3,8) = y(4:6)
        jac(4:6,1:3) = y(8) * J
        jac(4:6,8) = F
        jac(7,8) = 1._qp
        !! Assign Hessian
        hes(4:6,1:3,1:3) = y(8)*H
        hes(4:6,8,1:3) = J
        hes(4:6,1:3,8) = J
        hes(1:3,8,4:6) = eye
        hes(1:3,4:6,8) = eye
    end subroutine
    subroutine allderivs_thirdbody(me, mu, y, rbods, vbods, abods, &
                                 & acc, jac, hes)
        ! acc_kepler: method to compute Keplerian derivatives
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! y              real (:)       Augmented dynamical state vector
        !                               field point in first 3 elements,
        !                               could be state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8)       Augmented acceleration vector
        class(dynamicsModel), intent(in)  :: me
        real(qp),             intent(in)  :: mu
        real(qp),             intent(in)  :: y(:), &
                                           & rbods(3), vbods(3), abods(3)
        real(qp),             intent(out) :: acc(size(y)), &
                                          & jac(size(y), size(y)), &
                                          & hes(size(y),size(y),size(y))
        real(qp)                          :: pos(3), rsp(3), &
                                             rspdotrsp, rspnorm, &
                                           & onebyrsp, onebyrsp2, &
                                           & onebyrsp3, onebyrsp5, &
                                           & onebyrsp7, &
                                           & rpdotrp, rpnorm, &
                                           & onebyrp, onebyrp2, &
                                           & onebyrp3, onebyrp5, &
                                           & onebyrp7, &
                                           & F(3), J(3,3), H(3,3,3), &
                                           & rsprsp(3,3), eye(3,3), &
                                           & rsprsprsp(3,3,3), &
                                           & dijrk(3,3,3), dkirj(3,3,3), &
                                           & djkri(3,3,3), abyt(3), &
                                           & abytbyt(3), Jbyt(3,3), &
                                           & rspdotvp, rpdotvp, &
                                           & vpdotvp, rspdotap, &
                                           & rpdotap, vprsp(3,3), rspvp(3,3)
        integer i
        acc = 0._qp
        jac = 0._qp
        hes = 0._qp
        eye = 0._qp
        do i=1,3
            eye(i,i) = 1._qp
        end do
        pos = y(:3)
        rsp = rbods - pos
        rspdotrsp = sum(rsp**2)
        rspnorm = sqrt(rspdotrsp)
        onebyrsp = 1._qp/rspnorm
        onebyrsp2 = onebyrsp*onebyrsp
        onebyrsp3 = onebyrsp2*onebyrsp
        onebyrsp5 = onebyrsp3*onebyrsp2
        onebyrsp7 = onebyrsp5*onebyrsp2

        rpdotrp = sum(rbods**2)
        rpnorm = sqrt(rpdotrp)
        onebyrp = 1._qp/rpnorm
        onebyrp2 = onebyrp*onebyrp
        onebyrp3 = onebyrp2*onebyrp
        onebyrp5 = onebyrp3*onebyrp2
        onebyrp7 = onebyrp5*onebyrp2

        rsprsp = outer_2vec(rsp,rsp)
        rsprsprsp = outer_3vec(rsp,rsp,rsp)
        djkri = outer_kron(rsp,1)
        dkirj = outer_kron(rsp,2)
        dijrk = outer_kron(rsp,3)

        F = mu*(rsp*onebyrsp3 - rbods*onebyrp3)
        J = mu*(3._qp*onebyrsp5*rsprsp - eye*onebyrsp3)
        H = 3._qp*mu*(5*onebyrsp7*rsprsprsp &
                    - onebyrsp5*(dijrk + dkirj + djkri))

        rspdotvp = sum(rsp*vbods)
        rpdotvp = sum(rbods*vbods)
        vpdotvp = sum(vbods**2)
        rspdotap = sum(rsp*abods)
        rpdotap = sum(rbods*abods)


        abyt = mu*(- 3._qp*rspdotvp*onebyrsp5*rsp &
                   + 3._qp*rpdotvp*onebyrp5*rbods &
                   + vbods*onebyrsp3 &
                   - vbods*onebyrp3 &
                  )
        abytbyt = mu*( &
                    3._qp*(5._qp*onebyrsp7*rspdotvp**2 &
                         - onebyrsp5*(rspdotap + vpdotvp))*rsp &
                   -3._qp*(5._qp*onebyrp7*rpdotvp**2 &
                         - onebyrp5*(rpdotap + vpdotvp))*rbods &
                   -6._qp*(onebyrsp5*(rspdotvp) &
                         - onebyrp5*rpdotvp)*vbods &
                     + (onebyrsp3 - onebyrp3)*abods &
                     )
        vprsp = outer_2vec(vbods,rsp)
        rspvp = outer_2vec(rsp,vbods)
        Jbyt = 3._qp*mu*( &
                        - 5*onebyrsp7*rspdotvp*rsprsp &
                        + onebyrsp5*(vprsp + rspvp + rspdotvp*eye) &
                        )
        !! Assign accel
        acc(:3) = y(4:6)
        acc(4:6) = F
        acc(7) = 1._qp
        acc(8) = 0._qp
        !! Regularization
        acc = acc*y(8)
        !! Jacobian assignment
        jac(1:3,4:6) = y(8) * eye
        jac(1:3,8) = y(4:6)
        jac(4:6,1:3) = y(8) * J
        jac(4:6,7) = y(8) * abyt
        jac(4:6,8) = F
        jac(7,8) = 1._qp
        !! Assign Hessian
        hes(4:6,1:3,1:3) = y(8)*H
        hes(4:6,7,1:3) = y(8) * Jbyt
        hes(4:6,1:3,7) = y(8) * Jbyt
        hes(4:6,8,1:3) = J
        hes(4:6,1:3,8) = J
        hes(1:3,8,4:6) = eye
        hes(1:3,4:6,8) = eye
        hes(4:6,7,7) = y(8) * abytbyt
        hes(4:6,8,7) = abyt
        hes(4:6,7,8) = abyt
    end subroutine
    subroutine allderivs_sh(me, time, y, acc, jac, hes)
        ! allderivs_sh: compute the dynamics, jacobian, and hessian
        !               due to the gravitation of an extended body.
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! time           real           time in seconds past J2000 to evaluate
        ! y              real (8)       state vector from central body to
        !                               field point
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! acc            float (8)      time derivative of 6-d state,
        !                               i.e. the dynamics for the system:
        !                               (xdot, ydot, zdot, xddot, yddot, zddot)
        ! jac            float (8,8)    jacobian of dynamics: 
        !                               jac(i,j) = df_i/dx_j
        ! hes            float (8,8,8)  hessian of dynamics: 
        !                               hes(i,j,k) = d^2f_i/(dx_j*dx_k)
        class(dynamicsModel), intent(inout) :: me
        real(qp),             intent(in)    :: time, y(:)
        real(qp),             intent(out)   :: acc(size(y)), &
                                             & jac(size(y),size(y)), &
                                             & hes(size(y),size(y),size(y))
        real(dp)                            :: RIU(3,3), RIUdot(3,3), RIUdotdot(3,3)
        real(dp)                            :: V_rot_sh, FU(3), JU(3,3), &
                                             & HU(3,3,3), &
                                             & FI(3), JI(3,3), HI(3,3,3), &
                                             & aIt(3), aItt(3), &
                                             & JIt(3,3), &
                                             & time_d, ri(3), vi(3), &
                                             & ru(3), rIt(3), rItt(3), inter(3,3,3)
        acc = 0._qp
        jac = 0._qp
        hes = 0._qp
        time_d = real(time,dp)
        RIU = me%rot%call(time_d)
        RIUdot = me%rot%calldot(time_d)
        RIUdotdot = me%rot%callddot(time_d)
        ri = real(y(:3),dp)
        vi = real(y(4:6),dp)
        ru = mmult(transpose(RIU),ri)
        call shpines(real(me%central_body_ref_radius,dp), &
                   & real(me%central_body_mu,dp), &
                   & me%pdat, &
                   & me%degord, &
                   & 3, &
                   & ru, &
                   & V_rot_sh, &
                   & FU, &
                   & JU, &
                   & HU &
                   &)
        ! This is not really the rotating frame velocity,
        ! just the partial of the rotating radius with respect
        ! to time
        rIt = mmult(transpose(RIUdot),ri)
        rItt = mmult(transpose(RIUdotdot),ri)
        FI = mmult(RIU,FU)
        acc = y(8) * [y(4:6),real(FI,qp), 1._qp, 0._qp]
        JI = mmult(RIU,mmult(JU,transpose(RIU)))
        HI = mattens(RIU,quad(transpose(RIU),HU,3),3)
        aIt = mmult(RIUdot,FU) &
            + mmult(mmult(RIU,JU),rIt)
        ! these lines implement the operation
        ! ia, abc, jb -> ijc with RIU, HU, RIU
        inter = mattens(RIU,HU,3)
        inter = reshape(inter,[3,3,3],order=[2,1,3])
        inter = mattens(RIU,inter,3)
        inter = reshape(inter,[3,3,3],order=[2,1,3])
        JIt = mmult(RIUdot,mmult(JU,transpose(RIU))) &
            + mmult(RIU,mmult(JU,transpose(RIUdot))) &
            + vectens3(rIt,inter,3)
        aItt = mmult(RIUdotdot,FU) &
             + 2*mmult(RIUdot,mmult(JU,rIt)) &
             + mmult(RIU,vectensquad(rIt,HU,3)) &
             + mmult(RIU,mmult(JU,rItt))
        ! Assign acceleration in quad
        ! Fill in Jacobian nonzero blocks
        jac(1:3,4:6) = y(8)*eyemat(3)
        jac(1:3,8)   = y(4:6)
        jac(4:6,1:3) = y(8)*real(JI,qp)
        jac(4:6,7)   = y(8)*real(aIt,qp)
        jac(4:6,8)   = real(FI,qp)
        jac(7,8)     = 1._qp
        ! Fill in Hessian nonzero blocks
        hes(4:6,1:3,1:3) = y(8)*real(HI,qp)
        hes(4:6,7,1:3)   = y(8)*real(JIt,qp)
        hes(4:6,1:3,7)   = y(8)*real(JIt,qp)
        hes(4:6,8,1:3)   = real(JI,qp)
        hes(4:6,1:3,8)   = real(JI,qp)
        hes(1:3,8,4:6)   = eyemat(3)
        hes(1:3,4:6,8)   = eyemat(3)
        hes(4:6,7,7)     = y(8)*real(aItt,qp)
        hes(4:6,8,7)     = real(aIt,qp)
        hes(4:6,7,8)     = real(aIt,qp)
    end subroutine
    subroutine accelonly_sh(me, time, y, acc)
        ! allderivs_sh: compute the dynamics, jacobian, and hessian
        !               due to the gravitation of an extended body.
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! time           real           time in seconds past J2000 to evaluate
        ! y              real (8)       state vector from central body to
        !                               field point
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! acc            float (8)      time derivative of 6-d state,
        !                               i.e. the dynamics for the system:
        !                               (xdot, ydot, zdot, xddot, yddot, zddot)
        class(dynamicsModel), intent(inout) :: me
        real(qp),             intent(in)    :: time, y(:)
        real(qp),             intent(out)   :: acc(size(y))
        real(dp)                            :: RIU(3,3), RIUdot(3,3)
        real(dp)                            :: V_rot_sh, FU(3), JU(3,3), &
                                             & HU(3,3,3), &
                                             & FI(3), &
                                             & time_d, ri(3), vi(3), &
                                             & ru(3), vu(3)
        acc = 0._qp
        time_d = real(time,dp)
        RIU = me%rot%call(time_d)
        ri = real(y(:3),dp)
        vi = real(y(4:6),dp)
        ru = mmult(transpose(RIU),ri)
        call shpines(real(me%central_body_ref_radius,dp), &
                   & real(me%central_body_mu,dp), &
                   & me%pdat, &
                   & me%degord, &
                   & 1, &
                   & ru, &
                   & V_rot_sh, &
                   & FU, &
                   & JU, &
                   & HU &
                   &)
        vu = mmult(transpose(RIUdot),ri) + mmult(transpose(RIU),vi)
        FI = mmult(RIU,FU)
        acc = y(8) * [y(4:6),real(FI,qp), 1._qp, 0._qp]
    end subroutine
    function eoms(me,t,y) result(res)
        ! eoms: method to compute dynamics function of extended state vector
        !       i.e. [xdot, stmdot, sttdot]
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! t              real           time in sec past J2000 to evaluate
        !                               the dynamics
        ! y              real (584)     Extended dyamics vector:
        !                               584 = 8 + 8 ** 2 + 8 ** 3
        !                        y = [x, reshape(stm,(64)), reshape(stt,(512))]
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (584)     big vector of dynamics
        class(dynamicsModel), intent(inout) :: me
        real(qp)            , intent(in)    :: t, y(:)
        real(qp)                            :: stm(m,m), &
                                               stt(m,m,m), &
                                               statedot(m), &
                                               stmdot(m,m), &
                                               sttdot(m,m,m), &
                                               acc(m), &
                                               jac(m,m), &
                                               hes(m,m,m)
        real(qp)                            :: res(size(y))
        stm = reshape(y((m + 1):(m + m**2)),[m,m])
        stt = reshape(y((m + m**2 + 1):),[m,m,m])
        me%state = y(:m)
        call me%get_derivs(y(7), acc, jac, hes)
        statedot = acc
        stmdot = matmul(jac,stm)
        sttdot = mattens(jac,stt,m) + quad(stm,hes,m)
        res = [statedot, reshape(stmdot,[m**2]), reshape(sttdot,[m**3])]
    end function eoms
    function eoms_rails(me,t,y) result(res)
        ! eoms: method to compute dynamics function of extended state vector
        !       i.e. [xdot, stmdot, sttdot]
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! t              real           time in sec past J2000 to evaluate
        !                               the dynamics
        ! y              real (577)     Extended dyamics vector:
        !                               577 = 1 + 8 ** 2 + 8 ** 3
        !                        y = [t, reshape(stm,(64)), reshape(stt,(512))]
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (577)     big vector of dynamics
        class(dynamicsModel), intent(inout) :: me
        real(qp)            , intent(in)    :: t, y(:)
        real(qp)                            :: state(6), &
                                               stm(m,m), &
                                               stt(m,m,m), &
                                               stmdot(m,m), &
                                               sttdot(m,m,m), &
                                               acc(m), &
                                               jac(m,m), &
                                               hes(m,m,m), &
                                               r
        real(qp)                            :: res(size(y))
        stm = reshape(y(2:(1 + m**2)),[m,m])
        stt = reshape(y((2 + m**2):),[m,m,m])
        call me%get_derivs(y(1), acc, jac, hes)
        stmdot = matmul(jac,stm)
        sttdot = mattens(jac,stt,m) + quad(stm,hes,m)
        if (me%regularize) then
            state = me%trajstate(y(1))
            r = sqrt(sum(state(:3)**2))
            res(1) = r**(1.5_qp)*me%tof
        else 
            res(1) = me%tof
        endif
        res(2:1 + 8**2) = reshape(stmdot,[8**2])
        res(2+8**2:) = reshape(sttdot,[8**3])
    end function eoms_rails
    function trajstate(me,time) result(res)
        ! trajstate: method to return reference state at a time
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! time           real           time in sec past J2000 to evaluate
        !                               the dynamics
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6)       the reference state in the SPICE
        !                               J2000 coordinate frame
        !                               centered on the central body
        class(dynamicsmodel), intent(inout) :: me
        real(qp),             intent(in)    :: time
        real(dp)                            :: traj_state(6)
        real(qp)                            :: res(6)
        traj_state(:3) = me%bod_db%call(real(time,dp),me%traj_id,'p')
        traj_state(4:) = me%bod_db%call(real(time,dp),me%traj_id,'v')
        res = real(traj_state,qp)
    end function trajstate
    function outer_2vec(veca, vecb) result(res)
        ! outer_2vec: convenience function to take the outer product
        !             of two vectors
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! veca           real (n)       first vector to compute outer product
        ! vecb           real (m)       second vector to compute outer product
        !                               
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (n,m)     outer product of veca and vecb
        real(qp), intent(in) :: veca(:), vecb(:)
        real(qp)             :: res(size(veca),size(vecb))
        integer j, n,m
        n = size(veca); m = size(vecb)
        do j = 1,m
            res(:,j) = veca*vecb(j)
        end do
    end function outer_2vec
    function outer_3vec(veca, vecb, vecc) result(res)
        ! outer_3vec: convenience function to take the outer product
        !             of three vectors that appears in some of the gravity
        !             Hessians. The output is equivalent to:
        !             res(i,j,k) = veca(i)*vecb(j)*vecc(k)
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! veca           real (n)       first vector to compute outer product
        ! vecb           real (m)       second vector to compute outer product
        ! vecc           real (l)       third vector to compute outer product
        !                               
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (n,m,l)   outer product of veca, vecb, and vecc
        real(qp), intent(in) :: veca(:), vecb(:), vecc(:)
        real(qp)             :: res(size(veca),size(vecb),size(vecc))
        integer i,j, k, n, m, l
        n = size(veca); m = size(vecb); l = size(vecc)
        do i=1,n
            do j= 1,m
                do k=1,l
                    res(i,j,k) = veca(i)*vecb(j)*vecc(k)
                end do
            end do
        end do
    end function outer_3vec
    function outer_kron(vec, ind) result(res)
        ! outer_kron: convenience function to take the three index kronecker 
        !             delta appearing in some gravitational Hessians. 
        !             The output changes based on the value of ind:
        !             ind = 1 => res(i,j,k) = vec(i)*delta(j,k)
        !             ind = 2 => res(i,j,k) = vec(j)*delta(i,k)
        !             ind = 3 => res(i,j,k) = vec(k)*delta(i,j)
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! vec            real (n)       vector for computing outer product
        ! ind            integer        index not to appear in the delta
        !                               
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (n,n,n)   three index kronecker delta
        integer,  intent(in) :: ind
        real(qp), intent(in) :: vec(:)
        real(qp)             :: res(size(vec),size(vec),size(vec))
        integer i, j, n
        n = size(vec)
        res = 0._qp
        select case(ind)
        case(1)
            do i=1,n
                do j= 1,n
                    res(i,j,j) = vec(i)
                end do
            end do
        case(2)
            do i=1,n
                do j= 1,n
                    res(j,i,j) = vec(i)
                end do
            end do
        case(3)
            do i=1,n
                do j= 1,n
                    res(j,j,i) = vec(i)
                end do
            end do
        end select
    end function outer_kron
end module makemodel
