module makemodel
    use, intrinsic :: iso_fortran_env, only: wp=>real64, dp=>real64, qp=>real128
    use astkindmodule
    use shdatamodule
    use pinesmodule
    use tensorops
    use cheby, only: spice_subset
    use quat, only: rothist
    implicit none
    integer :: m=8

    type :: dynamicsModel
        type(spice_subset)    :: bod_db
        type(rothist)         :: rot
        integer               :: num_bodies, &
                               & central_body, &
                               & traj_id
        integer,  allocatable :: bodylist(:) ! num_bodies
        logical               :: shgrav, &
                               & tgt_on_rails
        real(qp)              :: central_body_ref_radius, &
                               & central_body_rad(3), &
                               & state(8), &
                               & central_body_mu, &
                               & tof
        real(qp), allocatable :: nbody_mus(:)     ! (num_bodies)
        type(SHarmData)       :: shdat
        type(PinesData)       :: pdat
        contains
            procedure:: init => init_dm
            procedure :: eoms
            procedure :: eoms_rails
            procedure :: get_derivs
            procedure :: acc_kepler
            procedure :: jac_kepler
            procedure :: hes_kepler
            procedure :: acc_nbody
            procedure :: jac_nbody
            procedure :: hes_nbody
            procedure :: fd_acc_kepler
            procedure :: fd_jac_kepler
            procedure :: fd_acc_nbody
            procedure :: fd_jac_nbody
            procedure :: fd_hes_kepler
            procedure :: fd_hes_nbody
            procedure :: allderivs_sh
            procedure :: trajstate
            procedure :: new_bodies
            procedure :: new_gravstatus
    end type dynamicsModel

    contains

    subroutine init_dm(me, subspice, traj_id, central_body, bodylist, &
                     & central_body_mu, central_body_ref_radius, mu_list, &
                     & shgrav, Cbar, Sbar, rot)
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
        ! Cbar           real (:,:)     4 pi (Kaula) normalized cosine Stokes 
        !                               for central body
        ! Sbar           real (:,:)     4 pi (Kaula) normalized sine Stokes 
        !                               for central body
        ! tgt_on_rails   logical        TRUE if target trajectory should not
        !                               be integrated
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
        logical,              intent(in)              :: shgrav
        real(qp),             intent(in)              :: Cbar(:,:), &
                                                         Sbar(:,:)
        real(dp)                                      :: rotmat(3,3)
        me%num_bodies = size(bodylist)
        me%shgrav = shgrav
        me%traj_id = traj_id
        me%bod_db = subspice
        allocate(me%bodylist(me%num_bodies),  &
                 me%nbody_mus(me%num_bodies))
        me%bodylist = bodylist
        me%central_body = central_body
        me%central_body_ref_radius = real(central_body_ref_radius,qp)
        me%central_body_mu = real(central_body_mu,qp)
        me%nbody_mus = real(mu_list,qp)
        me%tgt_on_rails = tgt_on_rails

        if (me%shgrav) then
            if (present(rot)) then
                me%rot = rot
            else
                print *, "Error: rotation history required. None provided"
                error stop
            endif
            call shdat_from_table( real(me%central_body_ref_radius,dp), & 
                                 & real(me%central_body_mu,dp), &
                                 & real(Cbar,dp), real(Sbar,dp), &
                                 & size(Cbar,1),size(Cbar,2), me%shdat &
                                 )
            call pinesinit(size(Cbar,1),me%shdat%Cml,me%shdat%Sml,me%pdat)
        endif
    end subroutine init_dm
    subroutine new_bodies(me, bodylist, mulist)
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
    subroutine new_gravstatus(me, shgrav, Cbar, Sbar)
        class(dynamicsModel), intent(inout) :: me
        logical,              intent(in)    :: shgrav
        real(dp),             intent(in)    :: Cbar(:,:), Sbar(:,:)
        me%shgrav = shgrav
        if (me%shgrav) then
            call shdat_from_table( real(me%central_body_ref_radius,dp), & 
                                 & real(me%central_body_mu,dp), &
                                 & real(Cbar,dp), real(Sbar,dp), &
                                 & size(Cbar,1),size(Cbar,2), me%shdat &
                                 )
            call pinesinit(size(Cbar,1),me%shdat%Cml,me%shdat%Sml,me%pdat)
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
        real(qp)                            :: acc_2b(m), jac_2b(m,m), hes_2b(m,m,m)
        real(qp)                            :: acc_nb(m), jac_nb(m,m), hes_nb(m,m,m)
        real(qp)                            :: nbody_radii(3,me%num_bodies), &
                                               nbody_vels(3,me%num_bodies), &
                                               nbody_accs(3,me%num_bodies)
        real(qp)                            :: r_bod_up(3), v_bod_up(3), a_bod_up(3), &
                                               y(m)
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
            ! for bodies i = 2. . .
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
            ! PLACEHOLDER
            ! TODO: Put in extended body gravity
            acc_2b = acc_kepler(me, me%central_body_mu, y)
            jac_2b = jac_kepler(me, me%central_body_mu, y)
            hes_2b = hes_kepler(me, me%central_body_mu, y)
            ! END PLACEHOLDER
        else
            acc_2b = acc_kepler(me, me%central_body_mu, y)
            jac_2b = jac_kepler(me, me%central_body_mu, y)
            hes_2b = hes_kepler(me, me%central_body_mu, y)
        endif
        ! for bodies i = 2. . . N NOTE: CAN BE DONE IN ||
        ! Compute point mass accel, jac, hess for i=1
        acc_nb = 0._qp
        jac_nb = 0._qp
        hes_nb = 0._qp
        do i=1,me%num_bodies
            r_bod_up = nbody_radii(:,i)
            v_bod_up = nbody_vels(:,i)
            a_bod_up = nbody_accs(:,i)
            acc_nb   = acc_nb + acc_nbody(me, me%nbody_mus(i), y,r_bod_up)
            jac_nb   = jac_nb + jac_nbody(me, me%nbody_mus(i), y,r_bod_up,v_bod_up)
            hes_nb   = hes_nb + hes_nbody(me, me%nbody_mus(i), y,r_bod_up,v_bod_up,a_bod_up)
        end do
        acc_nb(:3) = 0._qp
        acc = acc_2b + acc_nb
        jac = jac_2b
        jac(4:,:) = jac(4:,:) + jac_nb(4:,:)
        hes = hes_2b
        hes(4:,:,:) = hes(4:,:,:) + hes_nb(4:,:,:)
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
    function acc_kepler(me, mu, y) result (res)
        ! acc_kepler: method to compute Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! y              real (:)       Augmented dynamical state vector
        !                               field point in first 3 elements,
        !                               could be state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8)       Augmented acceleration vector
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: mu
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y))
        res = xdot(y)
        contains 
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0

                real(qp), dimension(8) :: res

                x0  =  mu*y(8)/(y(1)**2 + y(2)**2 + y(3)**2)**(3.0_qp/2.0_qp)

                res(1) =  y(4)*y(8)
                res(2) =  y(5)*y(8)
                res(3) =  y(6)*y(8)
                res(4) =  -x0*y(1)
                res(5) =  -x0*y(2)
                res(6) =  -x0*y(3)
                res(7) =  y(8)
                res(8) =  0
            end function xdot
            ! END AUTOCODE OUTPUT FOR XDOT
    end function
    function jac_kepler(me, mu, y) result (res)
        ! jac_kepler: method to compute Keplerian jacobian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! y              real (:)       Augmented state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8)     Keplerian Jacobian matrix 
        !                                         df_i
        !                               J(i,j) = ------
        !                                         dx_j
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: mu
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y))
        res = reshape(jac(y),[m,m])
        contains
            ! BEGIN AUTOCODE OUTPUT FOR JAC
            function jac(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10
                real(qp), dimension(64) :: inter
                real(qp), dimension(8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  mu/x3**(3.0_qp/2.0_qp)
                x5  =  -x4*y(8)
                x6  =  3*mu*y(8)/x3**(5.0_qp/2.0_qp)
                x7  =  x6*y(1)
                x8  =  x7*y(2)
                x9  =  x7*y(3)
                x10  =  x6*y(2)*y(3)

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x0*x6 + x5
                inter(5) =  x8
                inter(6) =  x9
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x8
                inter(13) =  x1*x6 + x5
                inter(14) =  x10
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x9
                inter(21) =  x10
                inter(22) =  x2*x6 + x5
                inter(23) =  0
                inter(24) =  0
                inter(25) =  y(8)
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  y(8)
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  y(8)
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  0
                inter(53) =  0
                inter(54) =  0
                inter(55) =  0
                inter(56) =  0
                inter(57) =  y(4)
                inter(58) =  y(5)
                inter(59) =  y(6)
                inter(60) =  -x4*y(1)
                inter(61) =  -x4*y(2)
                inter(62) =  -x4*y(3)
                inter(63) =  1
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function hes_kepler(me, mu, y) result (res)
        ! hes_kepler: method to compute Keplerian hessian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! y              real (:)       Augmented state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8,8)   Keplerian Hessian tensor H where
        !                                              df_i 
        !                               H(i,j,k) = -----------
        !                                           dx_j dx_k
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: mu
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y),size(y))
        res = reshape(hes(y),[m,m,m])
        contains
            ! BEGIN AUTOCODE OUTPUT FOR HES
            function hes(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  15*mu*y(8)/x3**(7.0_qp/2.0_qp)
                x6  =  -3*mu*x4*y(2)*y(8)
                x7  =  x0*x5
                x8  =  -x6 - x7*y(2)
                x9  =  -3*mu*x4*y(3)*y(8)
                x10  =  -x7*y(3) - x9
                x11  =  -3*mu*x4*y(1)*y(8)
                x12  =  x5*y(1)
                x13  =  -x1*x12 - x11
                x14  =  y(2)*y(3)
                x15  =  -x12*x14
                x16  =  -x11 - x12*x2
                x17  =  -mu/x3**(3.0_qp/2.0_qp)
                x18  =  3*mu*x4
                x19  =  x0*x18 + x17
                x20  =  x18*y(1)
                x21  =  x20*y(2)
                x22  =  x20*y(3)
                x23  =  -x1*x5*y(3) - x9
                x24  =  -x2*x5*y(2) - x6
                x25  =  x1*x18 + x17
                x26  =  x14*x18
                x27  =  x17 + x18*x2

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  9*mu*x4*y(1)*y(8) - x5*y(1)**3
                inter(5) =  x8
                inter(6) =  x10
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x8
                inter(13) =  x13
                inter(14) =  x15
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x10
                inter(21) =  x15
                inter(22) =  x16
                inter(23) =  0
                inter(24) =  0
                inter(25) =  0
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  0
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  0
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  0
                inter(53) =  0
                inter(54) =  0
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  x19
                inter(61) =  x21
                inter(62) =  x22
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x8
                inter(69) =  x13
                inter(70) =  x15
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x13
                inter(77) =  9*mu*x4*y(2)*y(8) - x5*y(2)**3
                inter(78) =  x23
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x15
                inter(85) =  x23
                inter(86) =  x24
                inter(87) =  0
                inter(88) =  0
                inter(89) =  0
                inter(90) =  0
                inter(91) =  0
                inter(92) =  0
                inter(93) =  0
                inter(94) =  0
                inter(95) =  0
                inter(96) =  0
                inter(97) =  0
                inter(98) =  0
                inter(99) =  0
                inter(100) =  0
                inter(101) =  0
                inter(102) =  0
                inter(103) =  0
                inter(104) =  0
                inter(105) =  0
                inter(106) =  0
                inter(107) =  0
                inter(108) =  0
                inter(109) =  0
                inter(110) =  0
                inter(111) =  0
                inter(112) =  0
                inter(113) =  0
                inter(114) =  0
                inter(115) =  0
                inter(116) =  0
                inter(117) =  0
                inter(118) =  0
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  x21
                inter(125) =  x25
                inter(126) =  x26
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x10
                inter(133) =  x15
                inter(134) =  x16
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x15
                inter(141) =  x23
                inter(142) =  x24
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x16
                inter(149) =  x24
                inter(150) =  9*mu*x4*y(3)*y(8) - x5*y(3)**3
                inter(151) =  0
                inter(152) =  0
                inter(153) =  0
                inter(154) =  0
                inter(155) =  0
                inter(156) =  0
                inter(157) =  0
                inter(158) =  0
                inter(159) =  0
                inter(160) =  0
                inter(161) =  0
                inter(162) =  0
                inter(163) =  0
                inter(164) =  0
                inter(165) =  0
                inter(166) =  0
                inter(167) =  0
                inter(168) =  0
                inter(169) =  0
                inter(170) =  0
                inter(171) =  0
                inter(172) =  0
                inter(173) =  0
                inter(174) =  0
                inter(175) =  0
                inter(176) =  0
                inter(177) =  0
                inter(178) =  0
                inter(179) =  0
                inter(180) =  0
                inter(181) =  0
                inter(182) =  0
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  x22
                inter(189) =  x26
                inter(190) =  x27
                inter(191) =  0
                inter(192) =  0
                inter(193) =  0
                inter(194) =  0
                inter(195) =  0
                inter(196) =  0
                inter(197) =  0
                inter(198) =  0
                inter(199) =  0
                inter(200) =  0
                inter(201) =  0
                inter(202) =  0
                inter(203) =  0
                inter(204) =  0
                inter(205) =  0
                inter(206) =  0
                inter(207) =  0
                inter(208) =  0
                inter(209) =  0
                inter(210) =  0
                inter(211) =  0
                inter(212) =  0
                inter(213) =  0
                inter(214) =  0
                inter(215) =  0
                inter(216) =  0
                inter(217) =  0
                inter(218) =  0
                inter(219) =  0
                inter(220) =  0
                inter(221) =  0
                inter(222) =  0
                inter(223) =  0
                inter(224) =  0
                inter(225) =  0
                inter(226) =  0
                inter(227) =  0
                inter(228) =  0
                inter(229) =  0
                inter(230) =  0
                inter(231) =  0
                inter(232) =  0
                inter(233) =  0
                inter(234) =  0
                inter(235) =  0
                inter(236) =  0
                inter(237) =  0
                inter(238) =  0
                inter(239) =  0
                inter(240) =  0
                inter(241) =  0
                inter(242) =  0
                inter(243) =  0
                inter(244) =  0
                inter(245) =  0
                inter(246) =  0
                inter(247) =  0
                inter(248) =  0
                inter(249) =  1
                inter(250) =  0
                inter(251) =  0
                inter(252) =  0
                inter(253) =  0
                inter(254) =  0
                inter(255) =  0
                inter(256) =  0
                inter(257) =  0
                inter(258) =  0
                inter(259) =  0
                inter(260) =  0
                inter(261) =  0
                inter(262) =  0
                inter(263) =  0
                inter(264) =  0
                inter(265) =  0
                inter(266) =  0
                inter(267) =  0
                inter(268) =  0
                inter(269) =  0
                inter(270) =  0
                inter(271) =  0
                inter(272) =  0
                inter(273) =  0
                inter(274) =  0
                inter(275) =  0
                inter(276) =  0
                inter(277) =  0
                inter(278) =  0
                inter(279) =  0
                inter(280) =  0
                inter(281) =  0
                inter(282) =  0
                inter(283) =  0
                inter(284) =  0
                inter(285) =  0
                inter(286) =  0
                inter(287) =  0
                inter(288) =  0
                inter(289) =  0
                inter(290) =  0
                inter(291) =  0
                inter(292) =  0
                inter(293) =  0
                inter(294) =  0
                inter(295) =  0
                inter(296) =  0
                inter(297) =  0
                inter(298) =  0
                inter(299) =  0
                inter(300) =  0
                inter(301) =  0
                inter(302) =  0
                inter(303) =  0
                inter(304) =  0
                inter(305) =  0
                inter(306) =  0
                inter(307) =  0
                inter(308) =  0
                inter(309) =  0
                inter(310) =  0
                inter(311) =  0
                inter(312) =  0
                inter(313) =  0
                inter(314) =  1
                inter(315) =  0
                inter(316) =  0
                inter(317) =  0
                inter(318) =  0
                inter(319) =  0
                inter(320) =  0
                inter(321) =  0
                inter(322) =  0
                inter(323) =  0
                inter(324) =  0
                inter(325) =  0
                inter(326) =  0
                inter(327) =  0
                inter(328) =  0
                inter(329) =  0
                inter(330) =  0
                inter(331) =  0
                inter(332) =  0
                inter(333) =  0
                inter(334) =  0
                inter(335) =  0
                inter(336) =  0
                inter(337) =  0
                inter(338) =  0
                inter(339) =  0
                inter(340) =  0
                inter(341) =  0
                inter(342) =  0
                inter(343) =  0
                inter(344) =  0
                inter(345) =  0
                inter(346) =  0
                inter(347) =  0
                inter(348) =  0
                inter(349) =  0
                inter(350) =  0
                inter(351) =  0
                inter(352) =  0
                inter(353) =  0
                inter(354) =  0
                inter(355) =  0
                inter(356) =  0
                inter(357) =  0
                inter(358) =  0
                inter(359) =  0
                inter(360) =  0
                inter(361) =  0
                inter(362) =  0
                inter(363) =  0
                inter(364) =  0
                inter(365) =  0
                inter(366) =  0
                inter(367) =  0
                inter(368) =  0
                inter(369) =  0
                inter(370) =  0
                inter(371) =  0
                inter(372) =  0
                inter(373) =  0
                inter(374) =  0
                inter(375) =  0
                inter(376) =  0
                inter(377) =  0
                inter(378) =  0
                inter(379) =  1
                inter(380) =  0
                inter(381) =  0
                inter(382) =  0
                inter(383) =  0
                inter(384) =  0
                inter(385) =  0
                inter(386) =  0
                inter(387) =  0
                inter(388) =  0
                inter(389) =  0
                inter(390) =  0
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  0
                inter(397) =  0
                inter(398) =  0
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  0
                inter(405) =  0
                inter(406) =  0
                inter(407) =  0
                inter(408) =  0
                inter(409) =  0
                inter(410) =  0
                inter(411) =  0
                inter(412) =  0
                inter(413) =  0
                inter(414) =  0
                inter(415) =  0
                inter(416) =  0
                inter(417) =  0
                inter(418) =  0
                inter(419) =  0
                inter(420) =  0
                inter(421) =  0
                inter(422) =  0
                inter(423) =  0
                inter(424) =  0
                inter(425) =  0
                inter(426) =  0
                inter(427) =  0
                inter(428) =  0
                inter(429) =  0
                inter(430) =  0
                inter(431) =  0
                inter(432) =  0
                inter(433) =  0
                inter(434) =  0
                inter(435) =  0
                inter(436) =  0
                inter(437) =  0
                inter(438) =  0
                inter(439) =  0
                inter(440) =  0
                inter(441) =  0
                inter(442) =  0
                inter(443) =  0
                inter(444) =  0
                inter(445) =  0
                inter(446) =  0
                inter(447) =  0
                inter(448) =  0
                inter(449) =  0
                inter(450) =  0
                inter(451) =  0
                inter(452) =  x19
                inter(453) =  x21
                inter(454) =  x22
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  x21
                inter(461) =  x25
                inter(462) =  x26
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  x22
                inter(469) =  x26
                inter(470) =  x27
                inter(471) =  0
                inter(472) =  0
                inter(473) =  1
                inter(474) =  0
                inter(475) =  0
                inter(476) =  0
                inter(477) =  0
                inter(478) =  0
                inter(479) =  0
                inter(480) =  0
                inter(481) =  0
                inter(482) =  1
                inter(483) =  0
                inter(484) =  0
                inter(485) =  0
                inter(486) =  0
                inter(487) =  0
                inter(488) =  0
                inter(489) =  0
                inter(490) =  0
                inter(491) =  1
                inter(492) =  0
                inter(493) =  0
                inter(494) =  0
                inter(495) =  0
                inter(496) =  0
                inter(497) =  0
                inter(498) =  0
                inter(499) =  0
                inter(500) =  0
                inter(501) =  0
                inter(502) =  0
                inter(503) =  0
                inter(504) =  0
                inter(505) =  0
                inter(506) =  0
                inter(507) =  0
                inter(508) =  0
                inter(509) =  0
                inter(510) =  0
                inter(511) =  0
                inter(512) =  0
                res = reshape(inter,[8,8,8])
            end function hes
            ! END AUTOCODE OUTPUT FOR HES
    end function
    function acc_nbody(me, mu, y, rbods) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! y              real (:)       Augmented dynamical state vector
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8)       Dynamics
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: mu
        real(qp),             intent(in) :: y(:), &
                                            rbods(3)
        real(qp)                         :: res(size(y))
        res = xdot(y)
        contains
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5

                real(qp), dimension(8) :: res

                associate(rj => rbods, mu_j => mu)
                    x0  =  (rj(1)**2 + rj(2)**2 + rj(3)**2)**(-3.0_qp/2.0_qp)
                    x1  =  rj(1) - y(1)
                    x2  =  rj(2) - y(2)
                    x3  =  rj(3) - y(3)
                    x4  =  (x1**2 + x2**2 + x3**2)**(-3.0_qp/2.0_qp)
                    x5  =  mu_j*y(8)

                    res(1) =  y(4)*y(8)
                    res(2) =  y(5)*y(8)
                    res(3) =  y(6)*y(8)
                    res(4) =  x5*(-rj(1)*x0 + x1*x4)
                    res(5) =  x5*(-rj(2)*x0 + x2*x4)
                    res(6) =  x5*(-rj(3)*x0 + x3*x4)
                    res(7) =  0._qp
                    res(8) =  0._qp
                end associate
            end function xdot
            ! END AUTOCODE OUTPUT FOR XDOT
    end function
    function jac_nbody(me, mu, y, rbods,vbods) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! y              real (:)       Dynamical state vector
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! vbods          real (3)       Velocity vector of third body
        !                               relative to central body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8)     third-body Jacobian matrix 
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: mu
        real(qp),             intent(in) :: y(:), &
                                            rbods(3), &
                                            vbods(3)
        real(qp)                         :: res(size(y),size(y))
        real(qp)                         :: rj(3), vj(3), mu_j
        rj=rbods; vj=vbods; mu_j = mu
        res = jac(y)
        contains
            ! BEGIN AUTOCODE OUTPUT FOR JAC
            function jac(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23
                real(qp), dimension(64) :: inter
                real(qp), dimension(8,8) :: res

                x0  =  rj(1) - y(1)
                x1  =  rj(2) - y(2)
                x2  =  rj(3) - y(3)
                x3  =  x0**2 + x1**2 + x2**2
                x4  =  x3**(-3.0_qp/2.0_qp)
                x5  =  -x4
                x6  =  3*rj(1)
                x7  =  x6 - 3*y(1)
                x8  =  x3**(-5.0_qp/2.0_qp)
                x9  =  x0*x8
                x10  =  mu_j*y(8)
                x11  =  x1*x8
                x12  =  x10*x7
                x13  =  x2*x8
                x14  =  3*rj(2)
                x15  =  x14 - 3*y(2)
                x16  =  x10*x15
                x17  =  3*rj(3)
                x18  =  x17 - 3*y(3)
                x19  =  x10*x18
                x20  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x21  =  x20**(-3.0_qp/2.0_qp)
                x22  =  (-vj(1)*x6 - vj(2)*x14 - vj(3)*x17)/x20**(5.0_qp/2.0_qp)
                x23  =  -3*vj(1)*x0 - 3*vj(2)*x1 - 3*vj(3)*x2

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x10*(x5 + x7*x9)
                inter(5) =  x11*x12
                inter(6) =  x12*x13
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x16*x9
                inter(13) =  x10*(x11*x15 + x5)
                inter(14) =  x13*x16
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x19*x9
                inter(21) =  x11*x19
                inter(22) =  x10*(x13*x18 + x5)
                inter(23) =  0
                inter(24) =  0
                inter(25) =  y(8)
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  y(8)
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  y(8)
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  x10*(-rj(1)*x22 - vj(1)*x21 + vj(1)*x4 + x0*x23*x8)
                inter(53) =  x10*(-rj(2)*x22 - vj(2)*x21 + vj(2)*x4 + x1*x23*x8)
                inter(54) =  x10*(-rj(3)*x22 - vj(3)*x21 + vj(3)*x4 + x2*x23*x8)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  y(4)
                inter(58) =  y(5)
                inter(59) =  y(6)
                inter(60) =  mu_j*(-rj(1)*x21 + x0*x4)
                inter(61) =  mu_j*(-rj(2)*x21 + x1*x4)
                inter(62) =  mu_j*(-rj(3)*x21 + x2*x4)
                inter(63) =  0
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function hes_nbody(me, mu, y, rbods, vbods, abods) result (res)
        ! jac_kepler: method to compute third body jacobian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! y              real (:)       Dynamical state vector
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! vbods          real (3)       Velocity vector of third body
        !                               relative to central body
        ! abods          real (3)       Acceleration vector of third body
        !                               relative to central body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8,8)   third body Hessian tensor H where
        !                                              df_i 
        !                               H(i,j,k) = -----------
        !                                           dx_j dx_k
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: mu
        real(qp),             intent(in) :: y(:), &
                                            rbods(3), &
                                            vbods(3), &
                                            abods(3)
        real(qp)                         :: res(size(y),size(y),size(y))
        real(qp)                         :: rj(3), vj(3), aj(3), mu_j
        rj=rbods; vj = vbods; aj = abods; mu_j = mu
        res = hes(y)
        contains
            ! BEGIN AUTOCODE OUTPUT FOR HES
            function hes(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                      & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                      & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                      & x48, x49, x50, x51, x52, x53, x54, x55, & 
                                      & x56, x57, x58, x59, x60, x61, x62, x63, & 
                                      & x64, x65, x66, x67, x68, x69, x70, x71, & 
                                      & x72, x73, x74, x75, x76, x77, x78, x79, & 
                                      & x80, x81, x82, x83, x84, x85, x86, x87, & 
                                      & x88, x89, x90, x91, x92, x93, x94, x95, & 
                                      & x96, x97, x98, x99, x100, x101, x102
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  rj(1) - y(1)
                x1  =  rj(2) - y(2)
                x2  =  rj(3) - y(3)
                x3  =  x0**2 + x1**2 + x2**2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  3*x0
                x6  =  x4*x5
                x7  =  3*rj(1)
                x8  =  x7 - 3*y(1)
                x9  =  x4*x8
                x10  =  5*rj(1)
                x11  =  x10 - 5*y(1)
                x12  =  x3**(-7.0_qp/2.0_qp)
                x13  =  mu_j*y(8)
                x14  =  3*x1
                x15  =  x14*x4
                x16  =  -x13*x15
                x17  =  x11*x12
                x18  =  x1*x17
                x19  =  x13*x8
                x20  =  3*x2
                x21  =  x20*x4
                x22  =  -x13*x21
                x23  =  x17*x2
                x24  =  3*rj(2)
                x25  =  x24 - 3*y(2)
                x26  =  x25*x4
                x27  =  -x13*x26
                x28  =  x0*x17
                x29  =  x13*x25
                x30  =  3*rj(3)
                x31  =  x30 - 3*y(3)
                x32  =  x31*x4
                x33  =  -x13*x32
                x34  =  x13*x31
                x35  =  vj(2)*x14
                x36  =  vj(3)*x20
                x37  =  -vj(1)*x5 - x35 - x36
                x38  =  x37*x4
                x39  =  -x38
                x40  =  vj(1)*x6 + vj(1)*x9 + x39
                x41  =  vj(1)*x15
                x42  =  vj(2)*x9
                x43  =  vj(1)*x21
                x44  =  vj(3)*x9
                x45  =  x3**(-3.0_qp/2.0_qp)
                x46  =  -x45
                x47  =  mu_j*(x0*x9 + x46)
                x48  =  mu_j*x9
                x49  =  x1*x48
                x50  =  x2*x48
                x51  =  5*rj(2)
                x52  =  x51 - 5*y(2)
                x53  =  x12*x52
                x54  =  x0*x53
                x55  =  -x13*x9
                x56  =  x1*x53
                x57  =  x2*x53
                x58  =  -x13*x6
                x59  =  vj(2)*x6
                x60  =  vj(1)*x26
                x61  =  vj(2)*x26 + x35*x4 + x39
                x62  =  vj(2)*x21
                x63  =  vj(3)*x26
                x64  =  mu_j*x26
                x65  =  x0*x64
                x66  =  mu_j*(x1*x26 + x46)
                x67  =  x2*x64
                x68  =  -x32
                x69  =  5*rj(3)
                x70  =  x69 - 5*y(3)
                x71  =  x12*x70
                x72  =  x0*x71
                x73  =  x1*x71
                x74  =  x2*x71
                x75  =  vj(3)*x6
                x76  =  vj(1)*x32
                x77  =  vj(3)*x15
                x78  =  vj(2)*x32
                x79  =  vj(3)*x32 + x36*x4 + x39
                x80  =  mu_j*x32
                x81  =  x0*x80
                x82  =  x1*x80
                x83  =  mu_j*(x2*x32 + x46)
                x84  =  x12*(-5*vj(1)*x0 - 5*vj(2)*x1 - 5*vj(3)*x2)
                x85  =  x0*x84
                x86  =  x1*x84
                x87  =  x2*x84
                x88  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x89  =  x88**(-3.0_qp/2.0_qp)
                x90  =  x88**(-5.0_qp/2.0_qp)
                x91  =  -vj(1)*x7 - vj(2)*x24 - vj(3)*x30
                x92  =  x90*x91
                x93  =  2*vj(1)
                x94  =  x91*(-vj(1)*x10 - vj(2)*x51 - vj(3)*x69)/x88**(7.0_qp/2.0_qp)
                x95  =  3*vj(1)**2 + 3*vj(2)**2 + 3*vj(3)**2
                x96  =  x90*(-aj(1)*x7 - aj(2)*x24 - aj(3)*x30 - x95)
                x97  =  x4*(-aj(1)*x5 - aj(2)*x14 - aj(3)*x20 - x95)
                x98  =  2*vj(2)
                x99  =  2*vj(3)
                x100  =  mu_j*(-rj(1)*x92 + vj(1)*x45 - vj(1)*x89 + x0*x37*x4)
                x101  =  mu_j*(-rj(2)*x92 + vj(2)*x45 - vj(2)*x89 + x1*x37*x4)
                x102  =  mu_j*(-rj(3)*x92 + vj(3)*x45 - vj(3)*x89 + x2*x37*x4)

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x13*(x0*x11*x12*x8 - x6 - 2*x9)
                inter(5) =  x16 + x18*x19
                inter(6) =  x19*x23 + x22
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x27 + x28*x29
                inter(13) =  x13*(x1*x11*x12*x25 - x9)
                inter(14) =  x23*x29
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x28*x34 + x33
                inter(21) =  x18*x34
                inter(22) =  x13*(x11*x12*x2*x31 - x9)
                inter(23) =  0
                inter(24) =  0
                inter(25) =  0
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  0
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  0
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  x13*(x28*x37 + x40)
                inter(53) =  x13*(x18*x37 + x41 + x42)
                inter(54) =  x13*(x23*x37 + x43 + x44)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  x47
                inter(61) =  x49
                inter(62) =  x50
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x13*(-x26 + x54*x8)
                inter(69) =  x19*x56 + x55
                inter(70) =  x19*x57
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x29*x54 + x58
                inter(77) =  x13*(x1*x12*x25*x52 - x15 - 2*x26)
                inter(78) =  x22 + x29*x57
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x34*x54
                inter(85) =  x33 + x34*x56
                inter(86) =  x13*(x12*x2*x31*x52 - x26)
                inter(87) =  0
                inter(88) =  0
                inter(89) =  0
                inter(90) =  0
                inter(91) =  0
                inter(92) =  0
                inter(93) =  0
                inter(94) =  0
                inter(95) =  0
                inter(96) =  0
                inter(97) =  0
                inter(98) =  0
                inter(99) =  0
                inter(100) =  0
                inter(101) =  0
                inter(102) =  0
                inter(103) =  0
                inter(104) =  0
                inter(105) =  0
                inter(106) =  0
                inter(107) =  0
                inter(108) =  0
                inter(109) =  0
                inter(110) =  0
                inter(111) =  0
                inter(112) =  0
                inter(113) =  0
                inter(114) =  0
                inter(115) =  0
                inter(116) =  x13*(x37*x54 + x59 + x60)
                inter(117) =  x13*(x37*x56 + x61)
                inter(118) =  x13*(x37*x57 + x62 + x63)
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  x65
                inter(125) =  x66
                inter(126) =  x67
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x13*(x68 + x72*x8)
                inter(133) =  x19*x73
                inter(134) =  x19*x74 + x55
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x29*x72
                inter(141) =  x13*(x25*x73 + x68)
                inter(142) =  x27 + x29*x74
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x34*x72 + x58
                inter(149) =  x16 + x34*x73
                inter(150) =  x13*(x12*x2*x31*x70 - x21 - 2*x32)
                inter(151) =  0
                inter(152) =  0
                inter(153) =  0
                inter(154) =  0
                inter(155) =  0
                inter(156) =  0
                inter(157) =  0
                inter(158) =  0
                inter(159) =  0
                inter(160) =  0
                inter(161) =  0
                inter(162) =  0
                inter(163) =  0
                inter(164) =  0
                inter(165) =  0
                inter(166) =  0
                inter(167) =  0
                inter(168) =  0
                inter(169) =  0
                inter(170) =  0
                inter(171) =  0
                inter(172) =  0
                inter(173) =  0
                inter(174) =  0
                inter(175) =  0
                inter(176) =  0
                inter(177) =  0
                inter(178) =  0
                inter(179) =  0
                inter(180) =  x13*(x37*x72 + x75 + x76)
                inter(181) =  x13*(x37*x73 + x77 + x78)
                inter(182) =  x13*(x37*x74 + x79)
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  x81
                inter(189) =  x82
                inter(190) =  x83
                inter(191) =  0
                inter(192) =  0
                inter(193) =  0
                inter(194) =  0
                inter(195) =  0
                inter(196) =  0
                inter(197) =  0
                inter(198) =  0
                inter(199) =  0
                inter(200) =  0
                inter(201) =  0
                inter(202) =  0
                inter(203) =  0
                inter(204) =  0
                inter(205) =  0
                inter(206) =  0
                inter(207) =  0
                inter(208) =  0
                inter(209) =  0
                inter(210) =  0
                inter(211) =  0
                inter(212) =  0
                inter(213) =  0
                inter(214) =  0
                inter(215) =  0
                inter(216) =  0
                inter(217) =  0
                inter(218) =  0
                inter(219) =  0
                inter(220) =  0
                inter(221) =  0
                inter(222) =  0
                inter(223) =  0
                inter(224) =  0
                inter(225) =  0
                inter(226) =  0
                inter(227) =  0
                inter(228) =  0
                inter(229) =  0
                inter(230) =  0
                inter(231) =  0
                inter(232) =  0
                inter(233) =  0
                inter(234) =  0
                inter(235) =  0
                inter(236) =  0
                inter(237) =  0
                inter(238) =  0
                inter(239) =  0
                inter(240) =  0
                inter(241) =  0
                inter(242) =  0
                inter(243) =  0
                inter(244) =  0
                inter(245) =  0
                inter(246) =  0
                inter(247) =  0
                inter(248) =  0
                inter(249) =  1
                inter(250) =  0
                inter(251) =  0
                inter(252) =  0
                inter(253) =  0
                inter(254) =  0
                inter(255) =  0
                inter(256) =  0
                inter(257) =  0
                inter(258) =  0
                inter(259) =  0
                inter(260) =  0
                inter(261) =  0
                inter(262) =  0
                inter(263) =  0
                inter(264) =  0
                inter(265) =  0
                inter(266) =  0
                inter(267) =  0
                inter(268) =  0
                inter(269) =  0
                inter(270) =  0
                inter(271) =  0
                inter(272) =  0
                inter(273) =  0
                inter(274) =  0
                inter(275) =  0
                inter(276) =  0
                inter(277) =  0
                inter(278) =  0
                inter(279) =  0
                inter(280) =  0
                inter(281) =  0
                inter(282) =  0
                inter(283) =  0
                inter(284) =  0
                inter(285) =  0
                inter(286) =  0
                inter(287) =  0
                inter(288) =  0
                inter(289) =  0
                inter(290) =  0
                inter(291) =  0
                inter(292) =  0
                inter(293) =  0
                inter(294) =  0
                inter(295) =  0
                inter(296) =  0
                inter(297) =  0
                inter(298) =  0
                inter(299) =  0
                inter(300) =  0
                inter(301) =  0
                inter(302) =  0
                inter(303) =  0
                inter(304) =  0
                inter(305) =  0
                inter(306) =  0
                inter(307) =  0
                inter(308) =  0
                inter(309) =  0
                inter(310) =  0
                inter(311) =  0
                inter(312) =  0
                inter(313) =  0
                inter(314) =  1
                inter(315) =  0
                inter(316) =  0
                inter(317) =  0
                inter(318) =  0
                inter(319) =  0
                inter(320) =  0
                inter(321) =  0
                inter(322) =  0
                inter(323) =  0
                inter(324) =  0
                inter(325) =  0
                inter(326) =  0
                inter(327) =  0
                inter(328) =  0
                inter(329) =  0
                inter(330) =  0
                inter(331) =  0
                inter(332) =  0
                inter(333) =  0
                inter(334) =  0
                inter(335) =  0
                inter(336) =  0
                inter(337) =  0
                inter(338) =  0
                inter(339) =  0
                inter(340) =  0
                inter(341) =  0
                inter(342) =  0
                inter(343) =  0
                inter(344) =  0
                inter(345) =  0
                inter(346) =  0
                inter(347) =  0
                inter(348) =  0
                inter(349) =  0
                inter(350) =  0
                inter(351) =  0
                inter(352) =  0
                inter(353) =  0
                inter(354) =  0
                inter(355) =  0
                inter(356) =  0
                inter(357) =  0
                inter(358) =  0
                inter(359) =  0
                inter(360) =  0
                inter(361) =  0
                inter(362) =  0
                inter(363) =  0
                inter(364) =  0
                inter(365) =  0
                inter(366) =  0
                inter(367) =  0
                inter(368) =  0
                inter(369) =  0
                inter(370) =  0
                inter(371) =  0
                inter(372) =  0
                inter(373) =  0
                inter(374) =  0
                inter(375) =  0
                inter(376) =  0
                inter(377) =  0
                inter(378) =  0
                inter(379) =  1
                inter(380) =  0
                inter(381) =  0
                inter(382) =  0
                inter(383) =  0
                inter(384) =  0
                inter(385) =  0
                inter(386) =  0
                inter(387) =  0
                inter(388) =  x13*(x40 + x8*x85)
                inter(389) =  x13*x41 + x13*x42 + x19*x86
                inter(390) =  x13*x43 + x13*x44 + x19*x87
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  x13*x59 + x13*x60 + x29*x85
                inter(397) =  x13*(x25*x86 + x61)
                inter(398) =  x13*x62 + x13*x63 + x29*x87
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  x13*x75 + x13*x76 + x34*x85
                inter(405) =  x13*x77 + x13*x78 + x34*x86
                inter(406) =  x13*(x31*x87 + x79)
                inter(407) =  0
                inter(408) =  0
                inter(409) =  0
                inter(410) =  0
                inter(411) =  0
                inter(412) =  0
                inter(413) =  0
                inter(414) =  0
                inter(415) =  0
                inter(416) =  0
                inter(417) =  0
                inter(418) =  0
                inter(419) =  0
                inter(420) =  0
                inter(421) =  0
                inter(422) =  0
                inter(423) =  0
                inter(424) =  0
                inter(425) =  0
                inter(426) =  0
                inter(427) =  0
                inter(428) =  0
                inter(429) =  0
                inter(430) =  0
                inter(431) =  0
                inter(432) =  0
                inter(433) =  0
                inter(434) =  0
                inter(435) =  0
                inter(436) =  x13*(aj(1)*x45 - aj(1)*x89 - rj(1)*x94 - rj(1)*x96 + x0*x97 + x37*x85 + x38*x93 - x92*x93)
                inter(437) =  x13*(aj(2)*x45 - aj(2)*x89 - rj(2)*x94 - rj(2)*x96 + x1*x97 + x37*x86 + x38*x98 - x92*x98)
                inter(438) =  x13*(aj(3)*x45 - aj(3)*x89 - rj(3)*x94 - rj(3)*x96 + x2*x97 + x37*x87 + x38*x99 - x92*x99)
                inter(439) =  0
                inter(440) =  0
                inter(441) =  0
                inter(442) =  0
                inter(443) =  0
                inter(444) =  x100
                inter(445) =  x101
                inter(446) =  x102
                inter(447) =  0
                inter(448) =  0
                inter(449) =  0
                inter(450) =  0
                inter(451) =  0
                inter(452) =  x47
                inter(453) =  x49
                inter(454) =  x50
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  x65
                inter(461) =  x66
                inter(462) =  x67
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  x81
                inter(469) =  x82
                inter(470) =  x83
                inter(471) =  0
                inter(472) =  0
                inter(473) =  1
                inter(474) =  0
                inter(475) =  0
                inter(476) =  0
                inter(477) =  0
                inter(478) =  0
                inter(479) =  0
                inter(480) =  0
                inter(481) =  0
                inter(482) =  1
                inter(483) =  0
                inter(484) =  0
                inter(485) =  0
                inter(486) =  0
                inter(487) =  0
                inter(488) =  0
                inter(489) =  0
                inter(490) =  0
                inter(491) =  1
                inter(492) =  0
                inter(493) =  0
                inter(494) =  0
                inter(495) =  0
                inter(496) =  0
                inter(497) =  0
                inter(498) =  0
                inter(499) =  0
                inter(500) =  x100
                inter(501) =  x101
                inter(502) =  x102
                inter(503) =  0
                inter(504) =  0
                inter(505) =  0
                inter(506) =  0
                inter(507) =  0
                inter(508) =  0
                inter(509) =  0
                inter(510) =  0
                inter(511) =  0
                inter(512) =  0
                res = reshape(inter,[8,8,8])
            end function hes
            ! END AUTOCODE OUTPUT FOR HES
    end function
    subroutine allderivs_sh(me, time, r, acc, jac, hes)
        ! allderivs_sh: compute the dynamics, jacobian, and hessian
        !               due to the gravitation of an extended body.
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! time           real           time in seconds past J2000 to evaluate
        ! r              real (3)       position vector from central body to
        !                               field point
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! NAME           TYPE           DESCRIPTION
        ! acc            float (6)      time derivative of 6-d state,
        !                               i.e. the dynamics for the system:
        !                               (xdot, ydot, zdot, xddot, yddot, zddot)
        ! jac            float (6,6)    jacobian of dynamics: 
        !                               jac(i,j) = df_i/dx_j
        ! hes            float (6,6,6)  hessian of dynamics: 
        !                               hes(i,j,k) = d^2f_i/(dx_j*dx_k)
        class(dynamicsModel), intent(inout)  :: me
        real(qp),             intent(in)     :: time, &
                                              & r(3)
        real(dp)                             :: pot_d, &
                                              & acc_d(3), &
                                              & jac_d(3,3), &
                                              & hes_d(3,3,3), &
                                              & acc_qp(3), &
                                              & jac_qp(3,3), &
                                              & hes_qp(3,3,3)
        real(qp),             intent(out)    :: acc(6), &
                                              & jac(6,6), &
                                              & hes(6,6,6)
        real(dp)                             :: xform_mat(3,3), &
                                              & r_bf(3)

        ! TODO: Call transform to get transformation matrix 
        ! to/from body fixed frame in xform_mat
        r_bf = mmult(xform_mat,r_bf)

        call shpines( me%shdat%Rbody, me%shdat%GM, me%pdat, &
                    & me%shdat%maxdeg, me%shdat%maxorder, &
                    & real(r_bf,dp), pot_d, acc_d, jac_d, hes_d &
                    & )
        ! TODO: add correct velocity transformation
        acc(4:) = real(acc_d,qp)
        jac(4:,:3) = real(jac_d,qp)
        hes(4:,:3,:3) = real(hes_d,qp)
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
        real(qp)                            :: stm(m,m), &
                                               stt(m,m,m), &
                                               stmdot(m,m), &
                                               sttdot(m,m,m), &
                                               acc(m), &
                                               jac(m,m), &
                                               hes(m,m,m)
        real(qp)                            :: res(size(y))
        stm = reshape(y((1+1):(1 + m**2)),[m,m])
        stt = reshape(y((1 + m**2 + 1):),[m,m,m])
        call me%get_derivs(y(1), acc, jac, hes)
        res(1) = me%tof
        stmdot = matmul(jac,stm)
        sttdot = mattens(jac,stt,m) + quad(stm,hes,m)
        res = [me%tof, reshape(stmdot,[m**2]), reshape(sttdot,[m**3])]
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
    !~~~~~~~~~~~ TESTING ~~~~~~~~~~~
    function fd_acc_kepler(me, y) result (res)
        ! acc_kepler: method to compute Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! y              real (:)       Augmented dynamical state vector
        !                               field point in first 3 elements,
        !                               could be state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8)       Augmented acceleration vector
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y))
        real(qp)                         ::  mu
        mu = me%central_body_mu
        res = xdot(y)
        contains 
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0

                real(qp), dimension(8) :: res

                x0  =  mu*y(8)/(y(1)**2 + y(2)**2 + y(3)**2)**(3.0_qp/2.0_qp)

                res(1) =  y(4)*y(8)
                res(2) =  y(5)*y(8)
                res(3) =  y(6)*y(8)
                res(4) =  -x0*y(1)
                res(5) =  -x0*y(2)
                res(6) =  -x0*y(3)
                res(7) =  y(8)
                res(8) =  0
            end function xdot
            ! END AUTOCODE OUTPUT FOR XDOT
    end function
    function fd_jac_kepler(me, y) result (res)
        ! jac_kepler: method to compute Keplerian jacobian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! y              real (:)       Augmented state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8)     Keplerian Jacobian matrix 
        !                                         df_i
        !                               J(i,j) = ------
        !                                         dx_j
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y))
        real(qp)                         :: t, mu
        t = y(7)
        mu = me%central_body_mu
        res = reshape(jac(y),[m,m])
        contains
            ! BEGIN AUTOCODE OUTPUT FOR JAC
            function jac(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10
                real(qp), dimension(64) :: inter
                real(qp), dimension(8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  mu/x3**(3.0_qp/2.0_qp)
                x5  =  -x4*y(8)
                x6  =  3*mu*y(8)/x3**(5.0_qp/2.0_qp)
                x7  =  x6*y(1)
                x8  =  x7*y(2)
                x9  =  x7*y(3)
                x10  =  x6*y(2)*y(3)

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x0*x6 + x5
                inter(5) =  x8
                inter(6) =  x9
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x8
                inter(13) =  x1*x6 + x5
                inter(14) =  x10
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x9
                inter(21) =  x10
                inter(22) =  x2*x6 + x5
                inter(23) =  0
                inter(24) =  0
                inter(25) =  y(8)
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  y(8)
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  y(8)
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  0
                inter(53) =  0
                inter(54) =  0
                inter(55) =  0
                inter(56) =  0
                inter(57) =  y(4)
                inter(58) =  y(5)
                inter(59) =  y(6)
                inter(60) =  -x4*y(1)
                inter(61) =  -x4*y(2)
                inter(62) =  -x4*y(3)
                inter(63) =  1
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function fd_hes_kepler(me, y) result (res)
        ! hes_kepler: method to compute Keplerian hessian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! y              real (:)       Augmented state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8,8)   Keplerian Hessian tensor H where
        !                                              df_i 
        !                               H(i,j,k) = -----------
        !                                           dx_j dx_k
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y),size(y)), mu
        mu = me%central_body_mu
        res = reshape(hes(y),[m,m,m])
        contains
            ! BEGIN AUTOCODE OUTPUT FOR HES
            function hes(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  15*mu*y(8)/x3**(7.0_qp/2.0_qp)
                x6  =  -3*mu*x4*y(2)*y(8)
                x7  =  x0*x5
                x8  =  -x6 - x7*y(2)
                x9  =  -3*mu*x4*y(3)*y(8)
                x10  =  -x7*y(3) - x9
                x11  =  -3*mu*x4*y(1)*y(8)
                x12  =  x5*y(1)
                x13  =  -x1*x12 - x11
                x14  =  y(2)*y(3)
                x15  =  -x12*x14
                x16  =  -x11 - x12*x2
                x17  =  -mu/x3**(3.0_qp/2.0_qp)
                x18  =  3*mu*x4
                x19  =  x0*x18 + x17
                x20  =  x18*y(1)
                x21  =  x20*y(2)
                x22  =  x20*y(3)
                x23  =  -x1*x5*y(3) - x9
                x24  =  -x2*x5*y(2) - x6
                x25  =  x1*x18 + x17
                x26  =  x14*x18
                x27  =  x17 + x18*x2

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  9*mu*x4*y(1)*y(8) - x5*y(1)**3
                inter(5) =  x8
                inter(6) =  x10
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x8
                inter(13) =  x13
                inter(14) =  x15
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x10
                inter(21) =  x15
                inter(22) =  x16
                inter(23) =  0
                inter(24) =  0
                inter(25) =  0
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  0
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  0
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  0
                inter(53) =  0
                inter(54) =  0
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  x19
                inter(61) =  x21
                inter(62) =  x22
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x8
                inter(69) =  x13
                inter(70) =  x15
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x13
                inter(77) =  9*mu*x4*y(2)*y(8) - x5*y(2)**3
                inter(78) =  x23
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x15
                inter(85) =  x23
                inter(86) =  x24
                inter(87) =  0
                inter(88) =  0
                inter(89) =  0
                inter(90) =  0
                inter(91) =  0
                inter(92) =  0
                inter(93) =  0
                inter(94) =  0
                inter(95) =  0
                inter(96) =  0
                inter(97) =  0
                inter(98) =  0
                inter(99) =  0
                inter(100) =  0
                inter(101) =  0
                inter(102) =  0
                inter(103) =  0
                inter(104) =  0
                inter(105) =  0
                inter(106) =  0
                inter(107) =  0
                inter(108) =  0
                inter(109) =  0
                inter(110) =  0
                inter(111) =  0
                inter(112) =  0
                inter(113) =  0
                inter(114) =  0
                inter(115) =  0
                inter(116) =  0
                inter(117) =  0
                inter(118) =  0
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  x21
                inter(125) =  x25
                inter(126) =  x26
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x10
                inter(133) =  x15
                inter(134) =  x16
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x15
                inter(141) =  x23
                inter(142) =  x24
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x16
                inter(149) =  x24
                inter(150) =  9*mu*x4*y(3)*y(8) - x5*y(3)**3
                inter(151) =  0
                inter(152) =  0
                inter(153) =  0
                inter(154) =  0
                inter(155) =  0
                inter(156) =  0
                inter(157) =  0
                inter(158) =  0
                inter(159) =  0
                inter(160) =  0
                inter(161) =  0
                inter(162) =  0
                inter(163) =  0
                inter(164) =  0
                inter(165) =  0
                inter(166) =  0
                inter(167) =  0
                inter(168) =  0
                inter(169) =  0
                inter(170) =  0
                inter(171) =  0
                inter(172) =  0
                inter(173) =  0
                inter(174) =  0
                inter(175) =  0
                inter(176) =  0
                inter(177) =  0
                inter(178) =  0
                inter(179) =  0
                inter(180) =  0
                inter(181) =  0
                inter(182) =  0
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  x22
                inter(189) =  x26
                inter(190) =  x27
                inter(191) =  0
                inter(192) =  0
                inter(193) =  0
                inter(194) =  0
                inter(195) =  0
                inter(196) =  0
                inter(197) =  0
                inter(198) =  0
                inter(199) =  0
                inter(200) =  0
                inter(201) =  0
                inter(202) =  0
                inter(203) =  0
                inter(204) =  0
                inter(205) =  0
                inter(206) =  0
                inter(207) =  0
                inter(208) =  0
                inter(209) =  0
                inter(210) =  0
                inter(211) =  0
                inter(212) =  0
                inter(213) =  0
                inter(214) =  0
                inter(215) =  0
                inter(216) =  0
                inter(217) =  0
                inter(218) =  0
                inter(219) =  0
                inter(220) =  0
                inter(221) =  0
                inter(222) =  0
                inter(223) =  0
                inter(224) =  0
                inter(225) =  0
                inter(226) =  0
                inter(227) =  0
                inter(228) =  0
                inter(229) =  0
                inter(230) =  0
                inter(231) =  0
                inter(232) =  0
                inter(233) =  0
                inter(234) =  0
                inter(235) =  0
                inter(236) =  0
                inter(237) =  0
                inter(238) =  0
                inter(239) =  0
                inter(240) =  0
                inter(241) =  0
                inter(242) =  0
                inter(243) =  0
                inter(244) =  0
                inter(245) =  0
                inter(246) =  0
                inter(247) =  0
                inter(248) =  0
                inter(249) =  1
                inter(250) =  0
                inter(251) =  0
                inter(252) =  0
                inter(253) =  0
                inter(254) =  0
                inter(255) =  0
                inter(256) =  0
                inter(257) =  0
                inter(258) =  0
                inter(259) =  0
                inter(260) =  0
                inter(261) =  0
                inter(262) =  0
                inter(263) =  0
                inter(264) =  0
                inter(265) =  0
                inter(266) =  0
                inter(267) =  0
                inter(268) =  0
                inter(269) =  0
                inter(270) =  0
                inter(271) =  0
                inter(272) =  0
                inter(273) =  0
                inter(274) =  0
                inter(275) =  0
                inter(276) =  0
                inter(277) =  0
                inter(278) =  0
                inter(279) =  0
                inter(280) =  0
                inter(281) =  0
                inter(282) =  0
                inter(283) =  0
                inter(284) =  0
                inter(285) =  0
                inter(286) =  0
                inter(287) =  0
                inter(288) =  0
                inter(289) =  0
                inter(290) =  0
                inter(291) =  0
                inter(292) =  0
                inter(293) =  0
                inter(294) =  0
                inter(295) =  0
                inter(296) =  0
                inter(297) =  0
                inter(298) =  0
                inter(299) =  0
                inter(300) =  0
                inter(301) =  0
                inter(302) =  0
                inter(303) =  0
                inter(304) =  0
                inter(305) =  0
                inter(306) =  0
                inter(307) =  0
                inter(308) =  0
                inter(309) =  0
                inter(310) =  0
                inter(311) =  0
                inter(312) =  0
                inter(313) =  0
                inter(314) =  1
                inter(315) =  0
                inter(316) =  0
                inter(317) =  0
                inter(318) =  0
                inter(319) =  0
                inter(320) =  0
                inter(321) =  0
                inter(322) =  0
                inter(323) =  0
                inter(324) =  0
                inter(325) =  0
                inter(326) =  0
                inter(327) =  0
                inter(328) =  0
                inter(329) =  0
                inter(330) =  0
                inter(331) =  0
                inter(332) =  0
                inter(333) =  0
                inter(334) =  0
                inter(335) =  0
                inter(336) =  0
                inter(337) =  0
                inter(338) =  0
                inter(339) =  0
                inter(340) =  0
                inter(341) =  0
                inter(342) =  0
                inter(343) =  0
                inter(344) =  0
                inter(345) =  0
                inter(346) =  0
                inter(347) =  0
                inter(348) =  0
                inter(349) =  0
                inter(350) =  0
                inter(351) =  0
                inter(352) =  0
                inter(353) =  0
                inter(354) =  0
                inter(355) =  0
                inter(356) =  0
                inter(357) =  0
                inter(358) =  0
                inter(359) =  0
                inter(360) =  0
                inter(361) =  0
                inter(362) =  0
                inter(363) =  0
                inter(364) =  0
                inter(365) =  0
                inter(366) =  0
                inter(367) =  0
                inter(368) =  0
                inter(369) =  0
                inter(370) =  0
                inter(371) =  0
                inter(372) =  0
                inter(373) =  0
                inter(374) =  0
                inter(375) =  0
                inter(376) =  0
                inter(377) =  0
                inter(378) =  0
                inter(379) =  1
                inter(380) =  0
                inter(381) =  0
                inter(382) =  0
                inter(383) =  0
                inter(384) =  0
                inter(385) =  0
                inter(386) =  0
                inter(387) =  0
                inter(388) =  0
                inter(389) =  0
                inter(390) =  0
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  0
                inter(397) =  0
                inter(398) =  0
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  0
                inter(405) =  0
                inter(406) =  0
                inter(407) =  0
                inter(408) =  0
                inter(409) =  0
                inter(410) =  0
                inter(411) =  0
                inter(412) =  0
                inter(413) =  0
                inter(414) =  0
                inter(415) =  0
                inter(416) =  0
                inter(417) =  0
                inter(418) =  0
                inter(419) =  0
                inter(420) =  0
                inter(421) =  0
                inter(422) =  0
                inter(423) =  0
                inter(424) =  0
                inter(425) =  0
                inter(426) =  0
                inter(427) =  0
                inter(428) =  0
                inter(429) =  0
                inter(430) =  0
                inter(431) =  0
                inter(432) =  0
                inter(433) =  0
                inter(434) =  0
                inter(435) =  0
                inter(436) =  0
                inter(437) =  0
                inter(438) =  0
                inter(439) =  0
                inter(440) =  0
                inter(441) =  0
                inter(442) =  0
                inter(443) =  0
                inter(444) =  0
                inter(445) =  0
                inter(446) =  0
                inter(447) =  0
                inter(448) =  0
                inter(449) =  0
                inter(450) =  0
                inter(451) =  0
                inter(452) =  x19
                inter(453) =  x21
                inter(454) =  x22
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  x21
                inter(461) =  x25
                inter(462) =  x26
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  x22
                inter(469) =  x26
                inter(470) =  x27
                inter(471) =  0
                inter(472) =  0
                inter(473) =  1
                inter(474) =  0
                inter(475) =  0
                inter(476) =  0
                inter(477) =  0
                inter(478) =  0
                inter(479) =  0
                inter(480) =  0
                inter(481) =  0
                inter(482) =  1
                inter(483) =  0
                inter(484) =  0
                inter(485) =  0
                inter(486) =  0
                inter(487) =  0
                inter(488) =  0
                inter(489) =  0
                inter(490) =  0
                inter(491) =  1
                inter(492) =  0
                inter(493) =  0
                inter(494) =  0
                inter(495) =  0
                inter(496) =  0
                inter(497) =  0
                inter(498) =  0
                inter(499) =  0
                inter(500) =  0
                inter(501) =  0
                inter(502) =  0
                inter(503) =  0
                inter(504) =  0
                inter(505) =  0
                inter(506) =  0
                inter(507) =  0
                inter(508) =  0
                inter(509) =  0
                inter(510) =  0
                inter(511) =  0
                inter(512) =  0
                res = reshape(inter,[8,8,8])
            end function hes
            ! END AUTOCODE OUTPUT FOR HES
    end function
    function fd_acc_nbody(me, y) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! y              real (:)       Augmented dynamical state vector
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8)       Dynamics
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y))
        real(qp)                         :: thisr(3)
        real(qp)                         :: t
        integer                          :: i
        t = y(7)
        res = 0._qp
        do i=1,me%num_bodies
            thisr= me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'p')
            res = res + xdot(y, thisr,me%nbody_mus(i))
        end do
        contains
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(y, rj, mu_j) result(res)
                implicit none
                real(qp), intent(in) :: y (:), rj(:), mu_j
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5

                real(qp), dimension(8) :: res

                x0  =  (rj(1)**2 + rj(2)**2 + rj(3)**2)**(-3.0_qp/2.0_qp)
                x1  =  rj(1) - y(1)
                x2  =  rj(2) - y(2)
                x3  =  rj(3) - y(3)
                x4  =  (x1**2 + x2**2 + x3**2)**(-3.0_qp/2.0_qp)
                x5  =  mu_j*y(8)

                res(1) =  y(4)*y(8)
                res(2) =  y(5)*y(8)
                res(3) =  y(6)*y(8)
                res(4) =  x5*(-rj(1)*x0 + x1*x4)
                res(5) =  x5*(-rj(2)*x0 + x2*x4)
                res(6) =  x5*(-rj(3)*x0 + x3*x4)
                res(7) =  0
                res(8) =  0
            end function xdot
            ! END AUTOCODE OUTPUT FOR XDOT
    end function
    function fd_jac_nbody(me, y) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! y              real (:)       Dynamical state vector
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! vbods          real (3)       Velocity vector of third body
        !                               relative to central body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8)     third-body Jacobian matrix 
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y))
        real(qp)                         :: thisr(3), thisv(3)
        real(qp)                         :: t
        integer                          :: i
        t = y(7)
        res = 0._qp
        do i=1,me%num_bodies
            thisr= real(me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'p'),qp)
            thisv= real(me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'v'),qp)
            res = res + jac(y, thisr, thisv,me%nbody_mus(i))
        end do
        contains
            ! BEGIN AUTOCODE OUTPUT FOR JAC
            function jac(y,rj,vj,mu_j) result(res)
                implicit none
                real(qp), intent(in) :: y (:), rj(:), vj(:), mu_j
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23
                real(qp), dimension(64) :: inter
                real(qp), dimension(8,8) :: res

                x0  =  rj(1) - y(1)
                x1  =  rj(2) - y(2)
                x2  =  rj(3) - y(3)
                x3  =  x0**2 + x1**2 + x2**2
                x4  =  x3**(-3.0_qp/2.0_qp)
                x5  =  -x4
                x6  =  3*rj(1)
                x7  =  x6 - 3*y(1)
                x8  =  x3**(-5.0_qp/2.0_qp)
                x9  =  x0*x8
                x10  =  mu_j*y(8)
                x11  =  x1*x8
                x12  =  x10*x7
                x13  =  x2*x8
                x14  =  3*rj(2)
                x15  =  x14 - 3*y(2)
                x16  =  x10*x15
                x17  =  3*rj(3)
                x18  =  x17 - 3*y(3)
                x19  =  x10*x18
                x20  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x21  =  x20**(-3.0_qp/2.0_qp)
                x22  =  (-vj(1)*x6 - vj(2)*x14 - vj(3)*x17)/x20**(5.0_qp/2.0_qp)
                x23  =  -3*vj(1)*x0 - 3*vj(2)*x1 - 3*vj(3)*x2

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x10*(x5 + x7*x9)
                inter(5) =  x11*x12
                inter(6) =  x12*x13
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x16*x9
                inter(13) =  x10*(x11*x15 + x5)
                inter(14) =  x13*x16
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x19*x9
                inter(21) =  x11*x19
                inter(22) =  x10*(x13*x18 + x5)
                inter(23) =  0
                inter(24) =  0
                inter(25) =  y(8)
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  y(8)
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  y(8)
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  x10*(-rj(1)*x22 - vj(1)*x21 + vj(1)*x4 + x0*x23*x8)
                inter(53) =  x10*(-rj(2)*x22 - vj(2)*x21 + vj(2)*x4 + x1*x23*x8)
                inter(54) =  x10*(-rj(3)*x22 - vj(3)*x21 + vj(3)*x4 + x2*x23*x8)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  y(4)
                inter(58) =  y(5)
                inter(59) =  y(6)
                inter(60) =  mu_j*(-rj(1)*x21 + x0*x4)
                inter(61) =  mu_j*(-rj(2)*x21 + x1*x4)
                inter(62) =  mu_j*(-rj(3)*x21 + x2*x4)
                inter(63) =  0
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function fd_hes_nbody(me, y) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! y              real (:)       Dynamical state vector
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! vbods          real (3)       Velocity vector of third body
        !                               relative to central body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (8,8)     third-body Jacobian matrix 
        class(dynamicsModel), intent(in) :: me
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y),size(y))
        real(qp)                         :: thisr(3), thisv(3), thisa(3), &
                                          & thishes(8,8,8)
        real(qp)                         :: t
        integer                          :: i
        t = y(7)
        res = 0._qp
        do i=1,me%num_bodies
            thisr= real(me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'p'),qp)
            thisv= real(me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'v'),qp)
            thisa= real(me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'a'),qp)
            thishes =  hes(y, thisr, thisv, thisa, me%nbody_mus(i))
            res(4:,:,:) = res(4:,:,:) + thishes(4:,:,:)
        end do
        contains
            ! BEGIN AUTOCODE OUTPUT FOR HES
            function hes(y,rj,vj,aj,mu_j) result(res)
                implicit none
                real(qp), intent(in) :: y (:), rj(3), vj(3), aj(3), mu_j
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                      & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                      & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                      & x48, x49, x50, x51, x52, x53, x54, x55, & 
                                      & x56, x57, x58, x59, x60, x61, x62, x63, & 
                                      & x64, x65, x66, x67, x68, x69, x70, x71, & 
                                      & x72, x73, x74, x75, x76, x77, x78, x79, & 
                                      & x80, x81, x82, x83, x84, x85, x86, x87, & 
                                      & x88, x89, x90, x91, x92, x93, x94, x95, & 
                                      & x96, x97, x98, x99, x100, x101, x102
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  rj(1) - y(1)
                x1  =  rj(2) - y(2)
                x2  =  rj(3) - y(3)
                x3  =  x0**2 + x1**2 + x2**2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  3*x0
                x6  =  x4*x5
                x7  =  3*rj(1)
                x8  =  x7 - 3*y(1)
                x9  =  x4*x8
                x10  =  5*rj(1)
                x11  =  x10 - 5*y(1)
                x12  =  x3**(-7.0_qp/2.0_qp)
                x13  =  mu_j*y(8)
                x14  =  3*x1
                x15  =  x14*x4
                x16  =  -x13*x15
                x17  =  x11*x12
                x18  =  x1*x17
                x19  =  x13*x8
                x20  =  3*x2
                x21  =  x20*x4
                x22  =  -x13*x21
                x23  =  x17*x2
                x24  =  3*rj(2)
                x25  =  x24 - 3*y(2)
                x26  =  x25*x4
                x27  =  -x13*x26
                x28  =  x0*x17
                x29  =  x13*x25
                x30  =  3*rj(3)
                x31  =  x30 - 3*y(3)
                x32  =  x31*x4
                x33  =  -x13*x32
                x34  =  x13*x31
                x35  =  vj(2)*x14
                x36  =  vj(3)*x20
                x37  =  -vj(1)*x5 - x35 - x36
                x38  =  x37*x4
                x39  =  -x38
                x40  =  vj(1)*x6 + vj(1)*x9 + x39
                x41  =  vj(1)*x15
                x42  =  vj(2)*x9
                x43  =  vj(1)*x21
                x44  =  vj(3)*x9
                x45  =  x3**(-3.0_qp/2.0_qp)
                x46  =  -x45
                x47  =  mu_j*(x0*x9 + x46)
                x48  =  mu_j*x9
                x49  =  x1*x48
                x50  =  x2*x48
                x51  =  5*rj(2)
                x52  =  x51 - 5*y(2)
                x53  =  x12*x52
                x54  =  x0*x53
                x55  =  -x13*x9
                x56  =  x1*x53
                x57  =  x2*x53
                x58  =  -x13*x6
                x59  =  vj(2)*x6
                x60  =  vj(1)*x26
                x61  =  vj(2)*x26 + x35*x4 + x39
                x62  =  vj(2)*x21
                x63  =  vj(3)*x26
                x64  =  mu_j*x26
                x65  =  x0*x64
                x66  =  mu_j*(x1*x26 + x46)
                x67  =  x2*x64
                x68  =  -x32
                x69  =  5*rj(3)
                x70  =  x69 - 5*y(3)
                x71  =  x12*x70
                x72  =  x0*x71
                x73  =  x1*x71
                x74  =  x2*x71
                x75  =  vj(3)*x6
                x76  =  vj(1)*x32
                x77  =  vj(3)*x15
                x78  =  vj(2)*x32
                x79  =  vj(3)*x32 + x36*x4 + x39
                x80  =  mu_j*x32
                x81  =  x0*x80
                x82  =  x1*x80
                x83  =  mu_j*(x2*x32 + x46)
                x84  =  x12*(-5*vj(1)*x0 - 5*vj(2)*x1 - 5*vj(3)*x2)
                x85  =  x0*x84
                x86  =  x1*x84
                x87  =  x2*x84
                x88  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x89  =  x88**(-3.0_qp/2.0_qp)
                x90  =  x88**(-5.0_qp/2.0_qp)
                x91  =  -vj(1)*x7 - vj(2)*x24 - vj(3)*x30
                x92  =  x90*x91
                x93  =  2*vj(1)
                x94  =  x91*(-vj(1)*x10 - vj(2)*x51 - vj(3)*x69)/x88**(7.0_qp/2.0_qp)
                x95  =  3*vj(1)**2 + 3*vj(2)**2 + 3*vj(3)**2
                x96  =  x90*(-aj(1)*x7 - aj(2)*x24 - aj(3)*x30 - x95)
                x97  =  x4*(-aj(1)*x5 - aj(2)*x14 - aj(3)*x20 - x95)
                x98  =  2*vj(2)
                x99  =  2*vj(3)
                x100  =  mu_j*(-rj(1)*x92 + vj(1)*x45 - vj(1)*x89 + x0*x37*x4)
                x101  =  mu_j*(-rj(2)*x92 + vj(2)*x45 - vj(2)*x89 + x1*x37*x4)
                x102  =  mu_j*(-rj(3)*x92 + vj(3)*x45 - vj(3)*x89 + x2*x37*x4)

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x13*(x0*x11*x12*x8 - x6 - 2*x9)
                inter(5) =  x16 + x18*x19
                inter(6) =  x19*x23 + x22
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x27 + x28*x29
                inter(13) =  x13*(x1*x11*x12*x25 - x9)
                inter(14) =  x23*x29
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x28*x34 + x33
                inter(21) =  x18*x34
                inter(22) =  x13*(x11*x12*x2*x31 - x9)
                inter(23) =  0
                inter(24) =  0
                inter(25) =  0
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  0
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  0
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  0
                inter(50) =  0
                inter(51) =  0
                inter(52) =  x13*(x28*x37 + x40)
                inter(53) =  x13*(x18*x37 + x41 + x42)
                inter(54) =  x13*(x23*x37 + x43 + x44)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  x47
                inter(61) =  x49
                inter(62) =  x50
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x13*(-x26 + x54*x8)
                inter(69) =  x19*x56 + x55
                inter(70) =  x19*x57
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x29*x54 + x58
                inter(77) =  x13*(x1*x12*x25*x52 - x15 - 2*x26)
                inter(78) =  x22 + x29*x57
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x34*x54
                inter(85) =  x33 + x34*x56
                inter(86) =  x13*(x12*x2*x31*x52 - x26)
                inter(87) =  0
                inter(88) =  0
                inter(89) =  0
                inter(90) =  0
                inter(91) =  0
                inter(92) =  0
                inter(93) =  0
                inter(94) =  0
                inter(95) =  0
                inter(96) =  0
                inter(97) =  0
                inter(98) =  0
                inter(99) =  0
                inter(100) =  0
                inter(101) =  0
                inter(102) =  0
                inter(103) =  0
                inter(104) =  0
                inter(105) =  0
                inter(106) =  0
                inter(107) =  0
                inter(108) =  0
                inter(109) =  0
                inter(110) =  0
                inter(111) =  0
                inter(112) =  0
                inter(113) =  0
                inter(114) =  0
                inter(115) =  0
                inter(116) =  x13*(x37*x54 + x59 + x60)
                inter(117) =  x13*(x37*x56 + x61)
                inter(118) =  x13*(x37*x57 + x62 + x63)
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  x65
                inter(125) =  x66
                inter(126) =  x67
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x13*(x68 + x72*x8)
                inter(133) =  x19*x73
                inter(134) =  x19*x74 + x55
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x29*x72
                inter(141) =  x13*(x25*x73 + x68)
                inter(142) =  x27 + x29*x74
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x34*x72 + x58
                inter(149) =  x16 + x34*x73
                inter(150) =  x13*(x12*x2*x31*x70 - x21 - 2*x32)
                inter(151) =  0
                inter(152) =  0
                inter(153) =  0
                inter(154) =  0
                inter(155) =  0
                inter(156) =  0
                inter(157) =  0
                inter(158) =  0
                inter(159) =  0
                inter(160) =  0
                inter(161) =  0
                inter(162) =  0
                inter(163) =  0
                inter(164) =  0
                inter(165) =  0
                inter(166) =  0
                inter(167) =  0
                inter(168) =  0
                inter(169) =  0
                inter(170) =  0
                inter(171) =  0
                inter(172) =  0
                inter(173) =  0
                inter(174) =  0
                inter(175) =  0
                inter(176) =  0
                inter(177) =  0
                inter(178) =  0
                inter(179) =  0
                inter(180) =  x13*(x37*x72 + x75 + x76)
                inter(181) =  x13*(x37*x73 + x77 + x78)
                inter(182) =  x13*(x37*x74 + x79)
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  x81
                inter(189) =  x82
                inter(190) =  x83
                inter(191) =  0
                inter(192) =  0
                inter(193) =  0
                inter(194) =  0
                inter(195) =  0
                inter(196) =  0
                inter(197) =  0
                inter(198) =  0
                inter(199) =  0
                inter(200) =  0
                inter(201) =  0
                inter(202) =  0
                inter(203) =  0
                inter(204) =  0
                inter(205) =  0
                inter(206) =  0
                inter(207) =  0
                inter(208) =  0
                inter(209) =  0
                inter(210) =  0
                inter(211) =  0
                inter(212) =  0
                inter(213) =  0
                inter(214) =  0
                inter(215) =  0
                inter(216) =  0
                inter(217) =  0
                inter(218) =  0
                inter(219) =  0
                inter(220) =  0
                inter(221) =  0
                inter(222) =  0
                inter(223) =  0
                inter(224) =  0
                inter(225) =  0
                inter(226) =  0
                inter(227) =  0
                inter(228) =  0
                inter(229) =  0
                inter(230) =  0
                inter(231) =  0
                inter(232) =  0
                inter(233) =  0
                inter(234) =  0
                inter(235) =  0
                inter(236) =  0
                inter(237) =  0
                inter(238) =  0
                inter(239) =  0
                inter(240) =  0
                inter(241) =  0
                inter(242) =  0
                inter(243) =  0
                inter(244) =  0
                inter(245) =  0
                inter(246) =  0
                inter(247) =  0
                inter(248) =  0
                inter(249) =  1
                inter(250) =  0
                inter(251) =  0
                inter(252) =  0
                inter(253) =  0
                inter(254) =  0
                inter(255) =  0
                inter(256) =  0
                inter(257) =  0
                inter(258) =  0
                inter(259) =  0
                inter(260) =  0
                inter(261) =  0
                inter(262) =  0
                inter(263) =  0
                inter(264) =  0
                inter(265) =  0
                inter(266) =  0
                inter(267) =  0
                inter(268) =  0
                inter(269) =  0
                inter(270) =  0
                inter(271) =  0
                inter(272) =  0
                inter(273) =  0
                inter(274) =  0
                inter(275) =  0
                inter(276) =  0
                inter(277) =  0
                inter(278) =  0
                inter(279) =  0
                inter(280) =  0
                inter(281) =  0
                inter(282) =  0
                inter(283) =  0
                inter(284) =  0
                inter(285) =  0
                inter(286) =  0
                inter(287) =  0
                inter(288) =  0
                inter(289) =  0
                inter(290) =  0
                inter(291) =  0
                inter(292) =  0
                inter(293) =  0
                inter(294) =  0
                inter(295) =  0
                inter(296) =  0
                inter(297) =  0
                inter(298) =  0
                inter(299) =  0
                inter(300) =  0
                inter(301) =  0
                inter(302) =  0
                inter(303) =  0
                inter(304) =  0
                inter(305) =  0
                inter(306) =  0
                inter(307) =  0
                inter(308) =  0
                inter(309) =  0
                inter(310) =  0
                inter(311) =  0
                inter(312) =  0
                inter(313) =  0
                inter(314) =  1
                inter(315) =  0
                inter(316) =  0
                inter(317) =  0
                inter(318) =  0
                inter(319) =  0
                inter(320) =  0
                inter(321) =  0
                inter(322) =  0
                inter(323) =  0
                inter(324) =  0
                inter(325) =  0
                inter(326) =  0
                inter(327) =  0
                inter(328) =  0
                inter(329) =  0
                inter(330) =  0
                inter(331) =  0
                inter(332) =  0
                inter(333) =  0
                inter(334) =  0
                inter(335) =  0
                inter(336) =  0
                inter(337) =  0
                inter(338) =  0
                inter(339) =  0
                inter(340) =  0
                inter(341) =  0
                inter(342) =  0
                inter(343) =  0
                inter(344) =  0
                inter(345) =  0
                inter(346) =  0
                inter(347) =  0
                inter(348) =  0
                inter(349) =  0
                inter(350) =  0
                inter(351) =  0
                inter(352) =  0
                inter(353) =  0
                inter(354) =  0
                inter(355) =  0
                inter(356) =  0
                inter(357) =  0
                inter(358) =  0
                inter(359) =  0
                inter(360) =  0
                inter(361) =  0
                inter(362) =  0
                inter(363) =  0
                inter(364) =  0
                inter(365) =  0
                inter(366) =  0
                inter(367) =  0
                inter(368) =  0
                inter(369) =  0
                inter(370) =  0
                inter(371) =  0
                inter(372) =  0
                inter(373) =  0
                inter(374) =  0
                inter(375) =  0
                inter(376) =  0
                inter(377) =  0
                inter(378) =  0
                inter(379) =  1
                inter(380) =  0
                inter(381) =  0
                inter(382) =  0
                inter(383) =  0
                inter(384) =  0
                inter(385) =  0
                inter(386) =  0
                inter(387) =  0
                inter(388) =  x13*(x40 + x8*x85)
                inter(389) =  x13*x41 + x13*x42 + x19*x86
                inter(390) =  x13*x43 + x13*x44 + x19*x87
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  x13*x59 + x13*x60 + x29*x85
                inter(397) =  x13*(x25*x86 + x61)
                inter(398) =  x13*x62 + x13*x63 + x29*x87
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  x13*x75 + x13*x76 + x34*x85
                inter(405) =  x13*x77 + x13*x78 + x34*x86
                inter(406) =  x13*(x31*x87 + x79)
                inter(407) =  0
                inter(408) =  0
                inter(409) =  0
                inter(410) =  0
                inter(411) =  0
                inter(412) =  0
                inter(413) =  0
                inter(414) =  0
                inter(415) =  0
                inter(416) =  0
                inter(417) =  0
                inter(418) =  0
                inter(419) =  0
                inter(420) =  0
                inter(421) =  0
                inter(422) =  0
                inter(423) =  0
                inter(424) =  0
                inter(425) =  0
                inter(426) =  0
                inter(427) =  0
                inter(428) =  0
                inter(429) =  0
                inter(430) =  0
                inter(431) =  0
                inter(432) =  0
                inter(433) =  0
                inter(434) =  0
                inter(435) =  0
                inter(436) =  x13*(aj(1)*x45 - aj(1)*x89 - rj(1)*x94 - rj(1)*x96 + x0*x97 + x37*x85 + x38*x93 - x92*x93)
                inter(437) =  x13*(aj(2)*x45 - aj(2)*x89 - rj(2)*x94 - rj(2)*x96 + x1*x97 + x37*x86 + x38*x98 - x92*x98)
                inter(438) =  x13*(aj(3)*x45 - aj(3)*x89 - rj(3)*x94 - rj(3)*x96 + x2*x97 + x37*x87 + x38*x99 - x92*x99)
                inter(439) =  0
                inter(440) =  0
                inter(441) =  0
                inter(442) =  0
                inter(443) =  0
                inter(444) =  x100
                inter(445) =  x101
                inter(446) =  x102
                inter(447) =  0
                inter(448) =  0
                inter(449) =  0
                inter(450) =  0
                inter(451) =  0
                inter(452) =  x47
                inter(453) =  x49
                inter(454) =  x50
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  x65
                inter(461) =  x66
                inter(462) =  x67
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  x81
                inter(469) =  x82
                inter(470) =  x83
                inter(471) =  0
                inter(472) =  0
                inter(473) =  1
                inter(474) =  0
                inter(475) =  0
                inter(476) =  0
                inter(477) =  0
                inter(478) =  0
                inter(479) =  0
                inter(480) =  0
                inter(481) =  0
                inter(482) =  1
                inter(483) =  0
                inter(484) =  0
                inter(485) =  0
                inter(486) =  0
                inter(487) =  0
                inter(488) =  0
                inter(489) =  0
                inter(490) =  0
                inter(491) =  1
                inter(492) =  0
                inter(493) =  0
                inter(494) =  0
                inter(495) =  0
                inter(496) =  0
                inter(497) =  0
                inter(498) =  0
                inter(499) =  0
                inter(500) =  x100
                inter(501) =  x101
                inter(502) =  x102
                inter(503) =  0
                inter(504) =  0
                inter(505) =  0
                inter(506) =  0
                inter(507) =  0
                inter(508) =  0
                inter(509) =  0
                inter(510) =  0
                inter(511) =  0
                inter(512) =  0
                res = reshape(inter,[8,8,8])
            end function hes
            ! END AUTOCODE OUTPUT FOR HES
    end function
end module makemodel
