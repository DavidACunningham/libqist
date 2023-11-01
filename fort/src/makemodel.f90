module makemodel
    use, intrinsic :: iso_fortran_env, only: wp=>real64, dp=>real64, qp=>real128
    use astkindmodule
    use shdatamodule
    use pinesmodule
    use tensorops
    use cheby, only: spice_subset
    implicit none
    integer :: m=8

    type :: dynamicsModel
        type(spice_subset)    :: bod_db
        integer               :: num_bodies, &
                               & central_body, &
                               & traj_id
        integer,  allocatable :: bodylist(:) ! num_bodies
        logical               :: shgrav, &
                               & tgt_on_rails
        real(qp)              :: central_body_ref_radius, &
                               & central_body_rad(3), &
                               & state(8)
        real(qp), allocatable :: nbody_mus(:), &     ! (num_bodies)
                               & central_body_mu, &
                               & nbody_radii(:,:), & ! (3,num_bodies)
                               & nbody_vels(:,:), &  ! (3,num_bodies)
                               & nbody_accs(:,:)     ! (3,num_bodies)
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
            procedure :: fd_hes_kepler
            procedure :: fd_acc_nbody
            procedure :: fd_jac_nbody
            procedure :: fd_hes_nbody
            procedure :: allderivs_sh
            procedure :: trajstate
    end type dynamicsModel

    contains

    subroutine init_dm(me, subspice, traj_id, central_body, bodylist, &
                     & central_body_mu, central_body_ref_radius, mu_list, &
                     & shgrav, Cbar, Sbar, tgt_on_rails)
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
        class(dynamicsModel), intent(inout) :: me
        type(spice_subset),   intent(in)    :: subspice
        integer,              intent(in)    :: traj_id, & 
                                               central_body, &
                                               bodylist(:)
        real(dp),             intent(in)    :: central_body_ref_radius, &
                                               central_body_mu, &
                                               mu_list(:)
        logical,              intent(in)    :: shgrav, &
                                             & tgt_on_rails
        real(dp),             intent(in)    :: Cbar(:,:), &
                                               Sbar(:,:)
        me%num_bodies = size(bodylist)
        me%shgrav = shgrav
        me%traj_id = traj_id
        me%bod_db = subspice
        allocate(me%bodylist(me%num_bodies),  &
                 me%nbody_mus(me%num_bodies), &
                 me%nbody_radii(3,me%num_bodies), &
                 me%nbody_vels(3,me%num_bodies),  &
                 me%nbody_accs(3,me%num_bodies)   &
                )
        me%bodylist = bodylist
        me%central_body_ref_radius = real(central_body_ref_radius,qp)
        me%central_body_mu = real(central_body_mu,qp)
        me%nbody_mus = real(mu_list,qp)
        me%tgt_on_rails = tgt_on_rails

        if (me%shgrav) then
            call shdat_from_table( real(me%central_body_ref_radius,dp), & 
                                 & real(me%central_body_mu,dp), &
                                 & real(Cbar,dp), real(Sbar,dp), &
                                 & size(Cbar,1),size(Cbar,2), me%shdat &
                                 )
            call pinesinit(size(Cbar,1),me%shdat%Cml,me%shdat%Sml,me%pdat)
        endif
    end subroutine init_dm
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
        real(qp)                            :: r_up(6), &
                                               r_bod_up(3), v_bod_up(3), a_bod_up(3), &
                                               y(m)
        integer i
        ! at time t (in s past j2000)
            ! get r from central body to traj
        if (me%tgt_on_rails) then
            y(:3) = me%bod_db%call(real(time,dp),me%traj_id,'p')
            y(4:6) = me%bod_db%call(real(time,dp),me%traj_id,'v')
        else
            y = me%state
        endif
            ! for bodies i = 2. . .
        do i=1,me%num_bodies
            ! get r from central body to body i
            me%nbody_radii(:,i) = me%bod_db%call(real(time,dp), &
                                               & me%bodylist(i),'p')
            me%nbody_vels(:,i) = me%bod_db%call(real(time,dp), &
                                               & me%bodylist(i),'v')
            me%nbody_accs(:,i) = me%bod_db%call(real(time,dp), &
                                               & me%bodylist(i),'a')
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
        !$OMP PARALLEL DO
        do i=1,me%num_bodies
            r_bod_up = real(me%nbody_radii(:,i),qp)
            v_bod_up = real(me%nbody_vels(:,i),qp)
            a_bod_up = real(me%nbody_accs(:,i),qp)
            acc_nb   = acc_nb + acc_nbody(me, me%nbody_mus(i), r_up,r_bod_up)
            jac_nb   = jac_nb + jac_nbody(me, me%nbody_mus(i), r_up,r_bod_up,v_bod_up)
            hes_nb   = hes_nb + hes_nbody(me, me%nbody_mus(i), r_up,r_bod_up,v_bod_up,a_bod_up)
        end do
        !$OMP END PARALLEL DO
        acc = acc_2b + acc_nb
        jac = jac_2b + jac_nb
        hes = hes_2b + hes_nb
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

                x0  =  mu/(y(1)**2 + y(2)**2 + y(3)**2)**(3.0_qp/2.0_qp)

                res(1) =  y(4)
                res(2) =  y(5)
                res(3) =  y(6)
                res(4) =  -x0*y(1)
                res(5) =  -x0*y(2)
                res(6) =  -x0*y(3)
                res(7) =  1
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
                                      & x8, x9, x10, x11, x12, x13, x14
                real(qp), dimension(64) :: inter
                real(qp), dimension(8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  mu/x3**(3.0_qp/2.0_qp)
                x5  =  -x4
                x6  =  mu/x3**(5.0_qp/2.0_qp)
                x7  =  3*x6
                x8  =  3*y(1)
                x9  =  x6*x8
                x10  =  x9*y(2)
                x11  =  x9*y(3)
                x12  =  3*y(2)
                x13  =  x12*x6*y(3)
                x14  =  x6*(-x12*y(5) - x8*y(4) - 3*y(3)*y(6))

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x0*x7 + x5
                inter(5) =  x10
                inter(6) =  x11
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x10
                inter(13) =  x1*x7 + x5
                inter(14) =  x13
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x11
                inter(21) =  x13
                inter(22) =  x2*x7 + x5
                inter(23) =  0
                inter(24) =  0
                inter(25) =  1
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  1
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  1
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  -x4*y(1)
                inter(50) =  -x4*y(2)
                inter(51) =  -x4*y(3)
                inter(52) =  -x14*y(1) - x4*y(4)
                inter(53) =  -x14*y(2) - x4*y(5)
                inter(54) =  -x14*y(3) - x4*y(6)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  0
                inter(61) =  0
                inter(62) =  0
                inter(63) =  0
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function hes_kepler(me, mu, y) result (res)
        ! jac_kepler: method to compute Keplerian jacobian
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
                                      & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                      & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                      & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                      & x48, x49, x50, x51, x52, x53, x54, x55, & 
                                      & x56, x57, x58, x59, x60, x61, x62, x63, & 
                                      & x64, x65, x66, x67
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  mu/x3**(7.0_qp/2.0_qp)
                x6  =  15*x5
                x7  =  -3*mu*x4*y(2)
                x8  =  x0*x6
                x9  =  -x7 - x8*y(2)
                x10  =  -3*mu*x4*y(3)
                x11  =  -x10 - x8*y(3)
                x12  =  -3*mu*x4*y(1)
                x13  =  x6*y(1)
                x14  =  -x1*x13 - x12
                x15  =  y(2)*y(3)
                x16  =  -x13*x15
                x17  =  -x12 - x13*x2
                x18  =  mu/x3**(3.0_qp/2.0_qp)
                x19  =  -x18
                x20  =  mu*x4
                x21  =  3*x20
                x22  =  x0*x21 + x19
                x23  =  x20*y(1)
                x24  =  3*x23
                x25  =  x24*y(2)
                x26  =  x24*y(3)
                x27  =  x0*x5
                x28  =  y(1)*y(4)
                x29  =  y(2)*y(5)
                x30  =  y(3)*y(6)
                x31  =  -3*x28 - 3*x29 - 3*x30
                x32  =  5*x31
                x33  =  x20*x31
                x34  =  -x33
                x35  =  6*x20
                x36  =  x28*x35 + x34
                x37  =  x32*x5
                x38  =  x37*y(1)
                x39  =  x21*y(2)
                x40  =  x24*y(5) + x39*y(4)
                x41  =  x38*y(2) + x40
                x42  =  x21*y(3)
                x43  =  x24*y(6) + x42*y(4)
                x44  =  x38*y(3) + x43
                x45  =  -x1*x6*y(3) - x10
                x46  =  -x2*x6*y(2) - x7
                x47  =  x1*x21 + x19
                x48  =  x39*y(3)
                x49  =  x29*x35 + x34
                x50  =  x39*y(6) + x42*y(5)
                x51  =  x15*x37 + x50
                x52  =  x19 + x2*x21
                x53  =  x30*x35 + x34
                x54  =  -5*x28 - 5*x29 - 5*x30
                x55  =  3*x54
                x56  =  x5*x55
                x57  =  x56*y(1)
                x58  =  x40 + x57*y(2)
                x59  =  x43 + x57*y(3)
                x60  =  x15*x56 + x50
                x61  =  mu**2
                x62  =  x3**(-3)
                x63  =  2*x33
                x64  =  x31*x5*x54
                x65  =  3*x18
                x66  =  x0*x65 + x1*x65 + x2*x65 - 3*y(4)**2 - 3*y(5)**2 - 3*y(6)**2
                x67  =  x20*x66

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  9*mu*x4*y(1) - x6*y(1)**3
                inter(5) =  x9
                inter(6) =  x11
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x9
                inter(13) =  x14
                inter(14) =  x16
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x11
                inter(21) =  x16
                inter(22) =  x17
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
                inter(49) =  x22
                inter(50) =  x25
                inter(51) =  x26
                inter(52) =  x27*x32 + x36
                inter(53) =  x41
                inter(54) =  x44
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  0
                inter(61) =  0
                inter(62) =  0
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x9
                inter(69) =  x14
                inter(70) =  x16
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x14
                inter(77) =  9*mu*x4*y(2) - x6*y(2)**3
                inter(78) =  x45
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x16
                inter(85) =  x45
                inter(86) =  x46
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
                inter(113) =  x25
                inter(114) =  x47
                inter(115) =  x48
                inter(116) =  x41
                inter(117) =  x1*x37 + x49
                inter(118) =  x51
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  0
                inter(125) =  0
                inter(126) =  0
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x11
                inter(133) =  x16
                inter(134) =  x17
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x16
                inter(141) =  x45
                inter(142) =  x46
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x17
                inter(149) =  x46
                inter(150) =  9*mu*x4*y(3) - x6*y(3)**3
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
                inter(177) =  x26
                inter(178) =  x48
                inter(179) =  x52
                inter(180) =  x44
                inter(181) =  x51
                inter(182) =  x2*x37 + x53
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  0
                inter(189) =  0
                inter(190) =  0
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
                inter(244) =  x22
                inter(245) =  x25
                inter(246) =  x26
                inter(247) =  0
                inter(248) =  0
                inter(249) =  0
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
                inter(308) =  x25
                inter(309) =  x47
                inter(310) =  x48
                inter(311) =  0
                inter(312) =  0
                inter(313) =  0
                inter(314) =  0
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
                inter(372) =  x26
                inter(373) =  x48
                inter(374) =  x52
                inter(375) =  0
                inter(376) =  0
                inter(377) =  0
                inter(378) =  0
                inter(379) =  0
                inter(380) =  0
                inter(381) =  0
                inter(382) =  0
                inter(383) =  0
                inter(384) =  0
                inter(385) =  0
                inter(386) =  0
                inter(387) =  0
                inter(388) =  x27*x55 + x36
                inter(389) =  x58
                inter(390) =  x59
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  x58
                inter(397) =  x1*x56 + x49
                inter(398) =  x60
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  x59
                inter(405) =  x60
                inter(406) =  x2*x56 + x53
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
                inter(433) =  -x18*y(4) - x33*y(1)
                inter(434) =  -x18*y(5) - x33*y(2)
                inter(435) =  -x18*y(6) - x33*y(3)
                inter(436) =  -x23*x66 + x61*x62*y(1) - x63*y(4) - x64*y(1)
                inter(437) =  x61*x62*y(2) - x63*y(5) - x64*y(2) - x67*y(2)
                inter(438) =  x61*x62*y(3) - x63*y(6) - x64*y(3) - x67*y(3)
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
                inter(452) =  0
                inter(453) =  0
                inter(454) =  0
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  0
                inter(461) =  0
                inter(462) =  0
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  0
                inter(469) =  0
                inter(470) =  0
                inter(471) =  0
                inter(472) =  0
                inter(473) =  0
                inter(474) =  0
                inter(475) =  0
                inter(476) =  0
                inter(477) =  0
                inter(478) =  0
                inter(479) =  0
                inter(480) =  0
                inter(481) =  0
                inter(482) =  0
                inter(483) =  0
                inter(484) =  0
                inter(485) =  0
                inter(486) =  0
                inter(487) =  0
                inter(488) =  0
                inter(489) =  0
                inter(490) =  0
                inter(491) =  0
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
        real(qp)                         :: rj(3), mu_j
        rj = rbods; mu_j = mu;
        res = xdot(y)
        contains
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
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
                res(7) =  y(8)
                res(8) =  0
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
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27
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
                x20  =  y(8)**2
                x21  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x22  =  x21**(-3.0_qp/2.0_qp)
                x23  =  mu_j*(-rj(1)*x22 + x0*x4)
                x24  =  mu_j*(-rj(2)*x22 + x1*x4)
                x25  =  mu_j*(-rj(3)*x22 + x2*x4)
                x26  =  (-vj(1)*x6 - vj(2)*x14 - vj(3)*x17)/x21**(5.0_qp/2.0_qp)
                x27  =  -3.0_qp/2.0_qp*x0*(2*vj(1) - 2*y(4)) - 3.0_qp/2.0_qp*x1*(2*vj(2) - 2*y(5)) - &
                  3.0_qp/2.0_qp*x2*(2*vj(3) - 2*y(6))

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
                inter(49) =  x20*x23
                inter(50) =  x20*x24
                inter(51) =  x20*x25
                inter(52) =  x10*(-rj(1)*x26 - vj(1)*x22 + x0*x27*x8 + x4*(vj(1) - y(4)))
                inter(53) =  x10*(-rj(2)*x26 - vj(2)*x22 + x1*x27*x8 + x4*(vj(2) - y(5)))
                inter(54) =  x10*(-rj(3)*x26 - vj(3)*x22 + x2*x27*x8 + x4*(vj(3) - y(6)))
                inter(55) =  0
                inter(56) =  0
                inter(57) =  y(4)
                inter(58) =  y(5)
                inter(59) =  y(6)
                inter(60) =  x23
                inter(61) =  x24
                inter(62) =  x25
                inter(63) =  1
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
                                      & x96, x97, x98, x99, x100, x101, x102, x103, & 
                                      & x104, x105, x106, x107, x108, x109, x110, x111, & 
                                      & x112, x113, x114, x115, x116, x117, x118, x119, & 
                                      & x120, x121, x122
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  rj(1) - y(1)
                x1  =  rj(2) - y(2)
                x2  =  rj(3) - y(3)
                x3  =  x0**2 + x1**2 + x2**2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  x0*x4
                x6  =  3*x5
                x7  =  3*rj(1)
                x8  =  x7 - 3*y(1)
                x9  =  x4*x8
                x10  =  5*rj(1)
                x11  =  x10 - 5*y(1)
                x12  =  x3**(-7.0_qp/2.0_qp)
                x13  =  mu_j*y(8)
                x14  =  x1*x4
                x15  =  3*x14
                x16  =  -x13*x15
                x17  =  x11*x12
                x18  =  x1*x17
                x19  =  x13*x8
                x20  =  x2*x4
                x21  =  3*x20
                x22  =  -x13*x21
                x23  =  x17*x2
                x24  =  3*rj(2)
                x25  =  x24 - 3*y(2)
                x26  =  x25*x4
                x27  =  x13*x26
                x28  =  -x27
                x29  =  x0*x17
                x30  =  x13*x25
                x31  =  3*rj(3)
                x32  =  x31 - 3*y(3)
                x33  =  x32*x4
                x34  =  x13*x33
                x35  =  -x34
                x36  =  x13*x32
                x37  =  y(8)**2
                x38  =  x3**(-3.0_qp/2.0_qp)
                x39  =  -x38
                x40  =  mu_j*(x0*x9 + x39)
                x41  =  mu_j*x9
                x42  =  x1*x41
                x43  =  x2*x41
                x44  =  2*vj(1)
                x45  =  x44 - 2*y(4)
                x46  =  x0*x45
                x47  =  2*vj(2)
                x48  =  x47 - 2*y(5)
                x49  =  x1*x48
                x50  =  2*vj(3)
                x51  =  x50 - 2*y(6)
                x52  =  x2*x51
                x53  =  -3.0_qp/2.0_qp*x46 - 3.0_qp/2.0_qp*x49 - 3.0_qp/2.0_qp*x52
                x54  =  -x4*x53
                x55  =  3*vj(1) - 3*y(4)
                x56  =  vj(1) - y(4)
                x57  =  x5*x55 + x54 + x56*x9
                x58  =  vj(2) - y(5)
                x59  =  x58*x9
                x60  =  x4*x55
                x61  =  x1*x60
                x62  =  vj(3) - y(6)
                x63  =  x62*x9
                x64  =  x2*x60
                x65  =  5*rj(2)
                x66  =  x65 - 5*y(2)
                x67  =  x12*x66
                x68  =  x0*x67
                x69  =  x13*x9
                x70  =  -x69
                x71  =  x1*x67
                x72  =  x2*x67
                x73  =  -x13*x6
                x74  =  mu_j*x26
                x75  =  x0*x74
                x76  =  mu_j*(x1*x26 + x39)
                x77  =  x2*x74
                x78  =  3*vj(2) - 3*y(5)
                x79  =  x5*x78
                x80  =  x26*x56
                x81  =  x14*x78 + x26*x58 + x54
                x82  =  x26*x62
                x83  =  x20*x78
                x84  =  -x33
                x85  =  5*rj(3)
                x86  =  x85 - 5*y(3)
                x87  =  x12*x86
                x88  =  x0*x87
                x89  =  x1*x87
                x90  =  x2*x87
                x91  =  mu_j*x33
                x92  =  x0*x91
                x93  =  x1*x91
                x94  =  mu_j*(x2*x33 + x39)
                x95  =  3*vj(3) - 3*y(6)
                x96  =  x5*x95
                x97  =  x33*x56
                x98  =  x14*x95
                x99  =  x33*x58
                x100  =  x20*x95 + x33*x62 + x54
                x101  =  -5.0_qp/2.0_qp*x46 - 5.0_qp/2.0_qp*x49 - 5.0_qp/2.0_qp*x52
                x102  =  x101*x12
                x103  =  x0*x102
                x104  =  x1*x102
                x105  =  x102*x2
                x106  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x107  =  x106**(-3.0_qp/2.0_qp)
                x108  =  x106**(-5.0_qp/2.0_qp)
                x109  =  -vj(1)*x7 - vj(2)*x24 - vj(3)*x31
                x110  =  x108*x109
                x111  =  mu_j*(-rj(1)*x110 - vj(1)*x107 + x0*x4*x53 + x38*x56)
                x112  =  mu_j*(-rj(2)*x110 - vj(2)*x107 + x1*x4*x53 + x38*x58)
                x113  =  mu_j*(-rj(3)*x110 - vj(3)*x107 + x2*x4*x53 + x38*x62)
                x114  =  x109*(-vj(1)*x10 - vj(2)*x65 - vj(3)*x85)/x106**(7.0_qp/2.0_qp)
                x115  =  x108*(-aj(1)*x7 - aj(2)*x24 - aj(3)*x31 - 3*vj(1)**2 - 3*vj(2)**2 - 3*vj &
                  (3)**2)
                x116  =  x13*(-rj(1)*x107 + x0*x38)
                x117  =  2*x116
                x118  =  x13*(-rj(2)*x107 + x1*x38)
                x119  =  2*x118
                x120  =  x13*(-rj(3)*x107 + x2*x38)
                x121  =  2*x120
                x122  =  -3.0_qp/2.0_qp*x0*(2*aj(1) - x117) - 3.0_qp/2.0_qp*x1*(2*aj(2) - x119) - &
                  3.0_qp/2.0_qp*x2*(2*aj(3) - x121) - 3.0_qp/2.0_qp*x45*x56 - 3.0_qp/ &
                  2.0_qp*x48*x58 - 3.0_qp/2.0_qp*x51*x62

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
                inter(12) =  x28 + x29*x30
                inter(13) =  x13*(x1*x11*x12*x25 - x9)
                inter(14) =  x23*x30
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x29*x36 + x35
                inter(21) =  x18*x36
                inter(22) =  x13*(x11*x12*x2*x32 - x9)
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
                inter(49) =  x37*x40
                inter(50) =  x37*x42
                inter(51) =  x37*x43
                inter(52) =  x13*(x29*x53 + x57)
                inter(53) =  x13*(x18*x53 + x59 + x61)
                inter(54) =  x13*(x23*x53 + x63 + x64)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  x40
                inter(61) =  x42
                inter(62) =  x43
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x13*(-x26 + x68*x8)
                inter(69) =  x19*x71 + x70
                inter(70) =  x19*x72
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x30*x68 + x73
                inter(77) =  x13*(x1*x12*x25*x66 - x15 - 2*x26)
                inter(78) =  x22 + x30*x72
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x36*x68
                inter(85) =  x35 + x36*x71
                inter(86) =  x13*(x12*x2*x32*x66 - x26)
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
                inter(113) =  x37*x75
                inter(114) =  x37*x76
                inter(115) =  x37*x77
                inter(116) =  x13*(x53*x68 + x79 + x80)
                inter(117) =  x13*(x53*x71 + x81)
                inter(118) =  x13*(x53*x72 + x82 + x83)
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  x75
                inter(125) =  x76
                inter(126) =  x77
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x13*(x8*x88 + x84)
                inter(133) =  x19*x89
                inter(134) =  x19*x90 + x70
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x30*x88
                inter(141) =  x13*(x25*x89 + x84)
                inter(142) =  x28 + x30*x90
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x36*x88 + x73
                inter(149) =  x16 + x36*x89
                inter(150) =  x13*(x12*x2*x32*x86 - x21 - 2*x33)
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
                inter(177) =  x37*x92
                inter(178) =  x37*x93
                inter(179) =  x37*x94
                inter(180) =  x13*(x53*x88 + x96 + x97)
                inter(181) =  x13*(x53*x89 + x98 + x99)
                inter(182) =  x13*(x100 + x53*x90)
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  x92
                inter(189) =  x93
                inter(190) =  x94
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
                inter(244) =  x40*y(8)
                inter(245) =  x1*x69
                inter(246) =  x2*x69
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
                inter(308) =  x0*x27
                inter(309) =  x76*y(8)
                inter(310) =  x2*x27
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
                inter(372) =  x0*x34
                inter(373) =  x1*x34
                inter(374) =  x94*y(8)
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
                inter(388) =  x13*(x103*x8 + x57)
                inter(389) =  x104*x19 + x13*x59 + x13*x61
                inter(390) =  x105*x19 + x13*x63 + x13*x64
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  x103*x30 + x13*x79 + x13*x80
                inter(397) =  x13*(x104*x25 + x81)
                inter(398) =  x105*x30 + x13*x82 + x13*x83
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  x103*x36 + x13*x96 + x13*x97
                inter(405) =  x104*x36 + x13*x98 + x13*x99
                inter(406) =  x13*(x100 + x105*x32)
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
                inter(433) =  x111*x37
                inter(434) =  x112*x37
                inter(435) =  x113*x37
                inter(436) =  x13*(-aj(1)*x107 - rj(1)*x114 - rj(1)*x115 + x0*x101*x12*x53 + x0*x122*x4 - x110*x44 + x38*(aj(1) - x116) + 2*x4*x53*x56)
                inter(437) =  x13*(-aj(2)*x107 - rj(2)*x114 - rj(2)*x115 + x1*x101*x12*x53 + x1*x122*x4 - x110*x47 + x38*(aj(2) - x118) + 2*x4*x53*x58)
                inter(438) =  x13*(-aj(3)*x107 - rj(3)*x114 - rj(3)*x115 + x101*x12*x2*x53 - x110*x50 + x122*x2*x4 + x38*(aj(3) - x120) + 2*x4*x53*x62)
                inter(439) =  0
                inter(440) =  0
                inter(441) =  x116
                inter(442) =  x118
                inter(443) =  x120
                inter(444) =  x111
                inter(445) =  x112
                inter(446) =  x113
                inter(447) =  0
                inter(448) =  0
                inter(449) =  0
                inter(450) =  0
                inter(451) =  0
                inter(452) =  x40
                inter(453) =  x42
                inter(454) =  x43
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  x75
                inter(461) =  x76
                inter(462) =  x77
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  x92
                inter(469) =  x93
                inter(470) =  x94
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
                inter(497) =  x117
                inter(498) =  x119
                inter(499) =  x121
                inter(500) =  x111
                inter(501) =  x112
                inter(502) =  x113
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
        !                               258 = 8 + 8 ** 2 + 8 ** 3
        !                        y = [x, reshape(stm,(64)), reshape(stt,(584))]
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
        call me%get_derivs(t, acc, jac, hes)
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
        ! y              real (584)     Extended dyamics vector:
        !                               258 = 8 + 8 ** 2 + 8 ** 3
        !                        y = [x, reshape(stm,(64)), reshape(stt,(584))]
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
        call me%get_derivs(t, acc, jac, hes)
        statedot = acc
        stmdot = matmul(jac,stm)
        sttdot = mattens(jac,stt,m) + quad(stm,hes,m)
        res = [statedot, reshape(stmdot,[m**2]), reshape(sttdot,[m**3])]
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
        real(qp)                         :: t, mu
        mu = me%central_body_mu
        t = y(7)
        res = xdot(y)
        contains 
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(y) result(res)
                implicit none
                real(qp), intent(in) :: y (:)
                real(qp)             :: &
                                      & x0

                real(qp), dimension(8) :: res

                x0  =  mu/(y(1)**2 + y(2)**2 + y(3)**2)**(3.0_qp/2.0_qp)

                res(1) =  y(4)
                res(2) =  y(5)
                res(3) =  y(6)
                res(4) =  -x0*y(1)
                res(5) =  -x0*y(2)
                res(6) =  -x0*y(3)
                res(7) =  1
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
                                      & x8, x9, x10, x11, x12, x13, x14
                real(qp), dimension(64) :: inter
                real(qp), dimension(8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  mu/x3**(3.0_qp/2.0_qp)
                x5  =  -x4
                x6  =  mu/x3**(5.0_qp/2.0_qp)
                x7  =  3*x6
                x8  =  3*y(1)
                x9  =  x6*x8
                x10  =  x9*y(2)
                x11  =  x9*y(3)
                x12  =  3*y(2)
                x13  =  x12*x6*y(3)
                x14  =  x6*(-x12*y(5) - x8*y(4) - 3*y(3)*y(6))

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  x0*x7 + x5
                inter(5) =  x10
                inter(6) =  x11
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x10
                inter(13) =  x1*x7 + x5
                inter(14) =  x13
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x11
                inter(21) =  x13
                inter(22) =  x2*x7 + x5
                inter(23) =  0
                inter(24) =  0
                inter(25) =  1
                inter(26) =  0
                inter(27) =  0
                inter(28) =  0
                inter(29) =  0
                inter(30) =  0
                inter(31) =  0
                inter(32) =  0
                inter(33) =  0
                inter(34) =  1
                inter(35) =  0
                inter(36) =  0
                inter(37) =  0
                inter(38) =  0
                inter(39) =  0
                inter(40) =  0
                inter(41) =  0
                inter(42) =  0
                inter(43) =  1
                inter(44) =  0
                inter(45) =  0
                inter(46) =  0
                inter(47) =  0
                inter(48) =  0
                inter(49) =  -x4*y(1)
                inter(50) =  -x4*y(2)
                inter(51) =  -x4*y(3)
                inter(52) =  -x14*y(1) - x4*y(4)
                inter(53) =  -x14*y(2) - x4*y(5)
                inter(54) =  -x14*y(3) - x4*y(6)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  0
                inter(61) =  0
                inter(62) =  0
                inter(63) =  0
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function fd_hes_kepler(me, y) result (res)
        ! jac_kepler: method to compute Keplerian jacobian
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
        real(qp)                         :: res(size(y),size(y),size(y))
        real(qp)                         :: t, mu
        t = y(7)
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
                                      & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                      & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                      & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                      & x48, x49, x50, x51, x52, x53, x54, x55, & 
                                      & x56, x57, x58, x59, x60, x61, x62, x63, & 
                                      & x64, x65, x66, x67
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  y(1)**2
                x1  =  y(2)**2
                x2  =  y(3)**2
                x3  =  x0 + x1 + x2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  mu/x3**(7.0_qp/2.0_qp)
                x6  =  15*x5
                x7  =  -3*mu*x4*y(2)
                x8  =  x0*x6
                x9  =  -x7 - x8*y(2)
                x10  =  -3*mu*x4*y(3)
                x11  =  -x10 - x8*y(3)
                x12  =  -3*mu*x4*y(1)
                x13  =  x6*y(1)
                x14  =  -x1*x13 - x12
                x15  =  y(2)*y(3)
                x16  =  -x13*x15
                x17  =  -x12 - x13*x2
                x18  =  mu/x3**(3.0_qp/2.0_qp)
                x19  =  -x18
                x20  =  mu*x4
                x21  =  3*x20
                x22  =  x0*x21 + x19
                x23  =  x20*y(1)
                x24  =  3*x23
                x25  =  x24*y(2)
                x26  =  x24*y(3)
                x27  =  x0*x5
                x28  =  y(1)*y(4)
                x29  =  y(2)*y(5)
                x30  =  y(3)*y(6)
                x31  =  -3*x28 - 3*x29 - 3*x30
                x32  =  5*x31
                x33  =  x20*x31
                x34  =  -x33
                x35  =  6*x20
                x36  =  x28*x35 + x34
                x37  =  x32*x5
                x38  =  x37*y(1)
                x39  =  x21*y(2)
                x40  =  x24*y(5) + x39*y(4)
                x41  =  x38*y(2) + x40
                x42  =  x21*y(3)
                x43  =  x24*y(6) + x42*y(4)
                x44  =  x38*y(3) + x43
                x45  =  -x1*x6*y(3) - x10
                x46  =  -x2*x6*y(2) - x7
                x47  =  x1*x21 + x19
                x48  =  x39*y(3)
                x49  =  x29*x35 + x34
                x50  =  x39*y(6) + x42*y(5)
                x51  =  x15*x37 + x50
                x52  =  x19 + x2*x21
                x53  =  x30*x35 + x34
                x54  =  -5*x28 - 5*x29 - 5*x30
                x55  =  3*x54
                x56  =  x5*x55
                x57  =  x56*y(1)
                x58  =  x40 + x57*y(2)
                x59  =  x43 + x57*y(3)
                x60  =  x15*x56 + x50
                x61  =  mu**2
                x62  =  x3**(-3)
                x63  =  2*x33
                x64  =  x31*x5*x54
                x65  =  3*x18
                x66  =  x0*x65 + x1*x65 + x2*x65 - 3*y(4)**2 - 3*y(5)**2 - 3*y(6)**2
                x67  =  x20*x66

                inter(1) =  0
                inter(2) =  0
                inter(3) =  0
                inter(4) =  9*mu*x4*y(1) - x6*y(1)**3
                inter(5) =  x9
                inter(6) =  x11
                inter(7) =  0
                inter(8) =  0
                inter(9) =  0
                inter(10) =  0
                inter(11) =  0
                inter(12) =  x9
                inter(13) =  x14
                inter(14) =  x16
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x11
                inter(21) =  x16
                inter(22) =  x17
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
                inter(49) =  x22
                inter(50) =  x25
                inter(51) =  x26
                inter(52) =  x27*x32 + x36
                inter(53) =  x41
                inter(54) =  x44
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  0
                inter(61) =  0
                inter(62) =  0
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x9
                inter(69) =  x14
                inter(70) =  x16
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x14
                inter(77) =  9*mu*x4*y(2) - x6*y(2)**3
                inter(78) =  x45
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x16
                inter(85) =  x45
                inter(86) =  x46
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
                inter(113) =  x25
                inter(114) =  x47
                inter(115) =  x48
                inter(116) =  x41
                inter(117) =  x1*x37 + x49
                inter(118) =  x51
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  0
                inter(125) =  0
                inter(126) =  0
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x11
                inter(133) =  x16
                inter(134) =  x17
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x16
                inter(141) =  x45
                inter(142) =  x46
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x17
                inter(149) =  x46
                inter(150) =  9*mu*x4*y(3) - x6*y(3)**3
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
                inter(177) =  x26
                inter(178) =  x48
                inter(179) =  x52
                inter(180) =  x44
                inter(181) =  x51
                inter(182) =  x2*x37 + x53
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  0
                inter(189) =  0
                inter(190) =  0
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
                inter(244) =  x22
                inter(245) =  x25
                inter(246) =  x26
                inter(247) =  0
                inter(248) =  0
                inter(249) =  0
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
                inter(308) =  x25
                inter(309) =  x47
                inter(310) =  x48
                inter(311) =  0
                inter(312) =  0
                inter(313) =  0
                inter(314) =  0
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
                inter(372) =  x26
                inter(373) =  x48
                inter(374) =  x52
                inter(375) =  0
                inter(376) =  0
                inter(377) =  0
                inter(378) =  0
                inter(379) =  0
                inter(380) =  0
                inter(381) =  0
                inter(382) =  0
                inter(383) =  0
                inter(384) =  0
                inter(385) =  0
                inter(386) =  0
                inter(387) =  0
                inter(388) =  x27*x55 + x36
                inter(389) =  x58
                inter(390) =  x59
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  x58
                inter(397) =  x1*x56 + x49
                inter(398) =  x60
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  x59
                inter(405) =  x60
                inter(406) =  x2*x56 + x53
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
                inter(433) =  -x18*y(4) - x33*y(1)
                inter(434) =  -x18*y(5) - x33*y(2)
                inter(435) =  -x18*y(6) - x33*y(3)
                inter(436) =  -x23*x66 + x61*x62*y(1) - x63*y(4) - x64*y(1)
                inter(437) =  x61*x62*y(2) - x63*y(5) - x64*y(2) - x67*y(2)
                inter(438) =  x61*x62*y(3) - x63*y(6) - x64*y(3) - x67*y(3)
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
                inter(452) =  0
                inter(453) =  0
                inter(454) =  0
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  0
                inter(461) =  0
                inter(462) =  0
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  0
                inter(469) =  0
                inter(470) =  0
                inter(471) =  0
                inter(472) =  0
                inter(473) =  0
                inter(474) =  0
                inter(475) =  0
                inter(476) =  0
                inter(477) =  0
                inter(478) =  0
                inter(479) =  0
                inter(480) =  0
                inter(481) =  0
                inter(482) =  0
                inter(483) =  0
                inter(484) =  0
                inter(485) =  0
                inter(486) =  0
                inter(487) =  0
                inter(488) =  0
                inter(489) =  0
                inter(490) =  0
                inter(491) =  0
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
            function xdot(y,rj,mu_j) result(res)
                implicit none
                real(qp), intent(in) :: y (:), rj(3), mu_j
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
                res(7) =  y(8)
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
            thisr= me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'p')
            thisv= me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'v')
            res = res + jac(y, thisr, thisv,me%nbody_mus(i))
        end do
        contains
            ! BEGIN AUTOCODE OUTPUT FOR JAC
            function jac(y,rj,vj,mu_j) result(res)
                implicit none
                real(qp), intent(in) :: y (:), rj(3), vj(3), mu_j
                real(qp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27
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
                x20  =  y(8)**2
                x21  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x22  =  x21**(-3.0_qp/2.0_qp)
                x23  =  mu_j*(-rj(1)*x22 + x0*x4)
                x24  =  mu_j*(-rj(2)*x22 + x1*x4)
                x25  =  mu_j*(-rj(3)*x22 + x2*x4)
                x26  =  (-vj(1)*x6 - vj(2)*x14 - vj(3)*x17)/x21**(5.0_qp/2.0_qp)
                x27  =  -3.0_qp/2.0_qp*x0*(2*vj(1) - 2*y(4)) - 3.0_qp/2.0_qp*x1*(2*vj(2) - 2*y(5)) - &
                  3.0_qp/2.0_qp*x2*(2*vj(3) - 2*y(6))

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
                inter(49) =  x20*x23
                inter(50) =  x20*x24
                inter(51) =  x20*x25
                inter(52) =  x10*(-rj(1)*x26 - vj(1)*x22 + x0*x27*x8 + x4*(vj(1) - y(4)))
                inter(53) =  x10*(-rj(2)*x26 - vj(2)*x22 + x1*x27*x8 + x4*(vj(2) - y(5)))
                inter(54) =  x10*(-rj(3)*x26 - vj(3)*x22 + x2*x27*x8 + x4*(vj(3) - y(6)))
                inter(55) =  0
                inter(56) =  0
                inter(57) =  y(4)
                inter(58) =  y(5)
                inter(59) =  y(6)
                inter(60) =  x23
                inter(61) =  x24
                inter(62) =  x25
                inter(63) =  1
                inter(64) =  0
                res = reshape(inter,[8,8])
            end function jac
            ! END AUTOCODE OUTPUT FOR JAC
    end function
    function fd_hes_nbody(me, y) result (res)
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
        real(qp),             intent(in) :: y(:)
        real(qp)                         :: res(size(y),size(y),size(y))
        real(qp)                         :: thisr(3), thisv(3), thisa(3)
        real(qp)                         :: t
        integer                          :: i
        t = y(7)
        res = 0._qp
        do i=1,me%num_bodies
            thisr= me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'p')
            thisv= me%bod_db%call(real(t,dp), &
                               & me%bodylist(i),'v')
            res = res + hes(y, thisr, thisv, thisa, me%nbody_mus(i))
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
                                      & x96, x97, x98, x99, x100, x101, x102, x103, & 
                                      & x104, x105, x106, x107, x108, x109, x110, x111, & 
                                      & x112, x113, x114, x115, x116, x117, x118, x119, & 
                                      & x120, x121, x122
                real(qp), dimension(512) :: inter
                real(qp), dimension(8,8,8) :: res

                x0  =  rj(1) - y(1)
                x1  =  rj(2) - y(2)
                x2  =  rj(3) - y(3)
                x3  =  x0**2 + x1**2 + x2**2
                x4  =  x3**(-5.0_qp/2.0_qp)
                x5  =  x0*x4
                x6  =  3*x5
                x7  =  3*rj(1)
                x8  =  x7 - 3*y(1)
                x9  =  x4*x8
                x10  =  5*rj(1)
                x11  =  x10 - 5*y(1)
                x12  =  x3**(-7.0_qp/2.0_qp)
                x13  =  mu_j*y(8)
                x14  =  x1*x4
                x15  =  3*x14
                x16  =  -x13*x15
                x17  =  x11*x12
                x18  =  x1*x17
                x19  =  x13*x8
                x20  =  x2*x4
                x21  =  3*x20
                x22  =  -x13*x21
                x23  =  x17*x2
                x24  =  3*rj(2)
                x25  =  x24 - 3*y(2)
                x26  =  x25*x4
                x27  =  x13*x26
                x28  =  -x27
                x29  =  x0*x17
                x30  =  x13*x25
                x31  =  3*rj(3)
                x32  =  x31 - 3*y(3)
                x33  =  x32*x4
                x34  =  x13*x33
                x35  =  -x34
                x36  =  x13*x32
                x37  =  y(8)**2
                x38  =  x3**(-3.0_qp/2.0_qp)
                x39  =  -x38
                x40  =  mu_j*(x0*x9 + x39)
                x41  =  mu_j*x9
                x42  =  x1*x41
                x43  =  x2*x41
                x44  =  2*vj(1)
                x45  =  x44 - 2*y(4)
                x46  =  x0*x45
                x47  =  2*vj(2)
                x48  =  x47 - 2*y(5)
                x49  =  x1*x48
                x50  =  2*vj(3)
                x51  =  x50 - 2*y(6)
                x52  =  x2*x51
                x53  =  -3.0_qp/2.0_qp*x46 - 3.0_qp/2.0_qp*x49 - 3.0_qp/2.0_qp*x52
                x54  =  -x4*x53
                x55  =  3*vj(1) - 3*y(4)
                x56  =  vj(1) - y(4)
                x57  =  x5*x55 + x54 + x56*x9
                x58  =  vj(2) - y(5)
                x59  =  x58*x9
                x60  =  x4*x55
                x61  =  x1*x60
                x62  =  vj(3) - y(6)
                x63  =  x62*x9
                x64  =  x2*x60
                x65  =  5*rj(2)
                x66  =  x65 - 5*y(2)
                x67  =  x12*x66
                x68  =  x0*x67
                x69  =  x13*x9
                x70  =  -x69
                x71  =  x1*x67
                x72  =  x2*x67
                x73  =  -x13*x6
                x74  =  mu_j*x26
                x75  =  x0*x74
                x76  =  mu_j*(x1*x26 + x39)
                x77  =  x2*x74
                x78  =  3*vj(2) - 3*y(5)
                x79  =  x5*x78
                x80  =  x26*x56
                x81  =  x14*x78 + x26*x58 + x54
                x82  =  x26*x62
                x83  =  x20*x78
                x84  =  -x33
                x85  =  5*rj(3)
                x86  =  x85 - 5*y(3)
                x87  =  x12*x86
                x88  =  x0*x87
                x89  =  x1*x87
                x90  =  x2*x87
                x91  =  mu_j*x33
                x92  =  x0*x91
                x93  =  x1*x91
                x94  =  mu_j*(x2*x33 + x39)
                x95  =  3*vj(3) - 3*y(6)
                x96  =  x5*x95
                x97  =  x33*x56
                x98  =  x14*x95
                x99  =  x33*x58
                x100  =  x20*x95 + x33*x62 + x54
                x101  =  -5.0_qp/2.0_qp*x46 - 5.0_qp/2.0_qp*x49 - 5.0_qp/2.0_qp*x52
                x102  =  x101*x12
                x103  =  x0*x102
                x104  =  x1*x102
                x105  =  x102*x2
                x106  =  rj(1)**2 + rj(2)**2 + rj(3)**2
                x107  =  x106**(-3.0_qp/2.0_qp)
                x108  =  x106**(-5.0_qp/2.0_qp)
                x109  =  -vj(1)*x7 - vj(2)*x24 - vj(3)*x31
                x110  =  x108*x109
                x111  =  mu_j*(-rj(1)*x110 - vj(1)*x107 + x0*x4*x53 + x38*x56)
                x112  =  mu_j*(-rj(2)*x110 - vj(2)*x107 + x1*x4*x53 + x38*x58)
                x113  =  mu_j*(-rj(3)*x110 - vj(3)*x107 + x2*x4*x53 + x38*x62)
                x114  =  x109*(-vj(1)*x10 - vj(2)*x65 - vj(3)*x85)/x106**(7.0_qp/2.0_qp)
                x115  =  x108*(-aj(1)*x7 - aj(2)*x24 - aj(3)*x31 - 3*vj(1)**2 - 3*vj(2)**2 - 3*vj &
                  (3)**2)
                x116  =  x13*(-rj(1)*x107 + x0*x38)
                x117  =  2*x116
                x118  =  x13*(-rj(2)*x107 + x1*x38)
                x119  =  2*x118
                x120  =  x13*(-rj(3)*x107 + x2*x38)
                x121  =  2*x120
                x122  =  -3.0_qp/2.0_qp*x0*(2*aj(1) - x117) - 3.0_qp/2.0_qp*x1*(2*aj(2) - x119) - &
                  3.0_qp/2.0_qp*x2*(2*aj(3) - x121) - 3.0_qp/2.0_qp*x45*x56 - 3.0_qp/ &
                  2.0_qp*x48*x58 - 3.0_qp/2.0_qp*x51*x62

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
                inter(12) =  x28 + x29*x30
                inter(13) =  x13*(x1*x11*x12*x25 - x9)
                inter(14) =  x23*x30
                inter(15) =  0
                inter(16) =  0
                inter(17) =  0
                inter(18) =  0
                inter(19) =  0
                inter(20) =  x29*x36 + x35
                inter(21) =  x18*x36
                inter(22) =  x13*(x11*x12*x2*x32 - x9)
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
                inter(49) =  x37*x40
                inter(50) =  x37*x42
                inter(51) =  x37*x43
                inter(52) =  x13*(x29*x53 + x57)
                inter(53) =  x13*(x18*x53 + x59 + x61)
                inter(54) =  x13*(x23*x53 + x63 + x64)
                inter(55) =  0
                inter(56) =  0
                inter(57) =  0
                inter(58) =  0
                inter(59) =  0
                inter(60) =  x40
                inter(61) =  x42
                inter(62) =  x43
                inter(63) =  0
                inter(64) =  0
                inter(65) =  0
                inter(66) =  0
                inter(67) =  0
                inter(68) =  x13*(-x26 + x68*x8)
                inter(69) =  x19*x71 + x70
                inter(70) =  x19*x72
                inter(71) =  0
                inter(72) =  0
                inter(73) =  0
                inter(74) =  0
                inter(75) =  0
                inter(76) =  x30*x68 + x73
                inter(77) =  x13*(x1*x12*x25*x66 - x15 - 2*x26)
                inter(78) =  x22 + x30*x72
                inter(79) =  0
                inter(80) =  0
                inter(81) =  0
                inter(82) =  0
                inter(83) =  0
                inter(84) =  x36*x68
                inter(85) =  x35 + x36*x71
                inter(86) =  x13*(x12*x2*x32*x66 - x26)
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
                inter(113) =  x37*x75
                inter(114) =  x37*x76
                inter(115) =  x37*x77
                inter(116) =  x13*(x53*x68 + x79 + x80)
                inter(117) =  x13*(x53*x71 + x81)
                inter(118) =  x13*(x53*x72 + x82 + x83)
                inter(119) =  0
                inter(120) =  0
                inter(121) =  0
                inter(122) =  0
                inter(123) =  0
                inter(124) =  x75
                inter(125) =  x76
                inter(126) =  x77
                inter(127) =  0
                inter(128) =  0
                inter(129) =  0
                inter(130) =  0
                inter(131) =  0
                inter(132) =  x13*(x8*x88 + x84)
                inter(133) =  x19*x89
                inter(134) =  x19*x90 + x70
                inter(135) =  0
                inter(136) =  0
                inter(137) =  0
                inter(138) =  0
                inter(139) =  0
                inter(140) =  x30*x88
                inter(141) =  x13*(x25*x89 + x84)
                inter(142) =  x28 + x30*x90
                inter(143) =  0
                inter(144) =  0
                inter(145) =  0
                inter(146) =  0
                inter(147) =  0
                inter(148) =  x36*x88 + x73
                inter(149) =  x16 + x36*x89
                inter(150) =  x13*(x12*x2*x32*x86 - x21 - 2*x33)
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
                inter(177) =  x37*x92
                inter(178) =  x37*x93
                inter(179) =  x37*x94
                inter(180) =  x13*(x53*x88 + x96 + x97)
                inter(181) =  x13*(x53*x89 + x98 + x99)
                inter(182) =  x13*(x100 + x53*x90)
                inter(183) =  0
                inter(184) =  0
                inter(185) =  0
                inter(186) =  0
                inter(187) =  0
                inter(188) =  x92
                inter(189) =  x93
                inter(190) =  x94
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
                inter(244) =  x40*y(8)
                inter(245) =  x1*x69
                inter(246) =  x2*x69
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
                inter(308) =  x0*x27
                inter(309) =  x76*y(8)
                inter(310) =  x2*x27
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
                inter(372) =  x0*x34
                inter(373) =  x1*x34
                inter(374) =  x94*y(8)
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
                inter(388) =  x13*(x103*x8 + x57)
                inter(389) =  x104*x19 + x13*x59 + x13*x61
                inter(390) =  x105*x19 + x13*x63 + x13*x64
                inter(391) =  0
                inter(392) =  0
                inter(393) =  0
                inter(394) =  0
                inter(395) =  0
                inter(396) =  x103*x30 + x13*x79 + x13*x80
                inter(397) =  x13*(x104*x25 + x81)
                inter(398) =  x105*x30 + x13*x82 + x13*x83
                inter(399) =  0
                inter(400) =  0
                inter(401) =  0
                inter(402) =  0
                inter(403) =  0
                inter(404) =  x103*x36 + x13*x96 + x13*x97
                inter(405) =  x104*x36 + x13*x98 + x13*x99
                inter(406) =  x13*(x100 + x105*x32)
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
                inter(433) =  x111*x37
                inter(434) =  x112*x37
                inter(435) =  x113*x37
                inter(436) =  x13*(-aj(1)*x107 - rj(1)*x114 - rj(1)*x115 + x0*x101*x12*x53 + x0*x122*x4 - x110*x44 + x38*(aj(1) - x116) + 2*x4*x53*x56)
                inter(437) =  x13*(-aj(2)*x107 - rj(2)*x114 - rj(2)*x115 + x1*x101*x12*x53 + x1*x122*x4 - x110*x47 + x38*(aj(2) - x118) + 2*x4*x53*x58)
                inter(438) =  x13*(-aj(3)*x107 - rj(3)*x114 - rj(3)*x115 + x101*x12*x2*x53 - x110*x50 + x122*x2*x4 + x38*(aj(3) - x120) + 2*x4*x53*x62)
                inter(439) =  0
                inter(440) =  0
                inter(441) =  x116
                inter(442) =  x118
                inter(443) =  x120
                inter(444) =  x111
                inter(445) =  x112
                inter(446) =  x113
                inter(447) =  0
                inter(448) =  0
                inter(449) =  0
                inter(450) =  0
                inter(451) =  0
                inter(452) =  x40
                inter(453) =  x42
                inter(454) =  x43
                inter(455) =  0
                inter(456) =  0
                inter(457) =  0
                inter(458) =  0
                inter(459) =  0
                inter(460) =  x75
                inter(461) =  x76
                inter(462) =  x77
                inter(463) =  0
                inter(464) =  0
                inter(465) =  0
                inter(466) =  0
                inter(467) =  0
                inter(468) =  x92
                inter(469) =  x93
                inter(470) =  x94
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
                inter(497) =  x117
                inter(498) =  x119
                inter(499) =  x121
                inter(500) =  x111
                inter(501) =  x112
                inter(502) =  x113
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
