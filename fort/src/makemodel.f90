module makemodel
    use globals
    use astkindmodule
    use shdatamodule
    use pinesmodule
    use tensorops
    implicit none

    type :: dynamicsModel
        integer               :: num_bodies, &
                               & central_body, &
                               & traj_id
        integer,  allocatable :: bodylist(:) ! num_bodies
        logical               :: shgrav
        real(dp)              :: central_body_ref_radius, &
                               & central_body_rad(3)
        real(wp), allocatable :: nbody_mus(:), &      ! (num_bodies)
                               & central_body_mu, &
                               & nbody_radii(:,:)     ! (3,num_bodies)
        type(SHarmData)       :: shdat
        type(PinesData)       :: pdat
        character(len=256)    :: kernelfile
        contains
            procedure:: init => init_dm
            procedure :: eoms
            procedure :: get_derivs
            procedure :: acc_kepler
            procedure :: jac_kepler
            procedure :: hes_kepler
            procedure :: acc_nbody
            procedure :: jac_nbody
            procedure :: hes_nbody
            procedure :: allderivs_sh
            procedure :: trajstate
    end type dynamicsModel

    contains

    subroutine init_dm(me, kernelfile, traj_id, central_body, bodylist, &
                     & shgrav, Cbar, Sbar)
        ! init_dm: method to initialize dynamicsModel object
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! kernelfile     character      The absolute directory path for the 
        !                               SPICE metakernel
        ! traj_id        integer        the integer ID of the reference 
        !                               trajectory in the SPICE kernel
        ! central_body   integer        the integer ID of the selected central
        !                               body in the SPICE kernel
        ! bodylist       integer (:)    list of integer IDs of the gravitating
        !                               bodies in the SPICE kernel, other than
        !                               the central body
        ! shgrav         logical        TRUE if central body will be modeled
        !                               by spherical harmonics
        ! Cbar           real (:,:)     4 pi (Kaula) normalized cosine Stokes 
        !                               for central body
        ! Sbar           real (:,:)     4 pi (Kaula) normalized sine Stokes 
        !                               for central body
        ! OUTPUTS:
        ! NONE
        class(dynamicsModel), intent(inout) :: me
        integer,              intent(in)    :: traj_id, & 
                                               central_body, &
                                               bodylist(:)
        logical,              intent(in)    :: shgrav
        real(wp),             intent(in)    :: Cbar(:,:), &
                                               Sbar(:,:)
        character(len=256),   intent(in)    :: kernelfile
        real(dp)                            :: mudum(1)
        real(dp)                            :: radii(3)
        integer                             :: i,N
        me%num_bodies = size(bodylist)
        me%shgrav = shgrav
        me%traj_id = traj_id
        me%kernelfile = kernelfile
        call furnsh(me%kernelfile)
        allocate(me%bodylist(me%num_bodies),  &
                 me%nbody_mus(me%num_bodies), &
                 me%nbody_radii(3,me%num_bodies) &
                )
        me%bodylist = bodylist
        call bodvcd(central_body, 'RADII', 3, N, radii)
        call bodvcd(central_body, 'GM', 1, N, mudum)
        me%central_body_ref_radius = radii(1)
        me%central_body_mu = real(mudum(1), wp)
        do i=1,me%num_bodies
            call bodvcd(bodylist(i), 'GM', 1, N, mudum)
            me%nbody_mus(i) = real(mudum(1), wp)
        end do

        if (me%shgrav) then
            call shdat_from_table( me%central_body_ref_radius, & 
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
        real(wp),             intent(in)    :: time
        real(wp),             intent(out)   :: acc(6), jac(6,6), hes(6,6,6)
        real(wp)                            :: acc_2b(6), jac_2b(6,6), hes_2b(6,6,6)
        real(wp)                            :: acc_nb(6), jac_nb(6,6), hes_nb(6,6,6)
        real(wp)                            :: r_up(6), r_bod_up(3)
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
            call spkgps(me%traj_id, time, "J2000", me%central_body, &
                      & me%nbody_radii(:,i), &     ! (3,num_bodies)
                      & lt_dum)
            ! get r from central body to body i (SPKGPS)

        end do

        ! For central body
            ! Choose point mass or SH model
        r_up = real(traj_state,wp)
        if (me%shgrav) then
            ! PLACEHOLDER
            ! TODO: Put in extended body gravity
            acc_2b = acc_kepler(me, me%central_body_mu, r_up)
            jac_2b = jac_kepler(me, me%central_body_mu, r_up)
            hes_2b = hes_kepler(me, me%central_body_mu, r_up)
            ! END PLACEHOLDER
        else
            acc_2b = acc_kepler(me, me%central_body_mu, r_up)
            jac_2b = jac_kepler(me, me%central_body_mu, r_up)
            hes_2b = hes_kepler(me, me%central_body_mu, r_up)
        endif
        ! for bodies i = 2. . . N NOTE: CAN BE DONE IN ||
        ! Compute point mass accel, jac, hess for i=1
        !$OMP PARALLEL DO
        do i=1,me%num_bodies
            r_bod_up = real(me%nbody_radii(:,i),wp)
            acc_nb   = acc_nb + acc_nbody(me, me%nbody_mus(i), r_up,r_bod_up)
            jac_nb   = jac_nb + jac_nbody(me, me%nbody_mus(i), r_up,r_bod_up)
            hes_nb   = hes_nb + hes_nbody(me, me%nbody_mus(i), r_up,r_bod_up)
        end do
        !$OMP END PARALLEL DO
        acc = acc_2b + acc_nb
        jac = jac_2b + jac_nb
        hes = hes_2b + hes_nb
    end subroutine get_derivs
    function acc_kepler(me, mu, r) result (res)
        ! acc_kepler: method to compute Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! r              real (3)       Distance from gravitating body to
        !                               field point
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6)       Keplerian acceleration vector with
        !                               zero padding (0, 0, 0, ax, ay, az)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3)
        real(wp)                         :: res(6)
        res = xdot(r)
        contains 
            function xdot(x) result(res)
                real(wp), intent(in) :: x (:)
                real(wp)             :: &
                                      & x0

                real(wp), dimension(6) :: res

                x0  =  mu/(x(1)**2 + x(2)**2 + x(3)**2)**(3.0_wp/2.0_wp)

                res(1) =  x(4)
                res(2) =  x(5)
                res(3) =  x(6)
                res(4) =  -x(1)*x0
                res(5) =  -x(2)*x0
                res(6) =  -x(3)*x0
                end function xdot
    end function
    function jac_kepler(me, mu, r) result (res)
        ! jac_kepler: method to compute Keplerian jacobian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! r              real (3)       Distance from gravitating body to
        !                               field point
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6,6)     Keplerian Jacobian matrix 
        !                                         df_i
        !                               J(i,j) = ------
        !                                         dx_j
        !                               
        !                               [  dr     dv  ]
        !                               [ ----   ---- ]
        !                               [  dr     dr  ]
        !                               [             ]
        !                               [  dv     dv  ]
        !                               [ ----   ---- ]
        !                               [  dr     dv  ]
        !                               or, if A is the gradient of 
        !                               Keplerian acceleration f with respect
        !                               to position vector r:
        !                               [ 0  I ]
        !                               [ A  0 ]
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3)
        real(wp)                         :: res(6,6)
        res = reshape(jac(r),[6,6])
        contains
            function jac(x) result(res)
                implicit none
                real(wp), intent(in) :: x (:)
                real(wp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9

                real(wp), dimension(36) :: res

                x0  =  x(1)**2
                x1  =  x(2)**2
                x2  =  x(3)**2
                x3  =  x0 + x1 + x2
                x4  =  -mu/x3**(3.0_wp/2.0_wp)
                x5  =  3*mu/x3**(5.0_wp/2.0_wp)
                x6  =  x(1)*x5
                x7  =  x(2)*x6
                x8  =  x(3)*x6
                x9  =  x(2)*x(3)*x5

                res(1) =  0
                res(2) =  0
                res(3) =  0
                res(4) =  x0*x5 + x4
                res(5) =  x7
                res(6) =  x8
                res(7) =  0
                res(8) =  0
                res(9) =  0
                res(10) =  x7
                res(11) =  x1*x5 + x4
                res(12) =  x9
                res(13) =  0
                res(14) =  0
                res(15) =  0
                res(16) =  x8
                res(17) =  x9
                res(18) =  x2*x5 + x4
                res(19) =  1
                res(20) =  0
                res(21) =  0
                res(22) =  0
                res(23) =  0
                res(24) =  0
                res(25) =  0
                res(26) =  1
                res(27) =  0
                res(28) =  0
                res(29) =  0
                res(30) =  0
                res(31) =  0
                res(32) =  0
                res(33) =  1
                res(34) =  0
                res(35) =  0
                res(36) =  0
            end function jac
    end function
    function hes_kepler(me, mu, r) result (res)
        ! jac_kepler: method to compute Keplerian jacobian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of body
        ! r              real (3)       position vector from gravitating body to
        !                               field point
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6,6,6)   Keplerian Hessian tensor H where
        !                                              df_i 
        !                               H(i,j,k) = -----------
        !                                           dx_j dx_k
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3)
        real(wp)                         :: res(6,6,6)
        res = reshape(hes(r),[6,6,6])
        contains
            ! BEGIN AUTOCODE OUTPUT FOR HES
            function hes(x) result(res)
                implicit none
                real(wp), intent(in) :: x (:)
                real(wp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23

                real(wp), dimension(216) :: res

                x0  =  x(1)**2
                x1  =  x(2)**2
                x2  =  x(3)**2
                x3  =  x0 + x1 + x2
                x4  =  5/x3
                x5  =  x0*x4
                x6  =  mu*x(1)
                x7  =  3/x3**(5.0_wp/2.0_wp)
                x8  =  x6*x7
                x9  =  mu*x7
                x10  =  x9*(1 - x5)
                x11  =  x(2)*x10
                x12  =  x(3)*x10
                x13  =  x1*x4
                x14  =  1 - x13
                x15  =  x14*x8
                x16  =  -15*x(2)*x(3)*x6/x3**(7.0_wp/2.0_wp)
                x17  =  x2*x4
                x18  =  1 - x17
                x19  =  x18*x8
                x20  =  x(2)*x9
                x21  =  x(3)*x9
                x22  =  x14*x21
                x23  =  x18*x20

                res(1) =  0
                res(2) =  0
                res(3) =  0
                res(4) =  x8*(3 - x5)
                res(5) =  x11
                res(6) =  x12
                res(7) =  0
                res(8) =  0
                res(9) =  0
                res(10) =  x11
                res(11) =  x15
                res(12) =  x16
                res(13) =  0
                res(14) =  0
                res(15) =  0
                res(16) =  x12
                res(17) =  x16
                res(18) =  x19
                res(19) =  0
                res(20) =  0
                res(21) =  0
                res(22) =  0
                res(23) =  0
                res(24) =  0
                res(25) =  0
                res(26) =  0
                res(27) =  0
                res(28) =  0
                res(29) =  0
                res(30) =  0
                res(31) =  0
                res(32) =  0
                res(33) =  0
                res(34) =  0
                res(35) =  0
                res(36) =  0
                res(37) =  0
                res(38) =  0
                res(39) =  0
                res(40) =  x11
                res(41) =  x15
                res(42) =  x16
                res(43) =  0
                res(44) =  0
                res(45) =  0
                res(46) =  x15
                res(47) =  x20*(3 - x13)
                res(48) =  x22
                res(49) =  0
                res(50) =  0
                res(51) =  0
                res(52) =  x16
                res(53) =  x22
                res(54) =  x23
                res(55) =  0
                res(56) =  0
                res(57) =  0
                res(58) =  0
                res(59) =  0
                res(60) =  0
                res(61) =  0
                res(62) =  0
                res(63) =  0
                res(64) =  0
                res(65) =  0
                res(66) =  0
                res(67) =  0
                res(68) =  0
                res(69) =  0
                res(70) =  0
                res(71) =  0
                res(72) =  0
                res(73) =  0
                res(74) =  0
                res(75) =  0
                res(76) =  x12
                res(77) =  x16
                res(78) =  x19
                res(79) =  0
                res(80) =  0
                res(81) =  0
                res(82) =  x16
                res(83) =  x22
                res(84) =  x23
                res(85) =  0
                res(86) =  0
                res(87) =  0
                res(88) =  x19
                res(89) =  x23
                res(90) =  x21*(3 - x17)
                res(91) =  0
                res(92) =  0
                res(93) =  0
                res(94) =  0
                res(95) =  0
                res(96) =  0
                res(97) =  0
                res(98) =  0
                res(99) =  0
                res(100) =  0
                res(101) =  0
                res(102) =  0
                res(103) =  0
                res(104) =  0
                res(105) =  0
                res(106) =  0
                res(107) =  0
                res(108) =  0
                res(109) =  0
                res(110) =  0
                res(111) =  0
                res(112) =  0
                res(113) =  0
                res(114) =  0
                res(115) =  0
                res(116) =  0
                res(117) =  0
                res(118) =  0
                res(119) =  0
                res(120) =  0
                res(121) =  0
                res(122) =  0
                res(123) =  0
                res(124) =  0
                res(125) =  0
                res(126) =  0
                res(127) =  0
                res(128) =  0
                res(129) =  0
                res(130) =  0
                res(131) =  0
                res(132) =  0
                res(133) =  0
                res(134) =  0
                res(135) =  0
                res(136) =  0
                res(137) =  0
                res(138) =  0
                res(139) =  0
                res(140) =  0
                res(141) =  0
                res(142) =  0
                res(143) =  0
                res(144) =  0
                res(145) =  0
                res(146) =  0
                res(147) =  0
                res(148) =  0
                res(149) =  0
                res(150) =  0
                res(151) =  0
                res(152) =  0
                res(153) =  0
                res(154) =  0
                res(155) =  0
                res(156) =  0
                res(157) =  0
                res(158) =  0
                res(159) =  0
                res(160) =  0
                res(161) =  0
                res(162) =  0
                res(163) =  0
                res(164) =  0
                res(165) =  0
                res(166) =  0
                res(167) =  0
                res(168) =  0
                res(169) =  0
                res(170) =  0
                res(171) =  0
                res(172) =  0
                res(173) =  0
                res(174) =  0
                res(175) =  0
                res(176) =  0
                res(177) =  0
                res(178) =  0
                res(179) =  0
                res(180) =  0
                res(181) =  0
                res(182) =  0
                res(183) =  0
                res(184) =  0
                res(185) =  0
                res(186) =  0
                res(187) =  0
                res(188) =  0
                res(189) =  0
                res(190) =  0
                res(191) =  0
                res(192) =  0
                res(193) =  0
                res(194) =  0
                res(195) =  0
                res(196) =  0
                res(197) =  0
                res(198) =  0
                res(199) =  0
                res(200) =  0
                res(201) =  0
                res(202) =  0
                res(203) =  0
                res(204) =  0
                res(205) =  0
                res(206) =  0
                res(207) =  0
                res(208) =  0
                res(209) =  0
                res(210) =  0
                res(211) =  0
                res(212) =  0
                res(213) =  0
                res(214) =  0
                res(215) =  0
                res(216) =  0
            end function hes
            ! END AUTOCODE OUTPUT FOR HES
    end function
    function acc_nbody(me, mu, r, rbods) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! r              real (3)       Position vector from central body to
        !                               field point
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6)       Acceleration vector with
        !                               zero padding (0, 0, 0, ax, ay, az)
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3), &
                                            rbods(3)
        real(wp)                         :: res(6)
        real(wp)                         :: xj, yj, zj, mu_j
        xj = rbods(1); yj = rbods(2); zj = rbods(3); mu_j = mu;
        res = xdot(r)
        contains
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            function xdot(x) result(res)
                implicit none
                real(wp), intent(in) :: x (:)
                real(wp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11

                real(wp), dimension(6) :: res

                x0  =  xj**2
                x1  =  yj**2
                x2  =  zj**2
                x3  =  x0 + x1 + x2
                x4  =  1/sqrt(x3)
                x5  =  1_wp/x3
                x6  =  x0*x5 + x1*x5 + x2*x5
                x7  =  2*x4
                x8  =  x(1)*(x(1) - x7*xj) + x(2)*(x(2) - x7*yj) + x(3)*(x(3) - x7*zj)
                x9  =  x8/x6
                x10  =  x4*x9*(3*x9 + 3 + x8**2/x6**2)/((x9 + 1)**1.5_wp + 1)
                x11  =  mu_j/((-x(1) + xj)**2 + (-x(2) + yj)**2 + (-x(3) + zj)**2)**(3.0_wp/2.0_wp &
                  )

                res(1) =  x(4)
                res(2) =  x(5)
                res(3) =  x(6)
                res(4) =  -x11*(x(1) + x10*xj)
                res(5) =  -x11*(x(2) + x10*yj)
                res(6) =  -x11*(x(3) + x10*zj)
            end function xdot
            ! END AUTOCODE OUTPUT FOR XDOT
    end function
    function jac_nbody(me, mu, r, rbods) result (res)
        ! acc_nbody: method to compute third-body Keplerian acceleration
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! r              real (3)       Position vector from central body to
        !                               field point
        ! rbods          real (3)       Position vector from central body to
        !                               third body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6,6)     third-body Jacobian matrix 
        !                                         df_i
        !                               J(i,j) = ------
        !                                         dx_j
        !                               
        !                               [  dr     dv  ]
        !                               [ ----   ---- ]
        !                               [  dr     dr  ]
        !                               [             ]
        !                               [  dv     dv  ]
        !                               [ ----   ---- ]
        !                               [  dr     dv  ]
        !                               or, if A is the gradient of 
        !                               third-body acceleration with respect
        !                               to position vector r:
        !                               [ 0  I ]
        !                               [ A  0 ]
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3), &
                                            rbods(3)
        real(wp)                         :: res(6,6)
        real(wp)                         :: xj, yj, zj, mu_j
        xj = rbods(1); yj = rbods(2); zj = rbods(3); mu_j = mu
        res = reshape(jac(r),[6,6])
        contains
            function jac(x) result(res)
                implicit none
                real(wp), intent(in) :: x (:)
                real(wp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27, x28

                real(wp), dimension(36) :: res

                x0  =  xj**2
                x1  =  yj**2
                x2  =  zj**2
                x3  =  x0 + x1 + x2
                x4  =  1/sqrt(x3)
                x5  =  1_wp/x3
                x6  =  1_wp/(x0*x5 + x1*x5 + x2*x5)
                x7  =  2*x4
                x8  =  -x7*xj
                x9  =  -x7*yj
                x10  =  -x7*zj
                x11  =  x6*(x(1)*(x(1) + x8) + x(2)*(x(2) + x9) + x(3)*(x(3) + x10)) + 1
                x12  =  x4*(x11**1.5_wp - 1)
                x13  =  x(1) + x12*xj
                x14  =  (-x(1) + xj)**2 + (-x(2) + yj)**2 + (-x(3) + zj)**2
                x15  =  mu_j/x14**(5.0_wp/2.0_wp)
                x16  =  x15*(-3*x(1) + 3*xj)
                x17  =  1.5_wp*sqrt(x11)*x4*x6
                x18  =  x17*(2*x(1) + x8)
                x19  =  mu_j/x14**(3.0_wp/2.0_wp)
                x20  =  x(2) + x12*yj
                x21  =  x18*x19
                x22  =  x(3) + x12*zj
                x23  =  x15*(-3*x(2) + 3*yj)
                x24  =  x17*(2*x(2) + x9)
                x25  =  x19*x24
                x26  =  x15*(-3*x(3) + 3*zj)
                x27  =  x17*(2*x(3) + x10)
                x28  =  x19*x27

                res(1) =  0
                res(2) =  0
                res(3) =  0
                res(4) =  -x13*x16 - x19*(x18*xj + 1)
                res(5) =  -x16*x20 - x21*yj
                res(6) =  -x16*x22 - x21*zj
                res(7) =  0
                res(8) =  0
                res(9) =  0
                res(10) =  -x13*x23 - x25*xj
                res(11) =  -x19*(x24*yj + 1) - x20*x23
                res(12) =  -x22*x23 - x25*zj
                res(13) =  0
                res(14) =  0
                res(15) =  0
                res(16) =  -x13*x26 - x28*xj
                res(17) =  -x20*x26 - x28*yj
                res(18) =  -x19*(x27*zj + 1) - x22*x26
                res(19) =  1
                res(20) =  0
                res(21) =  0
                res(22) =  0
                res(23) =  0
                res(24) =  0
                res(25) =  0
                res(26) =  1
                res(27) =  0
                res(28) =  0
                res(29) =  0
                res(30) =  0
                res(31) =  0
                res(32) =  0
                res(33) =  1
                res(34) =  0
                res(35) =  0
                res(36) =  0
            end function jac
    end function
    function hes_nbody(me, mu, r, rbods) result (res)
        ! jac_kepler: method to compute third body jacobian
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! mu             real           Gravitational parameter of third body
        ! r              real (3)       position vector from central body to
        !                               field point
        ! r              real (3)       position vector from central body to
        !                               third body
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (6,6,6)   third body Hessian tensor H where
        !                                              df_i 
        !                               H(i,j,k) = -----------
        !                                           dx_j dx_k
        class(dynamicsModel), intent(in) :: me
        real(wp),             intent(in) :: mu
        real(wp),             intent(in) :: r(3), &
                                            rbods(3)
        real(wp)                         :: res(6,6,6)
        real(wp)                         :: xj, yj, zj, mu_j
        xj = rbods(1); yj = rbods(2); zj = rbods(3); mu_j = mu
        res = reshape(hes(r),[6,6,6])
        contains
            ! BEGIN AUTOCODE OUTPUT FOR HES
            function hes(x) result(res)
                implicit none
                real(wp), intent(in) :: x (:)
                real(wp)             :: &
                                      & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                      & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                      & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                      & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                      & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                      & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                      & x48, x49, x50, x51, x52, x53, x54, x55, & 
                                      & x56, x57, x58, x59, x60, x61, x62, x63, & 
                                      & x64, x65, x66, x67, x68, x69, x70, x71, & 
                                      & x72, x73, x74, x75

                real(wp), dimension(216) :: res

                x0  =  -x(1) + xj
                x1  =  x0**2
                x2  =  -x(2) + yj
                x3  =  x2**2
                x4  =  -x(3) + zj
                x5  =  x4**2
                x6  =  x1 + x3 + x5
                x7  =  1_wp/x6
                x8  =  1/sqrt(xj**2 + yj**2 + zj**2)
                x9  =  x(2) - 2*x8*yj
                x10  =  -x(2)*x9
                x11  =  x(1)*(-x(1) + 2*x8*xj)
                x12  =  x(3) - 2*x8*zj
                x13  =  -x(3)*x12 + x11 - 1
                x14  =  -x10 - x13
                x15  =  sqrt(x14)
                x16  =  3.0_wp*x15
                x17  =  -x(1) + x8*xj
                x18  =  x8*xj
                x19  =  x17*x18
                x20  =  x16*x19 - 1
                x21  =  x14**1.5_wp - 1
                x22  =  x(1) + x18*x21
                x23  =  5*x7
                x24  =  3*x7
                x25  =  x24*(x1*x23 - 1)
                x26  =  1_wp/x15
                x27  =  3.0_wp*x15 + 3.0_wp*x17**2*x26
                x28  =  mu_j/x6**(3.0_wp/2.0_wp)
                x29  =  x8*yj
                x30  =  x(2) + x21*x29
                x31  =  x8*zj
                x32  =  x(3) + x21*x31
                x33  =  -x(2) + x8*yj
                x34  =  3.0_wp*x26
                x35  =  x33*x34
                x36  =  x6**(-2)
                x37  =  15*x0*x36
                x38  =  x2*x37
                x39  =  x15*x7
                x40  =  9.0_wp*x0*x39
                x41  =  x18*x33
                x42  =  x20*x24
                x43  =  x28*(-x19*x35 + x2*x42 - x22*x38 + x40*x41)
                x44  =  x29*x33
                x45  =  x17*x34
                x46  =  x17*x29
                x47  =  9.0_wp*x39
                x48  =  x2*x47
                x49  =  x16*x44 - 1
                x50  =  x0*x24
                x51  =  x28*(-x30*x38 - x44*x45 + x46*x48 + x49*x50)
                x52  =  x17*x31
                x53  =  x31*x33
                x54  =  x28*(-x32*x38 - x35*x52 + x40*x53 + x48*x52)
                x55  =  -x(3) + x8*zj
                x56  =  x34*x55
                x57  =  x37*x4
                x58  =  x18*x55
                x59  =  x28*(-x19*x56 - x22*x57 + x4*x42 + x40*x58)
                x60  =  x29*x55
                x61  =  x4*x47
                x62  =  x28*(-x30*x57 + x40*x60 - x45*x60 + x46*x61)
                x63  =  x31*x55
                x64  =  x16*x63 - 1
                x65  =  x28*(-x32*x57 - x45*x63 + x50*x64 + x52*x61)
                x66  =  x24*(x23*x3 - 1)
                x67  =  sqrt(x(2)*x9 - x13)
                x68  =  3.0_wp*x33**2/x67 + 3.0_wp*x67
                x69  =  15*x2*x36*x4
                x70  =  x28*(-x22*x69 - x35*x58 + x41*x61 + x48*x58)
                x71  =  x28*(x24*x4*x49 - x30*x69 - x44*x56 + x48*x60)
                x72  =  x28*(x2*x24*x64 - x32*x69 - x35*x63 + x53*x61)
                x73  =  x24*(x23*x5 - 1)
                x74  =  sqrt(x(3)*x12 - x10 - x11 + 1)
                x75  =  3.0_wp*x55**2/x74 + 3.0_wp*x74

                res(1) =  0
                res(2) =  0
                res(3) =  0
                res(4) =  x28*(6*x0*x20*x7 - x18*x27 - x22*x25)
                res(5) =  x28*(18.0*x0*x15*x17*x7*x8*yj - x25*x30 - x27*x29)
                res(6) =  x28*(18.0*x0*x15*x17*x7*x8*zj - x25*x32 - x27*x31)
                res(7) =  0
                res(8) =  0
                res(9) =  0
                res(10) =  x43
                res(11) =  x51
                res(12) =  x54
                res(13) =  0
                res(14) =  0
                res(15) =  0
                res(16) =  x59
                res(17) =  x62
                res(18) =  x65
                res(19) =  0
                res(20) =  0
                res(21) =  0
                res(22) =  0
                res(23) =  0
                res(24) =  0
                res(25) =  0
                res(26) =  0
                res(27) =  0
                res(28) =  0
                res(29) =  0
                res(30) =  0
                res(31) =  0
                res(32) =  0
                res(33) =  0
                res(34) =  0
                res(35) =  0
                res(36) =  0
                res(37) =  0
                res(38) =  0
                res(39) =  0
                res(40) =  x43
                res(41) =  x51
                res(42) =  x54
                res(43) =  0
                res(44) =  0
                res(45) =  0
                res(46) =  x28*(18.0*x15*x2*x33*x7*x8*xj - x18*x68 - x22*x66)
                res(47) =  x28*(6*x2*x49*x7 - x29*x68 - x30*x66)
                res(48) =  x28*(18.0*x15*x2*x33*x7*x8*zj - x31*x68 - x32*x66)
                res(49) =  0
                res(50) =  0
                res(51) =  0
                res(52) =  x70
                res(53) =  x71
                res(54) =  x72
                res(55) =  0
                res(56) =  0
                res(57) =  0
                res(58) =  0
                res(59) =  0
                res(60) =  0
                res(61) =  0
                res(62) =  0
                res(63) =  0
                res(64) =  0
                res(65) =  0
                res(66) =  0
                res(67) =  0
                res(68) =  0
                res(69) =  0
                res(70) =  0
                res(71) =  0
                res(72) =  0
                res(73) =  0
                res(74) =  0
                res(75) =  0
                res(76) =  x59
                res(77) =  x62
                res(78) =  x65
                res(79) =  0
                res(80) =  0
                res(81) =  0
                res(82) =  x70
                res(83) =  x71
                res(84) =  x72
                res(85) =  0
                res(86) =  0
                res(87) =  0
                res(88) =  x28*(18.0*x15*x4*x55*x7*x8*xj - x18*x75 - x22*x73)
                res(89) =  x28*(18.0*x15*x4*x55*x7*x8*yj - x29*x75 - x30*x73)
                res(90) =  x28*(-x31*x75 - x32*x73 + 6*x4*x64*x7)
                res(91) =  0
                res(92) =  0
                res(93) =  0
                res(94) =  0
                res(95) =  0
                res(96) =  0
                res(97) =  0
                res(98) =  0
                res(99) =  0
                res(100) =  0
                res(101) =  0
                res(102) =  0
                res(103) =  0
                res(104) =  0
                res(105) =  0
                res(106) =  0
                res(107) =  0
                res(108) =  0
                res(109) =  0
                res(110) =  0
                res(111) =  0
                res(112) =  0
                res(113) =  0
                res(114) =  0
                res(115) =  0
                res(116) =  0
                res(117) =  0
                res(118) =  0
                res(119) =  0
                res(120) =  0
                res(121) =  0
                res(122) =  0
                res(123) =  0
                res(124) =  0
                res(125) =  0
                res(126) =  0
                res(127) =  0
                res(128) =  0
                res(129) =  0
                res(130) =  0
                res(131) =  0
                res(132) =  0
                res(133) =  0
                res(134) =  0
                res(135) =  0
                res(136) =  0
                res(137) =  0
                res(138) =  0
                res(139) =  0
                res(140) =  0
                res(141) =  0
                res(142) =  0
                res(143) =  0
                res(144) =  0
                res(145) =  0
                res(146) =  0
                res(147) =  0
                res(148) =  0
                res(149) =  0
                res(150) =  0
                res(151) =  0
                res(152) =  0
                res(153) =  0
                res(154) =  0
                res(155) =  0
                res(156) =  0
                res(157) =  0
                res(158) =  0
                res(159) =  0
                res(160) =  0
                res(161) =  0
                res(162) =  0
                res(163) =  0
                res(164) =  0
                res(165) =  0
                res(166) =  0
                res(167) =  0
                res(168) =  0
                res(169) =  0
                res(170) =  0
                res(171) =  0
                res(172) =  0
                res(173) =  0
                res(174) =  0
                res(175) =  0
                res(176) =  0
                res(177) =  0
                res(178) =  0
                res(179) =  0
                res(180) =  0
                res(181) =  0
                res(182) =  0
                res(183) =  0
                res(184) =  0
                res(185) =  0
                res(186) =  0
                res(187) =  0
                res(188) =  0
                res(189) =  0
                res(190) =  0
                res(191) =  0
                res(192) =  0
                res(193) =  0
                res(194) =  0
                res(195) =  0
                res(196) =  0
                res(197) =  0
                res(198) =  0
                res(199) =  0
                res(200) =  0
                res(201) =  0
                res(202) =  0
                res(203) =  0
                res(204) =  0
                res(205) =  0
                res(206) =  0
                res(207) =  0
                res(208) =  0
                res(209) =  0
                res(210) =  0
                res(211) =  0
                res(212) =  0
                res(213) =  0
                res(214) =  0
                res(215) =  0
                res(216) =  0
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
        real(wp),             intent(in)     :: time, &
                                              & r(3)
        real(dp)                             :: pot_d, &
                                              & acc_d(3), &
                                              & jac_d(3,3), &
                                              & hes_d(3,3,3), &
                                              & acc_wp(3), &
                                              & jac_wp(3,3), &
                                              & hes_wp(3,3,3)
        real(wp),             intent(out)    :: acc(6), &
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
        acc(4:) = real(acc_d,wp)
        jac(4:,:3) = real(jac_d,wp)
        hes(4:,:3,:3) = real(hes_d,wp)
    end subroutine
    function eoms(me,t,y) result(res)
        ! eoms: method to compute dynamics function of extended state vector
        !       i.e. [xdot, stmdot, sttdot]
        !TODO: extend with derivatives of time
        ! INPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! t              real           time in sec past J2000 to evaluate
        !                               the dynamics
        ! y              real (258)     Extended dyamics vector:
        !                               258 = 6 + 6 ** 2 + 6 ** 3
        !                        y = [x, reshape(stm,(36)), reshape(stt,(216))]
        ! OUTPUTS:
        ! NAME           TYPE           DESCRIPTION
        ! res            real (258)     third body Hessian tensor H where
        class(dynamicsModel), intent(inout) :: me
        real(wp)            , intent(in)    :: t, y(:)
        real(wp)                            :: state(n), &
                                               stm(n,n), &
                                               stt(n,n,n), &
                                               statedot(n), &
                                               stmdot(n,n), &
                                               sttdot(n,n,n), &
                                               acc(n), &
                                               jac(n,n), &
                                               hes(n,n,n)
        real(wp)                            :: res(size(y))
        state = y(:n)
        stm = reshape(y((n + 1):(n + n**2)),[n,n])
        stt = reshape(y((n + n**2 + 1):),[n,n,n])
        call me%get_derivs(t, acc, jac, hes)
        ! call amatrix(cc,muu,nn,x,amat)
        ! call hess(cc,muu,nn,x,hes)
        stmdot = matmul(jac,stm)
        sttdot = mattens(jac,stt,n) + quad(stm,hes,n)
        res = [state, reshape(stmdot,[n**2]), reshape(sttdot,[n**3])]
    end function eoms
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
        real(wp),             intent(in)    :: time
        real(wp)                            :: lt_dum
        real(dp)                            :: traj_state
        real(wp)                            :: res(6)
        call spkgeo( me%traj_id, real(time,dp), "J2000", me%central_body, &
                   & traj_state, lt_dum)
        res = real(traj_state,wp)
    end function trajstate
end module makemodel
