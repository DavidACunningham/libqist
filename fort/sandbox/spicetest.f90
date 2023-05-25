program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    implicit none
    integer bod_id, N, traj_id
    real(dp) :: mudum(1), time, lt_dum, &
               spkgeo_out(6), spkgps_out(3), &
               stdxdotd(6), xdotd(6), &
               stdjacd(6**2), jacd(6**2), &
               stdhesd(6**3), hesd(6**3), &
               trajstate(6)
    real(wp) :: stdxdotq(6), xdotq(6), &
               stdjacq(6**2), jacq(6**2), &
               stdhesq(6**3), hesq(6**3)
    character(len=12) :: arg

    time = 1000._8
    traj_id = -999
    call FURNSH("/home/david/wrk/nstgro/qist/kernels/mk.tf")
    call get_command_argument(1, arg)
    read(arg,*) bod_id
    call bodvcd(bod_id, 'GM', 1, N, mudum)
    call spkgeo( bod_id, time, "J2000", 399, &
               & spkgeo_out, lt_dum)
    call spkgps( bod_id, time, "J2000", 399, &
                       & spkgps_out, &
                       & lt_dum)
    call spkgeo( traj_id, time, "J2000", 399, &
                       & trajstate, &
                       & lt_dum)
    print *, "BODVCD GM OUTPUT: ", mudum
    print *, "SPKGEO OUTPUT: ", spkgeo_out
    print *, "SPKGPS OUTPUT: ", spkgps_out
    print *, "TRAJ STATE:  ", trajstate


    !~~~ double prec differencing
    print *, "~~~~~~~ DOUBLE PRECISION PERFORMANCE ~~~~~~~"
    stdxdotd = stdxdot_d(mudum(1), trajstate, spkgps_out)
    xdotd = xdot_d(mudum(1), trajstate, spkgps_out) 
    print *, "ACCEL STD"
    print *, stdxdotd
    print *, "ACCEL ENCKE BATTIN"
    print *, xdotd
    print *, "ACCEL MAX NORM ERR"
    print *, maxval(abs((stdxdotd-xdotd)/xdotd))
    stdjacd = stdjac_d(mudum(1), trajstate, spkgps_out)
    jacd = jac_d(mudum(1), trajstate, spkgps_out) 
    print *, "JAC MAX NORM ERR"
    print *, maxval(abs((stdjacd-jacd)/jacd))
    stdhesd = stdhes_d(mudum(1), trajstate, spkgps_out)
    hesd = hes_d(mudum(1), trajstate, spkgps_out) 
    print *, "HES MAX NORM ERR"
    print *, maxval(abs((stdhesd-hesd)/hesd))
                        
    print *, "~~~~~~~ QUAD PRECISION PERFORMANCE ~~~~~~~"
    stdxdotq = stdxdot_q(real(mudum(1),wp), &
                         real(trajstate,wp), &
                         real(spkgps_out,wp))
    xdotq = xdot_q(real(mudum(1),wp), &
                   real(trajstate,wp), &
                   real(spkgps_out,wp))
    print *, "ACCEL STD"
    print *, stdxdotq
    print *, "ACCEL ENCKE BATTIN"
    print *, xdotq
    print *, "ACCEL MAX NORM ERR"
    print *, maxval(abs((stdxdotq-xdotq)/xdotq))
    stdjacq = stdjac_q(real(mudum(1),wp), &
                         real(trajstate,wp), &
                         real(spkgps_out,wp))
    jacq = jac_q(real(mudum(1),wp), &
                   real(trajstate,wp), &
                   real(spkgps_out,wp))
    print *, "JAC MAX NORM ERR"
    print *, maxval(abs((stdjacq-jacq)/jacq))
    stdhesq = stdhes_q(real(mudum(1),wp), &
                         real(trajstate,wp), &
                         real(spkgps_out,wp))
    hesq = hes_q(real(mudum(1),wp), &
                   real(trajstate,wp), &
                   real(spkgps_out,wp))
    print *, "HES MAX NORM ERR"
    print *, maxval(abs((stdhesq-hesq)/hesq))
    contains
        function xdot_d(mu_j,x,rj) result(res)
            implicit none
            real(dp), intent(in) :: x (:), rj(:), mu_j
            real(dp)             :: &
                                  & x0, x1, x2, x3, x4, xj, yj,zj

            real(dp), dimension(6) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);


            x0  =  xj**2 + yj**2 + zj**2
            x1  =  x(1)*(x(1) - 2*xj) + x(2)*(x(2) - 2*yj) + x(3)*(x(3) - 2*zj)
            x2  =  x1/x0
            x3  =  x2*(3*x2 + 3 + x1**2/x0**2)/((x2 + 1)**1.5_dp + 1)
            x4  =  mu_j/((-x(1) + xj)**2 + (-x(2) + yj)**2 + (-x(3) + zj)**2)**(3.0_dp/2.0_dp &
              )

            res(1) =  x(4)
            res(2) =  x(5)
            res(3) =  x(6)
            res(4) =  -x4*(x(1) + x3*xj)
            res(5) =  -x4*(x(2) + x3*yj)
            res(6) =  -x4*(x(3) + x3*zj)
        end function
        ! END AUTOCODE OUTPUT FOR XDOT
        ! BEGIN AUTOCODE OUTPUT FOR JAC
        function jac_d(mu_j,x,rj) result(res)
            implicit none
            real(dp), intent(in) :: x (:), rj(:), mu_j
            real(dp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                  & x16, x17, x18, x19, x20, x21, xj, yj, zj

            real(dp), dimension(36) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  1_dp/(xj**2 + yj**2 + zj**2)
            x1  =  -2*xj
            x2  =  -2*yj
            x3  =  -2*zj
            x4  =  x0*(x(1)*(x(1) + x1) + x(2)*(x(2) + x2) + x(3)*(x(3) + x3)) + 1
            x5  =  x4**1.5_dp - 1
            x6  =  x(1) + x5*xj
            x7  =  (-x(1) + xj)**2 + (-x(2) + yj)**2 + (-x(3) + zj)**2
            x8  =  mu_j/x7**(5.0_dp/2.0_dp)
            x9  =  x8*(-3*x(1) + 3*xj)
            x10  =  1.5_dp*x0*sqrt(x4)
            x11  =  x10*(2*x(1) + x1)
            x12  =  mu_j/x7**(3.0_dp/2.0_dp)
            x13  =  x(2) + x5*yj
            x14  =  x11*x12
            x15  =  x(3) + x5*zj
            x16  =  x8*(-3*x(2) + 3*yj)
            x17  =  x10*(2*x(2) + x2)
            x18  =  x12*x17
            x19  =  x8*(-3*x(3) + 3*zj)
            x20  =  x10*(2*x(3) + x3)
            x21  =  x12*x20

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  -x12*(x11*xj + 1) - x6*x9
            res(5) =  -x13*x9 - x14*yj
            res(6) =  -x14*zj - x15*x9
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  -x16*x6 - x18*xj
            res(11) =  -x12*(x17*yj + 1) - x13*x16
            res(12) =  -x15*x16 - x18*zj
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  -x19*x6 - x21*xj
            res(17) =  -x13*x19 - x21*yj
            res(18) =  -x12*(x20*zj + 1) - x15*x19
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
        end function
        ! END AUTOCODE OUTPUT FOR JAC
        ! BEGIN AUTOCODE OUTPUT FOR HES
        function hes_d(mu_j,x,rj) result(res)
            implicit none
            real(dp), intent(in) :: x (:), rj(:), mu_j
            real(dp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                  & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                  & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                  & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                  & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                  & x48, x49, x50, x51, x52, x53, xj, yj, zj

            real(dp), dimension(216) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  -x(1) + xj
            x1  =  x0**2
            x2  =  -x(2) + yj
            x3  =  x2**2
            x4  =  -x(3) + zj
            x5  =  x4**2
            x6  =  x1 + x3 + x5
            x7  =  1_dp/x6
            x8  =  xj**2 + yj**2 + zj**2
            x9  =  1_dp/x8
            x10  =  -x9*(x(1)*(-x(1) + 2*xj) + x(2)*(-x(2) + 2*yj) + x(3)*(-x(3) + 2*zj)) + &
              1
            x11  =  sqrt(x10)
            x12  =  3.0_dp*x9
            x13  =  x11*x12
            x14  =  x0*x13*xj - 1
            x15  =  x10**1.5_dp - 1
            x16  =  x(1) + x15*xj
            x17  =  5*x7
            x18  =  3*x7
            x19  =  x18*(x1*x17 - 1)
            x20  =  1_dp/x11
            x21  =  x20*x9
            x22  =  x12*(x1*x21 + x11)
            x23  =  mu_j/x6**(3.0_dp/2.0_dp)
            x24  =  x(2) + x15*yj
            x25  =  x(3) + x15*zj
            x26  =  3.0_dp*x20/x8**2
            x27  =  x26*xj
            x28  =  15/x6**2
            x29  =  x16*x28
            x30  =  x23*(9.0_dp*x0*x11*x7*x9*xj - x0*x27 - x0*x29 + 3*x14*x7)
            x31  =  x2*x30
            x32  =  x26*yj
            x33  =  x24*x28
            x34  =  x13*x2*yj - 1
            x35  =  9.0_dp*x11*x2*x7*x9*yj - x2*x32 - x2*x33 + 3*x34*x7
            x36  =  x0*x23
            x37  =  x35*x36
            x38  =  x26*zj
            x39  =  x25*x28
            x40  =  x2*x36*(18.0_dp*x11*x7*x9*zj - x38 - x39)
            x41  =  x30*x4
            x42  =  x36*x4*(18.0_dp*x11*x7*x9*yj - x32 - x33)
            x43  =  x13*x4*zj - 1
            x44  =  9.0_dp*x11*x4*x7*x9*zj - x38*x4 - x39*x4 + 3*x43*x7
            x45  =  x36*x44
            x46  =  x18*(x17*x3 - 1)
            x47  =  x12*(x11 + x21*x3)
            x48  =  x23*x4
            x49  =  x2*x48*(18.0_dp*x11*x7*x9*xj - x27 - x29)
            x50  =  x35*x48
            x51  =  x2*x23*x44
            x52  =  x18*(x17*x5 - 1)
            x53  =  x12*(x11 + x21*x5)

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  x23*(6*x0*x14*x7 - x16*x19 - x22*xj)
            res(5) =  x23*(18.0*x1*x11*x7*x9*yj - x19*x24 - x22*yj)
            res(6) =  x23*(18.0*x1*x11*x7*x9*zj - x19*x25 - x22*zj)
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  x31
            res(11) =  x37
            res(12) =  x40
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  x41
            res(17) =  x42
            res(18) =  x45
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
            res(40) =  x31
            res(41) =  x37
            res(42) =  x40
            res(43) =  0
            res(44) =  0
            res(45) =  0
            res(46) =  x23*(18.0*x11*x3*x7*x9*xj - x16*x46 - x47*xj)
            res(47) =  x23*(6*x2*x34*x7 - x24*x46 - x47*yj)
            res(48) =  x23*(18.0*x11*x3*x7*x9*zj - x25*x46 - x47*zj)
            res(49) =  0
            res(50) =  0
            res(51) =  0
            res(52) =  x49
            res(53) =  x50
            res(54) =  x51
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
            res(76) =  x41
            res(77) =  x42
            res(78) =  x45
            res(79) =  0
            res(80) =  0
            res(81) =  0
            res(82) =  x49
            res(83) =  x50
            res(84) =  x51
            res(85) =  0
            res(86) =  0
            res(87) =  0
            res(88) =  x23*(18.0*x11*x5*x7*x9*xj - x16*x52 - x53*xj)
            res(89) =  x23*(18.0*x11*x5*x7*x9*yj - x24*x52 - x53*yj)
            res(90) =  x23*(-x25*x52 + 6*x4*x43*x7 - x53*zj)
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
        end function
        ! BEGIN AUTOCODE OUTPUT FOR XDOT
        function xdot_q(mu_j,x,rj) result(res)
            implicit none
            real(wp), intent(in) :: x (:), rj(:), mu_j
            real(wp)             :: &
                                  & x0, x1, x2, x3, x4, xj, yj,zj

            real(wp), dimension(6) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);


            x0  =  xj**2 + yj**2 + zj**2
            x1  =  x(1)*(x(1) - 2*xj) + x(2)*(x(2) - 2*yj) + x(3)*(x(3) - 2*zj)
            x2  =  x1/x0
            x3  =  x2*(3*x2 + 3 + x1**2/x0**2)/((x2 + 1)**1.5_wp + 1)
            x4  =  mu_j/((-x(1) + xj)**2 + (-x(2) + yj)**2 + (-x(3) + zj)**2)**(3.0_wp/2.0_wp &
              )

            res(1) =  x(4)
            res(2) =  x(5)
            res(3) =  x(6)
            res(4) =  -x4*(x(1) + x3*xj)
            res(5) =  -x4*(x(2) + x3*yj)
            res(6) =  -x4*(x(3) + x3*zj)
        end function
        ! END AUTOCODE OUTPUT FOR XDOT
        ! BEGIN AUTOCODE OUTPUT FOR JAC
        function jac_q(mu_j,x,rj) result(res)
            implicit none
            real(wp), intent(in) :: x (:), rj(:), mu_j
            real(wp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                  & x16, x17, x18, x19, x20, x21, xj, yj, zj

            real(wp), dimension(36) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  1_wp/(xj**2 + yj**2 + zj**2)
            x1  =  -2*xj
            x2  =  -2*yj
            x3  =  -2*zj
            x4  =  x0*(x(1)*(x(1) + x1) + x(2)*(x(2) + x2) + x(3)*(x(3) + x3)) + 1
            x5  =  x4**1.5_wp - 1
            x6  =  x(1) + x5*xj
            x7  =  (-x(1) + xj)**2 + (-x(2) + yj)**2 + (-x(3) + zj)**2
            x8  =  mu_j/x7**(5.0_wp/2.0_wp)
            x9  =  x8*(-3*x(1) + 3*xj)
            x10  =  1.5_wp*x0*sqrt(x4)
            x11  =  x10*(2*x(1) + x1)
            x12  =  mu_j/x7**(3.0_wp/2.0_wp)
            x13  =  x(2) + x5*yj
            x14  =  x11*x12
            x15  =  x(3) + x5*zj
            x16  =  x8*(-3*x(2) + 3*yj)
            x17  =  x10*(2*x(2) + x2)
            x18  =  x12*x17
            x19  =  x8*(-3*x(3) + 3*zj)
            x20  =  x10*(2*x(3) + x3)
            x21  =  x12*x20

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  -x12*(x11*xj + 1) - x6*x9
            res(5) =  -x13*x9 - x14*yj
            res(6) =  -x14*zj - x15*x9
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  -x16*x6 - x18*xj
            res(11) =  -x12*(x17*yj + 1) - x13*x16
            res(12) =  -x15*x16 - x18*zj
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  -x19*x6 - x21*xj
            res(17) =  -x13*x19 - x21*yj
            res(18) =  -x12*(x20*zj + 1) - x15*x19
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
        end function
        ! END AUTOCODE OUTPUT FOR JAC
        ! BEGIN AUTOCODE OUTPUT FOR HES
        function hes_q(mu_j,x,rj) result(res)
            implicit none
            real(wp), intent(in) :: x (:), rj(:), mu_j
            real(wp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                  & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                  & x24, x25, x26, x27, x28, x29, x30, x31, & 
                                  & x32, x33, x34, x35, x36, x37, x38, x39, & 
                                  & x40, x41, x42, x43, x44, x45, x46, x47, & 
                                  & x48, x49, x50, x51, x52, x53, xj, yj, zj
            real(wp), dimension(216) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  -x(1) + xj
            x1  =  x0**2
            x2  =  -x(2) + yj
            x3  =  x2**2
            x4  =  -x(3) + zj
            x5  =  x4**2
            x6  =  x1 + x3 + x5
            x7  =  1_wp/x6
            x8  =  xj**2 + yj**2 + zj**2
            x9  =  1_wp/x8
            x10  =  -x9*(x(1)*(-x(1) + 2*xj) + x(2)*(-x(2) + 2*yj) + x(3)*(-x(3) + 2*zj)) + &
              1
            x11  =  sqrt(x10)
            x12  =  3.0_wp*x9
            x13  =  x11*x12
            x14  =  x0*x13*xj - 1
            x15  =  x10**1.5_wp - 1
            x16  =  x(1) + x15*xj
            x17  =  5*x7
            x18  =  3*x7
            x19  =  x18*(x1*x17 - 1)
            x20  =  1_wp/x11
            x21  =  x20*x9
            x22  =  x12*(x1*x21 + x11)
            x23  =  mu_j/x6**(3.0_wp/2.0_wp)
            x24  =  x(2) + x15*yj
            x25  =  x(3) + x15*zj
            x26  =  3.0_wp*x20/x8**2
            x27  =  x26*xj
            x28  =  15/x6**2
            x29  =  x16*x28
            x30  =  x23*(9.0_wp*x0*x11*x7*x9*xj - x0*x27 - x0*x29 + 3*x14*x7)
            x31  =  x2*x30
            x32  =  x26*yj
            x33  =  x24*x28
            x34  =  x13*x2*yj - 1
            x35  =  9.0_wp*x11*x2*x7*x9*yj - x2*x32 - x2*x33 + 3*x34*x7
            x36  =  x0*x23
            x37  =  x35*x36
            x38  =  x26*zj
            x39  =  x25*x28
            x40  =  x2*x36*(18.0_wp*x11*x7*x9*zj - x38 - x39)
            x41  =  x30*x4
            x42  =  x36*x4*(18.0_wp*x11*x7*x9*yj - x32 - x33)
            x43  =  x13*x4*zj - 1
            x44  =  9.0_wp*x11*x4*x7*x9*zj - x38*x4 - x39*x4 + 3*x43*x7
            x45  =  x36*x44
            x46  =  x18*(x17*x3 - 1)
            x47  =  x12*(x11 + x21*x3)
            x48  =  x23*x4
            x49  =  x2*x48*(18.0_wp*x11*x7*x9*xj - x27 - x29)
            x50  =  x35*x48
            x51  =  x2*x23*x44
            x52  =  x18*(x17*x5 - 1)
            x53  =  x12*(x11 + x21*x5)

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  x23*(6*x0*x14*x7 - x16*x19 - x22*xj)
            res(5) =  x23*(18.0*x1*x11*x7*x9*yj - x19*x24 - x22*yj)
            res(6) =  x23*(18.0*x1*x11*x7*x9*zj - x19*x25 - x22*zj)
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  x31
            res(11) =  x37
            res(12) =  x40
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  x41
            res(17) =  x42
            res(18) =  x45
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
            res(40) =  x31
            res(41) =  x37
            res(42) =  x40
            res(43) =  0
            res(44) =  0
            res(45) =  0
            res(46) =  x23*(18.0*x11*x3*x7*x9*xj - x16*x46 - x47*xj)
            res(47) =  x23*(6*x2*x34*x7 - x24*x46 - x47*yj)
            res(48) =  x23*(18.0*x11*x3*x7*x9*zj - x25*x46 - x47*zj)
            res(49) =  0
            res(50) =  0
            res(51) =  0
            res(52) =  x49
            res(53) =  x50
            res(54) =  x51
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
            res(76) =  x41
            res(77) =  x42
            res(78) =  x45
            res(79) =  0
            res(80) =  0
            res(81) =  0
            res(82) =  x49
            res(83) =  x50
            res(84) =  x51
            res(85) =  0
            res(86) =  0
            res(87) =  0
            res(88) =  x23*(18.0*x11*x5*x7*x9*xj - x16*x52 - x53*xj)
            res(89) =  x23*(18.0*x11*x5*x7*x9*yj - x24*x52 - x53*yj)
            res(90) =  x23*(-x25*x52 + 6*x4*x43*x7 - x53*zj)
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
        end function
        ! END AUTOCODE OUTPUT FOR HES
        function stdxdot_d(mu_j,x,rj) result(res)
            implicit none
            real(dp), intent(in) :: x (:), rj(:), mu_j
            real(dp)             :: &
                                  & x0, x1, x2, x3, x4, xj, yj, zj

            real(dp), dimension(6) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  (xj**2 + yj**2 + zj**2)**(-3.0_dp/2.0_dp)
            x1  =  -x(1) + xj
            x2  =  -x(2) + yj
            x3  =  -x(3) + zj
            x4  =  (x1**2 + x2**2 + x3**2)**(-3.0_dp/2.0_dp)

            res(1) =  x(4)
            res(2) =  x(5)
            res(3) =  x(6)
            res(4) =  mu_j*(-x0*xj + x1*x4)
            res(5) =  mu_j*(-x0*yj + x2*x4)
            res(6) =  mu_j*(-x0*zj + x3*x4)
        end function 
        ! END AUTOCODE OUTPUT FOR XDOT

        ! BEGIN AUTOCODE OUTPUT FOR JAC
        function stdjac_d(mu_j,x,rj) result(res)
            implicit none
            real(dp), intent(in) :: x (:) , rj(:), mu_j
            real(dp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, xj, yj, zj

            real(dp), dimension(36) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  -x(1) + xj
            x1  =  -x(2) + yj
            x2  =  -x(3) + zj
            x3  =  x0**2 + x1**2 + x2**2
            x4  =  -1/x3**(3.0_dp/2.0_dp)
            x5  =  x3**(-5.0_dp/2.0_dp)
            x6  =  x5*(-3*x(1) + 3*xj)
            x7  =  mu_j*x6
            x8  =  x5*(-3*x(2) + 3*yj)
            x9  =  mu_j*x8
            x10  =  x5*(-3*x(3) + 3*zj)
            x11  =  mu_j*x10

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  mu_j*(x0*x6 + x4)
            res(5) =  x1*x7
            res(6) =  x2*x7
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  x0*x9
            res(11) =  mu_j*(x1*x8 + x4)
            res(12) =  x2*x9
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  x0*x11
            res(17) =  x1*x11
            res(18) =  mu_j*(x10*x2 + x4)
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
        end function 
        ! END AUTOCODE OUTPUT FOR JAC

        ! BEGIN AUTOCODE OUTPUT FOR HES
        function stdhes_d(mu_j,x,rj) result(res)
            implicit none
            real(dp), intent(in) :: x (:), rj(:), mu_j
            real(dp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                  & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                  & x24, x25, x26, xj, yj, zj

            real(dp), dimension(216) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  -x(1) + xj
            x1  =  x0**2
            x2  =  -x(2) + yj
            x3  =  x2**2
            x4  =  -x(3) + zj
            x5  =  x4**2
            x6  =  x1 + x3 + x5
            x7  =  5/x6
            x8  =  x1*x7
            x9  =  mu_j*x0
            x10  =  3/x6**(5.0_dp/2.0_dp)
            x11  =  x10*x9
            x12  =  mu_j*x10
            x13  =  x12*(x8 - 1)
            x14  =  x13*x2
            x15  =  x13*x4
            x16  =  x3*x7
            x17  =  x16 - 1
            x18  =  x11*x17
            x19  =  15*x2*x4*x9/x6**(7.0_dp/2.0_dp)
            x20  =  x5*x7
            x21  =  x20 - 1
            x22  =  x11*x21
            x23  =  x12*x2
            x24  =  x12*x4
            x25  =  x17*x24
            x26  =  x21*x23

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  x11*(x8 - 3)
            res(5) =  x14
            res(6) =  x15
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  x14
            res(11) =  x18
            res(12) =  x19
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  x15
            res(17) =  x19
            res(18) =  x22
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
            res(40) =  x14
            res(41) =  x18
            res(42) =  x19
            res(43) =  0
            res(44) =  0
            res(45) =  0
            res(46) =  x18
            res(47) =  x23*(x16 - 3)
            res(48) =  x25
            res(49) =  0
            res(50) =  0
            res(51) =  0
            res(52) =  x19
            res(53) =  x25
            res(54) =  x26
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
            res(76) =  x15
            res(77) =  x19
            res(78) =  x22
            res(79) =  0
            res(80) =  0
            res(81) =  0
            res(82) =  x19
            res(83) =  x25
            res(84) =  x26
            res(85) =  0
            res(86) =  0
            res(87) =  0
            res(88) =  x22
            res(89) =  x26
            res(90) =  x24*(x20 - 3)
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
        end function 
        ! END AUTOCODE OUTPUT FOR HES
        ! BEGIN AUTOCODE OUTPUT FOR XDOT
        function stdxdot_q(mu_j,x,rj) result(res)
            implicit none
            real(wp), intent(in) :: x (:), rj(:), mu_j
            real(wp)             :: &
                                  & x0, x1, x2, x3, x4, xj, yj, zj

            real(wp), dimension(6) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  (xj**2 + yj**2 + zj**2)**(-3.0_wp/2.0_wp)
            x1  =  -x(1) + xj
            x2  =  -x(2) + yj
            x3  =  -x(3) + zj
            x4  =  (x1**2 + x2**2 + x3**2)**(-3.0_wp/2.0_wp)

            res(1) =  x(4)
            res(2) =  x(5)
            res(3) =  x(6)
            res(4) =  mu_j*(-x0*xj + x1*x4)
            res(5) =  mu_j*(-x0*yj + x2*x4)
            res(6) =  mu_j*(-x0*zj + x3*x4)
        end function 
        ! END AUTOCODE OUTPUT FOR XDOT

        ! BEGIN AUTOCODE OUTPUT FOR JAC
        function stdjac_q(mu_j,x,rj) result(res)
            implicit none
            real(wp), intent(in) :: x (:) , rj(:), mu_j
            real(wp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, xj, yj, zj

            real(wp), dimension(36) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  -x(1) + xj
            x1  =  -x(2) + yj
            x2  =  -x(3) + zj
            x3  =  x0**2 + x1**2 + x2**2
            x4  =  -1/x3**(3.0_wp/2.0_wp)
            x5  =  x3**(-5.0_wp/2.0_wp)
            x6  =  x5*(-3*x(1) + 3*xj)
            x7  =  mu_j*x6
            x8  =  x5*(-3*x(2) + 3*yj)
            x9  =  mu_j*x8
            x10  =  x5*(-3*x(3) + 3*zj)
            x11  =  mu_j*x10

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  mu_j*(x0*x6 + x4)
            res(5) =  x1*x7
            res(6) =  x2*x7
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  x0*x9
            res(11) =  mu_j*(x1*x8 + x4)
            res(12) =  x2*x9
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  x0*x11
            res(17) =  x1*x11
            res(18) =  mu_j*(x10*x2 + x4)
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
        end function 
        ! END AUTOCODE OUTPUT FOR JAC

        ! BEGIN AUTOCODE OUTPUT FOR HES
        function stdhes_q(mu_j,x,rj) result(res)
            implicit none
            real(wp), intent(in) :: x (:), rj(:), mu_j
            real(wp)             :: &
                                  & x0, x1, x2, x3, x4, x5, x6, x7, & 
                                  & x8, x9, x10, x11, x12, x13, x14, x15, & 
                                  & x16, x17, x18, x19, x20, x21, x22, x23, & 
                                  & x24, x25, x26, xj, yj, zj

            real(wp), dimension(216) :: res
            xj = rj(1); yj = rj(2); zj = rj(3);

            x0  =  -x(1) + xj
            x1  =  x0**2
            x2  =  -x(2) + yj
            x3  =  x2**2
            x4  =  -x(3) + zj
            x5  =  x4**2
            x6  =  x1 + x3 + x5
            x7  =  5/x6
            x8  =  x1*x7
            x9  =  mu_j*x0
            x10  =  3/x6**(5.0_wp/2.0_wp)
            x11  =  x10*x9
            x12  =  mu_j*x10
            x13  =  x12*(x8 - 1)
            x14  =  x13*x2
            x15  =  x13*x4
            x16  =  x3*x7
            x17  =  x16 - 1
            x18  =  x11*x17
            x19  =  15*x2*x4*x9/x6**(7.0_wp/2.0_wp)
            x20  =  x5*x7
            x21  =  x20 - 1
            x22  =  x11*x21
            x23  =  x12*x2
            x24  =  x12*x4
            x25  =  x17*x24
            x26  =  x21*x23

            res(1) =  0
            res(2) =  0
            res(3) =  0
            res(4) =  x11*(x8 - 3)
            res(5) =  x14
            res(6) =  x15
            res(7) =  0
            res(8) =  0
            res(9) =  0
            res(10) =  x14
            res(11) =  x18
            res(12) =  x19
            res(13) =  0
            res(14) =  0
            res(15) =  0
            res(16) =  x15
            res(17) =  x19
            res(18) =  x22
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
            res(40) =  x14
            res(41) =  x18
            res(42) =  x19
            res(43) =  0
            res(44) =  0
            res(45) =  0
            res(46) =  x18
            res(47) =  x23*(x16 - 3)
            res(48) =  x25
            res(49) =  0
            res(50) =  0
            res(51) =  0
            res(52) =  x19
            res(53) =  x25
            res(54) =  x26
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
            res(76) =  x15
            res(77) =  x19
            res(78) =  x22
            res(79) =  0
            res(80) =  0
            res(81) =  0
            res(82) =  x19
            res(83) =  x25
            res(84) =  x26
            res(85) =  0
            res(86) =  0
            res(87) =  0
            res(88) =  x22
            res(89) =  x26
            res(90) =  x24*(x20 - 3)
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
        end function 
end program
