module kepdyn
    implicit none
    function acc_kepler(me, mu, r) result (res)
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


end module kepdyn
