module qist
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use globals
    use tensorops, only:  sttchain, sttchain_invert, vectensquad, &
                          stminvert, sttinvert
    use denseLight, only: lightSol
    implicit none
    private
    !! THIS IS DUMB AND FOR DEBUGGING AND ALSO A HACK
    !! FOR INTEGRATING INVERSES AT QUAD PRECISION
    ! integer, parameter, public :: dp=wp

    !! THIS IS NOT DUMB AND SHOULD GIVE THE DEFAULT BEHAVIOR
    ! integer, parameter :: dp=selected_real_kind(15)

    type, public :: Itraj
        !! All type variables accessed through getter/setter methods
        real(dp)               :: t0, tf 
        logical                :: initq
        type(lightSol)         :: reftraj
        character(len=1000)    :: densefilename

    contains
        !! Convenience procedures/check state
        procedure :: init 
        procedure call
        procedure :: state 
        procedure :: stm 
        procedure :: stm_i 
        procedure :: stt
        procedure :: stt_i 
        generic, public :: prop => prop_once, prop_many
        procedure, private :: prop_once
        procedure, private :: prop_many

        procedure stts_ab
        procedure zmap
    end type Itraj

    contains
    ! Initializer
    subroutine init(self, t0, tf, filepath, trajfile)
        character(len=*), intent(in) :: filepath, trajfile
        !! The filename where the chebyshev coefficients are located
        real(dp),            intent(in) :: t0, tf
        !! The initial and final simulation independent variable (t)
        class(Itraj),        intent(out):: self
        !! The function _returns_ an instance of the type. 
        ! print *, "Initializing reference model"
        self%initq=.true.
        ! The initializer is running
        self%t0               = t0
        self%tf               = tf
        open(unit=49, file=filepath//trajfile, status="old", access="stream")
        call self%reftraj%read(49)
        close(49)
        ! Make sure the chebyshev order in the file is good to go
        ! Set the coefficient filename
        ! Allocate the chebyshev coefficient array in memory
    end subroutine init
    !! Calling functions 
    function call(self,t,lind,uind) result(res)
        !! Return a (uind-lind)-dimensional at t
        !! If no lind (uind) is passed, beginning (end) of the vector is used
        class(Itraj),      intent(inout) :: self
        real(dp),          intent(in) :: t
        integer, optional, intent(in) :: uind, lind
        !! The value of the independent variable to generate
        real(dp), allocatable         :: res(:)
        integer l, u
        l = 1
        u = plen
        if (present(lind)) l=lind
        if (present(uind)) u=uind
        allocate(res(u-l+1))
        res = self%reftraj%call(t,l,u)
    end function
    function state(self,t) result(res)
        !! Return a regularized state at time t
        class(ITraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp), dimension(n)   :: res
        !! The returned state
        res = self%call(t,uind=n)
    end function state
    function stm(self,t) result(res)
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp)                 :: packlin(plen)
        real(dp), dimension(n,n) :: res
        !! The returned stm
        packlin = self%call(t)
        res=0._dp
        res(1,1) = packlin(9)
        res(1,2) = packlin(15)
        res(1,3) = packlin(21)
        res(1,4) = packlin(27)
        res(1,5) = packlin(33)
        res(1,6) = packlin(39)
        res(1,7) = packlin(45)
        res(1,8) = packlin(51)
        res(2,1) = packlin(10)
        res(2,2) = packlin(16)
        res(2,3) = packlin(22)
        res(2,4) = packlin(28)
        res(2,5) = packlin(34)
        res(2,6) = packlin(40)
        res(2,7) = packlin(46)
        res(2,8) = packlin(52)
        res(3,1) = packlin(11)
        res(3,2) = packlin(17)
        res(3,3) = packlin(23)
        res(3,4) = packlin(29)
        res(3,5) = packlin(35)
        res(3,6) = packlin(41)
        res(3,7) = packlin(47)
        res(3,8) = packlin(53)
        res(4,1) = packlin(12)
        res(4,2) = packlin(18)
        res(4,3) = packlin(24)
        res(4,4) = packlin(30)
        res(4,5) = packlin(36)
        res(4,6) = packlin(42)
        res(4,7) = packlin(48)
        res(4,8) = packlin(54)
        res(5,1) = packlin(13)
        res(5,2) = packlin(19)
        res(5,3) = packlin(25)
        res(5,4) = packlin(31)
        res(5,5) = packlin(37)
        res(5,6) = packlin(43)
        res(5,7) = packlin(49)
        res(5,8) = packlin(55)
        res(6,1) = packlin(14)
        res(6,2) = packlin(20)
        res(6,3) = packlin(26)
        res(6,4) = packlin(32)
        res(6,5) = packlin(38)
        res(6,6) = packlin(44)
        res(6,7) = packlin(50)
        res(6,8) = packlin(56)
        res(7,7)=1._dp
        res(7,8)=1._dp
        res(8,8)=1._dp
    end function stm
    function stm_i(self,t) result(res)
        !! Return a regularized stm at time tau
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        real(dp), dimension(n,n) :: res
        !! The returned stm
        res = stminvert(self%stm(t), n)
    end function stm_i
    function stt(self,t) result(res)
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp)                   :: packlin(plen)
        real(dp), dimension(n,n,n) :: res
        !! The returned stm
        packlin = self%call(t)
        res(1,1,1) = packlin(57)
        res(1,1,2) = packlin(63)
        res(1,1,3) = packlin(75)
        res(1,1,4) = packlin(93)
        res(1,1,5) = packlin(117)
        res(1,1,6) = packlin(147)
        res(1,1,7) = packlin(183)
        res(1,1,8) = packlin(225)
        res(1,2,1) = packlin(63)
        res(1,2,2) = packlin(69)
        res(1,2,3) = packlin(81)
        res(1,2,4) = packlin(99)
        res(1,2,5) = packlin(123)
        res(1,2,6) = packlin(153)
        res(1,2,7) = packlin(189)
        res(1,2,8) = packlin(231)
        res(1,3,1) = packlin(75)
        res(1,3,2) = packlin(81)
        res(1,3,3) = packlin(87)
        res(1,3,4) = packlin(105)
        res(1,3,5) = packlin(129)
        res(1,3,6) = packlin(159)
        res(1,3,7) = packlin(195)
        res(1,3,8) = packlin(237)
        res(1,4,1) = packlin(93)
        res(1,4,2) = packlin(99)
        res(1,4,3) = packlin(105)
        res(1,4,4) = packlin(111)
        res(1,4,5) = packlin(135)
        res(1,4,6) = packlin(165)
        res(1,4,7) = packlin(201)
        res(1,4,8) = packlin(243)
        res(1,5,1) = packlin(117)
        res(1,5,2) = packlin(123)
        res(1,5,3) = packlin(129)
        res(1,5,4) = packlin(135)
        res(1,5,5) = packlin(141)
        res(1,5,6) = packlin(171)
        res(1,5,7) = packlin(207)
        res(1,5,8) = packlin(249)
        res(1,6,1) = packlin(147)
        res(1,6,2) = packlin(153)
        res(1,6,3) = packlin(159)
        res(1,6,4) = packlin(165)
        res(1,6,5) = packlin(171)
        res(1,6,6) = packlin(177)
        res(1,6,7) = packlin(213)
        res(1,6,8) = packlin(255)
        res(1,7,1) = packlin(183)
        res(1,7,2) = packlin(189)
        res(1,7,3) = packlin(195)
        res(1,7,4) = packlin(201)
        res(1,7,5) = packlin(207)
        res(1,7,6) = packlin(213)
        res(1,7,7) = packlin(219)
        res(1,7,8) = packlin(261)
        res(1,8,1) = packlin(225)
        res(1,8,2) = packlin(231)
        res(1,8,3) = packlin(237)
        res(1,8,4) = packlin(243)
        res(1,8,5) = packlin(249)
        res(1,8,6) = packlin(255)
        res(1,8,7) = packlin(261)
        res(1,8,8) = packlin(267)
        res(2,1,1) = packlin(58)
        res(2,1,2) = packlin(64)
        res(2,1,3) = packlin(76)
        res(2,1,4) = packlin(94)
        res(2,1,5) = packlin(118)
        res(2,1,6) = packlin(148)
        res(2,1,7) = packlin(184)
        res(2,1,8) = packlin(226)
        res(2,2,1) = packlin(64)
        res(2,2,2) = packlin(70)
        res(2,2,3) = packlin(82)
        res(2,2,4) = packlin(100)
        res(2,2,5) = packlin(124)
        res(2,2,6) = packlin(154)
        res(2,2,7) = packlin(190)
        res(2,2,8) = packlin(232)
        res(2,3,1) = packlin(76)
        res(2,3,2) = packlin(82)
        res(2,3,3) = packlin(88)
        res(2,3,4) = packlin(106)
        res(2,3,5) = packlin(130)
        res(2,3,6) = packlin(160)
        res(2,3,7) = packlin(196)
        res(2,3,8) = packlin(238)
        res(2,4,1) = packlin(94)
        res(2,4,2) = packlin(100)
        res(2,4,3) = packlin(106)
        res(2,4,4) = packlin(112)
        res(2,4,5) = packlin(136)
        res(2,4,6) = packlin(166)
        res(2,4,7) = packlin(202)
        res(2,4,8) = packlin(244)
        res(2,5,1) = packlin(118)
        res(2,5,2) = packlin(124)
        res(2,5,3) = packlin(130)
        res(2,5,4) = packlin(136)
        res(2,5,5) = packlin(142)
        res(2,5,6) = packlin(172)
        res(2,5,7) = packlin(208)
        res(2,5,8) = packlin(250)
        res(2,6,1) = packlin(148)
        res(2,6,2) = packlin(154)
        res(2,6,3) = packlin(160)
        res(2,6,4) = packlin(166)
        res(2,6,5) = packlin(172)
        res(2,6,6) = packlin(178)
        res(2,6,7) = packlin(214)
        res(2,6,8) = packlin(256)
        res(2,7,1) = packlin(184)
        res(2,7,2) = packlin(190)
        res(2,7,3) = packlin(196)
        res(2,7,4) = packlin(202)
        res(2,7,5) = packlin(208)
        res(2,7,6) = packlin(214)
        res(2,7,7) = packlin(220)
        res(2,7,8) = packlin(262)
        res(2,8,1) = packlin(226)
        res(2,8,2) = packlin(232)
        res(2,8,3) = packlin(238)
        res(2,8,4) = packlin(244)
        res(2,8,5) = packlin(250)
        res(2,8,6) = packlin(256)
        res(2,8,7) = packlin(262)
        res(2,8,8) = packlin(268)
        res(3,1,1) = packlin(59)
        res(3,1,2) = packlin(65)
        res(3,1,3) = packlin(77)
        res(3,1,4) = packlin(95)
        res(3,1,5) = packlin(119)
        res(3,1,6) = packlin(149)
        res(3,1,7) = packlin(185)
        res(3,1,8) = packlin(227)
        res(3,2,1) = packlin(65)
        res(3,2,2) = packlin(71)
        res(3,2,3) = packlin(83)
        res(3,2,4) = packlin(101)
        res(3,2,5) = packlin(125)
        res(3,2,6) = packlin(155)
        res(3,2,7) = packlin(191)
        res(3,2,8) = packlin(233)
        res(3,3,1) = packlin(77)
        res(3,3,2) = packlin(83)
        res(3,3,3) = packlin(89)
        res(3,3,4) = packlin(107)
        res(3,3,5) = packlin(131)
        res(3,3,6) = packlin(161)
        res(3,3,7) = packlin(197)
        res(3,3,8) = packlin(239)
        res(3,4,1) = packlin(95)
        res(3,4,2) = packlin(101)
        res(3,4,3) = packlin(107)
        res(3,4,4) = packlin(113)
        res(3,4,5) = packlin(137)
        res(3,4,6) = packlin(167)
        res(3,4,7) = packlin(203)
        res(3,4,8) = packlin(245)
        res(3,5,1) = packlin(119)
        res(3,5,2) = packlin(125)
        res(3,5,3) = packlin(131)
        res(3,5,4) = packlin(137)
        res(3,5,5) = packlin(143)
        res(3,5,6) = packlin(173)
        res(3,5,7) = packlin(209)
        res(3,5,8) = packlin(251)
        res(3,6,1) = packlin(149)
        res(3,6,2) = packlin(155)
        res(3,6,3) = packlin(161)
        res(3,6,4) = packlin(167)
        res(3,6,5) = packlin(173)
        res(3,6,6) = packlin(179)
        res(3,6,7) = packlin(215)
        res(3,6,8) = packlin(257)
        res(3,7,1) = packlin(185)
        res(3,7,2) = packlin(191)
        res(3,7,3) = packlin(197)
        res(3,7,4) = packlin(203)
        res(3,7,5) = packlin(209)
        res(3,7,6) = packlin(215)
        res(3,7,7) = packlin(221)
        res(3,7,8) = packlin(263)
        res(3,8,1) = packlin(227)
        res(3,8,2) = packlin(233)
        res(3,8,3) = packlin(239)
        res(3,8,4) = packlin(245)
        res(3,8,5) = packlin(251)
        res(3,8,6) = packlin(257)
        res(3,8,7) = packlin(263)
        res(3,8,8) = packlin(269)
        res(4,1,1) = packlin(60)
        res(4,1,2) = packlin(66)
        res(4,1,3) = packlin(78)
        res(4,1,4) = packlin(96)
        res(4,1,5) = packlin(120)
        res(4,1,6) = packlin(150)
        res(4,1,7) = packlin(186)
        res(4,1,8) = packlin(228)
        res(4,2,1) = packlin(66)
        res(4,2,2) = packlin(72)
        res(4,2,3) = packlin(84)
        res(4,2,4) = packlin(102)
        res(4,2,5) = packlin(126)
        res(4,2,6) = packlin(156)
        res(4,2,7) = packlin(192)
        res(4,2,8) = packlin(234)
        res(4,3,1) = packlin(78)
        res(4,3,2) = packlin(84)
        res(4,3,3) = packlin(90)
        res(4,3,4) = packlin(108)
        res(4,3,5) = packlin(132)
        res(4,3,6) = packlin(162)
        res(4,3,7) = packlin(198)
        res(4,3,8) = packlin(240)
        res(4,4,1) = packlin(96)
        res(4,4,2) = packlin(102)
        res(4,4,3) = packlin(108)
        res(4,4,4) = packlin(114)
        res(4,4,5) = packlin(138)
        res(4,4,6) = packlin(168)
        res(4,4,7) = packlin(204)
        res(4,4,8) = packlin(246)
        res(4,5,1) = packlin(120)
        res(4,5,2) = packlin(126)
        res(4,5,3) = packlin(132)
        res(4,5,4) = packlin(138)
        res(4,5,5) = packlin(144)
        res(4,5,6) = packlin(174)
        res(4,5,7) = packlin(210)
        res(4,5,8) = packlin(252)
        res(4,6,1) = packlin(150)
        res(4,6,2) = packlin(156)
        res(4,6,3) = packlin(162)
        res(4,6,4) = packlin(168)
        res(4,6,5) = packlin(174)
        res(4,6,6) = packlin(180)
        res(4,6,7) = packlin(216)
        res(4,6,8) = packlin(258)
        res(4,7,1) = packlin(186)
        res(4,7,2) = packlin(192)
        res(4,7,3) = packlin(198)
        res(4,7,4) = packlin(204)
        res(4,7,5) = packlin(210)
        res(4,7,6) = packlin(216)
        res(4,7,7) = packlin(222)
        res(4,7,8) = packlin(264)
        res(4,8,1) = packlin(228)
        res(4,8,2) = packlin(234)
        res(4,8,3) = packlin(240)
        res(4,8,4) = packlin(246)
        res(4,8,5) = packlin(252)
        res(4,8,6) = packlin(258)
        res(4,8,7) = packlin(264)
        res(4,8,8) = packlin(270)
        res(5,1,1) = packlin(61)
        res(5,1,2) = packlin(67)
        res(5,1,3) = packlin(79)
        res(5,1,4) = packlin(97)
        res(5,1,5) = packlin(121)
        res(5,1,6) = packlin(151)
        res(5,1,7) = packlin(187)
        res(5,1,8) = packlin(229)
        res(5,2,1) = packlin(67)
        res(5,2,2) = packlin(73)
        res(5,2,3) = packlin(85)
        res(5,2,4) = packlin(103)
        res(5,2,5) = packlin(127)
        res(5,2,6) = packlin(157)
        res(5,2,7) = packlin(193)
        res(5,2,8) = packlin(235)
        res(5,3,1) = packlin(79)
        res(5,3,2) = packlin(85)
        res(5,3,3) = packlin(91)
        res(5,3,4) = packlin(109)
        res(5,3,5) = packlin(133)
        res(5,3,6) = packlin(163)
        res(5,3,7) = packlin(199)
        res(5,3,8) = packlin(241)
        res(5,4,1) = packlin(97)
        res(5,4,2) = packlin(103)
        res(5,4,3) = packlin(109)
        res(5,4,4) = packlin(115)
        res(5,4,5) = packlin(139)
        res(5,4,6) = packlin(169)
        res(5,4,7) = packlin(205)
        res(5,4,8) = packlin(247)
        res(5,5,1) = packlin(121)
        res(5,5,2) = packlin(127)
        res(5,5,3) = packlin(133)
        res(5,5,4) = packlin(139)
        res(5,5,5) = packlin(145)
        res(5,5,6) = packlin(175)
        res(5,5,7) = packlin(211)
        res(5,5,8) = packlin(253)
        res(5,6,1) = packlin(151)
        res(5,6,2) = packlin(157)
        res(5,6,3) = packlin(163)
        res(5,6,4) = packlin(169)
        res(5,6,5) = packlin(175)
        res(5,6,6) = packlin(181)
        res(5,6,7) = packlin(217)
        res(5,6,8) = packlin(259)
        res(5,7,1) = packlin(187)
        res(5,7,2) = packlin(193)
        res(5,7,3) = packlin(199)
        res(5,7,4) = packlin(205)
        res(5,7,5) = packlin(211)
        res(5,7,6) = packlin(217)
        res(5,7,7) = packlin(223)
        res(5,7,8) = packlin(265)
        res(5,8,1) = packlin(229)
        res(5,8,2) = packlin(235)
        res(5,8,3) = packlin(241)
        res(5,8,4) = packlin(247)
        res(5,8,5) = packlin(253)
        res(5,8,6) = packlin(259)
        res(5,8,7) = packlin(265)
        res(5,8,8) = packlin(271)
        res(6,1,1) = packlin(62)
        res(6,1,2) = packlin(68)
        res(6,1,3) = packlin(80)
        res(6,1,4) = packlin(98)
        res(6,1,5) = packlin(122)
        res(6,1,6) = packlin(152)
        res(6,1,7) = packlin(188)
        res(6,1,8) = packlin(230)
        res(6,2,1) = packlin(68)
        res(6,2,2) = packlin(74)
        res(6,2,3) = packlin(86)
        res(6,2,4) = packlin(104)
        res(6,2,5) = packlin(128)
        res(6,2,6) = packlin(158)
        res(6,2,7) = packlin(194)
        res(6,2,8) = packlin(236)
        res(6,3,1) = packlin(80)
        res(6,3,2) = packlin(86)
        res(6,3,3) = packlin(92)
        res(6,3,4) = packlin(110)
        res(6,3,5) = packlin(134)
        res(6,3,6) = packlin(164)
        res(6,3,7) = packlin(200)
        res(6,3,8) = packlin(242)
        res(6,4,1) = packlin(98)
        res(6,4,2) = packlin(104)
        res(6,4,3) = packlin(110)
        res(6,4,4) = packlin(116)
        res(6,4,5) = packlin(140)
        res(6,4,6) = packlin(170)
        res(6,4,7) = packlin(206)
        res(6,4,8) = packlin(248)
        res(6,5,1) = packlin(122)
        res(6,5,2) = packlin(128)
        res(6,5,3) = packlin(134)
        res(6,5,4) = packlin(140)
        res(6,5,5) = packlin(146)
        res(6,5,6) = packlin(176)
        res(6,5,7) = packlin(212)
        res(6,5,8) = packlin(254)
        res(6,6,1) = packlin(152)
        res(6,6,2) = packlin(158)
        res(6,6,3) = packlin(164)
        res(6,6,4) = packlin(170)
        res(6,6,5) = packlin(176)
        res(6,6,6) = packlin(182)
        res(6,6,7) = packlin(218)
        res(6,6,8) = packlin(260)
        res(6,7,1) = packlin(188)
        res(6,7,2) = packlin(194)
        res(6,7,3) = packlin(200)
        res(6,7,4) = packlin(206)
        res(6,7,5) = packlin(212)
        res(6,7,6) = packlin(218)
        res(6,7,7) = packlin(224)
        res(6,7,8) = packlin(266)
        res(6,8,1) = packlin(230)
        res(6,8,2) = packlin(236)
        res(6,8,3) = packlin(242)
        res(6,8,4) = packlin(248)
        res(6,8,5) = packlin(254)
        res(6,8,6) = packlin(260)
        res(6,8,7) = packlin(266)
        res(6,8,8) = packlin(272)
    end function stt
    function stt_i(self,t) result(res)
        !! Return a regularized stm at time t
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        real(dp), dimension(n,n,n) :: res
        !! This can be improved to get both the STM and STT
        !! with one binary search
        res = sttinvert(self%stm(t),self%stt(t), n)
    end function stt_i
    function prop_once(self,ta, tb, xa, order) result(res)
        !! Propagates the relative state xa at ta
        !! to tb
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: ta, tb, xa(n)
        integer, intent(in), optional :: order
        integer :: ord
        !! xa is the initial relative state 
        !! should be dimension 6
        real(dp), dimension(n)   :: res
        real(dp) :: stmab(n,n), sttab(n,n,n)
        ord = 2
        if (present(order)) ord=order
        select case (ord)
        case (1)
            stmab = matmul(self%stm(tb),self%stm_i(ta))
            res(:n) = matmul(stmab,xa)
        case default
            call sttchain(self%stm_i(ta),self%stt_i(ta), &
                           & self%stm(tb), self%stt(tb), &
                           & stmab, sttab, n)
            res(:n) = matmul(stmab,xa) + 0.5_dp*vectensquad(xa,sttab,n)
        end select 
    end function prop_once
    function prop_many(self,ta, tb, xa, order) result(res)
        !! Propagates the relative state xa at ta
        !! to tb
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: ta, tb, xa(:,:) ! n rows, : columns
        integer, intent(in), optional :: order
        integer :: ord, i
        !! xa is the initial relative state 
        !! should be dimension 8
        real(dp) :: res(n,size(xa,2))
        real(dp) :: stmab(n,n), sttab(n,n,n)
        ord = 2
        if (present(order)) ord=order
        call sttchain(self%stm_i(ta),self%stt_i(ta), &
                       & self%stm(tb), self%stt(tb), &
                       & stmab, sttab, n)
        select case (ord)
        case (1)
            res(:n,:) = matmul(stmab,xa)
        case default
            !$OMP PARALLEL DO
            do i=1,size(xa,2)
                res(:n,i) = matmul(stmab,xa(:,i)) + &
                           0.5_dp*vectensquad(xa(:,i),sttab,n)
            end do
            !$OMP END PARALLEL DO
        end select 
    end function prop_many
    subroutine stts_ab(self, ta, tb, stm, stt)
        !! Return the STM and STT from ta to tb
        class(Itraj), intent(inout)  :: self
        real(dp),     intent(in)  :: ta, tb
        real(dp),     intent(out) :: stm(n,n), stt(n,n,n)

        call sttchain_invert(self%stm(ta),self%stt(ta), &
                    & self%stm(tb), self%stt(tb), &
                    & stm, stt, n)
    end subroutine stts_ab
    function zmap(self, t,order) result(res)
        !! Testing procedure that chains an stm and stt with their own
        !! inverse--you ought to get I and 0, respectively
        class(Itraj),      intent(inout) :: self
        real(dp),          intent(in) :: t
        integer, optional, intent(in) :: order
        integer                       :: ord, i
        real(dp)                      :: stm(n,n), stt(n,n,n), eye(n,n)
        real(dp)                      :: res

        ord=2
        if (present(order)) ord=order

        res = 0._dp
        eye = 0._dp
        do i=1,n; eye(i,i) = 1._dp; end do
        call sttchain(self%stm_i(t),self%stt_i(t), &
                    & self%stm(t), self%stt(t), &
                    & stm, stt, n)

        res = norm2(eye - stm)
        if (ord.eq.2) res = res + norm2(stt)
    end function zmap
end module qist
