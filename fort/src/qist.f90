module qist
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use globals
    use tensorops, only:  sttchain, sttchain_invert, vectensquad, &
                          stminvert, sttinvert, vectens3, mattens, &
                          quad
    use denseLight, only: lightSol
    use frkmin,     only: odesolution
    implicit none
    private
    type, public :: Itraj
        real(dp)                 :: t0, tf, tof
        logical                  :: initq, regularized
        type(lightSol)           :: reftraj
        type(odesolution)        :: kvtau
        contains
            generic,   public    :: init => init_nml, init_var
            procedure            :: call_raw
            procedure            :: call 
            procedure            :: state 
            procedure            :: stm 
            procedure            :: stm_i 
            procedure            :: stt
            procedure            :: stt_i 
            generic,   public    :: prop => prop_once, prop_many
            procedure, private   :: prop_once
            procedure, private   :: prop_many
            procedure, private   :: init_nml
            procedure, private   :: init_var
            procedure            :: stts_ab
            procedure            :: stt_update
            procedure            :: tensor_change_of_basis
            procedure            :: zmap
    end type Itraj
    contains
    ! Initializer
    subroutine init_nml(self, namefile)
        character(len=*), intent(in) :: namefile
        real(dp)                     :: t0, tf
        character(len=1000)          :: qist_filename
        character(len=1000)          :: kvtau_filename
        integer                      :: stat, num
        class(Itraj),  intent(inout) :: self
        namelist /ITRAJ_CONFIG/ qist_filename, &
                                kvtau_filename, &
                                t0, &
                                tf
        t0 = 0._dp
        tf = 0._dp
        qist_filename = ""
        kvtau_filename = ""
        inquire(file=trim(adjustl(namefile)), iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: Bad ITRAJ config namelist filename"
            print *, "error code ", stat
            stop
        end if
        open(file=namefile, status="old", &
             iostat=stat,newunit=num)
        read(unit=num, nml=ITRAJ_CONFIG, iostat=stat)
        if (stat .ne. 0) then 
            print *, "ERROR: bad ITRAJ config namelist format"
            print *, "error code ", stat
            stop
        end if
        close(num)
        call self%init(t0,tf,qist_filename,kvtau_filename)
    end subroutine init_nml
    subroutine init_var(self, t0, tf, trajfile, kvtaufile)
        character(len=*), intent(in) :: trajfile, kvtaufile
        real(dp),         intent(in) :: t0, tf
        integer                      :: stat, num
        class(Itraj),  intent(inout) :: self
        self%initq=.true.
        self%t0               = t0
        self%tf               = tf
        self%tof              = tf-t0
        self%regularized      = .false.
        inquire(file=trim(adjustl(trajfile)), iostat=stat)
        if (stat.ne.0) then
            print *, "ERROR: Trajectory file not found."
            print *, "error code ", stat
            stop
        end if
        open(newunit=num, file=trim(adjustl(trajfile)), &
             status="old", access="stream",iostat=stat)
        call self%reftraj%read(num)
        close(num)
        if (kvtaufile.ne."") then
            self%regularized = .true.
            inquire(file=trim(adjustl(kvtaufile)), iostat=stat)
            if (stat.ne.0) then
                print *, "ERROR: Trajectory file not found."
                print *, "error code ", stat
                stop
            end if
            open(newunit=num, file=trim(adjustl(kvtaufile)), &
                 status="old", access="stream",iostat=stat)
            call self%kvtau%read(num)
            close(num)
        end if
    end subroutine init_var
    !! Calling functions 
    function call_raw(self,t,lind,uind) result(res)
        !! Return a (uind-lind)-dimensional at t
        !! for this function, t goes from 0 to 1.
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
    function call(self,t,lind,uind) result(res)
        !! Return a (uind-lind)-dimensional at t
        !! If no lind (uind) is passed, beginning (end) of the vector is used
        class(Itraj),      intent(inout) :: self
        real(dp),          intent(in) :: t
        integer, optional, intent(in) :: uind, lind
        !! The value of the independent variable to generate
        real(dp), allocatable         :: res(:)
        real(dp)                      :: t_elapsed, tau, k, dum(1)
        integer l, u
        l = 1
        u = plen
        t_elapsed = t - self%t0
        tau = t_elapsed/self%tof
        if (self%regularized) then
            dum = self%kvtau%call(tau)
            k = dum(1)
        else
            k = tau
        end if
        if (present(lind)) l=lind
        if (present(uind)) u=uind
        allocate(res(u-l+1))
        res = self%call_raw(k,l,u)
    end function
    function state(self,t) result(res)
        !! Return a state at time t
        class(ITraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp)                 :: dum(1)
        real(dp)                 :: res
        !! The returned state
        dum = self%call(t,uind=1)
        res = dum(1)
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
        res(1,1) = packlin(2)
        res(1,2) = packlin(8)
        res(1,3) = packlin(14)
        res(1,4) = packlin(20)
        res(1,5) = packlin(26)
        res(1,6) = packlin(32)
        res(1,7) = packlin(38)
        res(1,8) = packlin(44)
        res(2,1) = packlin(3)
        res(2,2) = packlin(9)
        res(2,3) = packlin(15)
        res(2,4) = packlin(21)
        res(2,5) = packlin(27)
        res(2,6) = packlin(33)
        res(2,7) = packlin(39)
        res(2,8) = packlin(45)
        res(3,1) = packlin(4)
        res(3,2) = packlin(10)
        res(3,3) = packlin(16)
        res(3,4) = packlin(22)
        res(3,5) = packlin(28)
        res(3,6) = packlin(34)
        res(3,7) = packlin(40)
        res(3,8) = packlin(46)
        res(4,1) = packlin(5)
        res(4,2) = packlin(11)
        res(4,3) = packlin(17)
        res(4,4) = packlin(23)
        res(4,5) = packlin(29)
        res(4,6) = packlin(35)
        res(4,7) = packlin(41)
        res(4,8) = packlin(47)
        res(5,1) = packlin(6)
        res(5,2) = packlin(12)
        res(5,3) = packlin(18)
        res(5,4) = packlin(24)
        res(5,5) = packlin(30)
        res(5,6) = packlin(36)
        res(5,7) = packlin(42)
        res(5,8) = packlin(48)
        res(6,1) = packlin(7)
        res(6,2) = packlin(13)
        res(6,3) = packlin(19)
        res(6,4) = packlin(25)
        res(6,5) = packlin(31)
        res(6,6) = packlin(37)
        res(6,7) = packlin(43)
        res(6,8) = packlin(49)
        res(7,7)=1._dp
        res(7,8)= (t-self%t0)/(self%tf-self%t0)
        res(8,8)=1._dp
    end function stm
    function stm_i(self,t) result(res)
        !! Return an stm at time t
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
        res = 0._dp
        packlin = self%call(t)
        res(1,1,1) = packlin(50)
        res(1,1,2) = packlin(56)
        res(1,1,3) = packlin(68)
        res(1,1,4) = packlin(86)
        res(1,1,5) = packlin(110)
        res(1,1,6) = packlin(140)
        res(1,1,7) = packlin(176)
        res(1,1,8) = packlin(218)
        res(1,2,1) = packlin(56)
        res(1,2,2) = packlin(62)
        res(1,2,3) = packlin(74)
        res(1,2,4) = packlin(92)
        res(1,2,5) = packlin(116)
        res(1,2,6) = packlin(146)
        res(1,2,7) = packlin(182)
        res(1,2,8) = packlin(224)
        res(1,3,1) = packlin(68)
        res(1,3,2) = packlin(74)
        res(1,3,3) = packlin(80)
        res(1,3,4) = packlin(98)
        res(1,3,5) = packlin(122)
        res(1,3,6) = packlin(152)
        res(1,3,7) = packlin(188)
        res(1,3,8) = packlin(230)
        res(1,4,1) = packlin(86)
        res(1,4,2) = packlin(92)
        res(1,4,3) = packlin(98)
        res(1,4,4) = packlin(104)
        res(1,4,5) = packlin(128)
        res(1,4,6) = packlin(158)
        res(1,4,7) = packlin(194)
        res(1,4,8) = packlin(236)
        res(1,5,1) = packlin(110)
        res(1,5,2) = packlin(116)
        res(1,5,3) = packlin(122)
        res(1,5,4) = packlin(128)
        res(1,5,5) = packlin(134)
        res(1,5,6) = packlin(164)
        res(1,5,7) = packlin(200)
        res(1,5,8) = packlin(242)
        res(1,6,1) = packlin(140)
        res(1,6,2) = packlin(146)
        res(1,6,3) = packlin(152)
        res(1,6,4) = packlin(158)
        res(1,6,5) = packlin(164)
        res(1,6,6) = packlin(170)
        res(1,6,7) = packlin(206)
        res(1,6,8) = packlin(248)
        res(1,7,1) = packlin(176)
        res(1,7,2) = packlin(182)
        res(1,7,3) = packlin(188)
        res(1,7,4) = packlin(194)
        res(1,7,5) = packlin(200)
        res(1,7,6) = packlin(206)
        res(1,7,7) = packlin(212)
        res(1,7,8) = packlin(254)
        res(1,8,1) = packlin(218)
        res(1,8,2) = packlin(224)
        res(1,8,3) = packlin(230)
        res(1,8,4) = packlin(236)
        res(1,8,5) = packlin(242)
        res(1,8,6) = packlin(248)
        res(1,8,7) = packlin(254)
        res(1,8,8) = packlin(260)
        res(2,1,1) = packlin(51)
        res(2,1,2) = packlin(57)
        res(2,1,3) = packlin(69)
        res(2,1,4) = packlin(87)
        res(2,1,5) = packlin(111)
        res(2,1,6) = packlin(141)
        res(2,1,7) = packlin(177)
        res(2,1,8) = packlin(219)
        res(2,2,1) = packlin(57)
        res(2,2,2) = packlin(63)
        res(2,2,3) = packlin(75)
        res(2,2,4) = packlin(93)
        res(2,2,5) = packlin(117)
        res(2,2,6) = packlin(147)
        res(2,2,7) = packlin(183)
        res(2,2,8) = packlin(225)
        res(2,3,1) = packlin(69)
        res(2,3,2) = packlin(75)
        res(2,3,3) = packlin(81)
        res(2,3,4) = packlin(99)
        res(2,3,5) = packlin(123)
        res(2,3,6) = packlin(153)
        res(2,3,7) = packlin(189)
        res(2,3,8) = packlin(231)
        res(2,4,1) = packlin(87)
        res(2,4,2) = packlin(93)
        res(2,4,3) = packlin(99)
        res(2,4,4) = packlin(105)
        res(2,4,5) = packlin(129)
        res(2,4,6) = packlin(159)
        res(2,4,7) = packlin(195)
        res(2,4,8) = packlin(237)
        res(2,5,1) = packlin(111)
        res(2,5,2) = packlin(117)
        res(2,5,3) = packlin(123)
        res(2,5,4) = packlin(129)
        res(2,5,5) = packlin(135)
        res(2,5,6) = packlin(165)
        res(2,5,7) = packlin(201)
        res(2,5,8) = packlin(243)
        res(2,6,1) = packlin(141)
        res(2,6,2) = packlin(147)
        res(2,6,3) = packlin(153)
        res(2,6,4) = packlin(159)
        res(2,6,5) = packlin(165)
        res(2,6,6) = packlin(171)
        res(2,6,7) = packlin(207)
        res(2,6,8) = packlin(249)
        res(2,7,1) = packlin(177)
        res(2,7,2) = packlin(183)
        res(2,7,3) = packlin(189)
        res(2,7,4) = packlin(195)
        res(2,7,5) = packlin(201)
        res(2,7,6) = packlin(207)
        res(2,7,7) = packlin(213)
        res(2,7,8) = packlin(255)
        res(2,8,1) = packlin(219)
        res(2,8,2) = packlin(225)
        res(2,8,3) = packlin(231)
        res(2,8,4) = packlin(237)
        res(2,8,5) = packlin(243)
        res(2,8,6) = packlin(249)
        res(2,8,7) = packlin(255)
        res(2,8,8) = packlin(261)
        res(3,1,1) = packlin(52)
        res(3,1,2) = packlin(58)
        res(3,1,3) = packlin(70)
        res(3,1,4) = packlin(88)
        res(3,1,5) = packlin(112)
        res(3,1,6) = packlin(142)
        res(3,1,7) = packlin(178)
        res(3,1,8) = packlin(220)
        res(3,2,1) = packlin(58)
        res(3,2,2) = packlin(64)
        res(3,2,3) = packlin(76)
        res(3,2,4) = packlin(94)
        res(3,2,5) = packlin(118)
        res(3,2,6) = packlin(148)
        res(3,2,7) = packlin(184)
        res(3,2,8) = packlin(226)
        res(3,3,1) = packlin(70)
        res(3,3,2) = packlin(76)
        res(3,3,3) = packlin(82)
        res(3,3,4) = packlin(100)
        res(3,3,5) = packlin(124)
        res(3,3,6) = packlin(154)
        res(3,3,7) = packlin(190)
        res(3,3,8) = packlin(232)
        res(3,4,1) = packlin(88)
        res(3,4,2) = packlin(94)
        res(3,4,3) = packlin(100)
        res(3,4,4) = packlin(106)
        res(3,4,5) = packlin(130)
        res(3,4,6) = packlin(160)
        res(3,4,7) = packlin(196)
        res(3,4,8) = packlin(238)
        res(3,5,1) = packlin(112)
        res(3,5,2) = packlin(118)
        res(3,5,3) = packlin(124)
        res(3,5,4) = packlin(130)
        res(3,5,5) = packlin(136)
        res(3,5,6) = packlin(166)
        res(3,5,7) = packlin(202)
        res(3,5,8) = packlin(244)
        res(3,6,1) = packlin(142)
        res(3,6,2) = packlin(148)
        res(3,6,3) = packlin(154)
        res(3,6,4) = packlin(160)
        res(3,6,5) = packlin(166)
        res(3,6,6) = packlin(172)
        res(3,6,7) = packlin(208)
        res(3,6,8) = packlin(250)
        res(3,7,1) = packlin(178)
        res(3,7,2) = packlin(184)
        res(3,7,3) = packlin(190)
        res(3,7,4) = packlin(196)
        res(3,7,5) = packlin(202)
        res(3,7,6) = packlin(208)
        res(3,7,7) = packlin(214)
        res(3,7,8) = packlin(256)
        res(3,8,1) = packlin(220)
        res(3,8,2) = packlin(226)
        res(3,8,3) = packlin(232)
        res(3,8,4) = packlin(238)
        res(3,8,5) = packlin(244)
        res(3,8,6) = packlin(250)
        res(3,8,7) = packlin(256)
        res(3,8,8) = packlin(262)
        res(4,1,1) = packlin(53)
        res(4,1,2) = packlin(59)
        res(4,1,3) = packlin(71)
        res(4,1,4) = packlin(89)
        res(4,1,5) = packlin(113)
        res(4,1,6) = packlin(143)
        res(4,1,7) = packlin(179)
        res(4,1,8) = packlin(221)
        res(4,2,1) = packlin(59)
        res(4,2,2) = packlin(65)
        res(4,2,3) = packlin(77)
        res(4,2,4) = packlin(95)
        res(4,2,5) = packlin(119)
        res(4,2,6) = packlin(149)
        res(4,2,7) = packlin(185)
        res(4,2,8) = packlin(227)
        res(4,3,1) = packlin(71)
        res(4,3,2) = packlin(77)
        res(4,3,3) = packlin(83)
        res(4,3,4) = packlin(101)
        res(4,3,5) = packlin(125)
        res(4,3,6) = packlin(155)
        res(4,3,7) = packlin(191)
        res(4,3,8) = packlin(233)
        res(4,4,1) = packlin(89)
        res(4,4,2) = packlin(95)
        res(4,4,3) = packlin(101)
        res(4,4,4) = packlin(107)
        res(4,4,5) = packlin(131)
        res(4,4,6) = packlin(161)
        res(4,4,7) = packlin(197)
        res(4,4,8) = packlin(239)
        res(4,5,1) = packlin(113)
        res(4,5,2) = packlin(119)
        res(4,5,3) = packlin(125)
        res(4,5,4) = packlin(131)
        res(4,5,5) = packlin(137)
        res(4,5,6) = packlin(167)
        res(4,5,7) = packlin(203)
        res(4,5,8) = packlin(245)
        res(4,6,1) = packlin(143)
        res(4,6,2) = packlin(149)
        res(4,6,3) = packlin(155)
        res(4,6,4) = packlin(161)
        res(4,6,5) = packlin(167)
        res(4,6,6) = packlin(173)
        res(4,6,7) = packlin(209)
        res(4,6,8) = packlin(251)
        res(4,7,1) = packlin(179)
        res(4,7,2) = packlin(185)
        res(4,7,3) = packlin(191)
        res(4,7,4) = packlin(197)
        res(4,7,5) = packlin(203)
        res(4,7,6) = packlin(209)
        res(4,7,7) = packlin(215)
        res(4,7,8) = packlin(257)
        res(4,8,1) = packlin(221)
        res(4,8,2) = packlin(227)
        res(4,8,3) = packlin(233)
        res(4,8,4) = packlin(239)
        res(4,8,5) = packlin(245)
        res(4,8,6) = packlin(251)
        res(4,8,7) = packlin(257)
        res(4,8,8) = packlin(263)
        res(5,1,1) = packlin(54)
        res(5,1,2) = packlin(60)
        res(5,1,3) = packlin(72)
        res(5,1,4) = packlin(90)
        res(5,1,5) = packlin(114)
        res(5,1,6) = packlin(144)
        res(5,1,7) = packlin(180)
        res(5,1,8) = packlin(222)
        res(5,2,1) = packlin(60)
        res(5,2,2) = packlin(66)
        res(5,2,3) = packlin(78)
        res(5,2,4) = packlin(96)
        res(5,2,5) = packlin(120)
        res(5,2,6) = packlin(150)
        res(5,2,7) = packlin(186)
        res(5,2,8) = packlin(228)
        res(5,3,1) = packlin(72)
        res(5,3,2) = packlin(78)
        res(5,3,3) = packlin(84)
        res(5,3,4) = packlin(102)
        res(5,3,5) = packlin(126)
        res(5,3,6) = packlin(156)
        res(5,3,7) = packlin(192)
        res(5,3,8) = packlin(234)
        res(5,4,1) = packlin(90)
        res(5,4,2) = packlin(96)
        res(5,4,3) = packlin(102)
        res(5,4,4) = packlin(108)
        res(5,4,5) = packlin(132)
        res(5,4,6) = packlin(162)
        res(5,4,7) = packlin(198)
        res(5,4,8) = packlin(240)
        res(5,5,1) = packlin(114)
        res(5,5,2) = packlin(120)
        res(5,5,3) = packlin(126)
        res(5,5,4) = packlin(132)
        res(5,5,5) = packlin(138)
        res(5,5,6) = packlin(168)
        res(5,5,7) = packlin(204)
        res(5,5,8) = packlin(246)
        res(5,6,1) = packlin(144)
        res(5,6,2) = packlin(150)
        res(5,6,3) = packlin(156)
        res(5,6,4) = packlin(162)
        res(5,6,5) = packlin(168)
        res(5,6,6) = packlin(174)
        res(5,6,7) = packlin(210)
        res(5,6,8) = packlin(252)
        res(5,7,1) = packlin(180)
        res(5,7,2) = packlin(186)
        res(5,7,3) = packlin(192)
        res(5,7,4) = packlin(198)
        res(5,7,5) = packlin(204)
        res(5,7,6) = packlin(210)
        res(5,7,7) = packlin(216)
        res(5,7,8) = packlin(258)
        res(5,8,1) = packlin(222)
        res(5,8,2) = packlin(228)
        res(5,8,3) = packlin(234)
        res(5,8,4) = packlin(240)
        res(5,8,5) = packlin(246)
        res(5,8,6) = packlin(252)
        res(5,8,7) = packlin(258)
        res(5,8,8) = packlin(264)
        res(6,1,1) = packlin(55)
        res(6,1,2) = packlin(61)
        res(6,1,3) = packlin(73)
        res(6,1,4) = packlin(91)
        res(6,1,5) = packlin(115)
        res(6,1,6) = packlin(145)
        res(6,1,7) = packlin(181)
        res(6,1,8) = packlin(223)
        res(6,2,1) = packlin(61)
        res(6,2,2) = packlin(67)
        res(6,2,3) = packlin(79)
        res(6,2,4) = packlin(97)
        res(6,2,5) = packlin(121)
        res(6,2,6) = packlin(151)
        res(6,2,7) = packlin(187)
        res(6,2,8) = packlin(229)
        res(6,3,1) = packlin(73)
        res(6,3,2) = packlin(79)
        res(6,3,3) = packlin(85)
        res(6,3,4) = packlin(103)
        res(6,3,5) = packlin(127)
        res(6,3,6) = packlin(157)
        res(6,3,7) = packlin(193)
        res(6,3,8) = packlin(235)
        res(6,4,1) = packlin(91)
        res(6,4,2) = packlin(97)
        res(6,4,3) = packlin(103)
        res(6,4,4) = packlin(109)
        res(6,4,5) = packlin(133)
        res(6,4,6) = packlin(163)
        res(6,4,7) = packlin(199)
        res(6,4,8) = packlin(241)
        res(6,5,1) = packlin(115)
        res(6,5,2) = packlin(121)
        res(6,5,3) = packlin(127)
        res(6,5,4) = packlin(133)
        res(6,5,5) = packlin(139)
        res(6,5,6) = packlin(169)
        res(6,5,7) = packlin(205)
        res(6,5,8) = packlin(247)
        res(6,6,1) = packlin(145)
        res(6,6,2) = packlin(151)
        res(6,6,3) = packlin(157)
        res(6,6,4) = packlin(163)
        res(6,6,5) = packlin(169)
        res(6,6,6) = packlin(175)
        res(6,6,7) = packlin(211)
        res(6,6,8) = packlin(253)
        res(6,7,1) = packlin(181)
        res(6,7,2) = packlin(187)
        res(6,7,3) = packlin(193)
        res(6,7,4) = packlin(199)
        res(6,7,5) = packlin(205)
        res(6,7,6) = packlin(211)
        res(6,7,7) = packlin(217)
        res(6,7,8) = packlin(259)
        res(6,8,1) = packlin(223)
        res(6,8,2) = packlin(229)
        res(6,8,3) = packlin(235)
        res(6,8,4) = packlin(241)
        res(6,8,5) = packlin(247)
        res(6,8,6) = packlin(253)
        res(6,8,7) = packlin(259)
        res(6,8,8) = packlin(265)
    end function stt
    function stt_i(self,t) result(res)
        !! Return an inverse stt at time t
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
    subroutine stt_update(self, ta, tb, xa,new_stm,new_stt)
        !! Return the STM and STT for a perturbation to the relative trajectory
        class(Itraj), intent(inout)  :: self
        real(dp),     intent(in)     :: ta, tb, xa(n)
        real(dp),     intent(out)    :: new_stm(n,n), new_stt(n,n,n)
        real(dp)                     :: ref_stm(n,n), ref_stt(n,n,n), &
                                      & stmupdate(n,n)
        call self%stts_ab(ta,tb, ref_stm, ref_stt)
        stmupdate = vectens3(xa,ref_stt,n)
        new_stm = ref_stm + stmupdate
        new_stt = ref_stt
    end subroutine stt_update
    subroutine tensor_change_of_basis(self, RNO, old_stm, old_stt, &
                                      new_stm, new_stt)
        !! Transform an STM and STT from an old coordinate basis to a new one
        !! RNO defined by vec_new = RNO@vec_old
        class(Itraj), intent(inout)  :: self
        real(dp),     intent(in)     :: RNO(n,n), old_stm(n,n), old_stt(n,n,n)
        real(dp),     intent(out)    :: new_stm(n,n), new_stt(n,n,n)
        new_stm = mmult(RNO, mmult(old_stm,transpose(RNO)))
        new_stt = mattens(RNO,quad(transpose(RNO),old_stt,n),n)
    end subroutine tensor_change_of_basis
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
