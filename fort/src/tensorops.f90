module tensorops
    use globals
    implicit none
    interface sttinvert
        module procedure q_sttinvert
        module procedure d_sttinvert
    end interface sttinvert
    interface sttchain_invert
        module procedure q_sttchain_invert
        module procedure d_sttchain_invert
    end interface sttchain_invert
    interface stminvert
        module procedure q_stminvert
        module procedure d_stminvert
    end interface stminvert
    interface stminv
        module procedure q_stminv8
        module procedure d_stminv8
    end interface stminv
    interface sttchain
        module procedure q_sttchain
        module procedure d_sttchain
    end interface sttchain
    interface vectensquad
        module procedure q_vectensquad
        module procedure d_vectensquad
    end interface vectensquad
    interface mattens
        module procedure q_mattens
        module procedure d_mattens
    end interface mattens
    interface quad
        module procedure q_quad
        module procedure d_quad
    end interface quad
    interface vectens1
        module procedure q_vectens1
        module procedure d_vectens1
    end interface vectens1
    interface vectens3
        module procedure q_vectens3
        module procedure d_vectens3
    end interface vectens3

    contains

    pure function eyemat(n) result(res)
        implicit none
        integer, intent(in) :: n
        integer :: i
        real(wp) :: res(n,n)
        res = 0._wp
        do i=1,n
            res(i,i) = 1._wp
        end do
    end function eyemat
    pure function d_eyemat(n) result(res)
        implicit none
        integer, intent(in) :: n
        integer :: i
        real(dp) :: res(n,n)
        res = 0._dp
        do i=1,n
            res(i,i) = 1._dp
        end do
    end function d_eyemat

    pure function zmat(n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp) :: res(n,n)
        res = 0._wp
        res(:n/2,n/2:) = eyemat(n/2)
        res(n/2:,:n/2) = -eyemat(n/2)
    end function zmat
    pure function d_zmat(n) result(res)
        implicit none
        integer, intent(in) :: n
        real(8) :: res(n,n)
        res = 0._8
        res(:n/2,n/2:) = d_eyemat(n/2)
        res(n/2:,:n/2) = -d_eyemat(n/2)
    end function d_zmat

    pure function q_mattens(m,t,n) result(res)
        ! ia, ajk -> ijk
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: m(n,n), t(n,n,n)
        real(wp) :: res(n,n,n)
        res=reshape(mmult(m,reshape(t,[n,n**2])),[n,n,n])
    end function q_mattens
    pure function d_mattens(m,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n,n), t(n,n,n)
        real(dp) :: res(n,n,n)
        res=reshape(mmult(m,reshape(t,[n,n**2])),[n,n,n])
    end function d_mattens

    pure function q_quad(m,t,n) result(res)
        ! iab, aj, bk -> ijk
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: m(n,n), t(n,n,n)
        real(wp) :: res(n,n,n)
        res = reshape(t,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
    end function q_quad
    pure function d_quad(m,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n,n), t(n,n,n)
        real(dp) :: res(n,n,n)
        res = reshape(t,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
    end function d_quad

    pure function q_vectens1(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n,n)
        res = reshape(mmult(v,reshape(t,[n,n**2])),[n,n])
    end function q_vectens1
    pure function d_vectens1(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v(n), t(n,n,n)
        real(dp) :: res(n,n)
        res = reshape(mmult(v,reshape(t,[n,n**2])),[n,n])
    end function d_vectens1

    pure function vectens2(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n,n)
        res = vectens1(v,reshape(t,[n,n,n],order=[2,1,3]),n)
    end function vectens2
    pure function q_vectens3(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n,n)
        res = vectens1(v,reshape(t,[n,n,n],order=[3,1,2]),n)
    end function q_vectens3
    pure function d_vectens3(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v(n), t(n,n,n)
        real(dp) :: res(n,n)
        res = vectens1(v,reshape(t,[n,n,n],order=[3,1,2]),n)
    end function d_vectens3
    
    pure function q_vectensquad(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n)
        res = mmult(v,transpose(vectens3(v,t,n)))
    end function q_vectensquad

    pure function d_vectensquad(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v(n), t(n,n,n)
        real(dp) :: res(n)
        res = mmult(v,transpose(vectens3(v,t,n)))
    end function d_vectensquad

    pure function q_stminvert(stm,n) result(res)
        implicit none
        real(wp), intent(in) :: stm(n,n)
        integer, intent(in) :: n 
        real(wp)             :: res(n,n)
        ! symplectic version
        select case(n)
        case(6)
            res = -mmult(mmult(zmat(n),transpose(stm)),zmat(n))
        case(8)
            res = stminv(stm)
        end select
    end function q_stminvert

    pure function d_stminvert(stm,n) result(res)
        implicit none
        real(8), intent(in) :: stm(n,n)
        integer, intent(in) :: n 
        real(8)             :: res(n,n)
        select case(n)
        case(6)
        ! symplectic version
            res = -mmult(mmult(d_zmat(n),transpose(stm)),d_zmat(n))
        case(8)
            res = stminv(stm)
        end select
    end function d_stminvert

    pure function d_sttinvert(stm,stt,n) result (res)
        implicit none
        real(8), intent(in) :: stm(n,n), stt(n,n,n)
        integer, intent(in) :: n
        real(8) :: istm(n,n), res(n,n,n), inter(n,n,n)
        istm = stminvert(stm,n)
        inter = -mattens(istm,stt,n)
        res = quad(istm,inter,n)
    end function d_sttinvert

    pure function q_sttinvert(stm,stt,n) result (res)
        implicit none
        real(wp), intent(in) :: stm(n,n), stt(n,n,n)
        integer, intent(in) :: n
        real(wp) :: istm(n,n), res(n,n,n), inter(n,n,n)
        istm = stminvert(stm,n)
        inter = -mattens(istm,stt,n)
        res = quad(istm,inter,n)
    end function q_sttinvert
        
    pure subroutine d_sttchain_invert(stm1,stt1,stm2,stt2, stmab, sttab,n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: stm1(n,n), stm2(n,n), stt1(n,n,n), stt2(n,n,n)
        real(8)             :: istm1(n,n), istt1(n,n,n)
        real(8)             :: inter(n,n,n) 
        real(8), intent(out) :: stmab(n,n), sttab(n,n,n)
        istm1 = stminvert(stm1,n)
        istt1 = sttinvert(stm1,stt1,n)

        stmab = mmult(stm2,istm1)
        inter = reshape(istt1,[n,n,n],order=[2,1,3])
        sttab = reshape(mattens(stm2,inter,n),[n,n,n],order=[2,1,3])
        inter = reshape(stt2,[n,n,n],order=[2,1,3])
        inter = mattens(transpose(istm1),inter,n)
        inter = reshape(inter,[n,n,n],order=[2,1,3])
        inter = reshape(inter,[n,n,n],order=[3,2,1])
        inter = mattens(transpose(istm1),inter,n)
        sttab = reshape(inter,[n,n,n],order=[3,2,1])
    end subroutine d_sttchain_invert

    pure subroutine q_sttchain_invert(stm1,stt1,stm2,stt2, stmab, sttab,n)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: stm1(n,n), stm2(n,n), stt1(n,n,n), stt2(n,n,n)
        real(wp)             :: istm1(n,n), istt1(n,n,n)
        real(wp)             :: inter(n,n,n) 
        real(wp), intent(out) :: stmab(n,n), sttab(n,n,n)
        istm1 = stminvert(stm1,n)
        istt1 = sttinvert(stm1,stt1,n)

        stmab = mmult(stm2,istm1)
        inter = reshape(istt1,[n,n,n],order=[2,1,3])
        sttab = reshape(mattens(stm2,inter,n),[n,n,n],order=[2,1,3])
        inter = reshape(stt2,[n,n,n],order=[2,1,3])
        inter = mattens(transpose(istm1),inter,n)
        inter = reshape(inter,[n,n,n],order=[2,1,3])
        inter = reshape(inter,[n,n,n],order=[3,2,1])
        inter = mattens(transpose(istm1),inter,n)
        sttab = reshape(inter,[n,n,n],order=[3,2,1])
    end subroutine q_sttchain_invert
    pure subroutine q_sttchain(istm,istt,stm,stt,stmab,sttab,n) 
        implicit none
        integer, intent(in) :: n
      real(wp), intent(in)  :: istm(n,n), stm(n,n), istt(n,n,n), stt(n,n,n)
        real(wp), intent(out) :: stmab(n,n), sttab(n,n,n)

        stmab = mmult(stm,istm)
        sttab = mattens(stm,istt,n)
        sttab = sttab + quad(istm,stt,n)
    end subroutine q_sttchain
    pure subroutine d_sttchain(istm,istt,stm,stt,stmab,sttab,n) 
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in)  :: istm(n,n), stm(n,n), istt(n,n,n), stt(n,n,n)
        real(dp), intent(out) :: stmab(n,n), sttab(n,n,n)

        stmab = mmult(stm,istm)
        sttab = mattens(stm,istt,n)
        sttab = sttab + quad(istm,stt,n)
    end subroutine d_sttchain
    pure function q_stminv8(m) result(res)
        implicit none
        real(wp), intent(in) :: m(8,8)
        real(wp)             :: x(64),&
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
                              & x120, x121, x122, x123, x124, x125, x126, x127, & 
                              & x128, x129, x130, x131, x132, x133, x134, x135, & 
                              & x136, x137, x138, x139, x140, x141, x142, x143, & 
                              & x144, x145, x146, x147, x148, x149, x150, x151, & 
                              & x152, x153, x154, x155, x156, x157, x158, x159, & 
                              & x160, x161, x162, x163, x164, x165
        real(wp), dimension(64) :: inter
        real(wp), dimension(8,8) :: res
        x = reshape(m,[64])
        x0  =  1._wp/x(1)
        x1  =  x(2)*x0
        x2  =  x(10) - x(9)*x1
        x3  =  1._wp/x2
        x4  =  x(17)*x0
        x5  =  x(18) - x(2)*x4
        x6  =  x(9)*x0
        x7  =  x3*x6
        x8  =  x4 - x5*x7
        x9  =  x(3)*x0
        x10  =  x(11) - x(9)*x9
        x11  =  x10*x3
        x12  =  -x(19) + x(3)*x4 + x11*x5
        x13  =  1._wp/x12
        x14  =  x13*(x1*x11 - x9)
        x15  =  x(25)*x0
        x16  =  x3*(x(2)*x15 - x(26))
        x17  =  -x(25)*x9 + x(27) + x10*x16
        x18  =  x13*x8
        x19  =  x15 + x16*x6 + x17*x18
        x20  =  x(4)*x0
        x21  =  x(12) - x(9)*x20
        x22  =  x21*x3
        x23  =  -x(20) + x(4)*x4 + x22*x5
        x24  =  x13*x23
        x25  =  -x(28) + x(4)*x15 - x16*x21 + x17*x24
        x26  =  1._wp/x25
        x27  =  x26*(-x1*x22 + x14*x23 + x20)
        x28  =  x(33)*x0
        x29  =  x3*(x(2)*x28 - x(34))
        x30  =  x13*(-x(33)*x9 + x(35) + x10*x29)
        x31  =  -x(36) + x(4)*x28 - x21*x29 + x23*x30
        x32  =  x19*x26
        x33  =  -x28 - x29*x6 - x30*x8 + x31*x32
        x34  =  x(5)*x0
        x35  =  x(13) - x(9)*x34
        x36  =  x3*x35
        x37  =  -x(21) + x(5)*x4 + x36*x5
        x38  =  x13*x37
        x39  =  -x(29) + x(5)*x15 - x16*x35 + x17*x38
        x40  =  x26*x39
        x41  =  x(37) - x(5)*x28 + x29*x35 - x30*x37 + x31*x40
        x42  =  1._wp/x41
        x43  =  x42*(x1*x36 - x14*x37 + x27*x39 - x34)
        x44  =  x(41)*x0
        x45  =  x3*(x(2)*x44 - x(42))
        x46  =  x13*(-x(41)*x9 + x(43) + x10*x45)
        x47  =  x26*(x(41)*x20 - x(44) - x21*x45 + x23*x46)
        x48  =  -x(41)*x34 + x(45) + x35*x45 - x37*x46 + x39*x47
        x49  =  x33*x42
        x50  =  -x19*x47 + x44 + x45*x6 + x46*x8 + x48*x49
        x51  =  x(6)*x0
        x52  =  x(14) - x(9)*x51
        x53  =  x3*x52
        x54  =  -x(22) + x(6)*x4 + x5*x53
        x55  =  x13*x54
        x56  =  -x(30) + x(6)*x15 - x16*x52 + x17*x55
        x57  =  x26*x56
        x58  =  x(38) - x(6)*x28 + x29*x52 - x30*x54 + x31*x57
        x59  =  x42*x58
        x60  =  -x(46) + x(6)*x44 - x45*x52 + x46*x54 - x47*x56 + x48*x59
        x61  =  1._wp/x60
        x62  =  x61*(-x1*x53 + x14*x54 - x27*x56 + x43*x58 + x51)
        x63  =  x(49)*x0
        x64  =  x3*(x(2)*x63 - x(50))
        x65  =  x13*(-x(49)*x9 + x(51) + x10*x64)
        x66  =  x26*(x(49)*x20 - x(52) - x21*x64 + x23*x65)
        x67  =  x42*(-x(49)*x34 + x(53) + x35*x64 - x37*x65 + x39*x66)
        x68  =  -x(54) + x(6)*x63 - x52*x64 + x54*x65 - x56*x66 + x58*x67
        x69  =  x50*x61
        x70  =  x19*x66 - x33*x67 - x6*x64 - x63 - x65*x8 + x68*x69
        x71  =  x(7)*x0
        x72  =  x(15) - x(9)*x71
        x73  =  x3*x72
        x74  =  -x(23) + x(7)*x4 + x5*x73
        x75  =  x13*x74
        x76  =  -x(31) + x(7)*x15 - x16*x72 + x17*x75
        x77  =  x26*x76
        x78  =  x(39) - x(7)*x28 + x29*x72 - x30*x74 + x31*x77
        x79  =  x42*x78
        x80  =  -x(47) + x(7)*x44 - x45*x72 + x46*x74 - x47*x76 + x48*x79
        x81  =  x61*x80
        x82  =  x(55) - x(7)*x63 + x64*x72 - x65*x74 + x66*x76 - x67*x78 + x68*x81
        x83  =  1._wp/x82
        x84  =  x83*(x1*x73 - x14*x74 + x27*x76 - x43*x78 + x62*x80 - x71)
        x85  =  x(57)*x0
        x86  =  x3*(x(2)*x85 - x(58))
        x87  =  x13*(-x(57)*x9 + x(59) + x10*x86)
        x88  =  x26*(x(57)*x20 - x(60) - x21*x86 + x23*x87)
        x89  =  x42*(-x(57)*x34 + x(61) + x35*x86 - x37*x87 + x39*x88)
        x90  =  x61*(x(57)*x51 - x(62) - x52*x86 + x54*x87 - x56*x88 + x58*x89)
        x91  =  -x(57)*x71 + x(63) + x72*x86 - x74*x87 + x76*x88 - x78*x89 + x80*x90
        x92  =  x70*x83
        x93  =  -x19*x88 + x33*x89 - x50*x90 + x6*x86 + x8*x87 + x85 + x91*x92
        x94  =  x(8)*x0
        x95  =  x(16) - x(9)*x94
        x96  =  x3*x95
        x97  =  -x(24) + x(8)*x4 + x5*x96
        x98  =  x13*x97
        x99  =  -x(32) + x(8)*x15 - x16*x95 + x17*x98
        x100  =  x26*x99
        x101  =  x(40) - x(8)*x28 + x100*x31 + x29*x95 - x30*x97
        x102  =  x101*x42
        x103  =  -x(48) + x(8)*x44 + x102*x48 - x45*x95 + x46*x97 - x47*x99
        x104  =  x103*x61
        x105  =  x(56) - x(8)*x63 - x101*x67 + x104*x68 + x64*x95 - x65*x97 + x66*x99
        x106  =  x105*x83
        x107  =  1._wp/(-x(64) + x(8)*x85 + x101*x89 - x103*x90 + x106*x91 - x86*x95 + x87* &
          x97 - x88*x99)
        x108  =  x107*(-x1*x96 + x101*x43 - x103*x62 + x105*x84 + x14*x97 - x27*x99 + x94 &
          )
        x109  =  x3*x5
        x110  =  x109*x13
        x111  =  x110*x17 - x16
        x112  =  x111*x26
        x113  =  -x109*x30 + x112*x31 + x29
        x114  =  x113*x42
        x115  =  x109*x46 - x111*x47 + x114*x48 - x45
        x116  =  x115*x61
        x117  =  -x109*x65 + x111*x66 - x113*x67 + x116*x68 + x64
        x118  =  x117*x83
        x119  =  x109*x87 - x111*x88 + x113*x89 - x115*x90 + x118*x91 - x86
        x120  =  x13*x17
        x121  =  x120*x26
        x122  =  x121*x31 - x30
        x123  =  x122*x42
        x124  =  -x120*x47 + x123*x48 + x46
        x125  =  x124*x61
        x126  =  x120*x66 - x122*x67 + x125*x68 - x65
        x127  =  x126*x83
        x128  =  -x120*x88 + x122*x89 - x124*x90 + x127*x91 + x87
        x129  =  x26*x31
        x130  =  x129*x42
        x131  =  x130*x48 - x47
        x132  =  x131*x61
        x133  =  -x129*x67 + x132*x68 + x66
        x134  =  x133*x83
        x135  =  x129*x89 - x131*x90 + x134*x91 - x88
        x136  =  x42*x48
        x137  =  x136*x61
        x138  =  x137*x68 - x67
        x139  =  x138*x83
        x140  =  -x136*x90 + x139*x91 + x89
        x141  =  x61*x68
        x142  =  x141*x83
        x143  =  x142*x91 - x90
        x144  =  x83*x91
        x145  =  x26*(-x11*x24 + x22)
        x146  =  x42*(x11*x38 + x145*x39 - x36)
        x147  =  x61*(-x11*x55 - x145*x56 + x146*x58 + x53)
        x148  =  x83*(x11*x75 + x145*x76 - x146*x78 + x147*x80 - x73)
        x149  =  x107*(x101*x146 - x103*x147 + x105*x148 - x11*x98 - x145*x99 + x96)
        x150  =  x42*(-x24*x40 + x38)
        x151  =  x61*(x150*x58 + x24*x57 - x55)
        x152  =  x83*(-x150*x78 + x151*x80 - x24*x77 + x75)
        x153  =  x107*(x100*x24 + x101*x150 - x103*x151 + x105*x152 - x98)
        x154  =  x61*(-x40*x59 + x57)
        x155  =  x83*(x154*x80 + x40*x79 - x77)
        x156  =  x107*(x100 - x102*x40 - x103*x154 + x105*x155)
        x157  =  x83*(-x59*x81 + x79)
        x158  =  x107*(-x102 + x104*x59 + x105*x157)
        x159  =  x107*(x104 - x106*x81)
        x160  =  x107*x93
        x161  =  x107*x119
        x162  =  x107*x128
        x163  =  x107*x135
        x164  =  x107*x140
        x165  =  x107*x143

        inter(1) =  x0 - x108*x93 + x14*x8 - x19*x27 + x33*x43 - x50*x62 + x70*x84 + x(2)*x(9)*x3/x(1)**2
        inter(2) =  -x1*x3 - x108*x119 + x109*x14 - x111*x27 + x113*x43 - x115*x62 + x117*x84
        inter(3) =  x108*x128 + x120*x27 - x122*x43 + x124*x62 - x126*x84 - x14
        inter(4) =  x108*x135 - x129*x43 + x131*x62 - x133*x84 + x27
        inter(5) =  -x108*x140 - x136*x62 + x138*x84 + x43
        inter(6) =  x108*x143 - x141*x84 + x62
        inter(7) =  -x108*x144 + x84
        inter(8) =  x108
        inter(9) =  -x11*x18 - x145*x19 + x146*x33 - x147*x50 + x148*x70 - x149*x93 - x7
        inter(10) =  -x10*x13*x5/x2**2 - x111*x145 + x113*x146 - x115*x147 + x117*x148 - x119*x149 + x3
        inter(11) =  x11*x13 + x120*x145 - x122*x146 + x124*x147 - x126*x148 + x128*x149
        inter(12) =  -x129*x146 + x131*x147 - x133*x148 + x135*x149 + x145
        inter(13) =  -x136*x147 + x138*x148 - x140*x149 + x146
        inter(14) =  -x141*x148 + x143*x149 + x147
        inter(15) =  -x144*x149 + x148
        inter(16) =  x149
        inter(17) =  -x150*x33 + x151*x50 - x152*x70 + x153*x93 + x18 - x24*x32
        inter(18) =  x110 - x112*x24 - x113*x150 + x115*x151 - x117*x152 + x119*x153
        inter(19) =  x122*x150 - x124*x151 + x126*x152 - x128*x153 - x13 + x17*x23*x26/x12**2
        inter(20) =  x129*x150 - x131*x151 + x133*x152 - x135*x153 + x24*x26
        inter(21) =  x136*x151 - x138*x152 + x140*x153 - x150
        inter(22) =  x141*x152 - x143*x153 - x151
        inter(23) =  x144*x153 - x152
        inter(24) =  -x153
        inter(25) =  -x154*x50 + x155*x70 - x156*x93 + x32 - x40*x49
        inter(26) =  x112 - x114*x40 - x115*x154 + x117*x155 - x119*x156
        inter(27) =  -x121 + x123*x40 + x124*x154 - x126*x155 + x128*x156
        inter(28) =  x131*x154 - x133*x155 + x135*x156 - x26 + x31*x39*x42/x25**2
        inter(29) =  -x136*x154 + x138*x155 - x140*x156 - x40*x42
        inter(30) =  -x141*x155 + x143*x156 + x154
        inter(31) =  -x144*x156 + x155
        inter(32) =  x156
        inter(33) =  -x157*x70 + x158*x93 + x49 - x59*x69
        inter(34) =  x114 - x116*x59 - x117*x157 + x119*x158
        inter(35) =  -x123 + x125*x59 + x126*x157 - x128*x158
        inter(36) =  -x130 + x132*x59 + x133*x157 - x135*x158
        inter(37) =  -x138*x157 + x140*x158 + x42 - x48*x58*x61/x41**2
        inter(38) =  x141*x157 - x143*x158 + x59*x61
        inter(39) =  x144*x158 - x157
        inter(40) =  -x158
        inter(41) =  -x159*x93 + x69 - x81*x92
        inter(42) =  x116 - x118*x81 - x119*x159
        inter(43) =  -x125 + x127*x81 + x128*x159
        inter(44) =  -x132 + x134*x81 + x135*x159
        inter(45) =  x137 - x139*x81 - x140*x159
        inter(46) =  x143*x159 - x61 + x68*x80*x83/x60**2
        inter(47) =  -x144*x159 - x81*x83
        inter(48) =  x159
        inter(49) =  -x106*x160 + x92
        inter(50) =  -x106*x161 + x118
        inter(51) =  x106*x162 - x127
        inter(52) =  x106*x163 - x134
        inter(53) =  -x106*x164 + x139
        inter(54) =  x106*x165 - x142
        inter(55) =  -x105*x107*x91/x82**2 + x83
        inter(56) =  x106*x107
        inter(57) =  x160
        inter(58) =  x161
        inter(59) =  -x162
        inter(60) =  -x163
        inter(61) =  x164
        inter(62) =  -x165
        inter(63) =  x107*x144
        inter(64) =  -x107
        res = reshape(inter,[8,8])
    end function q_stminv8
    pure function d_stminv8(m) result(res)
        implicit none
        real(dp), intent(in) :: m(:,:)
        real(dp)             :: x(64),&
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
                              & x120, x121, x122, x123, x124, x125, x126, x127, & 
                              & x128, x129, x130, x131, x132, x133, x134, x135, & 
                              & x136, x137, x138, x139, x140, x141, x142, x143, & 
                              & x144, x145, x146, x147, x148, x149, x150, x151, & 
                              & x152, x153, x154, x155, x156, x157, x158, x159, & 
                              & x160, x161, x162, x163, x164, x165
        real(dp), dimension(64) :: inter
        real(dp), dimension(8,8) :: res
        x = reshape(m,[64])
        x0  =  1._dp/x(1)
        x1  =  x(2)*x0
        x2  =  x(10) - x(9)*x1
        x3  =  1._dp/x2
        x4  =  x(17)*x0
        x5  =  x(18) - x(2)*x4
        x6  =  x(9)*x0
        x7  =  x3*x6
        x8  =  x4 - x5*x7
        x9  =  x(3)*x0
        x10  =  x(11) - x(9)*x9
        x11  =  x10*x3
        x12  =  -x(19) + x(3)*x4 + x11*x5
        x13  =  1._dp/x12
        x14  =  x13*(x1*x11 - x9)
        x15  =  x(25)*x0
        x16  =  x3*(x(2)*x15 - x(26))
        x17  =  -x(25)*x9 + x(27) + x10*x16
        x18  =  x13*x8
        x19  =  x15 + x16*x6 + x17*x18
        x20  =  x(4)*x0
        x21  =  x(12) - x(9)*x20
        x22  =  x21*x3
        x23  =  -x(20) + x(4)*x4 + x22*x5
        x24  =  x13*x23
        x25  =  -x(28) + x(4)*x15 - x16*x21 + x17*x24
        x26  =  1._dp/x25
        x27  =  x26*(-x1*x22 + x14*x23 + x20)
        x28  =  x(33)*x0
        x29  =  x3*(x(2)*x28 - x(34))
        x30  =  x13*(-x(33)*x9 + x(35) + x10*x29)
        x31  =  -x(36) + x(4)*x28 - x21*x29 + x23*x30
        x32  =  x19*x26
        x33  =  -x28 - x29*x6 - x30*x8 + x31*x32
        x34  =  x(5)*x0
        x35  =  x(13) - x(9)*x34
        x36  =  x3*x35
        x37  =  -x(21) + x(5)*x4 + x36*x5
        x38  =  x13*x37
        x39  =  -x(29) + x(5)*x15 - x16*x35 + x17*x38
        x40  =  x26*x39
        x41  =  x(37) - x(5)*x28 + x29*x35 - x30*x37 + x31*x40
        x42  =  1._dp/x41
        x43  =  x42*(x1*x36 - x14*x37 + x27*x39 - x34)
        x44  =  x(41)*x0
        x45  =  x3*(x(2)*x44 - x(42))
        x46  =  x13*(-x(41)*x9 + x(43) + x10*x45)
        x47  =  x26*(x(41)*x20 - x(44) - x21*x45 + x23*x46)
        x48  =  -x(41)*x34 + x(45) + x35*x45 - x37*x46 + x39*x47
        x49  =  x33*x42
        x50  =  -x19*x47 + x44 + x45*x6 + x46*x8 + x48*x49
        x51  =  x(6)*x0
        x52  =  x(14) - x(9)*x51
        x53  =  x3*x52
        x54  =  -x(22) + x(6)*x4 + x5*x53
        x55  =  x13*x54
        x56  =  -x(30) + x(6)*x15 - x16*x52 + x17*x55
        x57  =  x26*x56
        x58  =  x(38) - x(6)*x28 + x29*x52 - x30*x54 + x31*x57
        x59  =  x42*x58
        x60  =  -x(46) + x(6)*x44 - x45*x52 + x46*x54 - x47*x56 + x48*x59
        x61  =  1._dp/x60
        x62  =  x61*(-x1*x53 + x14*x54 - x27*x56 + x43*x58 + x51)
        x63  =  x(49)*x0
        x64  =  x3*(x(2)*x63 - x(50))
        x65  =  x13*(-x(49)*x9 + x(51) + x10*x64)
        x66  =  x26*(x(49)*x20 - x(52) - x21*x64 + x23*x65)
        x67  =  x42*(-x(49)*x34 + x(53) + x35*x64 - x37*x65 + x39*x66)
        x68  =  -x(54) + x(6)*x63 - x52*x64 + x54*x65 - x56*x66 + x58*x67
        x69  =  x50*x61
        x70  =  x19*x66 - x33*x67 - x6*x64 - x63 - x65*x8 + x68*x69
        x71  =  x(7)*x0
        x72  =  x(15) - x(9)*x71
        x73  =  x3*x72
        x74  =  -x(23) + x(7)*x4 + x5*x73
        x75  =  x13*x74
        x76  =  -x(31) + x(7)*x15 - x16*x72 + x17*x75
        x77  =  x26*x76
        x78  =  x(39) - x(7)*x28 + x29*x72 - x30*x74 + x31*x77
        x79  =  x42*x78
        x80  =  -x(47) + x(7)*x44 - x45*x72 + x46*x74 - x47*x76 + x48*x79
        x81  =  x61*x80
        x82  =  x(55) - x(7)*x63 + x64*x72 - x65*x74 + x66*x76 - x67*x78 + x68*x81
        x83  =  1._dp/x82
        x84  =  x83*(x1*x73 - x14*x74 + x27*x76 - x43*x78 + x62*x80 - x71)
        x85  =  x(57)*x0
        x86  =  x3*(x(2)*x85 - x(58))
        x87  =  x13*(-x(57)*x9 + x(59) + x10*x86)
        x88  =  x26*(x(57)*x20 - x(60) - x21*x86 + x23*x87)
        x89  =  x42*(-x(57)*x34 + x(61) + x35*x86 - x37*x87 + x39*x88)
        x90  =  x61*(x(57)*x51 - x(62) - x52*x86 + x54*x87 - x56*x88 + x58*x89)
        x91  =  -x(57)*x71 + x(63) + x72*x86 - x74*x87 + x76*x88 - x78*x89 + x80*x90
        x92  =  x70*x83
        x93  =  -x19*x88 + x33*x89 - x50*x90 + x6*x86 + x8*x87 + x85 + x91*x92
        x94  =  x(8)*x0
        x95  =  x(16) - x(9)*x94
        x96  =  x3*x95
        x97  =  -x(24) + x(8)*x4 + x5*x96
        x98  =  x13*x97
        x99  =  -x(32) + x(8)*x15 - x16*x95 + x17*x98
        x100  =  x26*x99
        x101  =  x(40) - x(8)*x28 + x100*x31 + x29*x95 - x30*x97
        x102  =  x101*x42
        x103  =  -x(48) + x(8)*x44 + x102*x48 - x45*x95 + x46*x97 - x47*x99
        x104  =  x103*x61
        x105  =  x(56) - x(8)*x63 - x101*x67 + x104*x68 + x64*x95 - x65*x97 + x66*x99
        x106  =  x105*x83
        x107  =  1._dp/(-x(64) + x(8)*x85 + x101*x89 - x103*x90 + x106*x91 - x86*x95 + x87* &
          x97 - x88*x99)
        x108  =  x107*(-x1*x96 + x101*x43 - x103*x62 + x105*x84 + x14*x97 - x27*x99 + x94 &
          )
        x109  =  x3*x5
        x110  =  x109*x13
        x111  =  x110*x17 - x16
        x112  =  x111*x26
        x113  =  -x109*x30 + x112*x31 + x29
        x114  =  x113*x42
        x115  =  x109*x46 - x111*x47 + x114*x48 - x45
        x116  =  x115*x61
        x117  =  -x109*x65 + x111*x66 - x113*x67 + x116*x68 + x64
        x118  =  x117*x83
        x119  =  x109*x87 - x111*x88 + x113*x89 - x115*x90 + x118*x91 - x86
        x120  =  x13*x17
        x121  =  x120*x26
        x122  =  x121*x31 - x30
        x123  =  x122*x42
        x124  =  -x120*x47 + x123*x48 + x46
        x125  =  x124*x61
        x126  =  x120*x66 - x122*x67 + x125*x68 - x65
        x127  =  x126*x83
        x128  =  -x120*x88 + x122*x89 - x124*x90 + x127*x91 + x87
        x129  =  x26*x31
        x130  =  x129*x42
        x131  =  x130*x48 - x47
        x132  =  x131*x61
        x133  =  -x129*x67 + x132*x68 + x66
        x134  =  x133*x83
        x135  =  x129*x89 - x131*x90 + x134*x91 - x88
        x136  =  x42*x48
        x137  =  x136*x61
        x138  =  x137*x68 - x67
        x139  =  x138*x83
        x140  =  -x136*x90 + x139*x91 + x89
        x141  =  x61*x68
        x142  =  x141*x83
        x143  =  x142*x91 - x90
        x144  =  x83*x91
        x145  =  x26*(-x11*x24 + x22)
        x146  =  x42*(x11*x38 + x145*x39 - x36)
        x147  =  x61*(-x11*x55 - x145*x56 + x146*x58 + x53)
        x148  =  x83*(x11*x75 + x145*x76 - x146*x78 + x147*x80 - x73)
        x149  =  x107*(x101*x146 - x103*x147 + x105*x148 - x11*x98 - x145*x99 + x96)
        x150  =  x42*(-x24*x40 + x38)
        x151  =  x61*(x150*x58 + x24*x57 - x55)
        x152  =  x83*(-x150*x78 + x151*x80 - x24*x77 + x75)
        x153  =  x107*(x100*x24 + x101*x150 - x103*x151 + x105*x152 - x98)
        x154  =  x61*(-x40*x59 + x57)
        x155  =  x83*(x154*x80 + x40*x79 - x77)
        x156  =  x107*(x100 - x102*x40 - x103*x154 + x105*x155)
        x157  =  x83*(-x59*x81 + x79)
        x158  =  x107*(-x102 + x104*x59 + x105*x157)
        x159  =  x107*(x104 - x106*x81)
        x160  =  x107*x93
        x161  =  x107*x119
        x162  =  x107*x128
        x163  =  x107*x135
        x164  =  x107*x140
        x165  =  x107*x143

        inter(1) =  x0 - x108*x93 + x14*x8 - x19*x27 + x33*x43 - x50*x62 + x70*x84 + x(2)*x(9)*x3/x(1)**2
        inter(2) =  -x1*x3 - x108*x119 + x109*x14 - x111*x27 + x113*x43 - x115*x62 + x117*x84
        inter(3) =  x108*x128 + x120*x27 - x122*x43 + x124*x62 - x126*x84 - x14
        inter(4) =  x108*x135 - x129*x43 + x131*x62 - x133*x84 + x27
        inter(5) =  -x108*x140 - x136*x62 + x138*x84 + x43
        inter(6) =  x108*x143 - x141*x84 + x62
        inter(7) =  -x108*x144 + x84
        inter(8) =  x108
        inter(9) =  -x11*x18 - x145*x19 + x146*x33 - x147*x50 + x148*x70 - x149*x93 - x7
        inter(10) =  -x10*x13*x5/x2**2 - x111*x145 + x113*x146 - x115*x147 + x117*x148 - x119*x149 + x3
        inter(11) =  x11*x13 + x120*x145 - x122*x146 + x124*x147 - x126*x148 + x128*x149
        inter(12) =  -x129*x146 + x131*x147 - x133*x148 + x135*x149 + x145
        inter(13) =  -x136*x147 + x138*x148 - x140*x149 + x146
        inter(14) =  -x141*x148 + x143*x149 + x147
        inter(15) =  -x144*x149 + x148
        inter(16) =  x149
        inter(17) =  -x150*x33 + x151*x50 - x152*x70 + x153*x93 + x18 - x24*x32
        inter(18) =  x110 - x112*x24 - x113*x150 + x115*x151 - x117*x152 + x119*x153
        inter(19) =  x122*x150 - x124*x151 + x126*x152 - x128*x153 - x13 + x17*x23*x26/x12**2
        inter(20) =  x129*x150 - x131*x151 + x133*x152 - x135*x153 + x24*x26
        inter(21) =  x136*x151 - x138*x152 + x140*x153 - x150
        inter(22) =  x141*x152 - x143*x153 - x151
        inter(23) =  x144*x153 - x152
        inter(24) =  -x153
        inter(25) =  -x154*x50 + x155*x70 - x156*x93 + x32 - x40*x49
        inter(26) =  x112 - x114*x40 - x115*x154 + x117*x155 - x119*x156
        inter(27) =  -x121 + x123*x40 + x124*x154 - x126*x155 + x128*x156
        inter(28) =  x131*x154 - x133*x155 + x135*x156 - x26 + x31*x39*x42/x25**2
        inter(29) =  -x136*x154 + x138*x155 - x140*x156 - x40*x42
        inter(30) =  -x141*x155 + x143*x156 + x154
        inter(31) =  -x144*x156 + x155
        inter(32) =  x156
        inter(33) =  -x157*x70 + x158*x93 + x49 - x59*x69
        inter(34) =  x114 - x116*x59 - x117*x157 + x119*x158
        inter(35) =  -x123 + x125*x59 + x126*x157 - x128*x158
        inter(36) =  -x130 + x132*x59 + x133*x157 - x135*x158
        inter(37) =  -x138*x157 + x140*x158 + x42 - x48*x58*x61/x41**2
        inter(38) =  x141*x157 - x143*x158 + x59*x61
        inter(39) =  x144*x158 - x157
        inter(40) =  -x158
        inter(41) =  -x159*x93 + x69 - x81*x92
        inter(42) =  x116 - x118*x81 - x119*x159
        inter(43) =  -x125 + x127*x81 + x128*x159
        inter(44) =  -x132 + x134*x81 + x135*x159
        inter(45) =  x137 - x139*x81 - x140*x159
        inter(46) =  x143*x159 - x61 + x68*x80*x83/x60**2
        inter(47) =  -x144*x159 - x81*x83
        inter(48) =  x159
        inter(49) =  -x106*x160 + x92
        inter(50) =  -x106*x161 + x118
        inter(51) =  x106*x162 - x127
        inter(52) =  x106*x163 - x134
        inter(53) =  -x106*x164 + x139
        inter(54) =  x106*x165 - x142
        inter(55) =  -x105*x107*x91/x82**2 + x83
        inter(56) =  x106*x107
        inter(57) =  x160
        inter(58) =  x161
        inter(59) =  -x162
        inter(60) =  -x163
        inter(61) =  x164
        inter(62) =  -x165
        inter(63) =  x107*x144
        inter(64) =  -x107
        res = reshape(inter,[8,8])
    end function d_stminv8
    pure subroutine sttchain_loop(istm,istt,stm,stt,stmab,sttab,n)
        implicit none
        integer i, a, b
        integer, intent(in) :: n
        real(wp), intent(in)  :: istm(n,n), stm(n,n), istt(n,n,n), stt(n,n,n)
        real(wp), intent(out) :: stmab(n,n), sttab(n,n,n)
        do a=1,n
        do i=1,n
            stmab(i,a) = dot_product(stm(i,:),istm(:,a))
        do b=1,n
            sttab(b,i,a) = dot_product(stm(i,:),istt(:,a,b)) &
                + dot_product(mmult(stt(i,:,:),istm(:,b)),istm(:,a))
        enddo; enddo; enddo
        sttab = reshape(sttab,shape=[n,n,n],order=[2,3,1])
    end subroutine sttchain_loop
end module tensorops
