
! BEGIN AUTOCODE OUTPUT FOR STMINV
function stminv(m) result(res)
    implicit none
    real(8), intent(in) :: m(8,8)
    real(8)             :: x(64)
    real(8)             :: &
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
                          & x112, x113, x114, x115, x116
    real(8), dimension(64) :: inter
    real(8), dimension(8,8) :: res

    x = reshape(m,[64])
    x0  =  1_8/x(1)
    x1  =  x(2)*x0
    x2  =  x(10) - x(9)*x1
    x3  =  1_8/x2
    x4  =  x(17)*x0
    x5  =  x(18) - x(2)*x4
    x6  =  x(9)*x0
    x7  =  x3*x6
    x8  =  x4 - x5*x7
    x9  =  x(3)*x0
    x10  =  x(11) - x(9)*x9
    x11  =  x10*x3
    x12  =  -x(19) + x(3)*x4 + x11*x5
    x13  =  1_8/x12
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
    x26  =  1_8/x25
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
    x42  =  1_8/x41
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
    x60  =  1_8/(-x(46) + x(6)*x44 - x45*x52 + x46*x54 - x47*x56 + x48*x59)
    x61  =  x60*(-x1*x53 + x14*x54 - x27*x56 + x43*x58 + x51)
    x62  =  x3*x5
    x63  =  x13*x62
    x64  =  -x16 + x17*x63
    x65  =  x26*x64
    x66  =  x29 - x30*x62 + x31*x65
    x67  =  x42*x66
    x68  =  -x45 + x46*x62 - x47*x64 + x48*x67
    x69  =  x13*x17
    x70  =  x26*x69
    x71  =  -x30 + x31*x70
    x72  =  x42*x71
    x73  =  x46 - x47*x69 + x48*x72
    x74  =  x26*x31
    x75  =  x42*x74
    x76  =  -x47 + x48*x75
    x77  =  x42*x48
    x78  =  x(7)*x0
    x79  =  x(15) - x(9)*x78
    x80  =  x3*x79
    x81  =  -x(23) + x(7)*x4 + x5*x80
    x82  =  x13*x81
    x83  =  -x(31) + x(7)*x15 - x16*x79 + x17*x82
    x84  =  x26*x83
    x85  =  x(39) - x(7)*x28 + x29*x79 - x30*x81 + x31*x84
    x86  =  x42*x85
    x87  =  -x(47) + x(7)*x44 - x45*x79 + x46*x81 - x47*x83 + x48*x86
    x88  =  x1*x80 - x14*x81 + x27*x83 - x43*x85 + x61*x87 - x78
    x89  =  x(8)*x0
    x90  =  x(16) - x(9)*x89
    x91  =  x3*x90
    x92  =  -x(24) + x(8)*x4 + x5*x91
    x93  =  x13*x92
    x94  =  -x(32) + x(8)*x15 - x16*x90 + x17*x93
    x95  =  x26*x94
    x96  =  x(40) - x(8)*x28 + x29*x90 - x30*x92 + x31*x95
    x97  =  x42*x96
    x98  =  -x(48) + x(8)*x44 - x45*x90 + x46*x92 - x47*x94 + x48*x97
    x99  =  x26*(-x11*x24 + x22)
    x100  =  x42*(x11*x38 - x36 + x39*x99)
    x101  =  x60*(x100*x58 - x11*x55 + x53 - x56*x99)
    x102  =  -x100*x85 + x101*x87 + x11*x82 - x80 + x83*x99
    x103  =  x42*(-x24*x40 + x38)
    x104  =  x60*(x103*x58 + x24*x57 - x55)
    x105  =  x24*x84
    x106  =  x103*x85
    x107  =  x104*x87
    x108  =  x60*(-x40*x59 + x57)
    x109  =  x108*x87 + x40*x86 - x84
    x110  =  x50*x60
    x111  =  x60*x68
    x112  =  x60*x73
    x113  =  x60*x76
    x114  =  x60*x87
    x115  =  x114*x59
    x116  =  x60*x98

    inter(1) =  x0 + x14*x8 - x19*x27 + x33*x43 - x50*x61 + x(2)*x(9)*x3/x(1)**2
    inter(2) =  -x1*x3 + x14*x62 - x27*x64 + x43*x66 - x61*x68
    inter(3) =  -x14 + x27*x69 - x43*x71 + x61*x73
    inter(4) =  x27 - x43*x74 + x61*x76
    inter(5) =  x43 - x61*x77
    inter(6) =  x61
    inter(7) =  x88
    inter(8) =  -x(56)*x88 + x1*x91 - x14*x92 + x27*x94 - x43*x96 + x61*x98 - x89
    inter(9) =  x100*x33 - x101*x50 - x11*x18 - x19*x99 - x7
    inter(10) =  -x10*x13*x5/x2**2 + x100*x66 - x101*x68 + x3 - x64*x99
    inter(11) =  -x100*x71 + x101*x73 + x11*x13 + x69*x99
    inter(12) =  -x100*x74 + x101*x76 + x99
    inter(13) =  x100 - x101*x77
    inter(14) =  x101
    inter(15) =  x102
    inter(16) =  -x(56)*x102 - x100*x96 + x101*x98 + x11*x93 - x91 + x94*x99
    inter(17) =  -x103*x33 + x104*x50 + x18 - x24*x32
    inter(18) =  -x103*x66 + x104*x68 - x24*x65 + x63
    inter(19) =  x103*x71 - x104*x73 - x13 + x17*x23*x26/x12**2
    inter(20) =  x103*x74 - x104*x76 + x24*x26
    inter(21) =  -x103 + x104*x77
    inter(22) =  -x104
    inter(23) =  x105 + x106 - x107 - x82
    inter(24) =  x(56)*(-x105 - x106 + x107 + x82) + x103*x96 - x104*x98 + x24*x95 - x93
    inter(25) =  -x108*x50 + x32 - x40*x49
    inter(26) =  -x108*x68 - x40*x67 + x65
    inter(27) =  x108*x73 + x40*x72 - x70
    inter(28) =  x108*x76 - x26 + x31*x39*x42/x25**2
    inter(29) =  -x108*x77 - x40*x42
    inter(30) =  x108
    inter(31) =  x109
    inter(32) =  -x(56)*x109 + x108*x98 + x40*x97 - x95
    inter(33) =  -x110*x59 + x49
    inter(34) =  -x111*x59 + x67
    inter(35) =  x112*x59 - x72
    inter(36) =  x113*x59 - x75
    inter(37) =  x42 - x48*x58*x60/x41**2
    inter(38) =  x59*x60
    inter(39) =  x115 - x86
    inter(40) =  x(56)*(-x115 + x86) + x116*x59 - x97
    inter(41) =  x110
    inter(42) =  x111
    inter(43) =  -x112
    inter(44) =  -x113
    inter(45) =  x60*x77
    inter(46) =  -x60
    inter(47) =  -x114
    inter(48) =  x(56)*x114 - x116
    inter(49) =  0
    inter(50) =  0
    inter(51) =  0
    inter(52) =  0
    inter(53) =  0
    inter(54) =  0
    inter(55) =  1
    inter(56) =  -x(56)
    inter(57) =  0
    inter(58) =  0
    inter(59) =  0
    inter(60) =  0
    inter(61) =  0
    inter(62) =  0
    inter(63) =  0
    inter(64) =  1
    res = reshape(inter,[8,8])
end function stminv
! END AUTOCODE OUTPUT FOR STMINV
