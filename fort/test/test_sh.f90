module test_sh
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use tinysh, only: PinesData, shpines, pinesinit
    use test_util, only: mprint, tprint
    implicit none
    contains
        subroutine sh_test(testpass)
            logical, intent(inout) :: testpass
            real(dp)               :: Vtrue, dVtrue(3), d2Vtrue(3,3), d3Vtrue(3,3,3)
            real(dp),allocatable   :: Clm(:,:), Slm(:,:), &
                                      Cml(:,:), Sml(:,:)
            type(PinesData) :: pines
            real(dp) :: V, dV(3), d2V(3,3), d3V(3,3,3), Rbody, GM, cart(3)
            real(dp), parameter :: dtol=1.D-22
            Rbody = 6378.137_dp
            GM = 3.986D+05
            cart = Rbody

            Vtrue          = 36.081331327754611_dp
            dVtrue         = [-1.8849858008825431D-003, -1.8850101627501258D-003, -1.8870305095350294D-003]
            d2Vtrue(1,:)   = [-2.5341422557527843D-010,   2.9529211670113480D-007,   2.9582059131367937D-007]
            d2Vtrue(2,:)   = [ 2.9529211670113480D-007,  -2.4190181074974302D-010,   2.9582632428968837D-007]
            d2Vtrue(3,:)   = [ 2.9582059131367937D-007,   2.9582632428968837D-007,   4.9531603632519351D-010]
            d3Vtrue(1,1,:) = [ 6.1833060939829652D-011,  -3.0761459311059287D-011,  -3.0872344351507036D-011]
            d3Vtrue(1,2,:) = [-3.0761459311059287D-011,  -3.0764946515452407D-011,  -7.7255090267888913D-011]
            d3Vtrue(1,3,:) = [-3.0872344351507036D-011,  -7.7255090267888913D-011,  -3.1068114424377219D-011]
            d3Vtrue(2,1,:) = [-3.0761459311059287D-011,  -3.0764946515452407D-011,  -7.7255090267888913D-011]
            d3Vtrue(2,2,:) = [-3.0764946515452407D-011,   6.1830340360358116D-011,  -3.0876347278308212D-011]
            d3Vtrue(2,3,:) = [-7.7255090267888913D-011,  -3.0876347278308212D-011,  -3.1068881049298791D-011]
            d3Vtrue(3,1,:) = [-3.0872344351507036D-011,  -7.7255090267888913D-011,  -3.1068114424377219D-011]
            d3Vtrue(3,2,:) = [-7.7255090267888913D-011,  -3.0876347278308212D-011,  -3.1068881049298791D-011]
            d3Vtrue(3,3,:) = [-3.1068114424377219D-011,  -3.1068881049298791D-011,   6.1748691629815287D-011]


            allocate(Clm(0:6,0:6), Slm(0:6,0:6), &
                        Cml(7,7), Sml(7,7))
            Clm = 0._dp
            Slm = 0._dp
            Clm(0,0) = 1._dp
            Clm(2,0) = -4.841694573200D-04
            Clm(2,1) = -3.103431067239D-10
            Clm(2,2) =  2.439373415940D-06
            Clm(3,0) =  9.571647583412D-07
            Clm(3,1) =  2.030446637169D-06
            Clm(3,2) =  9.047646744100D-07
            Clm(3,3) =  7.212852551704D-07
            Clm(4,0) =  5.399815392137D-07
            Clm(4,1) = -5.361808133703D-07
            Clm(4,2) =  3.504921442703D-07
            Clm(4,3) =  9.908610311151D-07
            Clm(4,4) = -1.884924225276D-07
            Clm(5,0) =  6.865032345839D-08
            Clm(5,1) = -6.291457940968D-08
            Clm(5,2) =  6.520586031691D-07
            Clm(5,3) = -4.518313784464D-07
            Clm(5,4) = -2.953234091704D-07
            Clm(5,5) =  1.748143504694D-07
            Clm(6,1) = -7.594326587940D-08

            Slm(2,0) =  0.000000000000D+00
            Slm(2,1) =  1.410757509442D-09
            Slm(2,2) = -1.400294011836D-06
            Slm(3,0) =  0.000000000000D+00
            Slm(3,1) =  2.482406346848D-07
            Slm(3,2) = -6.190066246333D-07
            Slm(3,3) =  1.414400065165D-06
            Slm(4,0) =  0.000000000000D+00
            Slm(4,1) = -4.735769769691D-07
            Slm(4,2) =  6.625051657439D-07
            Slm(4,3) = -2.009508998058D-07
            Slm(4,4) =  3.088185785570D-07
            Slm(5,0) =  0.000000000000D+00
            Slm(5,1) = -9.434259860005D-08
            Slm(5,2) = -3.233430798143D-07
            Slm(5,3) = -2.149423673602D-07
            Slm(5,4) =  4.981057884405D-08
            Slm(5,5) = -6.693546770160D-07
            Slm(6,1) =  2.652568324970D-08
            Cml = transpose(Clm)
            Sml = transpose(Slm)
            call pinesinit(6,Cml,Sml,pines)
            call shpines(Rbody, GM, pines,6,3,cart, V, dV, d2V, d3V)
            testpass = .true.
            if (abs(V- Vtrue).ge.dtol) then
                testpass = .false.
                write(*,*) "FAIL Spherical harmonics potential FAIL"
                write(*,*) "Error: ", V - Vtrue
            endif
            if (any(abs(dV- dVtrue).ge.dtol)) then
                testpass = .false.
                write(*,*) "FAIL Spherical harmonics gradient FAIL"
                write(*,*) "Error: ", dV - dVtrue
            endif
            if (any(abs(d2V- d2Vtrue).ge.dtol)) then
                testpass = .false.
                write(*,*) "FAIL Spherical harmonics Jacobian FAIL"
                write(*,*) "Error: "
                call mprint(d2V - d2Vtrue)
            endif
            if (any(abs(d3V- d3Vtrue).ge.dtol)) then
                testpass = .false.
                write(*,*) "FAIL Spherical harmonics Hessian FAIL"
                write(*,*) "Error: "
                call tprint(d3V - d3Vtrue)
            endif
            if (testpass) write(*,*) "PASS All spherical harmonics tests PASS"
        end subroutine sh_test

end module test_sh
