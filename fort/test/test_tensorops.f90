module test_tensorops
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use tensorops, only: mattens, &
                         vectensquad, &
                         sttchain, &
                         sttchain_invert, &
                         stminvert, &
                         sttinvert, &
                         vectens1, &
                         vectens2, &
                         vectens3, &
                         quad, &
                         eyemat
    implicit none
    contains
        subroutine test_eyemat(testpass)
            character(len=100) :: msg
            real(dp), parameter :: dtol = 1.e-14
            real(qp), parameter :: qtol = 1.e-28
            real(dp) :: eye_d(3,3)
            real(qp) :: eye_q(3,3)
            logical, intent(inout) :: testpass

            eye_d(1,:) = [1._dp, 0._dp, 0._dp]
            eye_d(2,:) = [0._dp, 1._dp, 0._dp]
            eye_d(3,:) = [0._dp, 0._dp, 1._dp]

            eye_q(1,:) = [1._qp, 0._qp, 0._qp]
            eye_q(2,:) = [0._qp, 1._qp, 0._qp]
            eye_q(3,:) = [0._qp, 0._qp, 1._qp]
            msg = ""
            testpass = .true.
            if (norm2(eye_d-eyemat(3)).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//" Double precision identity matrix fail."
            endif
            if (norm2(eye_q-eyemat(3)).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//" Quad precision identity matrix fail."
            endif
            if (testpass) msg = " Identity matrix test pass."
            write (*,*) msg
        end subroutine
        subroutine test_basic_tensorops_d(testpass)
            character(len=100)     :: msg
            logical, intent(inout) :: testpass
            real(dp), parameter :: dtol = 1.e-13_dp
            real(dp) :: vec(3), mat(3,3), tens(3,3,3), &
                        stm6(6,6), stm8(8,8), mt_truth(3,3,3), &
                        quad_truth(3,3,3), vectens1_truth(3,3), &
                        vectens2_truth(3,3), vectens3_truth(3,3), &
                        vectensquad_truth(3), stminvert6_truth(6,6), &
                        stminvert8_truth(8,8), &
                        mt_test(3,3,3), &
                        quad_test(3,3,3), vectens1_test(3,3), &
                        vectens2_test(3,3), vectens3_test(3,3), &
                        vectensquad_test(3), stminvert6_test(6,6), &
                        stminvert8_test(8,8)
            ! vec
            vec =  [10._dp, 30._dp, 93._dp]

            ! mat
                mat(1,:) =  [ 5._dp,  6._dp, 14._dp]
                mat(2,:) =  [86._dp, 51._dp, 15._dp]
                mat(3,:) =  [34._dp, 25._dp, 71._dp]

            ! tens
                tens(1,1,:) =  [ 1._dp,  0._dp, 64._dp]
                tens(1,2,:) =  [21._dp, 80._dp, 32._dp]
                tens(1,3,:) =  [43._dp, 90._dp, 83._dp]
                tens(2,1,:) =  [41._dp,  4._dp, 52._dp]
                tens(2,2,:) =  [70._dp, 26._dp, 44._dp]
                tens(2,3,:) =  [41._dp, 95._dp, 30._dp]
                tens(3,1,:) =  [ 5._dp, 98._dp, 44._dp]
                tens(3,2,:) =  [17._dp, 68._dp, 57._dp]
                tens(3,3,:) =  [39._dp, 21._dp, 63._dp]

            !stm6
                stm6(1,:) =  [92._dp, 68._dp, 73._dp, 85._dp, 40._dp, 51._dp]
                stm6(2,:) =  [63._dp, 38._dp, 47._dp, 71._dp, 88._dp, 13._dp]
                stm6(3,:) =  [33._dp, 53._dp, 82._dp,  6._dp, 45._dp,  3._dp]
                stm6(4,:) =  [ 11._dp, 100._dp,  36._dp,  87._dp,  22._dp,  31._dp]
                stm6(5,:) =  [34._dp, 63._dp, 78._dp, 17._dp, 51._dp, 13._dp]
                stm6(6,:) =  [62._dp, 45._dp, 93._dp,  2._dp, 77._dp, 72._dp]

            !stm8
                stm8(1,:) =  [25._dp, 25._dp, 85._dp, 14._dp, 46._dp, 13._dp, 74._dp, 54._dp]
                stm8(2,:) =  [46._dp, 81._dp, 50._dp, 13._dp, 72._dp, 46._dp, 95._dp, 65._dp]
                stm8(3,:) =  [30._dp, 44._dp,  4._dp, 85._dp, 38._dp, 34._dp, 15._dp, 15._dp]
                stm8(4,:) =  [19._dp, 34._dp, 98._dp, 99._dp, 53._dp, 98._dp, 34._dp, 52._dp]
                stm8(5,:) =  [48._dp, 97._dp, 12._dp, 84._dp, 73._dp, 12._dp, 68._dp, 83._dp]
                stm8(6,:) =  [76._dp, 31._dp, 96._dp, 94._dp, 39._dp, 80._dp,  3._dp, 54._dp]
                stm8(7,:) =  [18._dp, 49._dp, 96._dp, 81._dp, 50._dp,  6._dp, 67._dp, 54._dp]
                stm8(8,:) =  [95._dp, 77._dp, 70._dp, 35._dp, 55._dp, 49._dp, 87._dp, 11._dp]

            ! mattens
                mt_truth(1,1,:) =  [ 321._dp, 1396._dp, 1248._dp]
                mt_truth(1,2,:) =  [ 763._dp, 1508._dp, 1222._dp]
                mt_truth(1,3,:) =  [1007._dp, 1314._dp, 1477._dp]
                mt_truth(2,1,:) =  [2252._dp, 1674._dp, 8816._dp]
                mt_truth(2,2,:) =  [5631._dp, 9226._dp, 5851._dp]
                mt_truth(2,3,:) =  [ 6374._dp, 12900._dp,  9613._dp]
                mt_truth(3,1,:) =  [1414._dp, 7058._dp, 6600._dp]
                mt_truth(3,2,:) =  [3671._dp, 8198._dp, 6235._dp]
                mt_truth(3,3,:) =  [5256._dp, 6926._dp, 8045._dp]

            ! quad
                quad_truth(1,1,:) =  [1071601._dp,  673928._dp,  613396._dp]
                quad_truth(1,2,:) =  [694234._dp, 438017._dp, 415539._dp]
                quad_truth(1,3,:) =  [916796._dp, 589107._dp, 677297._dp]
                quad_truth(2,1,:) =  [682067._dp, 452100._dp, 548500._dp]
                quad_truth(2,2,:) =  [456959._dp, 301671._dp, 358375._dp]
                quad_truth(2,3,:) =  [760713._dp, 481901._dp, 421133._dp]
                quad_truth(3,1,:) =  [867513._dp, 558130._dp, 660906._dp]
                quad_truth(3,2,:) =  [564690._dp, 363513._dp, 431889._dp]
                quad_truth(3,3,:) =  [551504._dp, 365197._dp, 523585._dp]

            ! vectens1
                vectens1_truth(1,:) =  [1705._dp, 9234._dp, 6292._dp]
                vectens1_truth(2,:) =  [3891._dp, 7904._dp, 6941._dp]
                vectens1_truth(3,:) =  [5287._dp, 5703._dp, 7589._dp]

            ! vectens2
                vectens2_truth(1,:) =  [ 4639._dp, 10770._dp,  9319._dp]
                vectens2_truth(2,:) =  [6323._dp, 9655._dp, 4630._dp]
                vectens2_truth(3,:) =  [4187._dp, 4973._dp, 8009._dp]

            ! vectens3
                vectens3_truth(1,:) =  [ 5962._dp,  5586._dp, 10849._dp]
                vectens3_truth(2,:) =  [5366._dp, 5572._dp, 6050._dp]
                vectens3_truth(3,:) =  [7082._dp, 7511._dp, 6879._dp]

            ! vectensquad
                vectensquad_truth =  [1236157._dp,  783470._dp,  935897._dp]

            ! STMinvert6
                stminvert6_truth(1,:) =  [ 87._dp,  17._dp,   2._dp, -85._dp, -71._dp,  -6._dp]
                stminvert6_truth(2,:) =  [ 22._dp,  51._dp,  77._dp, -40._dp, -88._dp, -45._dp]
                stminvert6_truth(3,:) =  [ 31._dp,  13._dp,  72._dp, -51._dp, -13._dp,  -3._dp]
                stminvert6_truth(4,:) =  [-11._dp, -34._dp, -62._dp,  92._dp,  63._dp,  33._dp]
                stminvert6_truth(5,:) =  [-100._dp,  -63._dp,  -45._dp,   68._dp,   38._dp,   53._dp]
                stminvert6_truth(6,:) =  [-36._dp, -78._dp, -93._dp,  73._dp,  47._dp,  82._dp]

            ! STMinvert8
                stminvert8_truth(1_dp,:) =  [ 0.01684607_dp, -0.02049215_dp, -0.00492154_dp,  0.00023436_dp,  0.01364028_dp,&
                0.00309208_dp, -0.0174303 _dp,  0.01145994_dp]
                stminvert8_truth(2_dp,:) =  [-0.05263277_dp,  0.03873808_dp, -0.01003245_dp, -0.01360167_dp, -0.01584197_dp,&
                0.01094427_dp,  0.03666172_dp, -0.00671567_dp]
                stminvert8_truth(3_dp,:) =  [-0.01830201_dp,  0.03185471_dp,  0.01231027_dp, -0.02074384_dp, -0.02911729_dp,&
                0.01512645_dp,  0.02883093_dp, -0.01319872_dp]
                stminvert8_truth(4_dp,:) =  [ 0.0121908 _dp, -0.03885463_dp, -0.01532077_dp,  0.02289227_dp,  0.02725751_dp,&
                -0.01398507_dp, -0.01438121_dp,  0.0160061 _dp]
                stminvert8_truth(5_dp,:) =  [ 0.00310332_dp,  0.11799515_dp,  0.1345172 _dp, -0.09211388_dp, -0.10954499_dp,&
                0.04305073_dp,  0.04573723_dp, -0.06976503_dp]
                stminvert8_truth(6_dp,:) =  [ 0.00150244_dp, -0.02015711_dp, -0.02413202_dp,  0.02898049_dp,  0.02012894_dp,&
                -0.01155453_dp, -0.02071762_dp,  0.0141881 _dp]
                stminvert8_truth(7_dp,:) =  [ 0.03374979_dp, -0.08094552_dp, -0.05577335_dp,  0.05802861_dp,  0.06445587_dp,&
                -0.03950086_dp, -0.04521933_dp,  0.04392145_dp]
                stminvert8_truth(8_dp,:) =  [ 0.01147961_dp, -0.03325185_dp, -0.04083109_dp,  0.02487578_dp,  0.03992714_dp,&
                -0.00644309_dp, -0.02256353_dp,  0.01025478_dp]

            mt_test = mattens(mat,tens,3)
            quad_test = quad(mat,tens,3)
            vectens1_test = vectens1(vec,tens,3)
            vectens2_test = vectens2(vec,tens,3)
            vectens3_test = vectens3(vec,tens,3)
            vectensquad_test = vectensquad(vec,tens,3)
            stminvert6_test = stminvert(stm6,6)
            stminvert8_test = stminvert(stm8,8)

            testpass = .true.
            msg = ""
            if (norm2(mt_test - mt_truth).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision mattens fail."
            endif
            if (norm2(quad_test - quad_test).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision matrix quad fail."
            endif
            if (norm2(vectens1_test - vectens1_truth).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision vectens1 fail."
            endif
            if (norm2(vectens2_test - vectens2_truth).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision vectens2 fail."
            endif
            if (norm2(vectens3_test - vectens3_truth).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision vectens3 fail."
            endif
            if (norm2(vectensquad_test - vectensquad_truth).ge.dtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision vectensquad fail."
            endif
            if (norm2(stminvert6_test - stminvert6_truth).ge.dtol) then
                print *, norm2(stminvert6_test - stminvert6_truth)
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision symplectic invert fail."
            endif
            if (norm2(stminvert8_test - stminvert8_truth).ge.1.e-7_dp) then
                print *, norm2(stminvert8_test - stminvert8_truth)
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Double precision 8x8 invert fail."
            endif

            if (testpass) msg = trim("All basic Double precision tensor operations pass.")
            write (*,*) msg
        end subroutine

        subroutine test_basic_tensorops_q(testpass)
            character(len=100)     :: msg
            logical, intent(inout) :: testpass
            real(qp), parameter :: qtol = 1.e-26_qp
            real(qp) :: vec(3), mat(3,3), tens(3,3,3), &
                        stm6(6,6), stm8(8,8), mt_truth(3,3,3), &
                        quad_truth(3,3,3), vectens1_truth(3,3), &
                        vectens2_truth(3,3), vectens3_truth(3,3), &
                        vectensquad_truth(3), stminvert6_truth(6,6), &
                        stminvert8_truth(8,8), &
                        mt_test(3,3,3), &
                        quad_test(3,3,3), vectens1_test(3,3), &
                        vectens2_test(3,3), vectens3_test(3,3), &
                        vectensquad_test(3), stminvert6_test(6,6), &
                        stminvert8_test(8,8)
            ! vec
            vec =  [10._qp, 30._qp, 93._qp]

            ! mat
                mat(1,:) =  [ 5._qp,  6._qp, 14._qp]
                mat(2,:) =  [86._qp, 51._qp, 15._qp]
                mat(3,:) =  [34._qp, 25._qp, 71._qp]

            ! tens
                tens(1,1,:) =  [ 1._qp,  0._qp, 64._qp]
                tens(1,2,:) =  [21._qp, 80._qp, 32._qp]
                tens(1,3,:) =  [43._qp, 90._qp, 83._qp]
                tens(2,1,:) =  [41._qp,  4._qp, 52._qp]
                tens(2,2,:) =  [70._qp, 26._qp, 44._qp]
                tens(2,3,:) =  [41._qp, 95._qp, 30._qp]
                tens(3,1,:) =  [ 5._qp, 98._qp, 44._qp]
                tens(3,2,:) =  [17._qp, 68._qp, 57._qp]
                tens(3,3,:) =  [39._qp, 21._qp, 63._qp]

            !stm6
                stm6(1,:) =  [92._qp, 68._qp, 73._qp, 85._qp, 40._qp, 51._qp]
                stm6(2,:) =  [63._qp, 38._qp, 47._qp, 71._qp, 88._qp, 13._qp]
                stm6(3,:) =  [33._qp, 53._qp, 82._qp,  6._qp, 45._qp,  3._qp]
                stm6(4,:) =  [ 11._qp, 100._qp,  36._qp,  87._qp,  22._qp,  31._qp]
                stm6(5,:) =  [34._qp, 63._qp, 78._qp, 17._qp, 51._qp, 13._qp]
                stm6(6,:) =  [62._qp, 45._qp, 93._qp,  2._qp, 77._qp, 72._qp]

            !stm8
                stm8(1,:) =  [25._qp, 25._qp, 85._qp, 14._qp, 46._qp, 13._qp, 74._qp, 54._qp]
                stm8(2,:) =  [46._qp, 81._qp, 50._qp, 13._qp, 72._qp, 46._qp, 95._qp, 65._qp]
                stm8(3,:) =  [30._qp, 44._qp,  4._qp, 85._qp, 38._qp, 34._qp, 15._qp, 15._qp]
                stm8(4,:) =  [19._qp, 34._qp, 98._qp, 99._qp, 53._qp, 98._qp, 34._qp, 52._qp]
                stm8(5,:) =  [48._qp, 97._qp, 12._qp, 84._qp, 73._qp, 12._qp, 68._qp, 83._qp]
                stm8(6,:) =  [76._qp, 31._qp, 96._qp, 94._qp, 39._qp, 80._qp,  3._qp, 54._qp]
                stm8(7,:) =  [18._qp, 49._qp, 96._qp, 81._qp, 50._qp,  6._qp, 67._qp, 54._qp]
                stm8(8,:) =  [95._qp, 77._qp, 70._qp, 35._qp, 55._qp, 49._qp, 87._qp, 11._qp]

            ! mattens
                mt_truth(1,1,:) =  [ 321._qp, 1396._qp, 1248._qp]
                mt_truth(1,2,:) =  [ 763._qp, 1508._qp, 1222._qp]
                mt_truth(1,3,:) =  [1007._qp, 1314._qp, 1477._qp]
                mt_truth(2,1,:) =  [2252._qp, 1674._qp, 8816._qp]
                mt_truth(2,2,:) =  [5631._qp, 9226._qp, 5851._qp]
                mt_truth(2,3,:) =  [ 6374._qp, 12900._qp,  9613._qp]
                mt_truth(3,1,:) =  [1414._qp, 7058._qp, 6600._qp]
                mt_truth(3,2,:) =  [3671._qp, 8198._qp, 6235._qp]
                mt_truth(3,3,:) =  [5256._qp, 6926._qp, 8045._qp]

            ! quad
                quad_truth(1,1,:) =  [1071601._qp,  673928._qp,  613396._qp]
                quad_truth(1,2,:) =  [694234._qp, 438017._qp, 415539._qp]
                quad_truth(1,3,:) =  [916796._qp, 589107._qp, 677297._qp]
                quad_truth(2,1,:) =  [682067._qp, 452100._qp, 548500._qp]
                quad_truth(2,2,:) =  [456959._qp, 301671._qp, 358375._qp]
                quad_truth(2,3,:) =  [760713._qp, 481901._qp, 421133._qp]
                quad_truth(3,1,:) =  [867513._qp, 558130._qp, 660906._qp]
                quad_truth(3,2,:) =  [564690._qp, 363513._qp, 431889._qp]
                quad_truth(3,3,:) =  [551504._qp, 365197._qp, 523585._qp]

            ! vectens1
                vectens1_truth(1,:) =  [1705._qp, 9234._qp, 6292._qp]
                vectens1_truth(2,:) =  [3891._qp, 7904._qp, 6941._qp]
                vectens1_truth(3,:) =  [5287._qp, 5703._qp, 7589._qp]

            ! vectens2
                vectens2_truth(1,:) =  [ 4639._qp, 10770._qp,  9319._qp]
                vectens2_truth(2,:) =  [6323._qp, 9655._qp, 4630._qp]
                vectens2_truth(3,:) =  [4187._qp, 4973._qp, 8009._qp]

            ! vectens3
                vectens3_truth(1,:) =  [ 5962._qp,  5586._qp, 10849._qp]
                vectens3_truth(2,:) =  [5366._qp, 5572._qp, 6050._qp]
                vectens3_truth(3,:) =  [7082._qp, 7511._qp, 6879._qp]

            ! vectensquad
                vectensquad_truth =  [1236157._qp,  783470._qp,  935897._qp]

            ! STMinvert6
                stminvert6_truth(1,:) =  [ 87._qp,  17._qp,   2._qp, -85._qp, -71._qp,  -6._qp]
                stminvert6_truth(2,:) =  [ 22._qp,  51._qp,  77._qp, -40._qp, -88._qp, -45._qp]
                stminvert6_truth(3,:) =  [ 31._qp,  13._qp,  72._qp, -51._qp, -13._qp,  -3._qp]
                stminvert6_truth(4,:) =  [-11._qp, -34._qp, -62._qp,  92._qp,  63._qp,  33._qp]
                stminvert6_truth(5,:) =  [-100._qp,  -63._qp,  -45._qp,   68._qp,   38._qp,   53._qp]
                stminvert6_truth(6,:) =  [-36._qp, -78._qp, -93._qp,  73._qp,  47._qp,  82._qp]

            ! STMinvert8
                stminvert8_truth(1_qp,:) =  [ 0.01684607_qp, -0.02049215_qp, -0.00492154_qp,  0.00023436_qp,  0.01364028_qp,&
                0.00309208_qp, -0.0174303 _qp,  0.01145994_qp]
                stminvert8_truth(2_qp,:) =  [-0.05263277_qp,  0.03873808_qp, -0.01003245_qp, -0.01360167_qp, -0.01584197_qp,&
                0.01094427_qp,  0.03666172_qp, -0.00671567_qp]
                stminvert8_truth(3_qp,:) =  [-0.01830201_qp,  0.03185471_qp,  0.01231027_qp, -0.02074384_qp, -0.02911729_qp,&
                0.01512645_qp,  0.02883093_qp, -0.01319872_qp]
                stminvert8_truth(4_qp,:) =  [ 0.0121908 _qp, -0.03885463_qp, -0.01532077_qp,  0.02289227_qp,  0.02725751_qp,&
                -0.01398507_qp, -0.01438121_qp,  0.0160061 _qp]
                stminvert8_truth(5_qp,:) =  [ 0.00310332_qp,  0.11799515_qp,  0.1345172 _qp, -0.09211388_qp, -0.10954499_qp,&
                0.04305073_qp,  0.04573723_qp, -0.06976503_qp]
                stminvert8_truth(6_qp,:) =  [ 0.00150244_qp, -0.02015711_qp, -0.02413202_qp,  0.02898049_qp,  0.02012894_qp,&
                -0.01155453_qp, -0.02071762_qp,  0.0141881 _qp]
                stminvert8_truth(7_qp,:) =  [ 0.03374979_qp, -0.08094552_qp, -0.05577335_qp,  0.05802861_qp,  0.06445587_qp,&
                -0.03950086_qp, -0.04521933_qp,  0.04392145_qp]
                stminvert8_truth(8_qp,:) =  [ 0.01147961_qp, -0.03325185_qp, -0.04083109_qp,  0.02487578_qp,  0.03992714_qp,&
                -0.00644309_qp, -0.02256353_qp,  0.01025478_qp]

            mt_test = mattens(mat,tens,3)
            quad_test = quad(mat,tens,3)
            vectens1_test = vectens1(vec,tens,3)
            vectens2_test = vectens2(vec,tens,3)
            vectens3_test = vectens3(vec,tens,3)
            vectensquad_test = vectensquad(vec,tens,3)
            stminvert6_test = stminvert(stm6,6)
            stminvert8_test = stminvert(stm8,8)

            testpass = .true.
            msg = ""
            if (norm2(mt_test - mt_truth).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision mattens fail."
            endif
            if (norm2(quad_test - quad_test).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision matrix quad fail."
            endif
            if (norm2(vectens1_test - vectens1_truth).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision vectens1 fail."
            endif
            if (norm2(vectens2_test - vectens2_truth).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision vectens2 fail."
            endif
            if (norm2(vectens3_test - vectens3_truth).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision vectens3 fail."
            endif
            if (norm2(vectensquad_test - vectensquad_truth).ge.qtol) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision vectensquad fail."
            endif
            if (norm2(stminvert6_test - stminvert6_truth).ge.qtol) then
                print *, norm2(stminvert6_test - stminvert6_truth)
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision symplectic invert fail."
            endif
            if (norm2(stminvert8_test - stminvert8_truth).ge.1.e-7_qp) then
                print *, norm2(stminvert8_test - stminvert8_truth)
                testpass = .false.
                msg = trim(msg)//new_line("a")//"Quad precision 8x8 invert fail."
            endif

            if (testpass) msg = trim("All basic Quad precision tensor operations pass.")
            write (*,*) msg
        end subroutine
end module test_tensorops
