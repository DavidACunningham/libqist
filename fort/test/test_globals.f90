module test_globals
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use globals, only: mmult
    implicit none
    contains
         subroutine test_mmult(testpass)
            logical, intent(inout) :: testpass
            character(len=300) :: msg
            real(dp) :: mata_d(8,8), &
                        matb_d(8,6), &
                        veca_d(8), &
                        vecb_d(8)
            real(qp) :: mata_q(8,8), &
                        matb_q(8,6), &
                        veca_q(8), &
                        vecb_q(8)
            real(dp), parameter :: dtol = 1.e-13_dp
            real(qp), parameter :: qtol = 1.e-25_qp

            call random_number(mata_d)
            call random_number(matb_d)
            call random_number(veca_d)
            call random_number(vecb_d)
            call random_number(mata_q)
            call random_number(matb_q)
            call random_number(veca_q)
            call random_number(vecb_q)
            testpass = .true.
            msg = ""
            if (norm2(mmult(mata_d,matb_d)-matmul(mata_d,matb_d)).ge.dtol) then
                testpass=.false.
                msg = trim(msg)//"Double precision matrix-matrix multiply fail."
            endif
            if (norm2(mmult(mata_q,matb_q)-matmul(mata_q,matb_q)).ge.qtol) then
                testpass=.false.
                msg = trim(msg)//new_line("a")//" Quad precision matrix-matrix multiply fail."
            endif
            if (norm2(mmult(mata_d,veca_d)-matmul(mata_d,veca_d)).ge.dtol) then
                testpass=.false.
                msg = trim(msg)//new_line("a")//" Double precision matrix-vector multiply fail."
            endif
            if (norm2(mmult(mata_q,veca_q)-matmul(mata_q,veca_q)).ge.qtol) then
                testpass=.false.
                msg = trim(msg)//new_line("a")//" Quad precision matrix-vector multiply fail."
            endif
            if (norm2(mmult(veca_d,mata_d)-matmul(veca_d,mata_d)).ge.dtol) then
                testpass=.false.
                msg = trim(msg)//new_line("a")//" Double precision vector-matrix multiply fail."
            endif
            if (norm2(mmult(veca_q,mata_q)-matmul(veca_q,mata_q)).ge.qtol) then
                testpass=.false.
                msg = trim(msg)//new_line("a")//" Quad precision vector-matrix multiply fail."
            endif
            if (testpass) msg = "Matrix multiplication tests passed."
            write (*,*) trim(msg)
        end subroutine test_mmult
end module test_globals
