module test_frkmin_q
    use, intrinsic :: iso_fortran_env, only: qp=>real128
    use frkmin_q, only: solve_ivp, ODESolution, RungeKutta
    implicit none
    real(qp), parameter :: dtol = 9.e-26_qp
    real(qp), parameter :: &
                         & eps = 3.e0_qp, &
                           kepy0(6)  = [0._qp, 1._qp, 0._qp,  &
                                  sqrt(2._qp/((1._qp+eps))), 0._qp, &
                                  sqrt(2._qp/((1._qp+eps)))], &
                          pi = 2._qp*asin(1._qp), &
                          r0 = sqrt(sum(kepy0(:3)**2)), &
                          v0 = sqrt(sum(kepy0(4:6)**2)), &
                          E0 = v0**2/2._qp - 1._qp/r0, &
                          a  = -1._qp/(2._qp*E0), &
                          P  = 2._qp*pi*sqrt(a**3._qp), &
                          rtol = 3.e-28_qp, &
                          atol = 1.e-30_qp
    contains 
        subroutine run_frk_tests_q(testpass)
            logical, intent(inout) :: testpass
            character(len=102)     :: msg
            logical testpassvec(3)
            testpass = .true.
            write (*,*) "Testing sparse Quad integration. . ."
            call frk_sparse_test(testpassvec(1))
            msg = ""
            if (.not.testpassvec(1)) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"FAIL Sparse Quad integration FAIL"
            end if
            write (*,*) "Testing dense Quad integration. . ."
            call frk_dense_test(testpassvec(2))
            if (.not.testpassvec(2)) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"FAIL Dense Quad integration FAIL"
            end if
            write (*,*) "Testing Quad integration read/write. . ."
            call frk_readwrite_test(testpassvec(3))
            if (.not.testpassvec(1)) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//"FAIL Quad integration read/write FAIL"
            end if
            if (testpass) msg = "PASS Quad precision integration tests PASS"
            write (*,*) trim(msg)
        end subroutine
        subroutine frk_sparse_test(testpass)
            logical, intent(inout) :: testpass
            type(ODESolution)   :: solution
            integer, parameter :: ntest = 1000
            real(qp)              :: t0f(2), xf(6), vf, rf, Ef, dt, dE, dx(6)
            t0f = [0._qp, P]
            solution = solve_ivp(eom, t0f, kepy0, dense_output=.false., &
                                 rtol=rtol, atol=atol)
            xf = solution%ys(:,size(solution%ts))
            dx = xf - kepy0
            dt = P - solution%ts(size(solution%ts))
            rf = sqrt(sum(xf(:3)**2))
            vf = sqrt(sum(xf(4:6)**2))
            Ef = vf**2/2._qp - 1._qp/rf
            dE = Ef - E0
            testpass = .true.
            if (any(abs([dE, dx, dt]).ge.dtol)) then
                testpass = .false.
            endif
        end subroutine frk_sparse_test
        subroutine frk_dense_test(testpass)
            type(ODESolution)   :: solution, dummy_solution
            logical, intent(inout) :: testpass
            real(qp)              :: t0f(2), &
                                     testtimes(50), &
                                     teststates(6,50), &
                                     truestates(6,50), &
                                     diffs(6,50)
            integer i
            t0f = [0._qp, P]
            testtimes = [(t0f(1) + i*(t0f(2)-t0f(1))/49._qp, i=0,49)]
            solution = solve_ivp(eom, t0f, kepy0, dense_output=.true., &
                                 rtol=rtol, atol=atol)
            do i=1,50
                teststates(:,i) = solution%call(testtimes(i))
                dummy_solution  = solve_ivp(eom, [t0f(1), testtimes(i)], &
                                            kepy0, dense_output=.false., &
                                            rtol=rtol, atol=atol)
                truestates(:,i) = dummy_solution%ys(:,size(dummy_solution%ts))
            end do
            diffs = teststates - truestates
            testpass =.true.
            if (any(abs(diffs).ge.dtol)) then
                testpass = .false.
            endif
        end subroutine frk_dense_test
        subroutine frk_readwrite_test(testpass)
            type(ODESolution)   :: solution, solution_after_save
            logical, intent(inout) :: testpass
            real(qp)              :: t0f(2), testtimes(50), deltas(6,50)
            integer i, num
            t0f = [0._qp, P]
            testtimes = [(t0f(1) + i*(t0f(2)-t0f(1))/49._qp, i=1,50)]
            solution = solve_ivp(eom, t0f, kepy0, dense_output=.true., &
                                 rtol=rtol, atol=atol)
            write(*,*) "Writing solution to temp file"
            open(newunit=num,file="./tempsolution.odesolution", &
                 status='replace',access='stream')
            call solution%write(num)
            close(num)
            write(*,*) "Reading solution from temp file"
            open(newunit=num,file="./tempsolution.odesolution", &
                 status='old',access='stream')
            call solution_after_save%read(num)
            close(num)
            do i=1,4
                deltas(:,i) = abs( solution%call(testtimes(i)) &
                            - solution_after_save%call(testtimes(i)))
            end do
            write(*,*) "Cleaning up solution temp file. . ."
            open(newunit=num,file="./tempsolution.odesolution", &
                 status='old',access='stream')
            close(num,status="delete")
            write(*,*) "Done."
            testpass = .true.
            if (any(deltas.gt.dtol)) then
                testpass = .false.
            end if
        end subroutine frk_readwrite_test
        function eom(self,t,y) result(res)
            class(RungeKutta), intent(inout) :: self
            real(qp),            intent(in)    :: t
            real(qp),            intent(in)    :: y(:)
            real(qp)                           :: res(size(y))
            res = [y(4:6), -y(:3)/norm2(y(:3))**3]
        end function eom
end module test_frkmin_q
