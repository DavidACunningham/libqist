module test_frkmin
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use frkmin, only: solve_ivp, ODESolution, RungeKutta
    implicit none
    real(dp), parameter :: dtol = 9.e-14_dp
    real(dp), parameter :: &
                         & eps = 3.e0_dp, &
                           kepy0(6)  = [0._dp, 1._dp, 0._dp,  &
                                  sqrt(2._dp/((1._dp+eps))), 0._dp, &
                                  sqrt(2._dp/((1._dp+eps)))], &
                          pi = 2._dp*asin(1._dp), &
                          r0 = sqrt(sum(kepy0(:3)**2)), &
                          v0 = sqrt(sum(kepy0(4:6)**2)), &
                          E0 = v0**2/2._dp - 1._dp/r0, &
                          a  = -1._dp/(2._dp*E0), &
                          P  = 2._dp*pi*sqrt(a**3._dp), &
                          rtol = 3.e-14_dp, &
                          atol = 1.e-20_dp
    contains 
        subroutine run_frk_tests(testpass)
            logical, intent(inout) :: testpass
            character(len=102)     :: msg
            logical testpassvec(3)
            testpass = .true.
            call frk_sparse_test(testpassvec(1))
            msg = ""
            if (.not.testpassvec(1)) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//" Sparse Double integration fail."
            end if
            call frk_dense_test(testpassvec(2))
            if (.not.testpassvec(2)) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//" Dense Double integration fail."
            end if
            call frk_readwrite_test(testpassvec(3))
            if (.not.testpassvec(1)) then
                testpass = .false.
                msg = trim(msg)//new_line("a")//" Double integration read/write fail."
            end if
            if (testpass) msg = "Double precision integration tests pass."
            write (*,*) trim(msg)
        end subroutine
        subroutine frk_sparse_test(testpass)
            logical, intent(inout) :: testpass
            type(ODESolution)   :: solution
            integer, parameter :: ntest = 1000
            real(dp)              :: t0f(2), xf(6), vf, rf, Ef, dt, dE, dx(6)
            t0f = [0._dp, P]
            solution = solve_ivp(eom, t0f, kepy0, dense_output=.false., &
                                 rtol=rtol, atol=atol)
            xf = solution%ys(:,size(solution%ts))
            dx = xf - kepy0
            dt = P - solution%ts(size(solution%ts))
            rf = sqrt(sum(xf(:3)**2))
            vf = sqrt(sum(xf(4:6)**2))
            Ef = vf**2/2._dp - 1._dp/rf
            dE = Ef - E0
            testpass = .true.
            if (any(abs([dE, dx, dt]).ge.dtol)) then
                testpass = .false.
            endif
        end subroutine frk_sparse_test
        subroutine frk_dense_test(testpass)
            type(ODESolution)   :: solution, dummy_solution
            logical, intent(inout) :: testpass
            real(dp)              :: t0f(2), &
                                     testtimes(50), &
                                     teststates(6,50), &
                                     truestates(6,50), &
                                     diffs(6,50)
            integer i
            t0f = [0._dp, P]
            testtimes = [(t0f(1) + i*(t0f(2)-t0f(1))/49._dp, i=0,49)]
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
            real(dp)              :: t0f(2), testtimes(50), deltas(6,50)
            integer i, num
            t0f = [0._dp, P]
            testtimes = [(t0f(1) + i*(t0f(2)-t0f(1))/49._dp, i=1,50)]
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
            real(dp),            intent(in)    :: t
            real(dp),            intent(in)    :: y(:)
            real(dp)                           :: res(size(y))
            res = [y(4:6), -y(:3)/norm2(y(:3))**3]
        end function eom
end module test_frkmin
