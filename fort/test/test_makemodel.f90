module test_makemodel
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use genqist, only: gqist
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use findiffmod
    use tensorops, only: mattens, quad
    use test_util, only: mprint, tprint
    implicit none
    real(qp), parameter :: qtol = 1.e-18_qp, qtol_tb = 1.e-16_qp, qtol_sh=2.e-13_qp
    real(dp), parameter :: dtol_fd = 2.e-4_dp
    contains
        subroutine test_kep_grav(testpass)
            logical, intent(inout) :: testpass
            character(len=200)     :: msg
            logical                :: testvec(2)
            integer, parameter     :: nnodes=1000
            real(qp)               :: t0, tf
            real(qp)               :: y(8), acc(8), jac(8,8), hes(8,8,8), &
                                      jac_fd(8,8), hes_fd(8,8,8), &
                                      fdsteps(8), rpert, vpert, tpert, &
                                      calltime, jacdiff(8,8), hesdiff(8,8,8)
            type(gqist)            :: gq
            
            acc = 0._qp
            jac = 0._qp
            hes = 0._qp
            rpert = 1._qp ! km
            vpert = .01_qp ! km/s
            tpert =  1000._qp ! s
            fdsteps(:3) = rpert
            fdsteps(4:6) = vpert
            fdsteps(7:8) = tpert
            t0 =  772509669.184287_qp  ! 2024 Jun 24 14:00 UTC
            tf =  772682469.184231_qp  ! 2024 Jun 26 14:00 UTC
            calltime =  (tf-t0)/2 + t0
            call gq%init("./test_config_namelists.nml")
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%regularize = .false.
            gq%dynmod%tof = gq%tf - gq%t0
            y = [5000._qp, 1000._qp, 500._qp, 0._qp, 0._qp, 0._qp, calltime, tf-t0]
            call gq%dynmod%allderivs_kepler(gq%dynmod%central_body_mu, y, &
                                         & acc, jac, hes)
            jac_fd = 0._qp
            jac_fd = findiff_multiscale(fd_acc,y, fdsteps, 9)
            hes_fd = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
            jacdiff = jac_fd - jac
            hesdiff = hes_fd - hes
            testvec(1) = all(abs(jac-jac_fd).lt.qtol)
            testvec(2) = all(abs(hes-hes_fd).lt.qtol)
            testpass = .true.
            msg = ""
            if (.not.testvec(1)) msg = trim(msg)//new_line("a")//"FAIL Keplerian Jacobian test FAIL"
            if (.not.testvec(2)) msg = trim(msg)//new_line("a")//"FAIL Keplerian Hessian test FAIL"
            testpass = all(testvec)
            if (testpass) msg = trim("PASS Keplerian derivatives test PASS")
            write (*,*) msg
            contains
            function fd_acc(x) result(res)
                real(qp), intent(in) :: x(:)
                real(qp)             :: res(size(x)), &
                                      & jac(size(x),size(x)), &
                                      & hes(size(x),size(x),size(x))
                res = 0._qp
                call gq%dynmod%allderivs_kepler(gq%dynmod%central_body_mu, x, res,jac,hes)
            end function fd_acc
            function fd_jac(x) result(res)
                real(qp), intent(in) :: x(:)
                real(qp)             :: res(size(x),size(x)), &
                                      & acc(size(x)), &
                                      & hes(size(x),size(x),size(x))
                res = 0._qp
                call gq%dynmod%allderivs_kepler(gq%dynmod%central_body_mu, x, acc,res,hes)
            end function fd_jac
        end subroutine
        subroutine test_threebody_grav(testpass)
            logical, intent(inout) :: testpass
            character(len=200)     :: msg
            logical                :: testvec(2)
            integer, parameter     :: nnodes=1000
            real(qp)               :: t0, tf
            real(qp)               :: y(8), acc(8), jac(8,8), hes(8,8,8), &
                                      jac_fd(8,8), hes_fd(8,8,8), &
                                      pos_tb(3), vel_tb(3), acc_tb(3), &
                                      fdsteps(8), rpert, vpert, tpert, &
                                      calltime, jacdiff(8,8), hesdiff(8,8,8)
            type(gqist)            :: gq
            
            acc = 0._qp
            jac = 0._qp
            hes = 0._qp
            rpert = .01_qp ! km
            vpert = .001_qp ! km/s
            tpert =  1000._qp ! s
            fdsteps(:3) = rpert
            fdsteps(4:6) = vpert
            fdsteps(7:8) = tpert
            t0 =  772509669.184287_qp  ! 2024 Jun 24 14:00 UTC
            tf =  772682469.184231_qp  ! 2024 Jun 26 14:00 UTC
            calltime =  (tf-t0)/2 + t0
            call gq%init("./test_config_namelists.nml")
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%regularize = .false.
            gq%dynmod%tof = gq%tf - gq%t0
            y = [5000._qp, 1000._qp, 500._qp, 0._qp, 0._qp, 0._qp, calltime, tf-t0]
            pos_tb = real(gq%dynmod%bod_db%call(real(calltime,dp), &
                                               & gq%dynmod%bodylist(1),'p'), qp)
            vel_tb = real(gq%dynmod%bod_db%call(real(calltime,dp), &
                                               & gq%dynmod%bodylist(1),'v'), qp)
            acc_tb = real(gq%dynmod%bod_db%call(real(calltime,dp), &
                                               & gq%dynmod%bodylist(1),'a'), qp)
            call gq%dynmod%allderivs_thirdbody(gq%dynmod%nbody_mus(1), y, &
                                         & pos_tb, vel_tb, acc_tb, acc, jac, hes)
            jac_fd = 0._qp
            jac_fd = findiff_multiscale(fd_acc,y, fdsteps, 9)
            hes_fd = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
            jacdiff = (jac- jac_fd)
            hesdiff = (hes- hes_fd)
            testvec(1) = all(abs(jac-jac_fd).lt.qtol_tb)
            testvec(2) = all(abs(hes-hes_fd).lt.qtol)
            testpass = .true.
            msg = ""
            if (.not.testvec(1)) msg = trim(msg)//new_line("a")//"FAIL Third body Jacobian test FAIL"
            if (.not.testvec(2)) msg = trim(msg)//new_line("a")//"FAIL Third body Hessian test FAIL"
            testpass = all(testvec)
            if (testpass) msg = trim("PASS Third body derivatives test PASS")
            write (*,*) msg
            contains
            function fd_acc(x) result(res)
                real(qp), intent(in) :: x(:)
                real(qp)             :: res(size(x)), &
                                      & jac(size(x),size(x)), &
                                      & hes(size(x),size(x),size(x)), &
                                      & thispos(3), thisvel(3), thisacc(3)
                res = 0._qp
                thispos = real(gq%dynmod%bod_db%call(real(x(7),dp), &
                                                   & gq%dynmod%bodylist(1),'p'), qp)
                thisvel = real(gq%dynmod%bod_db%call(real(x(7),dp), &
                                                   & gq%dynmod%bodylist(1),'v'), qp)
                thisacc = real(gq%dynmod%bod_db%call(real(x(7),dp), &
                                                   & gq%dynmod%bodylist(1),'a'), qp)
                call gq%dynmod%allderivs_thirdbody(gq%dynmod%nbody_mus(1), &
                                x, thispos, thisvel, thisacc, res,jac,hes)
            end function fd_acc
            function fd_jac(x) result(res)
                real(qp), intent(in) :: x(:)
                real(qp)             :: res(size(x),size(x)), &
                                      & acc(size(x)), &
                                      & hes(size(x),size(x),size(x)), &
                                      & thispos(3), thisvel(3), thisacc(3)
                res = 0._qp
                thispos = real(gq%dynmod%bod_db%call(real(x(7),dp), &
                                                   & gq%dynmod%bodylist(1),'p'), qp)
                thisvel = real(gq%dynmod%bod_db%call(real(x(7),dp), &
                                                   & gq%dynmod%bodylist(1),'v'), qp)
                thisacc = real(gq%dynmod%bod_db%call(real(x(7),dp), &
                                                   & gq%dynmod%bodylist(1),'a'), qp)
                call gq%dynmod%allderivs_thirdbody(gq%dynmod%nbody_mus(1), &
                                x, thispos, thisvel, thisacc, acc,res,hes)
            end function fd_jac
        end subroutine
        subroutine test_SH_grav(testpass)
            logical, intent(inout) :: testpass
            character(len=200)     :: msg
            logical                :: testvec(2)
            integer, parameter     :: nnodes=1000
            real(qp)               :: t0, tf
            real(qp)               :: y(8), acc(8), jac(8,8), hes(8,8,8), &
                                      jac_fd(8,8), hes_fd(8,8,8), &
                                      fdsteps(8), rpert, vpert, tpert, &
                                      calltime, jacdiff(8,8), hesdiff(8,8,8)
            type(gqist)            :: gq
            acc = 0._qp
            jac = 0._qp
            hes = 0._qp
            rpert = 10._qp ! km
            vpert = .001_qp ! km/s
            tpert =  1000._qp ! s
            fdsteps(:3) = rpert
            fdsteps(4:6) = vpert
            fdsteps(7:8) = tpert
            t0 =  772509669.184287_qp  ! 2024 Jun 24 14:00 UTC
            tf =  772682469.184231_qp  ! 2024 Jun 26 14:00 UTC
            calltime =  (tf-t0)/2 + t0
            call gq%init("./test_config_namelists.nml")
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%regularize = .false.
            y = [5000._qp, 1000._qp, 500._qp, 0._qp, 0._qp, 0._qp, calltime, tf-t0]
            call gq%dynmod%allderivs_sh(calltime, y, acc, jac, hes)
            jac_fd = 0._qp
            jac_fd = findiff_multiscale(fd_acc,y, fdsteps, 9)
            hes_fd = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
            jacdiff = jac - jac_fd
            hesdiff = hes - hes_fd
            testvec(1) = all(abs(jac-jac_fd).lt.qtol_sh)
            testvec(2) = all(abs(hes-hes_fd).lt.qtol_tb)
            testpass = .true.
            msg = ""
            if (.not.testvec(1)) msg = trim(msg)//new_line("a")//"FAIL Spherical harmonics gravity Jacobian test FAIL"
            if (.not.testvec(2)) msg = trim(msg)//new_line("a")//"FAIL Spherical harmonics gravity Hessian test FAIL"
            testpass = all(testvec)
            if (testpass) msg = trim("PASS Spherical harmonics derivatives test PASS")
            write (*,*) msg
            contains
            function fd_acc(x) result(res)
                real(qp), intent(in) :: x(:)
                real(qp)             :: res(size(x)), &
                                      & jac(size(x),size(x)), &
                                      & hes(size(x),size(x),size(x))
                res = 0._qp
                call gq%dynmod%allderivs_sh(x(7), x, res,jac,hes)
            end function fd_acc
            function fd_jac(x) result(res)
                real(qp), intent(in) :: x(:)
                real(qp)             :: res(size(x),size(x)), &
                                      & acc(size(x)), &
                                      & hes(size(x),size(x),size(x))
                res = 0._qp
                call gq%dynmod%allderivs_sh(x(7), x, acc,res,hes)
            end function fd_jac
        end subroutine
        subroutine end_to_end_integration_test(testpass)
            logical, intent(inout) :: testpass
            logical                :: testvec(6)
            type(gqist)            :: gq
            type(odesolution)      :: base_sol 
            character(len=1000)    :: qist_config_file, metakernel_filepath
            real(qp)               :: t0, tf, tof, &
                                      rtol, atol, x0(6), fdrstep, fdvstep, fdtstep, &
                                      epsvec(8)
            integer                :: fdord
            integer                :: stat, num
            real(qp)               :: init_state(8), eye(8,8), fd_stm(8,8), &
                                  & analytic_stm(8,8), analytic_stt(8,8,8), & 
                                  & analytic_final_state(8), &
                                  & jacdiff(8,8), hesdiff(8,8,8), stmdiff(8,8), &
                                  & sttdiff(8,8,8), &
                                  & divisor(8,8), init_stt(8**3), fd_stt(8,8,8), &
                                  & stt_divisor(8,8,8), jac(8,8), fd_jac(8,8), &
                                  & hes(8,8,8), fd_hes(8,8,8), acc(8), thistime

            real(dp), allocatable :: spice_state(:,:)
            real(dp)              :: lt_dum, state_dum(6)
            integer i, j
            namelist /FD_CONFIG/   metakernel_filepath, &
                                   qist_config_file, &
                                   x0, &
                                   t0, &
                                   tf, &
                                   rtol, &
                                   atol, &
                                   fdrstep, &
                                   fdvstep, &
                                   fdtstep, &
                                   fdord

            ! Read in namelist
            inquire(file="./test_fd_config.nml", iostat=stat)
            if (stat .ne. 0) then 
                print *, "ERROR: Bad FD config namelist filename"
                stop
            end if
            open(file="./test_fd_config.nml", status="old", &
                 iostat=stat,newunit=num)
            read(unit=num, nml=FD_CONFIG, iostat=stat)
            if (stat .ne. 0) then 
                print *, "ERROR: bad FD config namelist format"
                print *, "./test_fd_config.nml"
                print *, stat
                stop
            end if
            close(num)
            epsvec(1:3) = fdrstep
            epsvec(4:6) = fdvstep
            epsvec(7:8) = fdtstep
            call gq%init(qist_config_file,.true.)
            eye = 0._qp; do i=1,8; eye(i,i) = 1._qp; end do
            init_stt = 0._qp
            tof = tf-t0
            ! To integrate in real time, set tof to 1.
            call furnsh(trim(adjustl(metakernel_filepath)))
            call spkgeo(gq%dynmod%traj_id, &
                        real(t0,dp), &
                        "J2000", &
                        gq%dynmod%central_body, &
                        state_dum, &
                        lt_dum &
                       )
            gq%t0 = t0
            gq%tf = tf
            gq%dynmod%tof = tof
            init_state = [gq%dynmod%trajstate(t0), t0, tof]
            gq%dynmod%regularize = .false.
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = init_state

            testvec(1) = all(abs(state_dum - init_state(:6)).lt.1.e-6)
            if (.not. testvec(1)) then
                print *, "FAIL Kernel resample FAIL"
                print *, "Error: "
                print *, abs(state_dum - init_state(:6))
            endif
            call gq%dynmod%get_derivs(init_state(7), acc, jac, hes)
            fd_jac = findiff_multiscale(fd_wrap_acc, init_state, epsvec, fdord)
            fd_hes = findiffhes_multiscale(fd_wrap_jac, init_state, epsvec, fdord)

            divisor = 1._qp
            stt_divisor = 1._qp
            where (abs(jac).ge.1.e-14_qp); divisor = jac; end where
            where (abs(hes).ge.1.e-10_qp); stt_divisor = hes; end where
            jacdiff = (jac - fd_jac)/divisor
            hesdiff = (hes - fd_hes)/stt_divisor
            testvec(2) = all(jacdiff.lt.5.e-5_qp)
            testvec(3) = all(hesdiff.lt.5.e-5_qp)
            if (.not.testvec(2)) then
                print *, "FAIL Dynamics Model Jacobian test FAIL"
                print *, "FD JACOBIAN"
                call mprint(fd_jac)
                print *, "ANALYTIC JACOBIAN"
                call mprint(jac)
                print *, "NORMALIZED JACOBIAN ERROR"
                call mprint(jacdiff)
            endif
            if (.not.testvec(3)) then
                print *, "FAIL Dynamics Model Hessian test FAIL"
                print *, "FD HESSIAN"
                call tprint(fd_hes)
                print *, "ANALYTIC HESSIAN"
                call tprint(hes)
                print *, "NORMALIZED HESSIAN ERROR"
                call tprint(hesdiff)
            endif
            print *, "Integrating base case"
            base_sol = solve_ivp(fd_eoms,&
                               & [0._qp, 1._qp], &
                               & [init_state, &
                                  reshape(eye,[8**2]), &
                                  init_stt], &
                               & dense_output=.false.,&
                               & rtol=rtol, &
                               & atol=atol, &
                               & istep=0.5_qp &
                              & )
            print *, "Done."
            allocate(spice_state(6,size(base_sol%ts)))
            do i=1,size(base_sol%ts)
                thistime = base_sol%ts(i)*tof+ t0
                call spkgeo(gq%dynmod%traj_id, &
                            real(thistime,dp), &
                            "J2000", &
                            gq%dynmod%central_body, &
                            spice_state(:,i), &
                            lt_dum &
                           )

            end do
            analytic_final_state = base_sol%ys(:8,size(base_sol%ts))
            analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
            analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8]) 
            testvec(4) = all( &
                             (&
                               analytic_final_state(:6) - &
                               real(spice_state(:,size(base_sol%ts)),qp) &
                             )/analytic_final_state(:6).lt.2.e-3_qp &
                            )
            if (.not.testvec(4)) then
                print *, "FAIL Dynamics integration test FAIL"
                print *, "FINAL STATE ERROR"
                print *, real(&
                              (&
                               analytic_final_state(:6) - &
                               real(spice_state(:,size(base_sol%ts)),qp) &
                              )/analytic_final_state(:6),8 &
                             )
                print *, "FINAL POSITION ERROR"
            end if

            print *, "Integrating finite diff first order with step ", real(epsvec,4)
            fd_stm = findiff_multiscale(fd_integrate, init_state, epsvec, fdord)
            print *, "Done."

            print *, "Integrating finite diff second order with step ", real(epsvec,4)
            fd_stt = findiffhes_multiscale(fd_integrate_jac, init_state, epsvec, fdord)
            print *, "Done."

            divisor = 1._qp
            stt_divisor = 1._qp
            where (abs(analytic_stm).ge.1.e-14_qp); divisor = jac; end where
            where (abs(analytic_stt).ge.1.e-10_qp); stt_divisor = hes; end where
            stmdiff = (analytic_stm - fd_stm)/divisor
            sttdiff = (analytic_stt - fd_stt)/stt_divisor
            testvec(5) = all(jacdiff.lt.5.e-4_qp)
            testvec(6) = all(hesdiff.lt.5.e-4_qp)
            if (.not.testvec(5)) then
                print *, "FAIL Dynamics Model STM test FAIL"
                print *, "FD STM"
                call mprint(fd_stm)
                print *, "ANALYTIC STM"
                call mprint(analytic_stm)
                print *, "NORMALIZED STM ERROR"
                call mprint(stmdiff)
            endif
            if (.not.testvec(6)) then
                print *, "FAIL Dynamics Model STT test FAIL"
                print *, "FD STT"
                call tprint(fd_stt)
                print *, "ANALYTIC STT"
                call tprint(analytic_stt)
                print *, "NORMALIZED STT ERROR"
                call tprint(sttdiff)
            endif
            testpass = all(testvec)
            if (.not.testpass) then
                print *, "FAIL End to end dynamics test FAIL"
            else
                print *, "PASS End to end dynamics test PASS"
            endif
            contains
                function fd_wrap_acc(x) result(res)
                    real(qp), intent(in) :: x(:)
                    real(qp)             :: res(size(x))
                    real(qp)             :: jac(8,8), hes(8,8,8), &
                                            acc(8)

                    gq%dynmod%tgt_on_rails = .false.
                    gq%dynmod%state = x
                    call gq%dynmod%get_derivs(x(7), acc, jac, hes)
                    res = acc
                end function
                function fd_wrap_jac(x) result(res)
                    real(qp), intent(in) :: x(:)
                    real(qp)             :: res(size(x),size(x))
                    real(qp)             :: jac(8,8), hes(8,8,8), &
                                            acc(8)

                    gq%dynmod%tgt_on_rails = .false.
                    gq%dynmod%state = x
                    call gq%dynmod%get_derivs(x(7), acc, jac, hes)
                    res = jac
                end function
                function fd_integrate(x) result(res)
                    type(odesolution) :: fd_sol
                    real(qp), intent(in) :: x(:)
                    real(qp)             :: res(size(x))
                    fd_sol = solve_ivp(fd_eoms,&
                                     & [0._qp, 1._qp], &
                                     & x(:8), &
                                     & dense_output=.false., &
                                     & rtol=rtol, &
                                     & atol=atol, &
                                     & istep=0.5_qp &
                                    & )
                    res = fd_sol%ys(:,size(fd_sol%ts))
                end function fd_integrate
                function fd_integrate_jac(x) result(res)
                    type(odesolution) :: fd_sol
                    real(qp), intent(in) :: x(:)
                    real(qp)             :: res(size(x),size(x)), &
                                          & finalstate(8+8**2)
                    fd_sol = solve_ivp(fd_eoms,&
                                     & [0._qp, 1._qp], &
                                     & [x(:8), reshape(eye,[8**2])],&
                                     & dense_output=.false.,&
                                     & rtol=rtol, &
                                     & atol=atol, &
                                     & istep=0.5_qp &
                                    & )
                    finalstate = fd_sol%ys(:,size(fd_sol%ts))
                    res = reshape(finalstate(9:), [8,8])
                end function fd_integrate_jac
                function fd_eoms(me, x, y) result(res)
                    class(RungeKutta), intent(inout) :: me
                    real(qp),          intent(in)    :: x, y(:)
                    real(qp)                         :: res(size(y))
                    real(qp)                         :: stm(8,8), stmdot(8,8), &
                                                        jac(8,8), hes(8,8,8), &
                                                        stt(8,8,8), sttdot(8,8,8), &
                                                        acc(8)

                    gq%dynmod%tgt_on_rails = .false.
                    gq%dynmod%state = y(:8)
                    call gq%dynmod%get_derivs(y(7), acc, jac, hes)
                    if (size(y)==8) then
                        res = acc
                    else if (size(y)==8 + 8**2) then
                        res(:8) = acc
                        stm = reshape(y(9:),[8,8])
                        stmdot = matmul(jac,stm)
                        res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
                    else
                        res(:8) = acc
                        stm = reshape(y(9:8+8**2),[8,8])
                        stt = reshape(y(9+8**2:), [8,8,8])
                        stmdot = matmul(jac,stm)
                        sttdot = mattens(jac,stt,8) + quad(stm,hes,8)
                        res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
                        res(9+8**2:) = reshape(sttdot, [8**3])
                    endif
                end function fd_eoms
        end subroutine
end module test_makemodel
