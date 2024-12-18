module test_genqist
    use genqist, only: gqist, &
                       make_spice_subset, &
                       model_accuracy_check, &
                       make_rotation, &
                       generate_kernel, &
                       make_qist_model, &
                       configdata
    use subspice, only: spice_subset
    use tensorops, only: quad, mattens
    use quat, only: rothist
    use frkmin_q, only: odesolution, rungekutta, solve_ivp
    use qist, only: itraj
    use test_util, only: mprint, tprint
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    implicit none
    real(dp), parameter :: dtol = 1.e-12_dp
    real(qp), parameter :: qtol = 1.e-22_qp

    contains
        subroutine test_make_spice_subset(testpass)
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(spice_subset)     :: subspice
            logical                :: testvec(3)
            integer num, stat
            testpass = .true.
            call make_spice_subset("./test_config_namelists.nml")
            open(file="./test_resample.subspice", iostat=stat, &
                 status="old", access="stream", newunit=num)
             call subspice%read(num)
             close(num)
             testvec(1) = abs(subspice%a-772509669.184287_dp).le.dtol
             testvec(2) = abs(subspice%b-772682469.184231_dp).le.dtol
             testvec(3) = subspice%central_body.eq.301
             msg = "PASS Make SPICE subset test PASS"
             testpass = all(testvec)
             if (.not.testpass) msg = "FAIL Make SPICE subset test FAIL"
             write (*,*) trim(msg)
        end subroutine
        subroutine test_make_rotation(testpass)
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(rothist)          :: rot
            logical                :: testvec(4)
            integer num, stat
            real(dp)               :: spice_dcm(3,3), test_time, test_dcm(3,3)
            testpass = .true.
            call make_rotation("./test_config_namelists.nml")
            open(file="./test_rot.rot", iostat=stat, &
                 status="old", access="stream", newunit=num)
             call rot%read(num)
             close(num)
             testvec(1) = abs(rot%t0-772509669.184287_dp).le.dtol
             testvec(2) = abs(rot%tf-772682469.184231_dp).le.dtol
             testvec(3) = rot%degree.eq.400
             testpass = all(testvec)
             test_time = (rot%tf -rot%t0)/2._dp + rot%t0
             call furnsh("/home/david/wrk/nstgro/qist/kernels/mk_test.tf")
             call pxfrm2("J2000", "MOON_PA", rot%t0, test_time, spice_dcm)
             test_dcm = rot%call(test_time)
             spice_dcm = transpose(spice_dcm)
             testvec(4) = all(abs(spice_dcm-test_dcm).le.dtol)
             msg = "PASS Make rotation history test PASS"
             testpass = all(testvec)
             if (.not.testpass) msg = "FAIL Make rotation history test FAIL"
             write (*,*) trim(msg)
        end subroutine
        subroutine test_grav_load(testpass)
            type(configdata)       :: cd
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(rothist)          :: rot
            logical                :: testvec(4)
            real(qp), allocatable  :: C(:,:), S(:,:)
            real(qp)               :: mu, ref_radius
            testpass = .true.
            call cd%init("./test_config_namelists.nml")
            testvec(1) = abs(cd%central_body_ref_radius- 0.1738E+04_qp).le.qtol
            testvec(2) = abs(cd%central_body_mu-0.4902799806931690E+04_qp).le.qtol
            testvec(3) = abs(cd%cbar(0,0)-1._qp).le.qtol
            testvec(4) = abs(cd%sbar(2,1)-9.7726994478962992E-10_qp).le.qtol
            testpass = all(testvec)
             msg = "PASS Gravity loading test PASS"
             testpass = all(testvec)
             if (.not.testpass) msg = "FAIL Gravity loading test FAIL"
             write (*,*) trim(msg)
        end subroutine
        subroutine test_gq_init(testpass)
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(gqist)            :: gq
            logical                :: testvec(2)
            testpass = .true.
            call gq%init("./test_config_namelists.nml")
            testvec(1) = abs(gq%dynmod%central_body_ref_radius- 0.1738E+04_qp).le.qtol
            testvec(2) = abs(gq%dynmod%central_body_mu-0.4902799806931690E+04_qp).le.qtol
            testpass = all(testvec)
            msg = "PASS Dynamics model loading test PASS"
            testpass = all(testvec)
            if (.not.testpass) msg = "FAIL Dynamics model loading test FAIL"
            write (*,*) trim(msg)
        end subroutine
        subroutine test_kernel_write(testpass)
            logical, intent(inout) :: testpass
            type(spice_subset) :: testspice
            real(dp)           :: time = 772509669.18429_dp, &  ! 2024 Jun 24 14:00 UTC
                                  trajstate(6), earthstate(6), lt_dum, test_trajstate(6), test_earthstate(6)
            integer stat, num
            testpass =.false.
            call generate_kernel("./test_config_namelists.nml")
            call make_spice_subset("./test_config_namelists.nml",.true.)
            open(file="./test_resample_withtraj.subspice", iostat=stat, &
                 status="old", access="stream", newunit=num)
             call testspice%read(num)
             close(num)
            call furnsh("/home/david/wrk/nstgro/qist/kernels/mk_test_withorbit.tf")
            test_trajstate(:3) = testspice%call(time,-999,"p")
            test_trajstate(4:) = testspice%call(time,-999,"v")
            test_earthstate(:3) = testspice%call(time,399,"p")
            test_earthstate(4:) = testspice%call(time,399,"v")

            call spkgeo(-999, &
                        time, &
                        "J2000", &
                        301, &
                        trajstate, &
                        lt_dum &
                       )
            call spkgeo(399, &
                        time, &
                        "J2000", &
                        301, &
                        earthstate, &
                        lt_dum &
                       )

            print *, test_trajstate - trajstate
            print *, test_earthstate - earthstate
            testpass = .true.
        end subroutine
        subroutine test_make_qist(testpass)
            type(odesolution)   :: base_sol !, qistsol
            type(itraj)         :: it
            type(gqist)         :: gq
            real(qp)            :: tof
            real(qp)            :: init_state(8), eye(8,8), &
                                 & analytic_stm(8,8), analytic_stt(8,8,8), & 
                                 & divisor(8,8), init_stt(8**3), &
                                 & stt_divisor(8,8,8), test_t0, test_tf, &
                                 & stmdiff(8,8), sttdiff(8,8,8)
            real(dp)            :: itraj_stm0f(8,8), itraj_stt0f(8,8,8), &
                                 & itraj_stm_f(8,8), itraj_stt_f(8,8,8), &
                                 & itraj_stm_0(8,8), itraj_stt_0(8,8,8)
            integer i
            logical, intent(inout) :: testpass
            logical                :: testvec(4)
            testpass = .true.
            call make_qist_model("./test_config_namelists.nml")
            call gq%init("./test_config_namelists.nml",.true.)
            call it%init("./test_config_namelists.nml")
            call it%stts_ab(it%t0,it%tf,itraj_stm0f, itraj_stt0f)
            itraj_stm_0 = it%stm(it%t0)
            itraj_stt_0 = it%stt(it%t0)
            itraj_stm_f = it%stm(it%tf)
            itraj_stt_f = it%stt(it%tf)
            eye = 0._qp; do i=1,8; eye(i,i) = 1._qp; end do
            init_stt = 0._qp
            test_t0 = gq%t0
            test_tf = gq%tf
            tof = test_tf-test_t0
            gq%dynmod%tof = tof
            init_state = [gq%dynmod%trajstate(gq%t0), test_t0, tof]
            gq%dynmod%tgt_on_rails = .false.
            gq%dynmod%state = init_state
            gq%dynmod%regularize = .false.
            print *, "Integrating base case"
            base_sol = solve_ivp(fd_eoms,&
                               & [0._qp, 1._qp], &
                               & [init_state, &
                                  reshape(eye,[8**2]), &
                                  init_stt], &
                               & dense_output=.false.,&
                               & rtol=gq%rtol, &
                               & atol=gq%atol, &
                               & istep=0.5_qp &
                              & )
            print *, "Done."
            analytic_stm = reshape(base_sol%ys(9:8+8**2,size(base_sol%ts)), [8,8])
            analytic_stt = reshape(base_sol%ys(9+8**2:,size(base_sol%ts)), [8,8,8])
            divisor = 1._qp; where (abs(analytic_stm).ge.1.e-14_qp); divisor = analytic_stm; end where
            stt_divisor = 1._qp; where (abs(analytic_stt).ge.1.e-14_qp); stt_divisor = analytic_stt; end where

            stmdiff = (analytic_stm- real(itraj_stm_f,qp))/divisor
            sttdiff = (analytic_stt- real(itraj_stt_f,qp))/stt_divisor
            testvec(1) = all(stmdiff.lt.5.e-4_qp)
            testvec(2) = all(sttdiff.lt.5.e-4_qp)
            if (.not.testvec(1)) then
                print *, "FAIL Itraj STM FAIL"
                print *, "Itraj STM"
                call mprint(itraj_stm_f)
                print *, "Analytic STM"
                call mprint(analytic_stm)
                print *, "STM difference"
                call mprint(stmdiff)
            end if
            if (.not.testvec(2)) then
                print *, "FAIL Itraj STT FAIL"
                print *, "Itraj STT"
                call tprint(itraj_stt_f)
                print *, "Analytic STT"
                call tprint(analytic_stt)
                print *, "STT difference"
                call tprint(sttdiff)
            end if
            stmdiff = (analytic_stm- real(itraj_stm0f,qp))/divisor
            sttdiff = (analytic_stt- real(itraj_stt0f,qp))/stt_divisor
            testvec(3) = all(stmdiff.lt.5.e-4_qp)
            testvec(4) = all(sttdiff.lt.5.e-4_qp)
            if (.not.testvec(3)) then
                print *, "FAIL Itraj chained STM FAIL"
                print *, "Itraj chained STM"
                call mprint(itraj_stm0f)
                print *, "Analytic Chained STM"
                call mprint(analytic_stm)
                print *, "STM difference"
                call mprint(stmdiff)
            end if
            if (.not.testvec(4)) then
                print *, "FAIL Itraj chained STT FAIL"
                print *, "Itraj chained STT"
                call tprint(itraj_stt0f)
                print *, "Analytic STT"
                call tprint(analytic_stt)
                print *, "STT difference"
                call tprint(sttdiff)
            end if
            testpass = all(testvec)
            if (.not.testpass) then
                print *, "FAIL Make QIST test FAIL"
            else
                print *, "PASS Make QIST test PASS"
                print *, "QIST regularization? ", it%regularized
            end if
            contains
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
                    res(:8) = acc
                    stm = reshape(y(9:8+8**2),[8,8])
                    stt = reshape(y(9+8**2:), [8,8,8])
                    stmdot = matmul(jac,stm)
                    sttdot = mattens(jac,stt,8) + quad(stm,hes,8)
                    res(9:8 + 8 ** 2) = reshape(stmdot,[8**2])
                    res(9+8**2:) = reshape(sttdot, [8**3])
                end function fd_eoms
        end subroutine
end module test_genqist
