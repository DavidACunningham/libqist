module test_genqist
    use genqist, only: gqist, &
                       make_spice_subset, &
                       load_gravity_model, &
                       model_accuracy_check, &
                       make_rotation, &
                       generate_kernel, &
                       make_qist_model
    use cheby, only: spice_subset
    use quat, only: rothist
    use frkmin_q, only: odesolution, rungekutta, solve_ivp
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
             call furnsh("../../../kernels/mk_test.tf")
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
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(rothist)          :: rot
            logical                :: testvec(4)
            real(qp), allocatable  :: C(:,:), S(:,:)
            real(qp)               :: mu, ref_radius
            testpass = .true.
            call load_gravity_model("./test_config_namelists.nml", &
                                    ref_radius,mu, C, S, rot)
            testvec(1) = abs(ref_radius- 0.1738E+04_qp).le.qtol
            testvec(2) = abs(mu-0.4902799806931690E+04_qp).le.qtol
            testvec(3) = abs(C(1,1)-1._qp).le.qtol
            testvec(4) = abs(S(3,2)-9.7726994478962992E-10_qp).le.qtol
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
            call make_spice_subset("./test_config_namelists_genqist.nml")
            open(file="./test_resample_withtraj.subspice", iostat=stat, &
                 status="old", access="stream", newunit=num)
             call testspice%read(num)
             close(num)
            call furnsh("../../../kernels/mk_test_withorbit.tf")
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
            logical, intent(inout) :: testpass
            testpass = .false.
            call make_qist_model("./test_config_namelists_genqist.nml")
            testpass = .true.
        end subroutine
end module test_genqist
