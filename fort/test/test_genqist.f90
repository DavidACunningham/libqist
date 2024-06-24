module test_genqist
    use genqist, only: gqist, &
                       make_spice_subset, &
                       load_gravity_model, &
                       model_accuracy_check, &
                       make_rotation
    use cheby, only: spice_subset
    use quat, only: rothist
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
            call make_spice_subset("./test_resamp_config.nml")
            open(file="./test_resample.subspice", iostat=stat, &
                 status="old", access="stream", newunit=num)
             call subspice%read(num)
             close(num)
             testvec(1) = abs(subspice%a-772509669.184287_dp).le.dtol
             testvec(2) = abs(subspice%b-775101669.183457_dp).le.dtol
             testvec(3) = subspice%central_body.eq.301
             msg = " Make SPICE subset test pass."
             testpass = all(testvec)
             if (.not.testpass) msg = " Make SPICE subset test fail."
             write (*,*) trim(msg)
        end subroutine
        subroutine test_make_rotation(testpass)
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(rothist)          :: rot
            logical                :: testvec(3)
            integer num, stat
            testpass = .true.
            call make_rotation("./test_rot_config.nml")
            open(file="./test_rot.rot", iostat=stat, &
                 status="old", access="stream", newunit=num)
             call rot%read(num)
             close(num)
             testvec(1) = abs(rot%t0-772509669.184287_dp).le.dtol
             testvec(2) = abs(rot%tf-775101669.183457_dp).le.dtol
             testvec(3) = rot%degree.eq.400
             testpass = all(testvec)
             msg = " Make rotation history test pass."
             testpass = all(testvec)
             if (.not.testpass) msg = " Make rotation history test fail."
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
            call load_gravity_model("./test_grav_config.nml", &
                                    ref_radius,mu, C, S, rot)
            testvec(1) = abs(ref_radius- 0.1738E+04_qp).le.qtol
            testvec(2) = abs(mu-0.4902799806931690E+04_qp).le.qtol
            testvec(3) = abs(C(1,1)-1._qp).le.qtol
            testvec(4) = abs(S(3,2)-9.7726994478962992E-10_qp).le.qtol
            testpass = all(testvec)
             msg = " Gravity loading test pass."
             testpass = all(testvec)
             if (.not.testpass) msg = " Gravity loading test fail."
             write (*,*) trim(msg)
        end subroutine

        subroutine test_gq_init(testpass)
            character(len=175)     :: msg
            logical, intent(inout) :: testpass
            type(gqist)            :: gq
            logical                :: testvec(2)
            testpass = .true.
            call gq%init("./test_dyn_config_notraj.nml")
            testvec(1) = abs(gq%dynmod%central_body_ref_radius- 0.1738E+04_qp).le.qtol
            testvec(2) = abs(gq%dynmod%central_body_mu-0.4902799806931690E+04_qp).le.qtol
            testpass = all(testvec)
            msg = " Dynamics model loading test pass."
            testpass = all(testvec)
            if (.not.testpass) msg = " Dynamics model loading test fail."
            write (*,*) trim(msg)
        end subroutine

end module test_genqist
