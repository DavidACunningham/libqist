module test_quat
    use quat
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use findiffmod
    use test_util, only: mprint
    implicit none
    contains

        subroutine quat_ops_test(testpass)
            logical, intent(inout) :: testpass
            real(dp), parameter    :: ang1=0.775406922180032_dp, &
                                      ang2=4.533519244025973_dp
            real(dp), parameter    :: ax1(3) = [0.841804704648762_dp, &
                                                0.190555373344514_dp, &
                                                0.505028206064516_dp], &
                             ax2(3) = [0.140882407647218_dp, &
                                       0.590274029950995_dp, &
                                       0.794813636509173_dp], &
                             quata_true(4) = [0.925779737552099_dp, &
                                              0.318255485448565_dp, &
                                              0.072041997999874_dp, &
                                              0.190932642688595_dp], &
                             quatc_true(4) =  [-0.777047634387741_dp,&
                                               -0.061399134234913_dp, &
                                                0.546670929859078_dp, &
                                                0.305905238030719_dp], &
                             dcma_true(3,3) = reshape([ 0.916709352960271_dp, &
                                                       -0.307667621584484_dp, &
                                                        0.254920765775777_dp, &
                                                        0.399378665769022_dp, &
                                                        0.724516343875695_dp, &
                                                       -0.561758621460812_dp, &
                                                       -0.011859322228432_dp, &
                                                        0.616779297711542_dp, &
                                                        0.787046793012169_dp],&
                                                        [3,3],order=[2,1]), &
                            dtol = 1.e-14_dp
            type(quaternion)       :: quata, quatb, quatc, quatd

            call quata%axang(ax1,ang1)
            call quatb%axang(ax2,ang2)
            call quatd%fromdcm(dcma_true)

            quatc = quatb*quata

            testpass = .true.
            if (norm2(quata%q-quata_true).ge.dtol) then
                testpass = .false.
                write (*,*) "FAIL Quaternion axang initialization FAIL"
                write (*,*) "Error: ", quata%q - quata_true
            endif
            if (norm2(quatd%q-quata_true).ge.dtol) then
                write(*,*) "FAIL Quaternion DCM initialization FAIL"
                write (*,*) "Error: ", quatd%q - quata_true
            endif
            if (norm2(quatc%q-quatc_true).ge.dtol) then
                write(*,*) "FAIL Quaternion composition FAIL"
                write (*,*) "Error: ", quatc%q - quatc_true
            endif
            if (norm2(quata%asdcm()-dcma_true).ge.dtol) then
                write(*,*) "FAIL Quaternion to DCM conversion FAIL"
                write (*,*) "Error: "
                call mprint(quata%asdcm())
                print *, ""
                call mprint(dcma_true)
                print *, ""
                call mprint(quata%asdcm() - dcma_true)
            endif
            if (testpass) write (*,*) "PASS Basic quaternion operations tests PASS"
        end subroutine quat_ops_test

        subroutine rot_hist_test(testpass)
            logical, intent(inout) :: testpass
            integer, parameter :: nnodes=80
            real(dp), parameter :: pi = asin(1._dp)*2._dp, &
                                   dtol = 1.e-13_dp
            real(dp), dimension(3,3) :: rotmat
            real(dp)                 :: t0, tf
            real(dp)                 :: qdot(4), qdot_fd(4), &
                                        dcmdot(3,3), dcmddot(3,3), &
                                        dcmdot_fd(3,3), dcmddot_fd(3,3), &
                                        qdum(4,1), midpoint_time
            type(rothist)            :: rot

            t0 = 0.12_dp
            tf = 4._dp*pi
            midpoint_time = (tf-t0)/2._dp + t0
            rotmat(1,:) = [1._dp,         0._dp,         0._dp]
            rotmat(2,:) = [0._dp, cos(pi/4._dp), sin(pi/4._dp)]
            rotmat(3,:) = [0._dp,-sin(pi/4._dp), cos(pi/4._dp)]
            call rot%init(fitfun, t0, tf, nnodes, rotmat)

            dcmdot     = rot%calldot(midpoint_time)
            dcmddot    = rot%callddot(midpoint_time)
            qdot       = rot%callqdot(midpoint_time)
            qdum       = real(findiffmat_td(fd_quat, real(midpoint_time,qp), &
                              (tf-t0)*1.e-3_qp, 9,real(qdum,qp) ),dp)
            qdot_fd    = qdum(:,1)
            dcmdot     =  rot%calldot(real(midpoint_time,dp))
            dcmddot    = rot%callddot(real(midpoint_time,dp))
            dcmdot_fd  = real(findiffmat_td(fd_rot, real(midpoint_time,qp), &
                              (tf-t0)*1.e-3_qp, 9, real(dcmdot,qp)),dp)
            dcmddot_fd = real(findiffmat_td(fd_rotdot, real(midpoint_time,qp), &
                              (tf-t0)*1.e-3_qp, 9, real(dcmdot,qp)),dp)

            testpass = .true.
            if (any(abs(qdot-qdot_fd).ge.dtol)) then
                testpass = .false.
                write (*,*) "FAIL Rothist quaternion derivative FAIL"
                write (*,*) "Error: ", (qdot-qdot_fd)
            endif
            if (any(abs(dcmdot_fd-dcmdot).ge.dtol)) then
                testpass = .false.
                write(*,*) "FAIL Rothist DCM derivative FAIL"
                write (*,*) "Error: "
                call mprint(dcmdot_fd - dcmdot)
            endif
            if (any(abs(dcmddot_fd-dcmddot).ge.dtol)) then
                testpass = .false.
                write(*,*) "FAIL Rothist DCM second derivative FAIL"
                write (*,*) "Error: "
                call mprint(dcmddot_fd - dcmddot)
            endif
            if (testpass) write (*,*) "PASS Rothist tests PASS"
            contains
            function fitfun(me, ta,tb) result(res)
                class(rothist), intent(inout) :: me
                real(dp), intent(in)          :: ta, tb
                real(dp)                      :: res(4), mat(3,3), x
                type(quaternion)              :: qclass
                x = tb
                mat(1,:) = [ cos(x*pi), sin(x*pi), 0._dp]
                mat(2,:) = [-sin(x*pi), cos(x*pi), 0._dp]
                mat(3,:) = [     0._dp,     0._dp, 1._dp]
                call qclass%fromdcm(mat)
                res = qclass%q
            end function
            function fd_quat(x) result(res)
                real(qp), intent(in) :: x
                real(qp), allocatable :: res(:,:)
                allocate(res(4,1))
                res(:,1) = real(rot%callq(real(x,dp)),qp)
            end function fd_quat
            function fd_rot(x) result(res)
                real(qp), intent(in) :: x
                real(qp), allocatable :: res(:,:)
                allocate(res(3,3))
                res = real(rot%call(real(x,dp)),qp)
            end function fd_rot
            function fd_rotdot(x) result(res)
                real(qp), intent(in) :: x
                real(qp), allocatable :: res(:,:)
                allocate(res(3,3))
                res = real(rot%calldot(real(x,dp)),qp)
            end function fd_rotdot
        end subroutine rot_hist_test
end module test_quat
