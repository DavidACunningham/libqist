module test_makemodel
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use genqist, only: gqist
    use findiffmod
    implicit none
    real(qp), parameter :: qtol = 1.e-21_qp, qtol_tb = 1.e-17_qp, qtol_sh=2.e-14_qp

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
                                      calltime
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
            tf =  775101669.183457_qp  ! 2024 Jul 24 14:00 UTC
            calltime =  (tf-t0)/2 + t0
            call gq%init("./test_config_namelists.nml")
            gq%dynmod%tgt_on_rails = .false.
            y = [5000._qp, 1000._qp, 500._qp, 0._qp, 0._qp, 0._qp, calltime, tf-t0]
            call gq%dynmod%allderivs_kepler(gq%dynmod%central_body_mu, y, &
                                         & acc, jac, hes)
            jac_fd = 0._qp
            jac_fd = findiff_multiscale(fd_acc,y, fdsteps, 9)
            hes_fd = findiffhes_multiscale(fd_jac,y, fdsteps, 9)
            testvec(1) = all(abs(jac-jac_fd).lt.qtol)
            testvec(2) = all(abs(hes-hes_fd).lt.qtol)
            testpass = .true.
            msg = ""
            if (.not.testvec(1)) msg = trim(msg)//new_line("a")//" Keplerian Jacobian test fail."
            if (.not.testvec(2)) msg = trim(msg)//new_line("a")//" Keplerian Hessian test fail."
            testpass = all(testvec)
            if (testpass) msg = trim(" Keplerian derivatives test pass")
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
            tf =  775101669.183457_qp  ! 2024 Jul 24 14:00 UTC
            calltime =  (tf-t0)/2 + t0
            call gq%init("./test_config_namelists.nml")
            gq%dynmod%tgt_on_rails = .false.
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
            if (.not.testvec(1)) msg = trim(msg)//new_line("a")//" Third body Jacobian test fail."
            if (.not.testvec(2)) msg = trim(msg)//new_line("a")//" Third body Hessian test fail."
            testpass = all(testvec)
            if (testpass) msg = trim(" Third body derivatives test pass")
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
            tf =  775101669.183457_qp  ! 2024 Jul 24 14:00 UTC
            calltime =  (tf-t0)/2 + t0
            call gq%init("./test_config_namelists.nml")
            gq%dynmod%tgt_on_rails = .false.
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
            if (.not.testvec(1)) msg = trim(msg)//new_line("a")//" Spherical harmonics gravity Jacobian test fail."
            if (.not.testvec(2)) msg = trim(msg)//new_line("a")//" Spherical harmonics gravity Hessian test fail."
            testpass = all(testvec)
            if (testpass) msg = trim(" Third body derivatives test pass")
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
end module test_makemodel
