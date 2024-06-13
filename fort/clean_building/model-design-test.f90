program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use genqist, only: model_accuracy_check, gqist, load_gravity_model
    use quat, only: rothist
    use util, only: print_to_file
    implicit none
    type(gqist) :: gq
    type(rothist) :: rot
    integer, parameter :: npoints=500
    real(qp), allocatable :: Cbar(:,:), Sbar(:,:)
    real(qp)              :: ref_radius, mu
    integer     :: num_bodies, current_bodylist(30)
    character(len=1000) :: arg, gravfile
    real(dp)    :: current_mulist(30), &
                 & spicepoints(6,npoints), &
                 & testpoints(6,npoints), testtimes(npoints)

    current_bodylist = 0
    current_mulist = 0._dp
    call get_command_argument(1,arg)
    call gq%init(trim(adjustl(arg)))

    gravfile = "/home/david/wrk/nstgro/qist/libqist/fort/data/moon8by8.nml"
    call load_gravity_model(gravfile, ref_radius, mu, Cbar, Sbar, rot)
    num_bodies = gq%dynmod%num_bodies
    current_bodylist(:num_bodies) = gq%dynmod%bodylist
    current_mulist(:num_bodies)   = real(gq%dynmod%nbody_mus,dp)

    call model_accuracy_check(gq, real(gq%t0,dp), real(gq%tf,dp), &
                              current_bodylist(:num_bodies), &
                              current_mulist(:num_bodies), &
                              gq%dynmod%shgrav, &
                              real(Cbar,dp), &
                              real(Sbar,dp), &
                              testpoints, &
                              spicepoints, &
                              testtimes &
                             )

    call print_to_file("testtimes",testtimes)
    call print_to_file("spicex", spicepoints(1,:))
    call print_to_file("spicey", spicepoints(2,:))
    call print_to_file("spicez", spicepoints(3,:))
    call print_to_file("spicexdot", spicepoints(4,:))
    call print_to_file("spiceydot", spicepoints(5,:))
    call print_to_file("spicezdot", spicepoints(6,:))
    call print_to_file("testx", testpoints(1,:))
    call print_to_file("testy", testpoints(2,:))
    call print_to_file("testz", testpoints(3,:))
    call print_to_file("testxdot", testpoints(4,:))
    call print_to_file("testydot", testpoints(5,:))
    call print_to_file("testzdot", testpoints(6,:))
    print *, "RMS POSITION DIFFERENCE (KM): "
    print *, real(rms(norm2(testpoints(:3,:)-spicepoints(:3,:),dim=1)),4)
    print *, "RMS VELOCITY DIFFERENCE (KM/S): "
    print *, real(rms(norm2(testpoints(4:6,:)-spicepoints(4:6,:),dim=1)),4)

    print *, "RMS POSITION DIFFERENCE (ND): "
    print *, real(rms(norm2(testpoints(:3,:)-spicepoints(:3,:),dim=1)/norm2(spicepoints(:3,:),dim=1)),4)
    print *, "RMS VELOCITY DIFFERENCE (ND): "
    print *, real(rms(norm2(testpoints(4:6,:)-spicepoints(4:6,:),dim=1)/norm2(spicepoints(4:6,:),dim=1)),4)
    print *, ""
    print *, ""
    print *, "MEAN POSITION DIFFERENCE (KM): "
    print *, real(mean(norm2(testpoints(:3,:)-spicepoints(:3,:),dim=1)),4)
    print *, "MEAN VELOCITY DIFFERENCE (KM/S): "
    print *, real(mean(norm2(testpoints(4:6,:)-spicepoints(4:6,:),dim=1)),4)

    print *, "MEAN POSITION DIFFERENCE (ND): "
    print *, real(mean(norm2(testpoints(:3,:)-spicepoints(:3,:),dim=1)/norm2(spicepoints(:3,:),dim=1)),4)
    print *, "MEAN VELOCITY DIFFERENCE (ND): "
    print *, real(mean(norm2(testpoints(4:6,:)-spicepoints(4:6,:),dim=1)/norm2(spicepoints(4:6,:),dim=1)),4)

    contains
        function rms(val) result(res)
            real(dp), intent(in) :: val(:)
            integer              :: sz
            real(dp)             :: res
            sz = size(val)
            res = sqrt(sum(val**2)/sz)
        end function rms

        function mean(val) result(res)
            real(dp), intent(in) :: val(:)
            integer              :: sz
            real(dp)             :: res
            sz = size(val)
            res = sum(val)/sz
        end function mean

end program main
