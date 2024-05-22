program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use genqist, only: model_accuracy_check, gqist
    implicit none
    type(gqist) :: gq
    integer, parameter :: npoints=500
    real(dp), parameter :: C(2,2) = 0._dp, S(2,2) = 0._dp
    integer     :: num_bodies, current_bodylist(30)
    real(dp)    :: current_mulist(30), &
                 & spicepoints(6,npoints), &
                 & testpoints(6,npoints), testtimes(npoints)

    current_bodylist = 0
    current_mulist = 0._dp
    call gq%init("./qist_config.nml")

    num_bodies = gq%dynmod%num_bodies
    current_bodylist(:num_bodies) = gq%dynmod%bodylist
    current_mulist(:num_bodies)   = gq%dynmod%nbody_mus

    call model_accuracy_check(gq, real(gq%t0,dp), real(gq%tf,dp), &
                              current_bodylist(:num_bodies), &
                              current_mulist(:num_bodies), &
                              gq%dynmod%shgrav, &
                              C, &
                              S, &
                              testpoints, &
                              spicepoints, &
                              testtimes &
                             )

    print *, "RMS POSITION DIFFERENCE (KM): "
    print *, real(rms(norm2(testpoints(:3,:)-spicepoints(:3,:),dim=1)),4)
    print *, "RMS VELOCITY DIFFERENCE (KM/S): "
    print *, real(rms(norm2(testpoints(3:6,:)-spicepoints(3:6,:),dim=1)),4)


    contains
        function rms(val) result(res)
            real(dp), intent(in) :: val(:)
            integer              :: sz
            real(dp)             :: res
            sz = size(val)
            res = sqrt(sum(val)/sz)
        end function rms


end program main
