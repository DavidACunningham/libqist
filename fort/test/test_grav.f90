program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use makemodel, only: dynamicsModel
    use quat, only: rothist, quaternion
    use cheby, only: spice_subset
    implicit none
    real(dp), parameter :: epoch=10000._dp, epoch2=336000._dp
    integer, parameter :: nnodes=300
    real(dp), dimension(3,3) :: rotmat_comp
    real(qp)                 :: rtol, atol, t0, tf
    real(qp)                 :: y(8), acc(8), jac(8,8), hes(8,8,8)
    type(rothist)            :: rot
    type(dynamicsModel)      :: dyn
    type(spice_subset)       :: subspice
    integer num, i, j
    !! ALL THAT'S NEEDED TO SETUP ROTS
    open(file="/home/david/wrk/nstgro/qist/libqist/fort/data/20240524_gw_resample.subspice", &
         newunit=num)
    call subspice%read(num)
    close(num)
    rtol = 1.e-13
    atol = 1.e-20
    t0 = 769269009.185
    tf = 769845939.185
    call furnsh("/home/david/wrk/nstgro/qist/kernels/mk.tf")
    call pxform('MOON_PA','J2000',epoch,rotmat_comp)
    call rot%init(fitfun, real(t0,dp), real(tf,dp), nnodes, rotmat_comp)

    call dyn%init(subspice, &
                  -60000, &
                  301, &
                  [399, 10], &
                  4902.800118_qp, &
                  1740._qp, &
                  [398600.5_qp, 1.32712440041279419e11_qp], &
                  .true., &
                  .true., &
                  rot, &
                  reshape([1._qp, 0._qp, 0._qp, 0._qp],[2,2]), &
                  reshape([0._qp, 0._qp, 0._qp, 0._qp],[2,2]) &
                 )
    !! ALL THAT'S NEEDED TO SETUP ROTS

    y = [dyn%trajstate(t0), t0, 1._qp]
    call dyn%allderivs_sh(t0, y, acc, jac, hes)
    do i=1,8
        print *, real(jac(i,:),4)
    end do
    contains
    function fitfun(me, ta,tb) result(res)
        class(rothist), intent(inout) :: me
        real(dp), intent(in)          :: ta, tb
        real(dp)                      :: res(4), mat(3,3)
        type(quaternion)              :: qclass
        call pxfrm2('MOON_PA','MOON_PA',ta,tb,mat)
        call qclass%fromdcm(mat)
        res = qclass%q
    end function
end program main
