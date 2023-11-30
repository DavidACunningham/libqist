program main
    use, intrinsic :: iso_fortran_env, only: wp=>real64
    use findiffmod
    use makemodel
    implicit none
    type(dynamicsmodel) :: dummy
    real(wp)            :: mu, r(6), rbods(3), &
                         & true_jac(6,6), fd_jac(6,6), &
                         & true_hes(6,6,6), fd_hes(6,6,6)
    mu = 1._wp
    call random_number(r)
    call random_number(rbods)
    rbods = 100*rbods
    true_jac = dummy%jac_nbody(mu,r,rbods)
    true_hes = dummy%hes_nbody(mu,r,rbods)
    fd_jac = findiff(accwrap, r, 1.e-3_wp,9)
    fd_hes = findiffhes(jacwrap, r, 1.e-3_wp,9)
    print *, "___JACOBIAN___"
    print *, norm2(true_jac - fd_jac)
    print *, "___HESSIAN___"
    print *, norm2(true_hes - fd_hes)

    contains
        function accwrap(rprime) result(res)
            real(wp), intent(in) :: rprime(:)
            real(wp)             :: res(size(rprime))

            res = dummy%acc_nbody(mu,rprime,rbods)

        end function accwrap
        function jacwrap(rprime) result(res)
            real(wp), intent(in) :: rprime(:)
            real(wp)             :: res(size(rprime),size(rprime))

            res = dummy%jac_nbody(mu,rprime,rbods)

        end function jacwrap
end program main

