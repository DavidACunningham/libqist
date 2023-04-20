module globals
     ! use, intrinsic :: iso_c_binding
     use, intrinsic :: ieee_arithmetic, only: ieee_value, &
                       ieee_is_finite, ieee_is_nan, ieee_negative_inf, &
                       ieee_positive_inf, ieee_quiet_nan, ieee_copy_sign,&
                       ieee_is_normal
    use, intrinsic :: iso_fortran_env, only: wp =>real128, dp =>real64
    implicit none
    integer, parameter            :: n=6, statesize=n+n**2+n**3,&
                                      plen = n+n*(n-1) + n**2*(n-1)/2, &
                                      stmlp = n+1, stmup = n+n*(n-1), &
                                      sttlp = n+n*(n-1)+1, &
                                      sttup = plen, &
                                      stml = n+1, stmu = n+n**2, &
                                      sttl = n+n**2+1, &
                                      sttu = n+n**2+n**3
    contains

    pure function mmult(matA, matB) result(res)
        real(wp), intent(in) :: matA(:,:), matB(:,:)
        real(wp)             :: res(size(matA,1),size(matB,2))
        res = matmul(matA,matB)
    end function mmult
    pure function qesolve(a,b,c) result(res)
        !! Returns roots of a quadratic equation
        !! Using stable solution method
        !! returns result as [re(x1), im(x1), re(x2), im(x2)]
        real(wp), intent(in) :: a,b,c
        complex(wp)          :: x1, x2, d
        real(wp)             :: res(4)

        d = sqrt(cmplx(b*b-4._wp*a*c))
        if (b.ge.0) then
            ! use b positive method
            x1 = -(b+d)/(2._wp*a)
            x2 = -2._wp*c/(b+d)
        else
            x1 = 2._wp*c/(d-b)
            x2 = (d-b)/(2._wp*a)
        end if
        res = [real(x1), aimag(x1), real(x2), aimag(x2)]
    end function qesolve

end module globals
