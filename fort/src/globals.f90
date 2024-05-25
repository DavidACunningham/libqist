module globals
     ! use, intrinsic :: iso_c_binding
     use, intrinsic :: ieee_arithmetic, only: ieee_value, &
                       ieee_is_finite, ieee_is_nan, ieee_negative_inf, &
                       ieee_positive_inf, ieee_quiet_nan, ieee_copy_sign, &
                       ieee_is_normal
    use, intrinsic :: iso_fortran_env, only: wp => real128, dp => real64
    implicit none
    ! Convenience variables for indexing dynamical state and packed state
    integer, parameter :: n=8, &
                          statesize=1+n**2+n**3,&
                          ! plen = Packed state LENgth
                          plen = 1 + n*(n-2) + (n-2)*(n*(n+1)/2), &
                          ! stmlp = STM Lower index, Packed
                          stmlp = 1 + 1, &
                          ! stmup = STM Upper index, Packed
                          stmup = 1 + n*(n-2), &
                          ! sttlp = STT Lower index, Packed
                          sttlp = 1 + n*(n-2)+ 1, &
                          ! sttup = STT Upper index, Packed
                          sttup = plen, &
                          ! stml = STM Lower index, unpacked
                          stml = 1+1, & 
                          ! stmu = STM Upper index, unpacked
                          stmu = 1+n**2, &
                          ! sttl = STT Lower index, unpacked
                          sttl = 1+n**2+1, &
                          sttu = 1+n**2+n**3
    interface mmult
        ! Wraps interface to locally-defined matrix multiplication functions
        module procedure mmult_matmat_q
        module procedure mmult_matvec_q
        module procedure mmult_matmat_d
        module procedure mmult_matvec_d
        module procedure mmult_matvec_vecfirst_q
        module procedure mmult_matvec_vecfirst_d
    end interface 
    contains
    pure function mmult_matmat_q(matA, matB) result(res)
        ! mmult_matmat_q: function handling matrix-matrix multiplication
        ! INPUTS:
        ! NAME    TYPE          COMMENTS
        ! matA    real128 (:,:) left matrix
        ! matB    real128 (:,:) right matrix
        ! OUTPUTS:
        ! NAME    TYPE          COMMENTS
        ! res     real128 (:,:) output matrix
        real(wp), intent(in) :: matA(:,:), matB(:,:)
        real(wp)             :: res(size(matA,1),size(matB,2))
        res = matmul(matA,matB)
    end function mmult_matmat_q
    pure function mmult_matvec_q(matA, matB) result(res)
        ! mmult_matvec_q: function handling matrix-vector multiplication
        ! INPUTS:
        ! NAME    TYPE          COMMENTS
        ! matA    real128 (:,:) left matrix
        ! matB    real128 (:)   vector
        ! OUTPUTS:
        ! NAME    TYPE          COMMENTS
        ! res     real128 (:)   output matrix
        real(wp), intent(in) :: matA(:,:), matB(:)
        real(wp)             :: res(size(matB))
        res = matmul(matA,matB)
    end function mmult_matvec_q
    pure function mmult_matvec_vecfirst_q(matA, matB) result(res)
        ! mmult_matvec_vecfirst_q: function handling matrix-vector multiplication
        ! INPUTS:
        ! NAME    TYPE          COMMENTS
        ! matA    real128 (:) vector
        ! matB    real128 (:,:)   right matrix
        ! OUTPUTS:
        ! NAME    TYPE          COMMENTS
        ! res     real128 (:)   output matrix
        real(wp), intent(in) :: matA(:), matB(:,:)
        real(wp)             :: res(size(matB,2))
        res = matmul(matA,matB)
    end function mmult_matvec_vecfirst_q
    pure function mmult_matmat_d(matA, matB) result(res)
        ! mmult_matmat_d: function handling matrix-matrix multiplication
        ! INPUTS:
        ! NAME    TYPE          COMMENTS
        ! matA    real64 (:,:)  left matrix
        ! matB    real64 (:,:)  vector
        ! OUTPUTS:
        ! NAME    TYPE          COMMENTS
        ! res     real64 (:,:)  output matrix
        real(dp), intent(in) :: matA(:,:), matB(:,:)
        real(dp)             :: res(size(matA,1),size(matB,2))
        res = matmul(matA,matB)
    end function mmult_matmat_d
    pure function mmult_matvec_d(matA, matB) result(res)
        ! mmult_matvec_d: function handling matrix-vector multiplication
        ! INPUTS:
        ! NAME    TYPE          COMMENTS
        ! matA    real64 (:,:)  left matrix
        ! matB    real64 (:)    vector
        ! OUTPUTS:
        ! NAME    TYPE          COMMENTS
        ! res     real64 (:)    output matrix
        real(dp), intent(in) :: matA(:,:), matB(:)
        real(dp)             :: res(size(matB))
        res = matmul(matA,matB)
    end function mmult_matvec_d
    pure function mmult_matvec_vecfirst_d(matA, matB) result(res)
        ! mmult_matvec_vecfirst_d: function handling matrix-vector multiplication
        ! INPUTS:
        ! NAME    TYPE          COMMENTS
        ! matA    real64 (:) vector
        ! matB    real64 (:,:)   right matrix
        ! OUTPUTS:
        ! NAME    TYPE          COMMENTS
        ! res     real64 (:)   output matrix
        real(dp), intent(in) :: matA(:), matB(:,:)
        real(dp)             :: res(size(matB,2))
        res = matmul(matA,matB)
    end function mmult_matvec_vecfirst_d
    pure function qesolve(a,b,c) result(res)
        ! qesolve: Returns roots of a quadratic equation
        !          using stable solution method
        ! INPUTS:
        ! NAME      TYPE     COMMENTS
        ! a         real     quadratic coefficient
        ! b         real     linear coefficient
        ! c         real     constant coefficient
        ! OUTPUTS:
        ! NAME      TYPE     COMMENTS
        ! res       real (4) Solution vector, roots are r1 and r2
        !                    [Re(r1), Im(r1), Re(r2), Im(r2)]
        real(wp), intent(in) :: a,b,c
        complex(wp)          :: x1, x2, d
        real(wp)             :: res(4)

        d = sqrt(cmplx(b*b-4._wp*a*c))
        if (b.ge.0) then
            x1 = -(b+d)/(2._wp*a)
            x2 = -2._wp*c/(b+d)
        else
            x1 = 2._wp*c/(d-b)
            x2 = (d-b)/(2._wp*a)
        end if
        res = [real(x1), aimag(x1), real(x2), aimag(x2)]
    end function qesolve
end module globals
