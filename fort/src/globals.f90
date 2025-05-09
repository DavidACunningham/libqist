! Title: globals.f90 
! Description:
!     Module containing global functions, constants, 
!     and subroutines for the QIST package
!
! References:
!     None
! 
! author: David Cunningham
! Last edited: See git log
module globals
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
                          ! sttu = STT Upper index, unpacked
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
    subroutine writeLog(message)
        character(len=1000), intent(in) :: message
        integer                         :: io, date_time(8)
        logical                         :: dasein
        character(len=50)               :: logfile
        character(len=12)               :: real_clock(3)
        logfile = 'QISTLog.txt'
        call date_and_time( &
                           real_clock(1), &
                           real_clock(2), &
                           real_clock(3), &
                           date_time      &
                          )
        inquire(file='QISTLog.txt', exist=dasein)
        if (dasein) then
            open( &
               & newunit=io, &
               & file=trim(adjustl(logfile)),  &
               & position='append', &
               & status = 'old', &
               & action= 'write' &
               & )
        else
            open( &
               & newunit=io, &
               & file=trim(adjustl(logfile)),  &
               & status = 'replace', &
               & action= 'write' &
               & )
        endif
        write(io,*) real_clock(1)//":"//real_clock(2)//":"//real_clock(3)//":: "//trim(adjustl(message))
        close(io)
    end subroutine writeLog
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
end module globals
