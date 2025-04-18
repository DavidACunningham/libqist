        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  8 12:22:29 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHCALL_V__genmod
          INTERFACE 
            FUNCTION CHCALL_V(A,B,COEFFS,X)
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B
              REAL(KIND=8), INTENT(IN) :: COEFFS(:)
              REAL(KIND=8), INTENT(IN) :: X(:)
              REAL(KIND=8) :: CHCALL_V(SIZE(X))
            END FUNCTION CHCALL_V
          END INTERFACE 
        END MODULE CHCALL_V__genmod
