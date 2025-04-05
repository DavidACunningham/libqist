        !COMPILER-GENERATED INTERFACE MODULE: Sat Apr  5 16:10:14 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHDERIV_S__genmod
          INTERFACE 
            FUNCTION CHDERIV_S(COEFFS,A,B) RESULT(RES)
              REAL(KIND=8), INTENT(IN) :: COEFFS(:)
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B
              REAL(KIND=8) :: RES(SIZE(COEFFS)-1)
            END FUNCTION CHDERIV_S
          END INTERFACE 
        END MODULE CHDERIV_S__genmod
