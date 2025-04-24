        !COMPILER-GENERATED INTERFACE MODULE: Tue Apr  8 12:22:29 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHFIT_S__genmod
          INTERFACE 
            FUNCTION CHFIT_S(N,FI) RESULT(RES)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: FI(N)
              REAL(KIND=8) :: RES(N)
            END FUNCTION CHFIT_S
          END INTERFACE 
        END MODULE CHFIT_S__genmod
