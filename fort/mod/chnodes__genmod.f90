        !COMPILER-GENERATED INTERFACE MODULE: Sat Apr  5 16:10:14 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CHNODES__genmod
          INTERFACE 
            FUNCTION CHNODES(N,A,B) RESULT(RES)
              INTEGER(KIND=4), INTENT(IN) :: N
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: B
              REAL(KIND=8) :: RES(N)
            END FUNCTION CHNODES
          END INTERFACE 
        END MODULE CHNODES__genmod
