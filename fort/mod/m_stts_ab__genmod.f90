        !COMPILER-GENERATED INTERFACE MODULE: Sat Apr  5 16:10:26 2025
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE M_STTS_AB__genmod
          INTERFACE 
            SUBROUTINE M_STTS_AB(TAUA,TAUB,STM_O,STT_O)
              REAL(KIND=8), INTENT(IN) :: TAUA
              REAL(KIND=8), INTENT(IN) :: TAUB
              REAL(KIND=8), INTENT(OUT) :: STM_O(8,8)
              REAL(KIND=8), INTENT(OUT) :: STT_O(8,8,8)
            END SUBROUTINE M_STTS_AB
          END INTERFACE 
        END MODULE M_STTS_AB__genmod
