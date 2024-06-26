module solveke
    implicit none
    contains
        function solve(m, ecc) result(res)
            real(kind=8), intent(in) :: m, ecc
            real(kind=8)             :: res
            real(kind=16)            :: thisE, miss, nextE, mwork, eccwork
            mwork = real(m,16)
            eccwork = real(ecc,16)
            thisE = mwork
            miss = 1._16
            do while ( miss .gt. 1.e-20_16)
                nextE = thisE - (thisE - eccwork*sin(thisE) - mwork)/(1-eccwork*cos(thisE))
                miss = nextE - eccwork*sin(nextE) - mwork
                thisE = nextE
            end do
            res = real(thisE, 8)
        end function solve
end module solveke
