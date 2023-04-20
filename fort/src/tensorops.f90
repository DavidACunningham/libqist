module tensorops
    use globals
    implicit none
    interface sttinvert
        module procedure q_sttinvert
        module procedure d_sttinvert
    end interface sttinvert
    interface sttchain_invert
        module procedure q_sttchain_invert
        module procedure d_sttchain_invert
    end interface sttchain_invert
    interface stminvert
        module procedure q_stminvert
        module procedure d_stminvert
    end interface stminvert
    interface sttchain
        module procedure q_sttchain
        module procedure d_sttchain
    end interface sttchain
    interface vectensquad
        module procedure q_vectensquad
        module procedure d_vectensquad
    end interface vectensquad
    interface mattens
        module procedure q_mattens
        module procedure d_mattens
    end interface mattens
    interface quad
        module procedure q_quad
        module procedure d_quad
    end interface quad
    interface vectens1
        module procedure q_vectens1
        module procedure d_vectens1
    end interface vectens1
    interface vectens3
        module procedure q_vectens3
        module procedure d_vectens3
    end interface vectens3

    contains

    pure function eyemat(n) result(res)
        implicit none
        integer, intent(in) :: n
        integer :: i
        real(wp) :: res(n,n)
        res = 0._wp
        do i=1,n
            res(i,i) = 1._wp
        end do
    end function eyemat
    pure function d_eyemat(n) result(res)
        implicit none
        integer, intent(in) :: n
        integer :: i
        real(dp) :: res(n,n)
        res = 0._dp
        do i=1,n
            res(i,i) = 1._dp
        end do
    end function d_eyemat

    pure function zmat(n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp) :: res(n,n)
        res = 0._wp
        res(:n/2,n/2:) = eyemat(n/2)
        res(n/2:,:n/2) = -eyemat(n/2)
    end function zmat
    pure function d_zmat(n) result(res)
        implicit none
        integer, intent(in) :: n
        real(8) :: res(n,n)
        res = 0._8
        res(:n/2,n/2:) = d_eyemat(n/2)
        res(n/2:,:n/2) = -d_eyemat(n/2)
    end function d_zmat

    pure function q_mattens(m,t,n) result(res)
        ! ia, ajk -> ijk
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: m(n,n), t(n,n,n)
        real(wp) :: res(n,n,n)
        res=reshape(mmult(m,reshape(t,[n,n**2])),[n,n,n])
    end function q_mattens
    pure function d_mattens(m,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n,n), t(n,n,n)
        real(dp) :: res(n,n,n)
        res=reshape(mmult(m,reshape(t,[n,n**2])),[n,n,n])
    end function d_mattens

    pure function q_quad(m,t,n) result(res)
        ! iab, aj, bk -> ijk
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: m(n,n), t(n,n,n)
        real(wp) :: res(n,n,n)
        res = reshape(t,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
    end function q_quad
    pure function d_quad(m,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n,n), t(n,n,n)
        real(dp) :: res(n,n,n)
        res = reshape(t,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
        res = mattens(transpose(m),res,n)
        res = reshape(res,[n,n,n],order=[2,3,1])
    end function d_quad

    pure function q_vectens1(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n,n)
        res = reshape(mmult(v,reshape(t,[n,n**2])),[n,n])
    end function q_vectens1
    pure function d_vectens1(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v(n), t(n,n,n)
        real(dp) :: res(n,n)
        res = reshape(mmult(v,reshape(t,[n,n**2])),[n,n])
    end function d_vectens1

    pure function vectens2(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n,n)
        res = vectens1(v,reshape(t,[n,n,n],order=[2,1,3]),n)
    end function vectens2

    pure function q_vectens3(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n,n)
        res = vectens1(v,reshape(t,[n,n,n],order=[3,2,1]),n)
    end function q_vectens3
    pure function d_vectens3(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v(n), t(n,n,n)
        real(dp) :: res(n,n)
        res = vectens1(v,reshape(t,[n,n,n],order=[3,2,1]),n)
    end function d_vectens3
    
    pure function q_vectensquad(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: v(n), t(n,n,n)
        real(wp) :: res(n)
        res = mmult(v,vectens3(v,t,n))
    end function q_vectensquad

    pure function d_vectensquad(v,t,n) result(res)
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in) :: v(n), t(n,n,n)
        real(dp) :: res(n)
        res = mmult(v,vectens3(v,t,n))
    end function d_vectensquad

    pure function q_stminvert(stm,n) result(res)
        implicit none
        real(wp), intent(in) :: stm(n,n)
        integer, intent(in) :: n 
        real(wp)             :: res(n,n)
        ! symplectic version
        res = -mmult(mmult(zmat(n),transpose(stm)),zmat(n))
    end function q_stminvert

    pure function d_stminvert(stm,n) result(res)
        implicit none
        real(8), intent(in) :: stm(n,n)
        integer, intent(in) :: n 
        real(8)             :: res(n,n)
        ! symplectic version
        res = -mmult(mmult(d_zmat(n),transpose(stm)),d_zmat(n))
    end function d_stminvert

    pure function d_sttinvert(stm,stt,n) result (res)
        implicit none
        real(8), intent(in) :: stm(n,n), stt(n,n,n)
        integer, intent(in) :: n
        real(8) :: istm(n,n), res(n,n,n), inter(n,n,n)
        istm = stminvert(stm,n)
        inter = -mattens(istm,stt,n)
        res = quad(istm,inter,n)
    end function d_sttinvert

    pure function q_sttinvert(stm,stt,n) result (res)
        implicit none
        real(wp), intent(in) :: stm(n,n), stt(n,n,n)
        integer, intent(in) :: n
        real(wp) :: istm(n,n), res(n,n,n), inter(n,n,n)
        istm = stminvert(stm,n)
        inter = -mattens(istm,stt,n)
        res = quad(istm,inter,n)
    end function q_sttinvert
        
    pure subroutine d_sttchain_invert(stm1,stt1,stm2,stt2, stmab, sttab,n)
        implicit none
        integer, intent(in) :: n
        real(8), intent(in) :: stm1(n,n), stm2(n,n), stt1(n,n,n), stt2(n,n,n)
        real(8)             :: istm1(n,n), istt1(n,n,n)
        real(8)             :: inter(n,n,n) 
        real(8), intent(out) :: stmab(n,n), sttab(n,n,n)
        istm1 = stminvert(stm1,n)
        istt1 = sttinvert(stm1,stt1,n)

        stmab = mmult(stm2,istm1)
        inter = reshape(istt1,[n,n,n],order=[2,1,3])
        sttab = reshape(mattens(stm2,inter,n),[n,n,n],order=[2,1,3])
        inter = reshape(stt2,[n,n,n],order=[2,1,3])
        inter = mattens(transpose(istm1),inter,n)
        inter = reshape(inter,[n,n,n],order=[2,1,3])
        inter = reshape(inter,[n,n,n],order=[3,2,1])
        inter = mattens(transpose(istm1),inter,n)
        sttab = reshape(inter,[n,n,n],order=[3,2,1])
    end subroutine d_sttchain_invert

    pure subroutine q_sttchain_invert(stm1,stt1,stm2,stt2, stmab, sttab,n)
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in) :: stm1(n,n), stm2(n,n), stt1(n,n,n), stt2(n,n,n)
        real(wp)             :: istm1(n,n), istt1(n,n,n)
        real(wp)             :: inter(n,n,n) 
        real(wp), intent(out) :: stmab(n,n), sttab(n,n,n)
        istm1 = stminvert(stm1,n)
        istt1 = sttinvert(stm1,stt1,n)

        stmab = mmult(stm2,istm1)
        inter = reshape(istt1,[n,n,n],order=[2,1,3])
        sttab = reshape(mattens(stm2,inter,n),[n,n,n],order=[2,1,3])
        inter = reshape(stt2,[n,n,n],order=[2,1,3])
        inter = mattens(transpose(istm1),inter,n)
        inter = reshape(inter,[n,n,n],order=[2,1,3])
        inter = reshape(inter,[n,n,n],order=[3,2,1])
        inter = mattens(transpose(istm1),inter,n)
        sttab = reshape(inter,[n,n,n],order=[3,2,1])
    end subroutine q_sttchain_invert

    pure subroutine q_sttchain(istm,istt,stm,stt,stmab,sttab,n) 
        implicit none
        integer, intent(in) :: n
        real(wp), intent(in)  :: istm(n,n), stm(n,n), istt(n,n,n), stt(n,n,n)
        real(wp), intent(out) :: stmab(n,n), sttab(n,n,n)

        stmab = mmult(stm,istm)
        sttab = mattens(stm,istt,n)
        sttab = sttab + quad(istm,stt,n)
    end subroutine q_sttchain
    pure subroutine d_sttchain(istm,istt,stm,stt,stmab,sttab,n) 
        implicit none
        integer, intent(in) :: n
        real(dp), intent(in)  :: istm(n,n), stm(n,n), istt(n,n,n), stt(n,n,n)
        real(dp), intent(out) :: stmab(n,n), sttab(n,n,n)

        stmab = mmult(stm,istm)
        sttab = mattens(stm,istt,n)
        sttab = sttab + quad(istm,stt,n)
    end subroutine d_sttchain

    pure subroutine sttchain_loop(istm,istt,stm,stt,stmab,sttab,n)
        implicit none
        integer i, a, b
        integer, intent(in) :: n
        real(wp), intent(in)  :: istm(n,n), stm(n,n), istt(n,n,n), stt(n,n,n)
        real(wp), intent(out) :: stmab(n,n), sttab(n,n,n)
        do a=1,n
        do i=1,n
            stmab(i,a) = dot_product(stm(i,:),istm(:,a))
        do b=1,n
            sttab(b,i,a) = dot_product(stm(i,:),istt(:,a,b)) &
                + dot_product(mmult(stt(i,:,:),istm(:,b)),istm(:,a))
        enddo; enddo; enddo
        sttab = reshape(sttab,shape=[n,n,n],order=[2,3,1])
    end subroutine sttchain_loop

end module tensorops
