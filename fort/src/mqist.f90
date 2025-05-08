! Title: mqist.f90 
! Description:
!    A thin wrapper for QIST functions
!    to allow a pre-generated QIST library
!    to be called by matlab using the loadlib function.
!    needs to be kept in sync with mqist.h.
!
! References:
!   None
! 
! author: David Cunningham
! Last edited: See git log
subroutine m_init_n(namefile_c)
    ! Initialize a QIST library from a compatible namefile
    !DEC$ ATTRIBUTES DLLEXPORT :: m_init_n
    use q_inter
    use, intrinsic :: iso_c_binding
    implicit none
    character(kind=c_char), dimension(*), intent(in) :: namefile_c
    character(len=:), allocatable :: namefile
    integer strlen, i
    !! The initial and final simulation independent variable (tau)
    !! The function _returns_ an instance of the type. 
    strlen = 0
    do
        if(namefile_c(strlen+1) == c_null_char) exit
        strlen = strlen+1
    end do
    allocate(character(len=strlen) :: namefile)
    do i=1,strlen
        namefile(i:i) = namefile_c(i)
    end do
    call init_n(namefile)
end subroutine m_init_n

subroutine m_state(tau,state_o)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_state
    use q_inter
    implicit none
    !! Return a state at time tau
    integer, parameter :: n=8
    real(8),     intent(in) :: tau
    !! The value of tau at which to get the state
    real(8), intent(out), dimension(n)   :: state_o
    !! The returned state
    state_o = state(tau)
end subroutine m_state

subroutine m_stm(tau,stm_o)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_stm
    use q_inter
    implicit none
    !! Return an stm at time tau
    integer, parameter :: n=8
    real(8),     intent(in) :: tau
    !! The value of tau at which to get the state
    real(8), intent(out), dimension(n,n) :: stm_o
    !! The returned stm
    stm_o = stm(tau)
end subroutine m_stm

subroutine m_stt(tau,stt_o)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_stt
    use q_inter
    implicit none
    !! Return an stt at time tau
    integer, parameter :: n=8
    real(8),     intent(in) :: tau
    !! The value of tau at which to get the stt
    real(8), intent(out), dimension(n,n,n) :: stt_o
    !! The returned stt
    stt_o = stt(tau)
end subroutine m_stt

subroutine m_stm_i(tau,stm_i_o)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_stm_i
    use q_inter
    implicit none
    !! Return an inverse stm at time tau
    integer, parameter :: n=8
    real(8),     intent(in) :: tau
    !! The value of tau at which to get the state
    real(8), intent(out), dimension(n,n) :: stm_i_o
    !! The returned stm
    stm_i_o = stm_i(tau)
end subroutine m_stm_i

subroutine m_stt_i(tau,stt_i_o)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_stt_i
    use q_inter
    implicit none
    !! Return an inverse stt at time tau
    integer, parameter :: n=8
    real(8),     intent(in) :: tau
    !! The value of tau at which to get the state
    real(8), intent(out), dimension(n,n,n) :: stt_i_o
    !! The returned stt
    stt_i_o = stt_i(tau)
end subroutine m_stt_i

!! physical time propagation
subroutine m_prop_once(ta, tb, xa, order, xb)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_prop_once
    use q_inter
    implicit none
    !! Propagates the relative state xa at time ta
    !! to time tb.
    integer :: o 
    real(8), intent(in)  :: ta, tb, xa(6)
    integer, intent(in)  :: order
    !! xa is the initial relative state 
    !! should be dimension 6
    real(8), intent(out) :: xb(8)
    o = order
    xb = prop(ta,tb,[xa, 0._dp, 0._dp],o)
end subroutine m_prop_once

subroutine m_stts_ab(taua, taub, stm_o, stt_o)
    !DEC$ ATTRIBUTES DLLEXPORT :: m_stts_ab
    use q_inter
    implicit none
    !! Return the STM and STT from taua to taub
    integer, parameter :: n=8
    real(8),     intent(in)  :: taua, taub
    real(8),     intent(out) :: stm_o(n,n), stt_o(n,n,n)
    call stts_ab(taua,taub,stm_o,stt_o)
end subroutine m_stts_ab
