! Title: pyqist.f90 
! Description:
!    A thin wrapper for QIST functions
!    to allow a pre-generated QIST library
!    to be called by a python and compiled by 
!    f2py.
!
! References:
!   None
! 
! author: David Cunningham
! Last edited: See git log
module pq
    use q_inter
    implicit none
    contains
    !! Initializers
    subroutine pw_init_v(t0, tf, trajfile,kvtaufile)
        ! Initialize from filenames--namefile initializer is preferred
        character(len=*), intent(in) :: trajfile, kvtaufile
        real(8),            intent(in) ::t0, tf
        !! The initial and final simulation independent variable (tau)
        !! The function _returns_ an instance of the type. 
        call init_v(t0, tf, trajfile, kvtaufile)
    end subroutine pw_init_v
    subroutine pw_init_n(namefile)
        ! Initialize from namefile
        character(len=*), intent(in) :: namefile
        !! The function _returns_ an instance of the type. 
        call init_n(namefile)
    end subroutine pw_init_n
     function pw_state(tau) result(res)
        !! Return a state at time tau
        real(8),     intent(in) :: tau
        !! The value of tau at which to get the state
        real(8)            :: res
        !! The returned state
        res = state(tau)
    end function pw_state
     function pw_stm(tau) result(res)
        !! Return a stm at time tau
        integer, parameter :: n=8
        real(8),     intent(in) :: tau
        !! The value of tau at which to get the stm
        real(8), dimension(n,n) :: res
        !! The returned stm
        res = stm(tau)
    end function pw_stm
     function pw_stt(tau) result(res)
        !! Return a stt at time tau
        integer, parameter :: n=8
        real(8),     intent(in) :: tau
        !! The value of tau at which to get the stt
        real(8), dimension(n,n,n) :: res
        !! The returned stt
        res = stt(tau)
    end function pw_stt
     function pw_stm_i(tau) result(res)
        !! Return an inverse stm at time tau
        integer, parameter :: n=8
        real(8),     intent(in) :: tau
        !! The value of tau at which to get the stm
        real(8), dimension(n,n) :: res
        !! The returned inverse stm
        res = stm_i(tau)
    end function pw_stm_i
     function pw_stt_i(tau) result(res)
        !! Return an inverse stt at time tau
        integer, parameter :: n=8
        real(8),     intent(in) :: tau
        !! The value of tau at which to get the stt
        real(8), dimension(n,n,n) :: res
        !! The returned stt
        res = stt_i(tau)
    end function pw_stt_i

    !! physical time propagation
    function pw_prop_once(ta, tb, xa, order) result(res)
        !! Propagates the relative state xa at time ta
        !! to time tb.
        integer :: o 
        real(8),     intent(in) :: ta, tb, xa(8)
        integer, intent(in), optional :: order
        !! xa is the initial relative state 
        !! should be dimension 8
        real(8)                 :: res(8)
        o = 2
        if (present(order)) o=order
        res = prop(ta,tb,xa,o)
    end function pw_prop_once

    function pw_prop_many(ta, tb, xa, order) result(res)
        !! Propagates the relative state xa at time ta
        !! to time tb.
        integer :: o
        real(8),     intent(in) :: ta, tb, xa(:,:)
        integer, intent(in), optional :: order
        !! xa is the initial relative state 
        !! should be dimension 8,:
        real(8)                 :: res(8,size(xa,2))
        if (present(order)) o=order
        res = prop(ta,tb,xa, o)
    end function pw_prop_many

    function pw_prop_back(tb, ta, xb, order) result(res)
        !! Propagates the relative state xb at time tb
        !! to time ta.
        integer :: o 
        real(8),     intent(in) :: ta, tb, xb(8)
        integer, intent(in), optional :: order
        !! xa is the initial relative state 
        !! should be dimension 8
        real(8)                 :: res(8)
        o = 2
        if (present(order)) o=order
        res = prop_back(tb,ta,xb,o)
    end function pw_prop_back

    subroutine pw_stts_ab(taua, taub, stm, stt)
        !! Return the STM and STT from taua to taub
        integer, parameter :: n=8
        real(8),     intent(in)  :: taua, taub
        real(8),     intent(out) :: stm(n,n), stt(n,n,n)
        call stts_ab(taua,taub,stm,stt)
    end subroutine pw_stts_ab
    subroutine pw_stt_update(ta, tb, xa,new_stm,new_stt)
        !! Return the STM and STT for a perturbation to the relative trajectory
        integer, parameter :: n=8
        real(8),     intent(in)     :: ta, tb, xa(n)
        real(8),     intent(out)    :: new_stm(n,n), new_stt(n,n,n)
        call stt_update(ta,tb, xa, new_stm, new_stt)
    end subroutine pw_stt_update
    subroutine pw_tensor_change_of_basis(RNOf, RON0, old_stm, old_stt, &
                                      new_stm, new_stt)
        !! Transform an STM and STT from an old coordinate basis to a new one
        !! RNO defined by vec_new = RNO@vec_old
        integer, parameter :: n=8
        real(8),     intent(in)     :: RNOf(n,n), RON0(n,n), old_stm(n,n), old_stt(n,n,n)
        real(8),     intent(out)    :: new_stm(n,n), new_stt(n,n,n)
        call tensor_change_of_basis( RNOf, RON0, old_stm, old_stt, &
                                      new_stm, new_stt)
    end subroutine pw_tensor_change_of_basis

     function pw_zmap(tau,order) result(res)
        !! Testing procedure that chains an stm and stt with their own
        !! inverse--you ought to get I and 0, respectively
        integer, parameter :: n=7
        real(8),          intent(in) :: tau
        integer, optional, intent(in) :: order
        integer                       :: ord
        real(8)                      :: res

        ord=2
        if (present(order)) ord=order
        res = zmap(tau,ord)
    end function pw_zmap
end module pq
