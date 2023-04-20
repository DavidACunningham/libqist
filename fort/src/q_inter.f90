module q_inter
    use qist, only: Itraj
    use globals, only: wp
    implicit none
    integer, parameter :: dp=selected_real_kind(15)
    integer, parameter, private :: n=6
    type(Itraj) :: it
    interface prop
        module procedure prop_once
        module procedure prop_many
    end interface prop
    contains

    !! Initializers
    subroutine init_i(t0, tf,filepath, trajfile)
        character(len=*), intent(in) :: filepath, trajfile
        real(dp),            intent(in) ::t0, tf
        !! The initial and final simulation independent variable (t)
        !! The function _returns_ an instance of the type. 
        call it%init(t0, tf,filepath, trajfile)
    end subroutine init_i
    !! Calling functions 
    !! Physical time calling functions
    function state(t) result(res)
        !! Return a nonregularized state at time t
        real(dp),       intent(in) :: t
        !! The value of the independent variable to generate
        real(dp)                   :: res(n)
        res = it%state(t)
    end function state
     function stm(t) result(res)
        !! Return a nonregularized stm at time t
        real(dp),     intent(in) :: t
        !! The value of t at which to get the stm
        real(dp), dimension(n,n) :: res
        !! The returned stm
        res = it%stm(t)
    end function stm
     function stm_i(t) result(res)
        !! Return the ISTT from a time t to t0
        real(dp),     intent(in) :: t
        real(dp), dimension(n,n) :: res
        res = it%stm_i(t)
    end function stm_i
     function stt(t) result(res)
        real(dp),     intent(in) :: t
        !! The value of t at which to get the stm
        real(dp), dimension(n,n,n) :: res
        res = it%stt(t)
    end function stt
     function stt_i(t) result(res)
        !! Return the ISTT from a time t to t0
        real(dp),     intent(in) :: t
        real(dp), dimension(n,n,n) :: res
        res = it%stt_i(t)
    end function stt_i
     function prop_once(ta, tb, xa,order) result(res)
        !! Propagates the relative state xa at ta
        !! to tb
        real(dp),     intent(in) :: ta, tb, xa(7)
        integer, intent(in), optional :: order
        integer o
        !! xa is the initial relative state 
        !! should be dimension 7
        real(dp), dimension(n)   :: res
        o=2
        if (present(order)) o = order
        res = it%prop(ta,tb,xa,order)
    end function prop_once
    function prop_many(ta, tb, xa,order) result(res)
        !! Propagates the relative state xa at ta
        !! to tb
        real(dp),     intent(in) :: ta, tb, xa(:,:)
        integer, intent(in), optional :: order
        integer o
        !! xa is the initial relative state 
        !! should be dimension 7
        real(dp), dimension(7,size(xa,2))   :: res, xapad(7,size(xa,2))
        xapad(:6,:) = xa
        xapad(7,:) = 0._8
        o=2
        if (present(order)) o = order
        res = it%prop(ta,tb,xapad,o)
    end function prop_many

    !! physical time propagation

    subroutine stts_ab(ta, tb, stm, stt)
        !! Return the STM and STT from ta to tb
        real(dp),     intent(in)  :: ta, tb
        real(dp),     intent(out) :: stm(n,n), stt(n,n,n)
        call it%stts_ab(ta,tb,stm,stt)
    end subroutine stts_ab

     function zmap(t,order) result(res)
        !! Testing procedure that chains an stm and stt with their own
        !! inverse--you ought to get I and 0, respectively
        real(dp),          intent(in) :: t
        integer, optional, intent(in) :: order
        integer                       :: ord
        real(dp)                      :: res

        ord=2
        if (present(order)) ord=order
        res = it%zmap(t,ord)
    end function zmap
end module q_inter
