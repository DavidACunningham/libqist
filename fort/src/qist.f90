module qist
    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use globals
    use tensorops, only:  sttchain, sttchain_invert, vectensquad, &
                          stminvert, sttinvert
    use denseLight, only: lightSol
    implicit none
    private
    !! THIS IS DUMB AND FOR DEBUGGING AND ALSO A HACK
    !! FOR INTEGRATING INVERSES AT QUAD PRECISION
    ! integer, parameter, public :: dp=wp

    !! THIS IS NOT DUMB AND SHOULD GIVE THE DEFAULT BEHAVIOR
    ! integer, parameter :: dp=selected_real_kind(15)

    type, public :: Itraj
        !! All type variables accessed through getter/setter methods
        real(dp)               :: t0, tf 
        logical                :: initq
        type(lightSol)         :: gw
        character(len=1000)    :: densefilename

    contains
        !! Convenience procedures/check state
        procedure :: init 
        procedure call
        procedure :: state 
        procedure :: stm 
        procedure :: stm_i 
        procedure :: stt
        procedure :: stt_i 
        generic, public :: prop => prop_once, prop_many
        procedure, private :: prop_once
        procedure, private :: prop_many

        procedure stts_ab
        procedure zmap
    end type Itraj

    contains
    ! Initializer
    subroutine init(self, t0, tf, filepath, trajfile)
        character(len=*), intent(in) :: filepath, trajfile
        !! The filename where the chebyshev coefficients are located
        real(dp),            intent(in) :: t0, tf
        !! The initial and final simulation independent variable (t)
        class(Itraj),        intent(out):: self
        !! The function _returns_ an instance of the type. 
        ! print *, "Initializing reference model"
        self%initq=.true.
        ! The initializer is running
        self%t0               = t0
        self%tf               = tf
        open(unit=49, file=filepath//trajfile, status="old", access="stream")
        call self%gw%read(49)
        close(49)
        ! Make sure the chebyshev order in the file is good to go
        ! Set the coefficient filename
        ! Allocate the chebyshev coefficient array in memory
    end subroutine init
    !! Calling functions 
    function call(self,t,lind,uind) result(res)
        !! Return a (uind-lind)-dimensional at t
        !! If no lind (uind) is passed, beginning (end) of the vector is used
        class(Itraj),      intent(inout) :: self
        real(dp),          intent(in) :: t
        integer, optional, intent(in) :: uind, lind
        !! The value of the independent variable to generate
        real(dp), allocatable         :: res(:)
        integer l, u
        l = 1
        u = plen
        if (present(lind)) l=lind
        if (present(uind)) u=uind
        allocate(res(u-l+1))
        res = self%gw%call(t,l,u)
    end function
    function state(self,t) result(res)
        !! Return a regularized state at time t
        class(ITraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp), dimension(n)   :: res
        !! The returned state
        res = self%call(t,uind=n)
    end function state
    function stm(self,t) result(res)
        !! Return a regularized stm at time t
        !! Return a regularized stm at time t
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp)                 :: packlin(n*n)
        real(dp), dimension(n,n) :: res
        !! The returned stm
        packlin = self%call(t,lind=stml,uind=stmu)
        res = reshape(packlin,[n,n])
    end function stm
    function stm_i(self,t) result(res)
        !! Return a regularized stm at time tau
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        real(dp), dimension(n,n) :: res
        !! The returned stm
        res = stminvert(self%stm(t), n)
    end function stm_i
    function stt(self,t) result(res)
        !! Return a regularized stm at time t
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        !! The value of t at which to get the state
        real(dp)                 :: packlin(n*n*n)
        real(dp), dimension(n,n,n) :: res
        !! The returned stm
        packlin = self%call(t,lind=sttl,uind=sttu)
        res = reshape(packlin,[n,n,n])
    end function stt
    function stt_i(self,t) result(res)
        !! Return a regularized stm at time t
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: t
        real(dp), dimension(n,n,n) :: res
        !! This can be improved to get both the STM and STT
        !! with one binary search
        res = sttinvert(self%stm(t),self%stt(t), n)
    end function stt_i
    function prop_once(self,ta, tb, xa, order) result(res)
        !! Propagates the relative state xa at ta
        !! to tb
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: ta, tb, xa(6)
        integer, intent(in), optional :: order
        integer :: ord
        !! xa is the initial relative state 
        !! should be dimension 6
        real(dp), dimension(n)   :: res
        real(dp) :: stmab(n,n), sttab(n,n,n)
        ord = 2
        if (present(order)) ord=order
        select case (ord)
        case (1)
            stmab = matmul(self%stm(tb),self%stm_i(ta))
            res(:n) = matmul(stmab,xa)
        case default
            call sttchain(self%stm_i(ta),self%stt_i(ta), &
                           & self%stm(tb), self%stt(tb), &
                           & stmab, sttab, n)
            res(:n) = matmul(stmab,xa) + 0.5_dp*vectensquad(xa,sttab,n)
        end select 
    end function prop_once
    function prop_many(self,ta, tb, xa, order) result(res)
        !! Propagates the relative state xa at ta
        !! to tb
        class(Itraj), intent(inout) :: self
        real(dp),     intent(in) :: ta, tb, xa(:,:) ! n rows, : columns
        integer, intent(in), optional :: order
        integer :: ord, i
        !! xa is the initial relative state 
        !! should be dimension 7
        real(dp) :: res(7,size(xa,2))
        real(dp) :: stmab(n,n), sttab(n,n,n)
        ord = 2
        if (present(order)) ord=order
        call sttchain(self%stm_i(ta),self%stt_i(ta), &
                       & self%stm(tb), self%stt(tb), &
                       & stmab, sttab, n)
        select case (ord)
        case (1)
            res(:n,:) = matmul(stmab,xa)
        case default
            !$OMP PARALLEL DO
            do i=1,size(xa,2)
                res(:n,i) = matmul(stmab,xa(:,i)) + &
                           0.5_dp*vectensquad(xa(:,i),sttab,n)
            end do
            !$OMP END PARALLEL DO
        end select 
    end function prop_many
    subroutine stts_ab(self, ta, tb, stm, stt)
        !! Return the STM and STT from ta to tb
        class(Itraj), intent(inout)  :: self
        real(dp),     intent(in)  :: ta, tb
        real(dp),     intent(out) :: stm(n,n), stt(n,n,n)

        call sttchain_invert(self%stm(ta),self%stt(ta), &
                    & self%stm(tb), self%stt(tb), &
                    & stm, stt, n)
    end subroutine stts_ab
    function zmap(self, t,order) result(res)
        !! Testing procedure that chains an stm and stt with their own
        !! inverse--you ought to get I and 0, respectively
        class(Itraj),      intent(inout) :: self
        real(dp),          intent(in) :: t
        integer, optional, intent(in) :: order
        integer                       :: ord, i
        real(dp)                      :: stm(n,n), stt(n,n,n), eye(n,n)
        real(dp)                      :: res

        ord=2
        if (present(order)) ord=order

        res = 0._dp
        eye = 0._dp
        do i=1,n; eye(i,i) = 1._dp; end do
        call sttchain(self%stm_i(t),self%stt_i(t), &
                    & self%stm(t), self%stt(t), &
                    & stm, stt, n)

        res = norm2(eye - stm)
        if (ord.eq.2) res = res + norm2(stt)
    end function zmap
end module qist
