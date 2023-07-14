module genqist
    use globals
    use makemodel
    use frkmain, only: solve_ivp, Odesolution, RungeKutta
    implicit none

    type gqist
        type(dynamicsmodel)   :: dynmod
        real(wp), allocatable :: initstate(:)
        real(wp)              :: rtol, atol
        contains
        procedure init
        procedure integrate
    end type gqist
    contains
        subroutine init(me,t0, tf, kernelfile, traj_id, central_body, bodylist, shgrav, Cbar, Sbar)
            class(gqist), intent(inout) :: me
            real(wp),     intent(in)    :: t0, tf
            integer,              intent(in)    :: traj_id, & 
                                                   central_body, &
                                                   bodylist(:)
            logical,              intent(in)    :: shgrav
            real(wp),             intent(in)    :: Cbar(:,:), &
                                                   Sbar(:,:)
            character(len=256),   intent(in)    :: kernelfile

            call me%dynmod%init(kernelfile, traj_id, central_body, bodylist, shgrav, Cbar, Sbar)

        end subroutine
        function integrate(meqist, t0, tf) result(res)
            class(gqist), intent(inout) :: meqist
            type(ODESolution)           :: res
            real(wp), intent(in)        :: t0, tf

            res = solve_ivp(myint_eoms,&
                                [t0, tf], &
                                meqist%dynmod%trajstate(t0),&
                                method="DOP853",&
                                dense_output=.true.,&
                                rtol=meqist%rtol*1.E-2_wp, atol=meqist%atol, istep=1._wp)
            
            contains
                function myint_eoms(me, x, y) result(res)
                    class(RungeKutta), intent(inout) :: me
                    real(wp),          intent(in)    :: x, y(:)
                    real(wp)                         :: res(size(y))
                    res = meqist%dynmod%eoms(x,y)
                end function
        end function integrate
end module genqist
