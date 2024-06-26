program main
    use util, only: print_to_file
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(odesolution)    :: base_sol
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, tof = tf,&
                            rtol = 1.e-25_qp, atol = 1.e-27_qp
    integer, parameter   :: nnodes = 1000
    real(qp), parameter  :: central_body_mu=398600.5_qp
    real(qp)             :: init_state(6)
    integer i
    real(qp)             :: times(nnodes), states(6,nnodes)

    times = [((tf-t0)/(nnodes-1) * i + t0, i=0,nnodes-1)]
    init_state = [ &
              4.898587196589413e-13_qp, &
              8000.0_qp, &
              0.0_qp, &
              -2.956795170934981_qp, &
              1.8105148709099377e-16_qp, &
              8.123727966096327_qp &
                 ]
    print *, "Integrating base case"
    base_sol = solve_ivp(xdot_kep,&
                  & [t0, tf], &
                  & init_state, &
                  & dense_output=.true.,&
                  & rtol=rtol, &
                  & atol=atol)
    print *, "Done."
    do i=1,nnodes
        states(:,i) = base_sol%call(times(i))
    end do
    call print_to_file("base_sol_x", real(states(1,:),dp))
    call print_to_file("base_sol_y", real(states(2,:),dp))
    call print_to_file("base_sol_z", real(states(3,:),dp))
    contains
        function xdot_kep(me, x,y) result(res)
            implicit none
            class(RungeKutta), intent(inout) :: me
            real(qp), intent(in) :: x, y (:)
            real(qp)             :: &
                                  & x0

            real(qp)             :: res(size(y))
            x0  =  central_body_mu/sqrt(sum(y(:3)**2))**3.0_qp
            res(1) =  y(4)
            res(2) =  y(5)
            res(3) =  y(6)
            res(4) =  -x0*y(1)
            res(5) =  -x0*y(2)
            res(6) =  -x0*y(3)
        end function xdot_kep
end program main
