program main
    use frkmin_q, only: solve_ivp, Odesolution, RungeKutta
    use util, only : print_to_file
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(odesolution)    :: base_sol
    real(qp), parameter  :: t0=0._qp, tf=4._qp*24._qp*3600._qp, tof = tf,&
                            rtol = 1.e-25_qp, atol = 1.e-27_qp
    integer, parameter   :: traj_id = -998, & 
                            central_body = 399, &
                            bodylist(3)= [10,301,5], &
                            nnodes = 5000
    real(qp), parameter  :: central_body_ref_radius=6378.137_qp, &
                            central_body_mu=398600.5_qp, &
                            mu_list(3)= [1.32712440041279419e11_qp, &
                                      &  4902.800118_qp, &
                                      &  126712764.1_qp]
    real(qp)             :: init_state(6)
    integer i, spkhand
    real(qp)             :: times(nnodes), states(nnodes,6)
    real(dp)             :: dtime(nnodes), x(6,nnodes), times_short(1000), states_short(6,1000), a, b

    call FURNSH("/home/david/wrk/nstgro/qist/kernels/mk_noorb.tf")
    times = [((tf-t0)/(nnodes-1) * i + t0, i=0,nnodes-1)]
    a = 0._dp
    b = 2._dp * 24._dp * 3600._dp
    times_short = [((b-a)/(1000-1)*i + a, i=0,1000-1)]
    init_state = [ &
              4.898587196589413e-13_qp, &
              8000.0_qp, &
              0.0_qp, &
              -2.956795170934981_qp, &
              1.8105148709099377e-16_qp, &
              8.123727966096327_qp &
                 ]

    print *, "Integrating base case"
    base_sol = solve_ivp(pert_eoms,&
                  & [t0, tf], &
                  & init_state, &
                  & dense_output=.true.,&
                  & rtol=rtol, &
                  & atol=atol)
    print *, "Done."
    print *, shape(base_sol%ts)
    print *, shape(base_sol%ys)
    do i=1,nnodes
        states(i,:) = base_sol%call(times(i))
        if (i.le.1000) then
            states_short(:,i) = real(base_sol%call(real(times_short(i),qp)),dp)
        endif
    end do
    dtime = real(times, dp)
    x = transpose(real(states, dp))
    call print_to_file("integrate_pre_resamp_x",states_short(1,:))
    call print_to_file("integrate_pre_resamp_y",states_short(2,:))
    call print_to_file("integrate_pre_resamp_z",states_short(3,:))
    call spkopn("perturbed_ker.bsp", "perturbedkernel",100,spkhand)
    call spkw13( &
           spkhand, & !handle
           -998, & ! body
           399, & ! center
           'J2000', & ! inframe
           dtime(1), & ! first
           dtime(nnodes), & ! last
           "test object ephemeris", & ! segid
           5, & ! degree
           nnodes, & ! nepocs
           x, & ! states
           dtime & ! epochs
          )
      call spkcls(spkhand)
    contains
        function pert_eoms(me,x,y) result(res)
            class(RungeKutta), intent(inout) :: me
            real(qp), intent(in) :: y(:), x
            real(qp)             :: res(size(y))
            real(dp)             :: radlist(size(bodylist),3), lt_dum
            integer i

            do i = 1, size(bodylist)
                call spkgeo( bodylist(i), real(x,dp), "J2000", central_body, &
                           & radlist(i,:), lt_dum)
            end do
            res = 0._qp
            res = xdot_kep(y)
            do i = 1, size(bodylist)
                res = res + xdot_nbody(y, mu_list(i),real(radlist(i,:),qp))
            end do
        end function pert_eoms
        function xdot_kep(y) result(res)
            implicit none
            real(qp), intent(in) :: y (:)
            real(qp)             :: &
                                  & x0

            real(qp), dimension(6) :: res
            x0  =  central_body_mu/sqrt(sum(y(:3)**2))**3.0_qp
            res(1) =  y(4)
            res(2) =  y(5)
            res(3) =  y(6)
            res(4) =  -x0*y(1)
            res(5) =  -x0*y(2)
            res(6) =  -x0*y(3)
        end function xdot_kep
        function xdot_nbody(y, mu, rj) result (res)
            ! acc_nbody: method to compute third-body Keplerian acceleration
            ! INPUTS:
            ! NAME           TYPE           DESCRIPTION
            ! mu             real           Gravitational parameter of third body
            ! y              real (:)       Augmented dynamical state vector
            ! rj          real (3)       Position vector from central body to
            !                               third body
            ! BEGIN AUTOCODE OUTPUT FOR XDOT
            implicit none
            real(qp), intent(in) :: y (:), mu, rj(3)
            real(qp)             :: &
                                  & x0, x1, x2, x3, x4
            real(qp), dimension(8) :: res
            x0  =  (rj(1)**2 + rj(2)**2 + rj(3)**2)**(-3.0_qp/2.0_qp)
            x1  =  rj(1) - y(1)
            x2  =  rj(2) - y(2)
            x3  =  rj(3) - y(3)
            x4  =  (x1**2 + x2**2 + x3**2)**(-3.0_qp/2.0_qp)
            ! velocities zeroed out 
            res(1) =  0._qp
            res(2) =  0._qp
            res(3) =  0._qp
            res(4) =  mu*(-rj(1)*x0 + x1*x4)
            res(5) =  mu*(-rj(2)*x0 + x2*x4)
            res(6) =  mu*(-rj(3)*x0 + x3*x4)
        end function xdot_nbody
end program main
