program main
    use genqist, only: gqist
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(gqist)         :: qist
    real(qp), parameter  :: t0=0._qp, tf=2._qp*24._qp*3600._qp, tof = tf,&
                            rtol = 1.e-17_qp, atol = 1.e-20_qp
    integer, parameter   :: traj_id = -998, & 
                            central_body = 399, &
                            bodylist(3)= [10,301,5]
    logical, parameter   :: shgrav = .false., rails=.true.
    integer, parameter   :: offset = 8
    real(qp), parameter  :: central_body_ref_radius=6381.137_qp, &
                            central_body_mu=398600.5_qp, &
                            mu_list(3)= [1.32712440041279419e11_qp, &
                                      &  4902.800118_qp, &
                                      &  126712764.1_qp]
    real(qp), parameter  :: Cbar(2,2) = 0._qp, &
                            Sbar(2,2) = 0._qp
    real(qp)             :: states(6,1000), time(1000)
    integer i

    call qist%init(t0, tf, &
                 & "./perturbed_reference.subspice", &
                 &  traj_id, &
                 &  399, &
                 &  bodylist, &
                 &  central_body_mu, &
                 &  central_body_ref_radius, &
                 &  mu_list, &
                 &  shgrav, &
                 &  Cbar, Sbar, rails, .true.)
    qist%dynmod%tof = tof
    time = [(i*(tof/1000._qp), i=1,1000)]
    do i = 1, 1000
        states(:3,i) = real(qist%dynmod%bod_db%call(real(time(i),dp),qist%dynmod%traj_id,'p'),qp)
        states(4:,i) = real(qist%dynmod%bod_db%call(real(time(i),dp),qist%dynmod%traj_id,'v'),qp)
    end do
    call print_to_file("testx",real(states(1,:),dp))
    call print_to_file("testy",real(states(2,:),dp))
    call print_to_file("testz",real(states(3,:),dp))
    contains
        subroutine print_to_file(fname, var)
            integer io,j
            real(dp), intent(in) :: var(:)
            character(len=*) :: fname
            open(newunit=io, file=trim(adjustl(fname))//".txt")
            do j = 1,size(var)
            write(io,*) var(j)
            end do
            close(io)
        end subroutine print_to_file
end program main
