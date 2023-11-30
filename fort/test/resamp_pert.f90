program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use cheby, only: chnodes, chderiv, chfit, chcall, pi, spice_subset
    implicit none
    type(spice_subset)    :: spiceb
    real(dp)            :: a, b, times(500), states(3,500)
    integer             ::  bodlist(4), deg, i

    a = 0._dp
    b = 2._dp * 24._dp * 3600._dp
    deg = 300
    times = [(i*b/500._dp, i=1,500)]
    bodlist = [10, 301, 5, -998]
    call spiceb%init("/home/david/wrk/nstgro/qist/kernels/mk.tf", 399, bodlist, a, b, deg)
    open(file="./perturbed_reference.subspice",unit=73,access="stream",status="replace")
    call spiceb%write(73)
    close(73)
    do i=1,500
        states(:,i) = spiceb%call(times(i), -998, 'p')
    end do
    call print_to_file('spice_resamp_x', states(1,:))
    call print_to_file('spice_resamp_y', states(2,:))
    call print_to_file('spice_resamp_z', states(3,:))
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
