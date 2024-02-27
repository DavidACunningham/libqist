program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use cheby, only: chnodes, chderiv, chfit, chcall, pi, spice_subset
    implicit none
    type(spice_subset)    :: spiceb, spicea
    real(dp)            :: a, b 
    real(dp),allocatable            :: times(:), states(:,:)
    integer             ::  bodlist(4), deg, i, ntest

    a = 0._dp
    b = 2._dp * 24._dp * 3600._dp
    deg = 2000
    ntest = 1000
    times = [((b-a)/(ntest-1)*i + a, i=0,ntest-1)]
    bodlist = [10, 301, 5, -998]
    call spicea%init("/home/david/wrk/nstgro/qist/kernels/mk.tf", 399, bodlist, a, b, deg)
    open(file="./perturbed_reference.subspice",unit=73,access="stream",status="replace")
    call spicea%write(73)
    close(73)
    open(file="./perturbed_reference.subspice",unit=73,access="stream",status="old")
    call spiceb%read(73)
    close(73)
    allocate(states(3,ntest))
    do i=1,ntest
        states(:,i) = spiceb%call(times(i), -998, 'p')
    end do
    call print_to_file('spice_resamp_x', states(1,:))
    call print_to_file('spice_resamp_y', states(2,:))
    call print_to_file('spice_resamp_z', states(3,:))
    call print_to_file('times_resamp', times)
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
        subroutine read_from_file(fname, var)
            integer io,j
            real(dp), intent(inout),allocatable :: var(:)
            real(dp)        :: dummy
            character(len=*) :: fname
            integer stat, nlines
            open(newunit=io, file=trim(adjustl(fname))//".txt")
            nlines = 0
            stat = 0
            do while (stat.eq.0)
                read(io,*,iostat=stat)
                nlines = nlines + 1
            end do
            rewind(io)
            if (allocated(var)) deallocate(var)
            allocate(var(nlines-1))
            do j = 1,nlines-1
                read(io,*) var(j)
                ! write(*,*) var(j)
            end do
            close(io)
        end subroutine read_from_file
end program main
