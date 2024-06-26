program main
    use cheby, only: spice_subset
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    type(spice_subset) :: subspice
    real(dp)           :: times(2000), spicepos(2000,3)
    integer i

    call subspice%init("/home/david/wrk/nstgro/qist/kernels/mk.tf", &
                     & 399, [10, 301, 5, -999], &
                     & 0._dp, 2._dp*24._dp*3600._dp, &
                     & 500)

    open(unit=75, file="./twbody_resample.subspice", &
       & access="stream", status="replace")
        call subspice%write(75)
    close(75)

    times = [(i*(2._dp*24._dp*3600._dp)/2000.,i = 1,2000)]

    spicepos = subspice%call(times,-999,"p")

    call print_to_file("times",times)
    call print_to_file("spiceposx",spicepos(:,1))
    call print_to_file("spiceposy",spicepos(:,2))
    call print_to_file("spiceposz",spicepos(:,3))
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
