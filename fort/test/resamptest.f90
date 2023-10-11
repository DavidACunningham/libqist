program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use cheby, only: chnodes, chderiv, chfit, chcall, pi, spice_subset
    implicit none
    character(len=12) :: arg
    integer, parameter  :: tlen=1000, indvar=1
    type(spice_subset)    :: spice
    real(dp)            :: a, b
    real(dp)            :: lt_dum, &
                           spkgeo_out(6)
    real(dp)            :: testpoints(tlen), &
                           truth(tlen,3), dtruth(tlen,3), &
                           cheb(tlen,3), dcheb(tlen,3), ddcheb(tlen,3)
    integer             :: i, bodlist(9), bod_id, deg

    call get_command_argument(1, arg)
    read(arg,*) bod_id
    call get_command_argument(2, arg)
    read(arg,*) a
    call get_command_argument(3, arg)
    read(arg,*) b
    call get_command_argument(4, arg)
    read(arg,*) deg
    bodlist = [1, 2, 301, 4, 5, 6, 7, 8, 9]
    call spice%init("/home/david/wrk/nstgro/qist/kernels/mk.tf", 399, bodlist, a, b, deg)
    testpoints = [(a + i*(b-a)/tlen, i=1,tlen)]
    do i=1,size(testpoints)
        call spkgeo(bod_id, testpoints(i), "J2000", 399, &
                          & spkgeo_out, lt_dum)
        truth(i,1) = spkgeo_out(1)
        truth(i,2) = spkgeo_out(2)
        truth(i,3) = spkgeo_out(3)
        dtruth(i,1) = spkgeo_out(4)
        dtruth(i,2) = spkgeo_out(5)
        dtruth(i,3) = spkgeo_out(6)
    end do
    cheb = spice%call(testpoints,bod_id,"p")
    dcheb = spice%call(testpoints,bod_id,"v")
    ddcheb = spice%call(testpoints,bod_id,"a")
    call print_to_file("truth", truth(:,indvar))
    call print_to_file("dtruth", dtruth(:,indvar))
    call print_to_file("cheb", cheb(:,indvar))
    call print_to_file("dcheb", dcheb(:,indvar))
    call print_to_file("ddcheb", ddcheb(:,indvar))
    call print_to_file("testpoints", testpoints)
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
