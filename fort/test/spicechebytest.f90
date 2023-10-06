program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use cheby, only: chnodes, chderiv, chfit, chcall, pi, vectorcheb
    implicit none
    character(len=12) :: arg
    integer, parameter  :: deg=100, tlen=1000, indvar=3
    type(vectorcheb)    :: pos, vel, acc
    real(dp)            :: a, b
    real(dp)            :: lt_dum, &
                           spkgeo_out(6), spkgps_out(3)
    real(dp)            :: nodes(deg), fnodes(deg,3), testpoints(tlen), &
                           truth(tlen,3), dtruth(tlen,3), ddtruth(tlen,3), &
                           cheb(tlen,3), dcheb(tlen,3), ddcheb(tlen,3)
    integer i, bod_id

    call FURNSH("/home/david/wrk/nstgro/qist/kernels/mk.tf")
    call get_command_argument(1, arg)
    read(arg,*) bod_id
    call get_command_argument(2, arg)
    read(arg,*) a
    call get_command_argument(3, arg)
    read(arg,*) b
    nodes = chnodes(deg, a, b)
    testpoints = [(a + i*b/tlen, i=1,tlen)]
    do i=1,size(nodes)
        call spkgps(bod_id, nodes(i), "J2000", 399, &
                          & spkgps_out, lt_dum)
        fnodes(i,1) = spkgps_out(1)
        fnodes(i,2) = spkgps_out(2)
        fnodes(i,3) = spkgps_out(3)
    end do
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
    call pos%fit(fnodes,a,b)
    vel = pos%deriv()
    acc = vel%deriv()
    cheb = pos%call(testpoints)
    dcheb = vel%call(testpoints)
    ddcheb = acc%call(testpoints)
    call print_to_file("truth", truth(:,indvar))
    call print_to_file("dtruth", dtruth(:,indvar))
    call print_to_file("cheb", cheb(:,indvar))
    call print_to_file("dcheb", dcheb(:,indvar))
    call print_to_file("ddcheb", ddcheb(:,indvar))
    call print_to_file("nodes", nodes)
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
