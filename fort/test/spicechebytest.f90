program main
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    use cheby, only: chnodes, chderiv, chfit, chcall, pi
    implicit none
    character(len=12) :: arg
    integer, parameter  :: deg=100, tlen=1000, elindex=3
    real(dp)            :: a, b
    real(dp)            :: lt_dum, &
                           spkgeo_out(6), spkgps_out(3)
    real(dp)            :: nodes(deg), fnodes(deg), testpoints(tlen), &
                           truth(tlen), dtruth(tlen), ddtruth(tlen), &
                           cheb(tlen), dcheb(tlen), &
                           coeffs(deg), dcoeffs(deg-1)
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
        fnodes(i) = spkgps_out(elindex)
    end do
    do i=1,size(testpoints)
        call spkgeo(bod_id, testpoints(i), "J2000", 399, &
                          & spkgeo_out, lt_dum)
        truth(i) = spkgeo_out(elindex)
        dtruth(i) = spkgeo_out(elindex+3)
    end do
    coeffs = chfit(deg,fnodes)
    dcoeffs = chderiv(coeffs,a,b)
    cheb = chcall(a,b,coeffs,testpoints)
    dcheb = chcall(a,b,dcoeffs,testpoints)
    call print_to_file("truth", truth)
    call print_to_file("dtruth", dtruth)
    call print_to_file("ddtruth", ddtruth)
    call print_to_file("coeffs", coeffs)
    call print_to_file("dcoeffs", dcoeffs)
    call print_to_file("cheb", cheb)
    call print_to_file("dcheb", dcheb)
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
