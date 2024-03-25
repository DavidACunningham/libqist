program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64
    use cheby, only: vectorcheb, chnodes, chfit, chcall
    use quat
    implicit none
    real(dp), parameter :: epoch=10000._dp, epoch2=336000._dp
    integer, parameter :: nnodes=30
    real(dp), dimension(3,3) :: rotmat, rotmat2, rotmat_comp
    real(dp)                 :: qquat(4), nodes(nnodes)
    real(dp)                 :: time_array, funcs(nnodes,4), quat2(4), &
                              & sfi(nnodes), scoeffs(nnodes), &
                              & time_test_points(1000), test_out(1000,4)
    integer i
    type(vectorcheb)         :: quatfit
    type(quaternion)         :: iquat

    call furnsh("./kernels/mk.tf")
    call pxform('J2000','IAU_EARTH',epoch,rotmat)
    call pxfrm2('J2000','IAU_EARTH',epoch,epoch2,rotmat2)
    call pxfrm2('IAU_EARTH','IAU_EARTH',epoch,epoch2,rotmat_comp)
    call m2q(rotmat,qquat)
    print * , real(qquat,4)
    print *, ''
    call pmat(rotmat)
    print*, ""
    call pmat(rotmat2)
    print*, ""
    call pmat(matmul(rotmat_comp,rotmat))
    print*, ""
    call pmat(rotmat2-matmul(rotmat_comp,rotmat))

    call pxfrm2('IAU_EARTH','IAU_EARTH',epoch,epoch+1.,rotmat_comp)
    call iquat%fromdcm(rotmat_comp)
    qquat = iquat%q
    nodes = chnodes(nnodes, epoch,epoch2)
    funcs(1,:) = qfit(epoch,nodes(1),qquat)
    do i = 2, nnodes
        funcs(i,:) = qfit(epoch,nodes(i),funcs(i-1,:))
    end do
    time_test_points = [(epoch + (epoch2-epoch)/1000 * i, i=1,1000)]
    call quatfit%fit(funcs,epoch,epoch2)

    do i = 1,1000
        test_out(i,:) = quatfit%call(time_test_points(i))
        test_out(i,:) = test_out(i,:)/norm2(test_out(i,:))
    end do
    call print_to_file("out1",test_out(:,1))
    call print_to_file("out2",test_out(:,2))
    call print_to_file("out3",test_out(:,3))
    call print_to_file("out4",test_out(:,4))
    call print_to_file('tpoints', time_test_points)
    
    
    call pxfrm2('IAU_EARTH','IAU_EARTH',epoch,5000._dp,rotmat_comp)
     
    ! call m2q( rotmat_comp,quat)

    contains
    subroutine pmat(mat)
        real(dp), intent(in) :: mat(3,3)
        integer i
        do i = 1,3
            print *, real(mat(i,:),4)
        end do
    end subroutine

    function qfit(a,b,pquat) result(res)
        real(dp), intent(in) :: a,b, pquat(4)
        real(dp)             :: res(4), mat(3,3), dot
        type(quaternion)         :: qclass
        call pxfrm2('IAU_EARTH','IAU_EARTH',a,b,mat)
        call qclass%fromdcm(mat)
        dot = dot_product(pquat,qclass%q)
        if (dot.lt.0) then
            res = -qclass%q
        else
            res = qclass%q
        endif
    end function qfit
    
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
