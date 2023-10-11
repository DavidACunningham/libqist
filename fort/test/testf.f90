program main
    implicit none
    integer :: a(27), i

    a = [(i, i=1,27)]

    print *, a(:3)
    do i=2,9
    print *, i
    print *, a(3*(i-1)+1:3*i)
    end do

end program main
