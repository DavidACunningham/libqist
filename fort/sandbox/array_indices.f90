program array_indices
    implicit none
    integer, parameter :: n=7
    character(len=5), dimension(n,n,n) :: cubearray, truth
    character(len=5), dimension(n**3) :: linarrayA, linarrayB
    character(len=5), dimension(n**2*(n-1)/2) :: packarray
    character(len=5), dimension(n*(n-1)*(n-1)) :: removerbarray
    character(len=5), dimension(n*n*(n+1)/2) :: removeudarray
    integer :: i,j,k,l

    l=1
    do k=1,n
        do j=1,n
            do i=1,n
                write(cubearray(i,j,k),"(I1,I1,I1)") i,j,k
                ! write(cubearray(i,j,k),"(I3)") l
                l=l+1
            enddo
        enddo
    enddo
    truth = cubearray

    do i=1,n
    print "(A,I1)", "Page ", i
    do j=1,n
        print *, cubearray(i,j,:)
    enddo; enddo;
    linarrayA = reshape(cubearray,[n**3])
    linarrayB = linarrayA
    print *, ""
    print *, "linearized"
    print *, linarrayA
    print *, "removeud"
    removeudarray = remove_ud(linarrayB,n,n)
    print *, removeudarray
    print *, "constructud"
    linarrayB =  construct_ud(removeudarray,n,n)
    print *, linarrayB
    cubearray = reshape(linarrayB,[n,n,n])
    print *, "reconstructed"
    do i=1,n
    print "(A,I1)", "Page ", i
    do j=1,n
        print *, cubearray(i,j,:)
    enddo; enddo;
    print *, "removerb"
    removerbarray = remove_rb(linarrayA,n)
    print *, removerbarray
    print *, "constructrb"
    linarrayA = construct_rb(removerbarray,n)
    print *, linarrayA
    cubearray = reshape(linarrayA,[n,n,n])
    print *, "reconstructed"
    do i=1,n
    print "(A,I1)", "Page ", i
    do j=1,n
        print *, cubearray(i,j,:)
    enddo; enddo;
    packarray =  pack_STT(linarrayA,n)
    print *, "packed"
    print *, packarray
    linarrayB = unpack_STT(packarray,n)
    print *, "unpacked"
    print *, linarrayB
    cubearray = reshape(linarrayB,[n,n,n])
    print *, "reconstructed"
    do i=1,n
    print "(A,I1)", "Page ", i
    do j=1,n
        print *, cubearray(i,j,:)
    enddo; enddo;
    contains

    ! FUNCTION TO REMOVE UPPER DIAGONAL ELEMENTS OF EACH PAGE OF A 3D ARRAY 
    !(FOR SYMMETRIC PARTIALS)
    function remove_UD(linin, n,m) result(linout)
        implicit none
        integer, intent(in) :: n,m
        integer :: i,j,k
        character(len=5), dimension(n*m**2), intent(in)  :: linin
        character(len=5), dimension(n*m*(m+1)/2) :: linout
        i=1; k=0
        do j=1,size(linout)
            if (mod(i,m*n) == 1) then 
                i=i+k*n
                k=k+1
            endif
            linout(j) = linin(i)
            i=i+1
        enddo
    end function remove_UD

    ! FUNCTION TO RECONSTRUCT UPPER DIAGONAL ELEMENTS OF EACH PAGE OF A 3D ARRAY 
    ! (FOR SYMMETRIC PARTIALS)
    function construct_UD(linin, n,m) result(linout)
        implicit none
        integer, intent(in) :: n,m
        integer :: i,j,k,l
        character(len=5), dimension(n*m*(m+1)/2), intent(in)  :: linin
        character(len=5), dimension(n*m**2) :: linout
        linout = "0"
        l=1; j=1
        do k=1,m
        do j=k,m
        do i=1,n
            linout(ijktol(i,j,k,n,m)) = linin(l)
            if (j.ne.k) linout(ijktol(i,k,j,n,m)) = linout(ijktol(i,j,k,n,m))
            l=l+1
        enddo; enddo; enddo;
    end function construct_UD

    ! FUNCTION TO REMOVE RIGHTMOST COLUMN OF MATRIX 
    ! (ROWS 1-6 ARE ZERO, ROW 7 IS UNITY)
    function remove_R(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        character(len=5), dimension(n*n), intent(in)  :: linin
        character(len=5), dimension(n*(n-1)) :: linout
        linout = linin(:size(linout))
    end function remove_R

    ! FUNCTION TO REPLACE RIGHTMOST COLUMN OF MATRIX 
    ! (ROWS 1-6 ARE ZERO, ROW 7 IS UNITY)
    function construct_R(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        character(len=5), dimension(n*(n-1)), intent(in) :: linin
        character(len=5), dimension(n**2)  :: linout
        linout(:size(linin)) = linin
        linout(size(linin)+1:size(linout)-1) = "000  "
        linout(size(linout)) = "111  "
    end function construct_R

    ! FUNCTION TO REMOVE RIGHTMOST COLUMN AND BOTTOM ROW OF 
    ! EACH PAGE OF A 3D ARRAY (FOR STT PAGES)OF MATRIX 
    ! (ROWS 1-6 ARE ZERO, ROW 7 IS UNITY)
    function remove_RB(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        integer :: i,j
        character(len=5), dimension(:), intent(in)  :: linin
        character(len=5), dimension(n*(n-1)**2) :: linout
        i=1; j=1
        do j=1,size(linout)
            if (mod(i,n**2) == (n**2-n+1)) then
                i=i+n
            endif
            linout(j) = linin(i)
            i=i+1
        enddo
    end function remove_RB

    ! FUNCTION TO RECONSTRUCT RIGHTMOST COLUMN AND BOTTOM ROW OF 
    ! EACH PAGE OF A 3D ARRAY (FOR STT PAGES)OF MATRIX 
    ! (ROWS 1-6 ARE ZERO, ROW 7 IS UNITY)
    function construct_RB(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        integer :: i
        character(len=5), dimension(:), intent(in) :: linin
        character(len=5), dimension(n**3)   :: linout
        do i=0,n-2
        linout(i*n**2+1:i*n**2+n*(n-1)) = linin(i*n*(n-1)+1:(i+1)*n*(n-1))
        linout(i*n**2+n*(n-1)+1:(i+1)*n**2) = "000  "
        enddo
        linout(n**2*(n-1)+1:) = "000  "
    end function construct_RB

    function pack_STM(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        character(len=5), dimension(:), intent(in) :: linin
        character(len=5), dimension(n*(n-1))  :: linout
        linout = remove_R(linin,n)
    end function pack_STM

    function pack_STT(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        character(len=5), dimension(:), intent(in) :: linin
        character(len=5), dimension(n*(n-1)**2) :: linint
        character(len=5), dimension(n**2*(n-1)/2)  :: linout
        linint =  remove_RB(linin,n)
        linout = remove_UD(linint,n,n-1)
    end function pack_STT

    function unpack_STM(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        character(len=5), dimension(:), intent(in) :: linin
        character(len=5), dimension(n**2)  :: linout
        linout = construct_R(linin,n)
    end function unpack_STM

    function unpack_STT(linin,n) result(linout)
        implicit none
        integer, intent(in) :: n
        character(len=5), dimension(:), intent(in) :: linin
        character(len=5), dimension(n*(n-1)**2) :: linint
        character(len=5), dimension(n**3)  :: linout
        linint = construct_UD(linin,n,n-1)
        linout = construct_RB(linint,n)
    end function unpack_STT

    elemental function ijktol(i,j,k,n,m) result(l)
        implicit none
        integer, intent(in) :: i,j,k,n,m
        integer :: l
        l = i+(j-1)*n+(k-1)*m*n
    end function ijktol

    ! convert an unpacked (i,j,k) tuple to a linear
    ! index l which will get the correct array element
    elemental function utop(i,j,k,n) result(l)
        implicit none
        integer, intent(in) :: i,j,k,n
        integer :: l, l_o, zptr, uptr, p
        if ((j.eq.n).neqv.(k.eq.n)) then 
            l = zptr
        else if ((j.eq.n).and.(k.eq.n)) then
            l = uptr
        else
            if (j.le.k) then
                ! do the lower triangular thing (i.e. leave j and k alone)
                ! See postscript for explanation
                l = ijktol(i,j,k,n,n-1)
            else
                ! do the upper triangular thing (i.e. swap k and j)
                l = ijktol(i,k,j,n,n-1)
            endif
                l_o = l - sum([(p-1, p=1,k)])
        endif
    end function utop

end program array_indices

! i just tells you which page you're on. think about it like
! { [1 2 3]_1      X            X      }
! { [4 5 6]_2 [10 11 12]_4      X      }
! { [7 8 9]_3 [13 14 15]_5 [16 17 18]_6}
! so each block contains ((l-1)/n)*n+mod(l,n)

! { 1 X_4 X_7 }
! { 2 4_5 X_8 }
! { 3 5_6 6_9 }
! l_offset = l_true-sum_1^k((l-1)) to get the factor to add to the mod value
! (l_offset-1)/n*n + mod(l,n)
