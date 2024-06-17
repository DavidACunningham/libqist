program main
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    use tensorops, only: vectens1, vectens3, mattens, vectensquad
    real(dp) :: mat(3,3), tens(3,3,3), vec(3), inter(3,3,3), &
                output(3,3)
    integer i

    ! these lines implement the operation
    ! ia, abc, jb -> ijc with RIU, HU, RIU
    mat(1,:) = [1., 2., 3.]
    mat(2,:) = [4., 5., 6.]
    mat(3,:) = [7., 8., 9.]

    tens(1,1,:) = [10., 19., 28.]
    tens(1,2,:) = [11., 20., 29.]
    tens(1,3,:) = [12., 21., 30.]
    tens(2,1,:) = [13., 22., 31.]
    tens(2,2,:) = [14., 23., 32.]
    tens(2,3,:) = [15., 24., 33.]
    tens(3,1,:) = [16., 25., 34.]
    tens(3,2,:) = [17., 26., 35.]
    tens(3,3,:) = [18., 27., 36.]

    vec = [37., 38., 39.]
    inter = mattens(mat,tens,3)
    inter = reshape(inter,[3,3,3],order=[2,1,3])
    inter = mattens(mat,inter,3)
    inter = reshape(inter,[3,3,3],order=[2,1,3])
    output =  vectens3(vec,inter,3)
    do i = 1,3
        print *, real(output(i,:),4)
    end do
    

end program main
