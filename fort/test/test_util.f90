module test_util
    use, intrinsic :: iso_fortran_env, only: dp=>real64, qp=>real128
    implicit none
    contains
        subroutine mprint(mat)
            class(*), intent(in) :: mat(:,:)
            integer i
            select type(mat)
            type is (real(dp))
                do i= 1,size(mat,1)
                    write (*,*) mat(i,:)
                end do
            type is (real(qp))
                do i= 1,size(mat,1)
                    write (*,*) mat(i,:)
                end do
            end select
        end subroutine

        subroutine tprint(tens)
            class(*), intent(in) :: tens(:,:,:)
            integer i,j
            select type(tens)
            type is (real(dp))
                do i= 1,size(tens,1)
                do j= 1,size(tens,2)
                    write (*,*) tens(i,j,:)
                end do
                end do
            type is (real(qp))
                do i= 1,size(tens,1)
                do j= 1,size(tens,2)
                    write (*,*) tens(i,j,:)
                end do
                end do
            end select
        end subroutine
end module test_util
