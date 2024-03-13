program main
    use frkmain_d, only: Odesolution
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp=>real128
    implicit none
    type(Odesolution) :: qist_sol
    open(unit=75, file="qist_sol_perturbed.odesolution", access="stream", status="old")
    call qist_sol%read(75)
    close(75)
    print *, qist_sol%call(1._dp)
end program main
