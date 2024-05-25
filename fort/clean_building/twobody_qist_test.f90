program main
    use genqist, only: make_qist_model
    implicit none
    character(len=1000) :: arg
    call get_command_argument(1,arg)
    call make_qist_model(trim(adjustl(arg)))
end program main
