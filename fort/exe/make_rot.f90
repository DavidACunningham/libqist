program main
    use genqist, only: make_rotation
    implicit none
    character(len=500) :: arg
    integer cacount
    cacount = command_argument_count()
    call get_command_argument(1,arg)
    call make_rotation(trim(adjustl(arg)))
end program main
