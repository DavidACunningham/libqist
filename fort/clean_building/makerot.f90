program main
    use genqist, only: make_rotation
    character(len=1000) :: arg
    call get_command_argument(1,arg)
    call make_rotation(arg)
end program
