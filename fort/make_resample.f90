program main
    use genqist, only: make_spice_subset
    implicit none
    character(len=500) :: arg
    integer cacount
    cacount = command_argument_count()
    call get_command_argument(1,arg)
    call make_spice_subset(trim(adjustl(arg)),.false.)
end program main
