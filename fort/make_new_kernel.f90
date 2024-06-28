program main
    use genqist, only: make_spice_subset, &
                       make_rotation, &
                       generate_kernel
    implicit none
    character(len=500) :: arg
    call get_command_argument(1,arg)
    call make_spice_subset(trim(adjustl(arg)))
    call make_rotation(trim(adjustl(arg)))
    call generate_kernel(trim(adjustl(arg)))
end program main
