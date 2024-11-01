program main
    use genqist, only: make_spice_subset, &
                       make_rotation, &
                       generate_kernel_cheby
    implicit none
    character(len=500) :: arg, arg2, arg3
    integer cacount
    cacount = command_argument_count()
    call get_command_argument(1,arg)
    if (cacount>1) then
        call get_command_argument(2,arg2)
        if (trim(adjustl(arg2)).eq."True") then
            call make_spice_subset(trim(adjustl(arg)),.false.)
        endif
    endif
    if (cacount>2) then
        call get_command_argument(3,arg3)
        if (trim(adjustl(arg3)).eq."True") then
            call make_rotation(trim(adjustl(arg)))
        endif
    endif
    call generate_kernel_cheby(trim(adjustl(arg)))
end program main
