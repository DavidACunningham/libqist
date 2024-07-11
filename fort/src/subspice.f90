module subspice
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use cheby
    implicit none
    type spice_subset
        integer               :: nbods, ndeg, central_body
        integer, allocatable  :: bodlist(:)
        real(dp)              :: a, b
        real(dp), allocatable :: pcoeffs(:,:), & !(ndeg,3*nbods)
                               & vcoeffs(:,:), & !(ndeg,3*nbods)
                               & acoeffs(:,:)    !(ndeg,3*nbods)
        contains
        procedure          :: init => fitspice
        procedure          :: write => spicewrite
        procedure          :: read => spiceread
        procedure, private :: call_spice_one
        procedure, private :: call_spice_sev
        generic, public    :: call => call_spice_one, call_spice_sev
    end type spice_subset
    contains
        subroutine spicewrite(me,unit_num)
            class(spice_subset), intent(in) :: me
            integer,             intent(in) :: unit_num
            write(unit_num) me%nbods
            write(unit_num) me%ndeg
            write(unit_num) me%central_body
            write(unit_num) me%bodlist
            write(unit_num) me%a
            write(unit_num) me%b
            write(unit_num) me%pcoeffs
            write(unit_num) me%vcoeffs
            write(unit_num) me%acoeffs
        end subroutine spicewrite
        subroutine spiceread(me,unit_num)
            class(spice_subset), intent(inout) :: me
            integer,             intent(in) :: unit_num
            read(unit_num) me%nbods
            read(unit_num) me%ndeg
            read(unit_num) me%central_body
            allocate(me%bodlist(me%nbods))
            read(unit_num) me%bodlist
            read(unit_num) me%a
            read(unit_num) me%b
            allocate(me%pcoeffs(me%ndeg,3*me%nbods), &
                   & me%vcoeffs(me%ndeg,3*me%nbods), &
                   & me%acoeffs(me%ndeg,3*me%nbods))
            read(unit_num) me%pcoeffs
            read(unit_num) me%vcoeffs
            read(unit_num) me%acoeffs
            me%acoeffs(me%ndeg-1:,:) = 0._dp
        end subroutine spiceread
        subroutine fitspice(me, kernelfile, central_body, bodlist, a, b, ndeg)
            class(spice_subset), intent(inout) :: me
            character(len=*),  intent(in)    :: kernelfile
            integer,             intent(in)    :: central_body, &
                                                & bodlist(:), &
                                                & ndeg
            real(dp),            intent(in)    :: a,b
            real(dp)                           :: nodes(ndeg), &
                                                & pfi(ndeg,3*(size(bodlist)+1)), &
                                                & vfi(ndeg,3*(size(bodlist)+1)), &
                                                & spkgeo_out(6), lt_dum
            integer i,j
            ! load kernel
            call FURNSH(trim(adjustl(kernelfile)))
            if (allocated(me%pcoeffs)) deallocate(me%pcoeffs)
            if (allocated(me%vcoeffs)) deallocate(me%vcoeffs)
            if (allocated(me%acoeffs)) deallocate(me%acoeffs)
            if (allocated(me%bodlist)) deallocate(me%bodlist)
            ! assign type properties
            me%nbods = size(bodlist) + 1
            me%ndeg  = ndeg
            me%a     = a
            me%b     = b
            me%central_body = central_body
            ! allocate bodlist, p/v/acoeffs
            allocate(me%pcoeffs(ndeg,3*me%nbods), &
                     me%vcoeffs(ndeg,3*me%nbods), &
                     me%acoeffs(ndeg,3*me%nbods), &
                     me%bodlist(me%nbods))
            ! assign bodies
            me%bodlist(1) = central_body
            me%bodlist(2:) = bodlist
            ! get chebyshev nodes
            nodes = chnodes(me%ndeg, me%a, me%b)
            !  for each node, do:
            do j = 1,me%ndeg
                ! SPKGEO call at time node(i) for pos/vel of central_body 
                ! with observer solar_system_barycenter
                call spkgeo(me%bodlist(1), nodes(j), "J2000", 0, &
                                  & spkgeo_out, lt_dum)
                pfi(j,:3) = spkgeo_out(:3)
                vfi(j,:3) = spkgeo_out(4:)
            end do
            ! for each body in bodlist, do:
            do i=2,me%nbods
            !  for each node, do:
                do j=1,me%ndeg
                    ! SPKGEO call at time node(i) for pos/vel of body(i) with
                    ! observer central_body
                    call spkgeo(me%bodlist(i), nodes(j), "J2000", me%bodlist(1), &
                                      & spkgeo_out, lt_dum)
                    pfi(j,3*(i-1)+1:3*i) = spkgeo_out(:3)
                    vfi(j,3*(i-1)+1:3*i) = spkgeo_out(4:)
                end do
            end do
            ! Fit all coeffs
            do i=1, 3*me%nbods
                me%pcoeffs(:,i) = chfit(me%ndeg,pfi(:,i))
                me%vcoeffs(:,i) = chfit(me%ndeg,vfi(:,i))
                me%acoeffs(:,i) = chderiv(me%vcoeffs(:,i),me%a,me%b)
            end do
        end subroutine fitspice
        function call_spice_one(me,x, bod_id,mode) result(res)
            class(spice_subset), intent(in) :: me
            real(dp),            intent(in) :: x
            integer,             intent(in) :: bod_id
            character(len=1),    intent(in) :: mode
            real(dp)                        :: res(3)
            integer                         :: i
            i = findloc(me%bodlist,bod_id,dim=1) - 1
            select case(mode)
            case ("p")
                res(1) = chcall(me%a,me%b,me%pcoeffs(:,3*i+1),x)
                res(2) = chcall(me%a,me%b,me%pcoeffs(:,3*i+2),x)
                res(3) = chcall(me%a,me%b,me%pcoeffs(:,3*i+3),x)
            case ("v")
                res(1) = chcall(me%a,me%b,me%vcoeffs(:,3*i+1),x)
                res(2) = chcall(me%a,me%b,me%vcoeffs(:,3*i+2),x)
                res(3) = chcall(me%a,me%b,me%vcoeffs(:,3*i+3),x)
            case ("a")
                res(1) = chcall(me%a,me%b,me%acoeffs(:,3*i+1),x)
                res(2) = chcall(me%a,me%b,me%acoeffs(:,3*i+2),x)
                res(3) = chcall(me%a,me%b,me%acoeffs(:,3*i+3),x)
            end select
        end function
        function call_spice_sev(me,x,bod_id,mode) result(res)
            class(spice_subset), intent(in) :: me
            real(dp),            intent(in) :: x(:)
            integer,             intent(in) :: bod_id
            character(len=1),    intent(in) :: mode
            real(dp)                        :: res(size(x),3)
            integer                         :: i
            i = findloc(me%bodlist,bod_id,dim=1) - 1
            select case(mode)
            case ("p")
                res(:,1) = chcall(me%a,me%b,me%pcoeffs(:,3*i+1),x)
                res(:,2) = chcall(me%a,me%b,me%pcoeffs(:,3*i+2),x)
                res(:,3) = chcall(me%a,me%b,me%pcoeffs(:,3*i+3),x)
            case ("v")
                res(:,1) = chcall(me%a,me%b,me%vcoeffs(:,3*i+1),x)
                res(:,2) = chcall(me%a,me%b,me%vcoeffs(:,3*i+2),x)
                res(:,3) = chcall(me%a,me%b,me%vcoeffs(:,3*i+3),x)
            case ("a")
                res(:,1) = chcall(me%a,me%b,me%acoeffs(:,3*i+1),x)
                res(:,2) = chcall(me%a,me%b,me%acoeffs(:,3*i+2),x)
                res(:,3) = chcall(me%a,me%b,me%acoeffs(:,3*i+3),x)
            end select
        end function call_spice_sev
end module subspice
