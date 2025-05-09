! Title: denselight.f90 
! Description:
!     Types, functions, and parameters for reading, writing, packing, unpacking,
!     and calling integrator results that include the STM and STT for QIST.
! 
! author: David Cunningham
! Last edited: See git log
module denseLight
    use frkmin, only: ODESolution, DOP853DenseOutput
    use, intrinsic :: iso_fortran_env, only: wp=>real64
    implicit none
    ! This is the integer array that maps the unpacked state
    ! to the packed state, removing symmetric second derivatives
    ! and constant values from the stored variables.
    integer, parameter :: qistpack(265) = [ &  
            1, 2, 3, 4, 5, 6, 7, 10, 11, 12, &
            13, 14, 15, 18, 19, 20, 21, 22, &
            23, 26, 27, 28, 29, 30, 31, 34, 35, &
            36, 37, 38, 39, 42, 43, 44, 45, 46, &
            47, 50, 51, 52, 53, 54, 55, 58, 59, &
            60, 61, 62, 63, 66, 67, 68, 69, 70, &
            71, 130, 131, 132, 133, 134, 135, 138, &
            139, 140, 141, 142, 143, 194, 195, 196, &
            197, 198, 199, 202, 203, 204, 205, 206, &
            207, 210, 211, 212, 213, 214, 215, 258, &
            259, 260, 261, 262, 263, 266, 267, 268, &
            269, 270, 271, 274, 275, 276, 277, 278, &
            279, 282, 283, 284, 285, 286, 287, 322, &
            323, 324, 325, 326, 327, 330, 331, 332, &
            333, 334, 335, 338, 339, 340, 341, 342, &
            343, 346, 347, 348, 349, 350, 351, 354, &
            355, 356, 357, 358, 359, 386, 387, 388, &
            389, 390, 391, 394, 395, 396, 397, 398, &
            399, 402, 403, 404, 405, 406, 407, 410, &
            411, 412, 413, 414, 415, 418, 419, 420, &
            421, 422, 423, 426, 427, 428, 429, 430, &
            431, 450, 451, 452, 453, 454, 455, 458, &
            459, 460, 461, 462, 463, 466, 467, 468, &
            469, 470, 471, 474, 475, 476, 477, 478, &
            479, 482, 483, 484, 485, 486, 487, 490, &
            491, 492, 493, 494, 495, 498, 499, 500, &
            501, 502, 503, 514, 515, 516, 517, 518, &
            519, 522, 523, 524, 525, 526, 527, 530, &
            531, 532, 533, 534, 535, 538, 539, 540, &
            541, 542, 543, 546, 547, 548, 549, 550, &
            551, 554, 555, 556, 557, 558, 559, 562, &
            563, 564, 565, 566, 567, 570, 571, 572, &
            573, 574, 575  &
                                          ], &
                          ncoeff = 7

    type lightSol
        ! A "light", i.e. trimmed down, implementation of the 
        ! ODESolution type, including read and write functions.
        integer               :: ndim, n_segments, fs(3)
        real(WP)              :: t_min, t_max 
        logical               :: ascending
        real(WP), allocatable :: t_olds(:), ts(:), ts_interp(:), &
                                 ts_sorted(:), &
                                 t_mins(:), t_maxs(:), hs(:), &
                                 !Note Frev in DOP853DenseOutput has shape (n_coeff, ndim)
                                 !For lightsol it is (ndim, n_coeff, n_segments)
                                 !to be in good column major order
                                 y_olds(:,:), Frevs(:,:,:), ys(:,:)
        contains
            procedure :: allocator
            generic :: convert => convert_from_file, &
                                  convert_from_object, &
                                  convert_from_object_and_pack, &
                                  convert_from_file_and_pack
            generic :: call => call_single, call_
            procedure :: init
            procedure :: convert_from_object
            procedure :: convert_from_file
            procedure :: convert_from_object_and_pack
            procedure :: convert_from_file_and_pack
            procedure :: write => LSwrite
            procedure :: read => LSread
            procedure :: call_single
            procedure :: call_
            procedure :: scall
            procedure :: call_impl
            procedure :: pack => selectElements
            procedure :: deallocator
    end type lightSol
    contains
    subroutine selectElements(self, elementlist)
        ! Keep only the elements in elementlist from the calling
        ! instance
        class(lightSol), intent(inout) :: self
        integer,         intent(in)    :: elementlist(:)
        integer                        :: newndim
        real(WP)                       :: new_y_olds(size(elementlist),self%n_segments), &
                                          new_ys(size(elementlist), self%n_segments+1), &
                                          newFrevs(0:size(elementlist)-1, &
                                                   0:self%fs(2)-1, &
                                                   self%fs(3))

        newndim = size(elementlist)
        new_y_olds = self%y_olds(elementlist, :)
        new_ys = self%ys(elementlist, :)
        newFrevs = self%Frevs(elementlist-1, :, :)
        self%ndim = size(elementlist)
        self%fs = [self%ndim, ncoeff, self%n_segments]
        deallocate(self%y_olds, self%ys, self%Frevs)
        allocate(self%y_olds, mold=new_y_olds)
        allocate(self%ys, mold=new_ys)
        allocate(self%Frevs(0:self%fs(1) - 1, &
                            0:self%fs(2) - 1, &
                            self%fs(3) &
                           ) &
                )
        self%y_olds = new_y_olds
        self%ys = new_ys
        self%Frevs = newFrevs
    end subroutine selectElements
    subroutine init(self, &
                    ndim, n_segments, &
                    t_min, t_max , &
                    ascending, &
                    t_olds, ts, ts_interp, &
                    ts_sorted, &
                    t_mins, t_maxs, hs, &
                    y_olds, Frevs, ys)

        ! Initialize the lightSol type
        class(lightSol), intent(inout) :: self
        integer        , intent(in)    :: ndim, n_segments
        real(WP)       , intent(in)    :: t_min, t_max 
        logical        , intent(in)    :: ascending
        real(WP),        intent(in)    :: t_olds(:), ts(:), ts_interp(:), &
                                          ts_sorted(:), &
                                          t_mins(:), t_maxs(:), hs(:), &
                                          y_olds(:,:), Frevs(:,:,:), ys(:,:)
        self%ndim= ndim
        self%n_segments = n_segments
        self%t_min = t_min
        self%t_max = t_max 
        self%ascending = ascending
        self%fs = [self%ndim, ncoeff, self%n_segments]

        call self%allocator()
        self%t_olds = t_olds
        self%ts     = ts
        self%ts_interp = ts_interp
        self%ts_sorted = ts_sorted
        self%t_mins = t_mins
        self%t_maxs = t_maxs
        self%hs     = hs
        self%y_olds = y_olds
        self%Frevs  = Frevs
        self%ys     = ys
    end subroutine init
    subroutine LSwrite(self,unit_num)
        ! Write the LightSol to disk
        class(lightSol), intent(inout) :: self
        integer,         intent(in)    :: unit_num
        write(unit_num) self%fs
        write(unit_num) self%ndim
        write(unit_num) self%n_segments
        write(unit_num) self%ascending
        write(unit_num) self%t_min
        write(unit_num) self%t_max
        write(unit_num) self%t_olds
        write(unit_num) self%ts
        write(unit_num) self%ts_interp
        write(unit_num) self%ts_sorted
        write(unit_num) self%t_mins
        write(unit_num) self%t_maxs
        write(unit_num) self%hs
        write(unit_num) self%y_olds
        write(unit_num) self%Frevs
        write(unit_num) self%ys
    end subroutine LSwrite
    subroutine LSread(self,unit_num)
        ! Read the LightSol from disk
        class(lightSol), intent(inout) :: self
        integer,         intent(in)    :: unit_num
        ! Read non-allocatable quantities
        read(unit_num) self%fs
        read(unit_num) self%ndim
        read(unit_num) self%n_segments
        read(unit_num) self%ascending
        read(unit_num) self%t_min
        read(unit_num) self%t_max

        ! Reallocate correct dimensions
        call self%allocator()

        ! read allocatable quantities
        read(unit_num) self%t_olds
        read(unit_num) self%ts
        read(unit_num) self%ts_interp
        read(unit_num) self%ts_sorted
        read(unit_num) self%t_mins
        read(unit_num) self%t_maxs
        read(unit_num) self%hs
        read(unit_num) self%y_olds
        read(unit_num) self%Frevs
        read(unit_num) self%ys
    end subroutine LSread
    subroutine convert_from_object(self, ODESol)
        ! Convert an ODESolution with dense output
        ! into a LightSol object
        class(lightSol), intent(inout) :: self
        type(ODESolution), intent(in) :: ODESol
        integer i
        if (.not.ODESol%dense) then
            print *, "Cannot convert sparse ODESolution to lightSol"
            error stop
        end if
        self%n_segments = size(ODESol%interpolants)
        self%ndim = ODESol%interpolants(1)%n
        self%ascending = ODESol%ascending
        self%t_min = real(ODESol%t_min,wp)
        self%t_max = real(ODESol%t_max,wp)

        self%fs = [self%ndim, ncoeff, self%n_segments]
        call self%allocator()
        self%ts_sorted = real(ODESol%ts_sorted,wp)
        self%ts =real(ODESol%ts,wp)
        self%ys = real(ODESol%ys, wp)
        do i=1, self%n_segments
            self%Frevs(:,:,i) = transpose(real(ODESol%interpolants(i)%Frev, wp))
            self%y_olds(:,i) =  real(ODESol%interpolants(i)%y_old, wp)
            self%t_olds(i) =    real(ODESol%interpolants(i)%t_old, wp)
            self%ts_interp(i) = real(ODESol%interpolants(i)%t, wp)
            self%t_mins(i) =    real(ODESol%interpolants(i)%t_min, wp)
            self%t_maxs(i) =    real(ODESol%interpolants(i)%t_max, wp)
            self%hs(i) =        real(ODESol%interpolants(i)%h, wp)
        end do
    end subroutine convert_from_object
    subroutine convert_from_object_and_pack(self, ODESol,elementlist)
        ! Convert an ODESolution with dense output
        ! into a LightSol object and discard some state elements
        class(lightSol),  intent(inout) :: self
        type(ODESolution),  intent(in) :: ODESol
        integer,            intent(in) :: elementlist(:)
        call self%convert_from_object(ODESol)
        call self%pack(elementlist)
    end subroutine convert_from_object_and_pack
    subroutine convert_from_file_and_pack(self,filename, elementlist)
        ! Convert an ODESolution with dense output from a file on disk
        ! into a LightSol object and discard some state elements
        class(lightSol),  intent(inout) :: self
        character(len=*),  intent(in) :: filename
        integer,            intent(in) :: elementlist(:)
        call self%convert_from_file(filename)
        call self%pack(elementlist)
    end subroutine convert_from_file_and_pack
    subroutine convert_from_file(self,filename)
        ! Convert an ODESolution with dense output from a file on disk
        ! into a LightSol object
        class(lightSol), intent(inout) :: self
        type(ODESolution) :: ODESol
        character(*), intent(in)    :: filename
        open(file=trim(adjustl(filename)), unit=73, access="stream", status="old")
        call ODESol%read(73)
        close(73)
        call self%convert_from_object(ODESol)
    end subroutine convert_from_file
    function call_single(self, t,lind,uind) result(res)
        ! Call some elements of the interpolant, at one value of t
        class(lightSol), intent(inout) :: self
        real(WP),           intent(in) :: t
        integer,            intent(in) :: lind,uind
        real(WP)                       :: res(self%ndim)
        integer                        :: ind, segment
        ! Here we preserve a certain symmetry that when t is in ts,
        ! then we prioritize a segment with a lower index.
        ind = binsearch(self%ts_sorted, t)

        segment = min(max(ind, 1), self%n_segments)
        if (.not. self%ascending) then
            segment = self%n_segments - 1 - segment
        end if

        res = self%scall(t,segment,lind,uind)
    end function
    function call_(self, t,lind,uind) result(res)
        ! Call some elements of the interpolant, at multiple values of t
        class(lightSol), intent(inout) :: self
        real(WP),           intent(in) :: t(:)
        integer,            intent(in) :: lind, uind
        real(WP)                       :: res(self%ndim,size(t))
        integer                        :: segments(size(t))
        integer iter

        !Evaluate the solution.

        !Parameters
        !----------
        !t : real (n_points)
        !    Points to evaluate at.
 
        !Returns
        !-------
        !y : real(n_states, n_points)
        !    Computed values.
        ! t must be sorted, unlike the scipy version

        !$OMP PARALLEL DO
        do iter=1,size(t)
            segments(iter) = binsearch(self%ts_sorted, t(iter))
        end do
        !$OMP END PARALLEL DO
        where (segments < 0._WP) segments = 1
        where (segments > self%n_segments - 1) segments  = self%n_segments
        if (.not. self%ascending) then
            segments = self%n_segments - 1 - segments
        endif

        !$OMP PARALLEL DO
        do iter=1,size(t)
            res(:,iter) = self%scall(t(iter),segments(iter),lind,uind)
        end do
        !$OMP END PARALLEL DO
    end function call_
    function scall(self, t,segment,lind,uind) result(res)
        class(lightSol), intent(inout) :: self
        real(WP),           intent(in) :: t
        integer,            intent(in) :: segment, lind, uind
        real(WP) :: resarray(self%ndim,1), res(self%ndim)
        resarray = self%call_impl([t],segment,lind,uind)
        res = resarray(:,1)
    end function scall
    function call_impl(self, t,segment,lind,uind) result(res)
        ! THe actual implementation of the interpolant calling
        class(lightSol),    intent(in) :: self
        integer,            intent(in) :: segment, lind, uind
        real(WP),           intent(in)    :: t(:)
        real(WP)                          :: res(uind-lind+1,size(t)), &
                                             x(size(t)), &
                                             y(uind-lind+1,size(t))
        integer i, j
        x = (t - self%t_olds(segment)) / self%hs(segment)

        y = 0._WP
        do i=0,size(self%Frevs,2)-1
            forall (j=1:size(t)) y(:,j) = y(:,j) + self%Frevs(lind-1:uind-1,i,segment)
            if (mod(i,2)== 0) then
                forall (j=1:size(t)) y(:,j) = y(:,j)*x(j)
            else
                forall (j=1:size(t)) y(:,j) = y(:,j)*(1._WP - x(j))
            endif
        end do
        !$OMP PARALLEL DO
        do i=1,size(t)
            y(:,i) = y(:,i) + self%y_olds(lind:uind,segment)
        end do
        !$OMP END PARALLEL DO
        res = y
    end function
    function binsearch(array,val) result(res)
        ! Binary search in a sorted array for the interval containing val
        real(WP), intent(in) :: array(:), val
        real(WP)             :: workarray(size(array))
        logical              :: flipped
        integer :: j, top, bot, n, res
        n = size(array)
        if (array(1).lt.array(n)) then
            workarray = array
            flipped = .false.
        else
            workarray = array(n:1:-1)
            flipped = .true.
        end if
        res = -1
        bot= 1
        top= n
        j = (bot+top)/2
        do while (bot.lt.top)
            if (val.lt.workarray(j)) then
                top= j
            else if (workarray(j+1).lt.val) then
                bot= j+1
            else  
                res = j 
                if (val.eq.workarray(j+1)) res = j+1;
                if (flipped) res = n-res + 1
                exit
            end if 
            j = (bot+ top)/2
        end do
        if (bot==n) res = n-1
    end function binsearch
    subroutine deallocator(self)
        ! Convenience function for deallocating all allocated 
        ! lightSol variables
        class(lightSol), intent(inout) :: self
        if(allocated(self%t_olds)) deallocate(self%t_olds)
        if(allocated(self%ts)) deallocate(self%ts)
        if(allocated(self%ts_interp)) deallocate(self%ts_interp)
        if(allocated(self%ts_sorted)) deallocate(self%ts_sorted)
        if(allocated(self%t_mins)) deallocate(self%t_mins)
        if(allocated(self%t_maxs)) deallocate(self%t_maxs)
        if(allocated(self%hs)) deallocate(self%hs)
        if(allocated(self%y_olds)) deallocate(self%y_olds)
        if(allocated(self%Frevs)) deallocate(self%Frevs)
        if(allocated(self%ys)) deallocate(self%ys)
    end subroutine deallocator
    subroutine allocator(self)
        ! Convenience subroutine for allocating all allocatable
        ! LightSol variables
        class(lightSol), intent(inout) :: self
        call self%deallocator()
        allocate( &
                 self%t_olds(self%n_segments), &
                 self%ts(self%n_segments+1), &
                 self%ts_interp(self%n_segments), &
                 self%ts_sorted(self%n_segments+1), &
                 self%t_mins(self%n_segments), &
                 self%t_maxs(self%n_segments),&
                 self%hs(self%n_segments), &
                 self%y_olds(self%ndim,self%n_segments), &
                 self%Frevs(0:self%fs(1)-1,0:self%fs(2)-1,self%fs(3)), &
                 self%ys(self%ndim,self%n_segments+1) &
                 )
     end subroutine allocator
end module denseLight
