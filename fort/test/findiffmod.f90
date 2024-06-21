module findiffmod
    use, intrinsic :: iso_fortran_env, only: wp=>real128
    implicit none
    contains
function findiffhes(j, xstar,eps,order) result(res)
    implicit none
    interface
        function j(x)
            import
            ! This function should accept a state vector x
            ! and emit the jacobian matrix of a vector-valued function
            real(wp), intent(in) :: x(:)
            real(wp)             :: j(size(x),size(x))
        end function j
    end interface
    real(wp), intent(in) :: xstar(:), eps
    integer,  intent(in) :: order
    real(wp)             :: res(size(xstar,1),size(xstar,1), size(xstar,1))
    integer              :: i, dims
    ! Perturb initial conditions for deriv calculation
    real(wp)             :: per(size(xstar), size(xstar)), &
                          & bigpertmatrix(size(xstar),size(xstar)), &
                          & pstep(size(xstar),size(xstar)), &
                          & mstep(size(xstar),size(xstar)), &
                          & pdf(size(xstar),size(xstar),size(xstar)), &
                          & mdf(size(xstar),size(xstar),size(xstar))
    real(wp), allocatable:: &
                          & twopstep(:,:), &
                          & twomstep(:,:), &
                          & threepstep(:,:), &
                          & threemstep(:,:), &
                          & fourpstep(:,:), &
                          & fourmstep(:,:), &
                          & twopdf(:,:,:), &
                          & twomdf(:,:,:), &
                          & threepdf(:,:,:), &
                          & threemdf(:,:,:), &
                          & fourpdf(:,:,:), &
                          & fourmdf(:,:,:)
    ! Perturb initial conditions for deriv calculation
    dims = size(xstar)
    per = 0._wp
    do i=1,dims
        bigpertmatrix(:,i) = xstar
        per(i,i) = eps
    end do
    pstep = bigpertmatrix + per
    mstep = bigpertmatrix - per

    if (order > 3 ) then
        allocate(twopstep, twomstep, mold=per)
        allocate(twopdf, twomdf, mold=pdf)
        twopstep = pstep + per
        twomstep = mstep - per
    end if
    if (order > 5) then
        allocate(threepstep, threemstep, mold=per)
        allocate(threepdf, threemdf, mold=pdf)
        threepstep = twopstep + per
        threemstep = twomstep - per
    end if
    if (order > 7) then
        allocate(fourpstep, fourmstep, mold=per)
        allocate(fourpdf, fourmdf, mold=pdf)
        fourpstep = threepstep + per
        fourmstep = threemstep - per
    end if
    do i=1,dims
        pdf(:,:,i) = j(pstep(:,i))
        mdf(:,:,i) = j(mstep(:,i))
    end do
    if (order > 3 ) then
        do i=1,dims
            twopdf(:,:,i) = j(twopstep(:,i))
            twomdf(:,:,i) = j(twomstep(:,i))
        end do
    end if
    if (order > 5) then
        do i=1,dims
            threepdf(:,:,i) = j(threepstep(:,i))
            threemdf(:,:,i) = j(threemstep(:,i))
        end do
    end if
    if (order > 7) then
        do i=1,dims
            fourpdf(:,:,i) = j(fourpstep(:,i))
            fourmdf(:,:,i) = j(fourmstep(:,i))
        end do
    end if

    select case (order)
        case(3)
            res = (pdf - mdf) / (2._wp*eps)
        case(5)
            res = (8._wp*pdf - twopdf - 8._wp*mdf + twomdf)/(12._wp*eps)
        case(7)
            res = (45._wp*pdf - 9._wp*twopdf + threepdf - 45._wp*mdf + &
                &  9._wp*twomdf -threemdf)/(60._wp*eps)
        case(9)
            res = (3._wp*(fourmdf - fourpdf) + 32._wp*(threepdf-threemdf) + &
                &  168._wp*(twomdf-twopdf) + 672._wp*(pdf-mdf))/(840._wp*eps)
        case default
            print *, 'Warning: order not implemented. Using 3 point stencil'
            res = (pdf - mdf) /(2._wp*eps)
    end select
end function findiffhes
function findiff(f, xstar,eps,order) result(res)
    implicit none
    interface
        function f(x)
            import
            real(wp), intent(in) :: x(:)
            real(wp)             :: f(size(x))
        end function f
    end interface
    real(wp), intent(in) :: xstar(:), eps
    integer,  intent(in) :: order
    real(wp)             :: res(size(xstar),size(xstar)), &
                          & per(size(xstar), size(xstar)), &
                          & bigpertmatrix(size(xstar),size(xstar)), &
                          & pstep(size(xstar),size(xstar)), &
                          & mstep(size(xstar),size(xstar)), &
                          & pdf(size(xstar),size(xstar)), &
                          & mdf(size(xstar),size(xstar))
    real(wp), allocatable:: &
                          & twopstep(:,:), &
                          & twomstep(:,:), &
                          & twopdf(:,:), &
                          & twomdf(:,:), &
                          & threepstep(:,:), &
                          & threemstep(:,:), &
                          & threepdf(:,:), &
                          & threemdf(:,:), &
                          & fourpstep(:,:), &
                          & fourmstep(:,:), &
                          & fourpdf(:,:), &
                          & fourmdf(:,:)
    integer              :: i, dims
    ! Perturb initial conditions for deriv calculation
    dims = size(xstar)
    per = 0._wp
    do i=1,dims
        bigpertmatrix(:,i) = xstar
        per(i,i) = eps
    end do
    pstep = bigpertmatrix + per
    mstep = bigpertmatrix - per
    if (order > 3 ) then
        allocate(twopstep, twomstep, twopdf, twomdf, mold=per)
        twopstep = pstep + per
        twomstep = mstep - per
    end if
    if (order > 5) then
        allocate(threepstep, threemstep, threepdf, threemdf, mold=per)
        threepstep = twopstep + per
        threemstep = twomstep - per
    end if
    if (order > 7) then
        allocate(fourpstep, fourmstep, fourpdf, fourmdf, mold=per)
        fourpstep = threepstep + per
        fourmstep = threemstep - per
    end if
    do i=1,dims
        pdf(:,i) = f(pstep(:,i))
        mdf(:,i) = f(mstep(:,i))
    end do
    if (order > 3 ) then
        do i=1,dims
            twopdf(:,i) = f(twopstep(:,i))
            twomdf(:,i) = f(twomstep(:,i))
        end do
    end if
    if (order > 5) then
        do i=1,dims
            threepdf(:,i) = f(threepstep(:,i))
            threemdf(:,i) = f(threemstep(:,i))
        end do
    end if
    if (order > 7) then
        do i=1,dims
            fourpdf(:,i) = f(fourpstep(:,i))
            fourmdf(:,i) = f(fourmstep(:,i))
        end do
    end if

    select case (order)
        case(3)
            res = (pdf - mdf) / (2._wp*eps)
        case(5)
            res = (8._wp*pdf - twopdf - 8._wp*mdf + twomdf)/(12._wp*eps)
        case(7)
            res = (45._wp*pdf - 9._wp*twopdf + threepdf - 45._wp*mdf + &
                &  9._wp*twomdf -threemdf)/(60._wp*eps)
        case(9)
            res = (3._wp*(fourmdf - fourpdf) + 32._wp*(threepdf-threemdf) + &
                &  168._wp*(twomdf-twopdf) + 672._wp*(pdf-mdf))/(840._wp*eps)
        case default
            print *, 'Warning: order not implemented. Using 3 point stencil'
            res = (pdf - mdf) /(2._wp*eps)
        end select
end function findiff
function findiffmat_td(f, xstar,eps,order,shape_template) result(res)
    implicit none
    interface
        function f(x)
            import
            real(wp), intent(in)  :: x
            real(wp), allocatable :: f(:,:)
        end function f
    end interface
    real(wp), intent(in) :: xstar, eps, shape_template(:,:)
    integer,  intent(in) :: order
    real(wp)             :: res(size(shape_template,1),size(shape_template,2)), &
                          & refpoint, &
                          & pstep, &
                          & mstep, &
                          & twopstep, &
                          & twomstep, &
                          & threepstep, &
                          & threemstep, &
                          & fourpstep, &
                          & fourmstep, &
                          & pdf(size(shape_template,1),size(shape_template,2)), &
                          & mdf(size(shape_template,1),size(shape_template,2))
    real(wp)             :: &
                          & twopdf(size(shape_template,1),size(shape_template,2)), &
                          & twomdf(size(shape_template,1),size(shape_template,2)), &
                          & threepdf(size(shape_template,1),size(shape_template,2)), &
                          & threemdf(size(shape_template,1),size(shape_template,2)), &
                          & fourpdf(size(shape_template,1),size(shape_template,2)), &
                          & fourmdf(size(shape_template,1),size(shape_template,2))
    ! Perturb initial conditions for deriv calculation
    refpoint = xstar
    pstep = refpoint + eps
    mstep = refpoint - eps
    if (order > 3 ) then
        twopstep = pstep + eps
        twomstep = mstep - eps
    end if
    if (order > 5) then
        threepstep = twopstep + eps
        threemstep = twomstep - eps
    end if
    if (order > 7) then
        fourpstep = threepstep + eps
        fourmstep = threemstep - eps
    end if
        pdf = f(pstep)
        mdf = f(mstep)
    if (order > 3 ) then
        twopdf = f(twopstep)
        twomdf = f(twomstep)
    end if
    if (order > 5) then
        threepdf = f(threepstep)
        threemdf = f(threemstep)
    end if
    if (order > 7) then
        fourpdf = f(fourpstep)
        fourmdf = f(fourmstep)
    end if

    select case (order)
        case(3)
            res = (pdf - mdf) / (2._wp*eps)
        case(5)
            res = (8._wp*pdf - twopdf - 8._wp*mdf + twomdf)/(12._wp*eps)
        case(7)
            res = (45._wp*pdf - 9._wp*twopdf + threepdf - 45._wp*mdf + &
                &  9._wp*twomdf -threemdf)/(60._wp*eps)
        case(9)
            res = (3._wp*(fourmdf - fourpdf) + 32._wp*(threepdf-threemdf) + &
                &  168._wp*(twomdf-twopdf) + 672._wp*(pdf-mdf))/(840._wp*eps)
        case default
            print *, 'Warning: order not implemented. Using 3 point stencil'
            res = (pdf - mdf) /(2._wp*eps)
        end select
end function findiffmat_td

function findiff_multiscale(f, xstar, epsvec, order) result(res)
    implicit none
    interface
        function f(x)
            import
            real(wp), intent(in) :: x(:)
            real(wp)             :: f(size(x))
        end function f
    end interface
    real(wp), intent(in) :: xstar(:), epsvec(size(xstar))
    integer,  intent(in) :: order
    real(wp)             :: res(size(xstar),size(xstar)), &
                          & per(size(xstar), size(xstar)), &
                          & bigpertmatrix(size(xstar),size(xstar)), &
                          & pstep(size(xstar),size(xstar)), &
                          & mstep(size(xstar),size(xstar)), &
                          & pdf(size(xstar),size(xstar)), &
                          & mdf(size(xstar),size(xstar))
    real(wp), allocatable:: &
                          & twopstep(:,:), &
                          & twomstep(:,:), &
                          & twopdf(:,:), &
                          & twomdf(:,:), &
                          & threepstep(:,:), &
                          & threemstep(:,:), &
                          & threepdf(:,:), &
                          & threemdf(:,:), &
                          & fourpstep(:,:), &
                          & fourmstep(:,:), &
                          & fourpdf(:,:), &
                          & fourmdf(:,:)
    integer              :: i, dims
    ! Perturb initial conditions for deriv calculation
    res = 0._wp
    dims = size(xstar)
    per = 0._wp
    do i=1,dims
        bigpertmatrix(:,i) = xstar
        per(i,i) = epsvec(i)
    end do
    pstep = bigpertmatrix + per
    mstep = bigpertmatrix - per
    if (order > 3 ) then
        allocate(twopstep, twomstep, twopdf, twomdf, mold=per)
        twopstep = pstep + per
        twomstep = mstep - per
    end if
    if (order > 5) then
        allocate(threepstep, threemstep, threepdf, threemdf, mold=per)
        threepstep = twopstep + per
        threemstep = twomstep - per
    end if
    if (order > 7) then
        allocate(fourpstep, fourmstep, fourpdf, fourmdf, mold=per)
        fourpstep = threepstep + per
        fourmstep = threemstep - per
    end if
    do i=1,dims
        pdf(:,i) = f(pstep(:,i))
        mdf(:,i) = f(mstep(:,i))
    end do
    if (order > 3 ) then
        do i=1,dims
            twopdf(:,i) = f(twopstep(:,i))
            twomdf(:,i) = f(twomstep(:,i))
        end do
    end if
    if (order > 5) then
        do i=1,dims
            threepdf(:,i) = f(threepstep(:,i))
            threemdf(:,i) = f(threemstep(:,i))
        end do
    end if
    if (order > 7) then
        do i=1,dims
            fourpdf(:,i) = f(fourpstep(:,i))
            fourmdf(:,i) = f(fourmstep(:,i))
        end do
    end if

    select case (order)
        case(3)
            do i=1,dims
                res(:,i) = (pdf(:,i) - mdf(:,i)) / (2._wp*epsvec(i))
            end do
        case(5)
            do i=1,dims
            res(:,i) = (8._wp*pdf(:,i) - twopdf(:,i) - 8._wp*mdf(:,i) + twomdf(:,i))/(12._wp*epsvec(i))
            end do
        case(7)
            do i=1,dims
            res(:,i) = (45._wp*pdf(:,i) - 9._wp*twopdf(:,i) + threepdf(:,i) - 45._wp*mdf(:,i) + &
                &  9._wp*twomdf(:,i) -threemdf(:,i))/(60._wp*epsvec(i))
            end do
        case(9)
            do i=1,dims
            res(:,i) = (3._wp*(fourmdf(:,i) - fourpdf(:,i)) + 32._wp*(threepdf(:,i)-threemdf(:,i)) + &
                &  168._wp*(twomdf(:,i)-twopdf(:,i)) + 672._wp*(pdf(:,i)-mdf(:,i)))/(840._wp*epsvec(i))
            end do
        case default
            print *, 'Warning: order not implemented. Using 3 point stencil'
            do i=1,dims
            res(:,i) = (pdf(:,i) - mdf(:,i)) /(2._wp*epsvec(i))
            end do
        end select
end function findiff_multiscale
function findiffhes_multiscale(j, xstar,epsvec,order) result(res)
    implicit none
    interface
        function j(x)
            import
            ! This function should accept a state vector x
            ! and emit the jacobian matrix of a vector-valued function
            real(wp), intent(in) :: x(:)
            real(wp)             :: j(size(x),size(x))
        end function j
    end interface
    real(wp), intent(in) :: xstar(:), epsvec(size(xstar))
    integer,  intent(in) :: order
    real(wp)             :: res(size(xstar,1),size(xstar,1), size(xstar,1))
    integer              :: i, dims
    ! Perturb initial conditions for deriv calculation
    real(wp)             :: per(size(xstar), size(xstar)), &
                          & bigpertmatrix(size(xstar),size(xstar)), &
                          & pstep(size(xstar),size(xstar)), &
                          & mstep(size(xstar),size(xstar)), &
                          & pdf(size(xstar),size(xstar),size(xstar)), &
                          & mdf(size(xstar),size(xstar),size(xstar))
    real(wp), allocatable:: &
                          & twopstep(:,:), &
                          & twomstep(:,:), &
                          & threepstep(:,:), &
                          & threemstep(:,:), &
                          & fourpstep(:,:), &
                          & fourmstep(:,:), &
                          & twopdf(:,:,:), &
                          & twomdf(:,:,:), &
                          & threepdf(:,:,:), &
                          & threemdf(:,:,:), &
                          & fourpdf(:,:,:), &
                          & fourmdf(:,:,:)
    ! Perturb initial conditions for deriv calculation
    dims = size(xstar)
    per = 0._wp
    res = 0._wp
    do i=1,dims
        bigpertmatrix(:,i) = xstar
        per(i,i) = epsvec(i)
    end do
    pstep = bigpertmatrix + per
    mstep = bigpertmatrix - per

    if (order > 3 ) then
        allocate(twopstep, twomstep, mold=per)
        allocate(twopdf, twomdf, mold=pdf)
        twopstep = pstep + per
        twomstep = mstep - per
    end if
    if (order > 5) then
        allocate(threepstep, threemstep, mold=per)
        allocate(threepdf, threemdf, mold=pdf)
        threepstep = twopstep + per
        threemstep = twomstep - per
    end if
    if (order > 7) then
        allocate(fourpstep, fourmstep, mold=per)
        allocate(fourpdf, fourmdf, mold=pdf)
        fourpstep = threepstep + per
        fourmstep = threemstep - per
    end if
    do i=1,dims
        pdf(:,:,i) = j(pstep(:,i))
        mdf(:,:,i) = j(mstep(:,i))
    end do
    if (order > 3 ) then
        do i=1,dims
            twopdf(:,:,i) = j(twopstep(:,i))
            twomdf(:,:,i) = j(twomstep(:,i))
        end do
    end if
    if (order > 5) then
        do i=1,dims
            threepdf(:,:,i) = j(threepstep(:,i))
            threemdf(:,:,i) = j(threemstep(:,i))
        end do
    end if
    if (order > 7) then
        do i=1,dims
            fourpdf(:,:,i) = j(fourpstep(:,i))
            fourmdf(:,:,i) = j(fourmstep(:,i))
        end do
    end if

    select case (order)
        case(3)
            do i=1,dims
                res(:,:,i) = (pdf(:,:,i) - mdf(:,:,i)) / (2._wp*epsvec(i))
            end do
        case(5)
            do i=1,dims
                res(:,:,i) = (8._wp*pdf(:,:,i) - twopdf(:,:,i) - 8._wp*mdf(:,:,i) + twomdf(:,:,i))/(12._wp*epsvec(i))
            end do
        case(7)
            do i=1,dims
                res(:,:,i) = (45._wp*pdf(:,:,i) - 9._wp*twopdf(:,:,i) + threepdf(:,:,i) - 45._wp*mdf(:,:,i) + &
                    &  9._wp*twomdf(:,:,i) -threemdf(:,:,i))/(60._wp*epsvec(i))
            end do
        case(9)
            do i=1,dims
                res(:,:,i) = (3._wp*(fourmdf(:,:,i) - fourpdf(:,:,i)) + 32._wp*(threepdf(:,:,i)-threemdf(:,:,i)) + &
                    &  168._wp*(twomdf(:,:,i)-twopdf(:,:,i)) + 672._wp*(pdf(:,:,i)-mdf(:,:,i)))/(840._wp*epsvec(i))
            end do
        case default
            print *, 'Warning: order not implemented. Using 3 point stencil'
            do i=1,dims
                res(:,:,i) = (pdf(:,:,i) - mdf(:,:,i)) / (2._wp*epsvec(i))
            end do
    end select
end function findiffhes_multiscale
end module
