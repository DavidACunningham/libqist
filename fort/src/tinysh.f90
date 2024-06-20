! pinesmodule.f90 - Pines spherical harmonics evaluation 
!
! references:
! Fantino, E., and Casotto, S. "Methods of harmonics synthesis for global geopotential models
! and their first-, second-, and third-order gradients." Journal of Geodesy, 2009.
! 
! author: Sean McArdle
! updated by : David Cunningham 5/10/2024
module tinysh
    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter, public :: ap=real64
    integer, parameter, public :: sp=kind(1e0)
    real(ap),parameter :: pi = acos(-1.0_ap)
    real(ap),parameter :: eye3x3(3,3)=reshape( &
        [[1.0_ap,0.0_ap,0.0_ap],&
        [0.0_ap,1.0_ap,0.0_ap],&
        [0.0_ap,0.0_ap,1.0_ap]],[3,3])

    type PinesData
        real(ap),allocatable,dimension(:) :: rhol
        real(ap),allocatable,dimension(:) :: cosmlam
        real(ap),allocatable,dimension(:) :: sinmlam
        real(ap),allocatable,dimension(:) :: cosphipow
        ! Precomputed constants
        real(ap),allocatable,dimension(:,:) :: gml_pc
        real(ap),allocatable,dimension(:,:) :: hml_pc
        ! Hemholtz polynomials and derivatives
        real(ap),allocatable,dimension(:,:) :: Hml
        real(ap),allocatable,dimension(:,:) :: dHml
        real(ap),allocatable,dimension(:,:) :: d2Hml
        real(ap),allocatable,dimension(:,:) :: d3Hml
        ! Auxiliary functions
        real(ap),allocatable,dimension(:,:) :: Lml
        real(ap),allocatable,dimension(:,:) :: dLml
        real(ap),allocatable,dimension(:,:) :: d2Lml
        real(ap),allocatable,dimension(:,:) :: Oml
        real(ap),allocatable,dimension(:,:) :: Qml
        real(ap),allocatable,dimension(:,:) :: Wml
        ! Vectorized working memory for speed
        real(ap),allocatable,dimension(:) :: Cvec
        real(ap),allocatable,dimension(:) :: Svec
        real(ap),allocatable,dimension(:) :: Cfactor
        real(ap),allocatable,dimension(:) :: Sfactor
        ! Lumped coefficients
        real(ap),allocatable,dimension(:) :: A1m, B1m
        real(ap),allocatable,dimension(:) :: A2m, B2m
        real(ap),allocatable,dimension(:) :: A3m, B3m
        real(ap),allocatable,dimension(:) :: A4m, B4m
        real(ap),allocatable,dimension(:) :: A5m, B5m
        real(ap),allocatable,dimension(:) :: A6m, B6m
        real(ap),allocatable,dimension(:) :: A7m, B7m
        real(ap),allocatable,dimension(:) :: A8m, B8m
        real(ap),allocatable,dimension(:) :: A9m, B9m
        real(ap),allocatable,dimension(:) :: A10m, B10m
        logical :: pinesallocated=.false.
    end type PinesData

    contains 

    ! pinesinit - allocate working memory for pines
    subroutine pinesinit(degmax, Cml, Sml, pdat)

        implicit none

        ! Input/output
        integer,intent(in) :: degmax
        real(ap),allocatable,dimension(:,:),intent(in) :: Cml ! C terms organized into order/degree
        real(ap),allocatable,dimension(:,:),intent(in) :: Sml ! S terms organized into order/degree
        real(ap),allocatable,dimension(:,:) :: Cml_buf, Sml_buf
        type(PinesData),intent(inout) :: pdat

        ! Local
        integer :: l,m
        integer :: idx

        call deallocatepines(pdat)

        if(.not.allocated(pdat%rhol)) allocate(pdat%rhol(degmax+1))
        if(.not.allocated(pdat%cosmlam)) allocate(pdat%cosmlam(degmax+2))
        if(.not.allocated(pdat%sinmlam)) allocate(pdat%sinmlam(degmax+2))
        if(.not.allocated(pdat%cosphipow)) allocate(pdat%cosphipow(degmax+2))
        if(.not.allocated(pdat%Hml)) allocate(pdat%Hml(degmax+4,degmax+4))
        if(.not.allocated(pdat%gml_pc)) allocate(pdat%gml_pc(degmax+4,degmax+4))
        if(.not.allocated(pdat%hml_pc)) allocate(pdat%hml_pc(degmax+4,degmax+4))
        if(.not.allocated(pdat%dHml)) allocate(pdat%dHml(degmax+1,degmax+1))
        if(.not.allocated(pdat%d2Hml)) allocate(pdat%d2Hml(degmax+1,degmax+1))
        if(.not.allocated(pdat%d3Hml)) allocate(pdat%d3Hml(degmax+1,degmax+1))
        if(.not.allocated(pdat%Lml)) allocate(pdat%Lml(degmax+1,degmax+1))
        if(.not.allocated(pdat%dLml)) allocate(pdat%dLml(degmax+1,degmax+1))
        if(.not.allocated(pdat%d2Lml)) allocate(pdat%d2Lml(degmax+1,degmax+1))
        if(.not.allocated(pdat%Oml)) allocate(pdat%Oml(degmax+1,degmax+1))
        if(.not.allocated(pdat%Qml)) allocate(pdat%Qml(degmax+1,degmax+1))
        if(.not.allocated(pdat%Wml)) allocate(pdat%Wml(degmax+1,degmax+1))
        if(.not.allocated(pdat%Cvec)) allocate(pdat%Cvec((degmax+1+2)*(degmax+1+1)/2))
        if(.not.allocated(pdat%Svec)) allocate(pdat%Svec((degmax+1+2)*(degmax+1+1)/2))
        if(.not.allocated(pdat%Cfactor)) allocate(pdat%Cfactor((degmax+1+2)*(degmax+1+1)/2))
        if(.not.allocated(pdat%Sfactor)) allocate(pdat%Sfactor((degmax+1+2)*(degmax+1+1)/2))
        if(.not.allocated(pdat%A1m)) allocate(pdat%A1m(degmax+2))
        if(.not.allocated(pdat%A2m)) allocate(pdat%A2m(degmax+2))
        if(.not.allocated(pdat%A3m)) allocate(pdat%A3m(degmax+2))
        if(.not.allocated(pdat%A4m)) allocate(pdat%A4m(degmax+2))
        if(.not.allocated(pdat%A5m)) allocate(pdat%A5m(degmax+2))
        if(.not.allocated(pdat%A6m)) allocate(pdat%A6m(degmax+2))
        if(.not.allocated(pdat%A7m)) allocate(pdat%A7m(degmax+2))
        if(.not.allocated(pdat%A8m)) allocate(pdat%A8m(degmax+2))
        if(.not.allocated(pdat%A9m)) allocate(pdat%A9m(degmax+2))
        if(.not.allocated(pdat%A10m)) allocate(pdat%A10m(degmax+2))
        if(.not.allocated(pdat%B1m)) allocate(pdat%B1m(degmax+2))
        if(.not.allocated(pdat%B2m)) allocate(pdat%B2m(degmax+2))
        if(.not.allocated(pdat%B3m)) allocate(pdat%B3m(degmax+2))
        if(.not.allocated(pdat%B4m)) allocate(pdat%B4m(degmax+2))
        if(.not.allocated(pdat%B5m)) allocate(pdat%B5m(degmax+2))
        if(.not.allocated(pdat%B6m)) allocate(pdat%B6m(degmax+2))
        if(.not.allocated(pdat%B7m)) allocate(pdat%B7m(degmax+2))
        if(.not.allocated(pdat%B8m)) allocate(pdat%B8m(degmax+2))
        if(.not.allocated(pdat%B9m)) allocate(pdat%B9m(degmax+2))
        if(.not.allocated(pdat%B10m)) allocate(pdat%B10m(degmax+2))

        ! Store Stokes coefficients in vectorized form for faster access
        allocate(Cml_buf(size(Cml,1)+1,size(Cml,2)+1))
        allocate(Sml_buf(size(Sml,1)+1,size(Sml,2)+1))
        Cml_buf = 0._ap
        Sml_buf = 0._ap
        Cml_buf(:size(Cml,1),:size(Cml,2)) = Cml
        Sml_buf(:size(Sml,1),:size(Sml,2)) = Sml
        idx = 0
        pdat%Cvec = 0._ap
        pdat%Svec = 0._ap
        do m = 0,degmax
            do l = m,degmax
                idx = idx + 1
                pdat%Cvec(idx) = Cml_buf(m+1,l+1)
                pdat%Svec(idx) = Sml_buf(m+1,l+1)
            end do
        end do

        ! Perform precomputations
        do l = 2,degmax+3
            do m = 0,l-2
                pdat%gml_pc(m+1,l+1) = sqrt((2.0_ap*l+1.0_ap)*(2.0_ap*l-1.0_ap)/(l+m)/(l-m))
                pdat%hml_pc(m+1,l+1) = sqrt((2.0_ap*l+1.0_ap)*(l-m-1.0_ap)*(l+m-1.0_ap)/(2.0_ap*l-3.0_ap)/(l+m)/(l-m))
            end do
        end do

        pdat%pinesallocated=.true.

    end subroutine pinesinit


    ! deallocatepines - deallocate working memory for pines
    subroutine deallocatepines(pdat)

        implicit none

        ! Input/output
        type(PinesData),intent(inout) :: pdat

        if (.not.pdat%pinesallocated) then
            print *, "WARNING: pines working memory not allocated, exiting (deallocatepines)"
            return
        end if

        if(allocated(pdat%Cvec)) deallocate(pdat%Cvec)
        if(allocated(pdat%Svec)) deallocate(pdat%Svec)
        if(allocated(pdat%Cfactor)) deallocate(pdat%Cfactor)
        if(allocated(pdat%Sfactor)) deallocate(pdat%Sfactor)
        if(allocated(pdat%rhol)) deallocate(pdat%rhol)
        if(allocated(pdat%cosmlam)) deallocate(pdat%cosmlam)
        if(allocated(pdat%sinmlam)) deallocate(pdat%sinmlam)
        if(allocated(pdat%cosphipow)) deallocate(pdat%cosphipow)
        if(allocated(pdat%Hml)) deallocate(pdat%Hml)
        if(allocated(pdat%gml_pc)) deallocate(pdat%gml_pc)
        if(allocated(pdat%hml_pc)) deallocate(pdat%hml_pc)
        if(allocated(pdat%dHml)) deallocate(pdat%dHml)
        if(allocated(pdat%d2Hml)) deallocate(pdat%d2Hml)
        if(allocated(pdat%d3Hml)) deallocate(pdat%d3Hml)
        if(allocated(pdat%Lml)) deallocate(pdat%Lml)
        if(allocated(pdat%dLml)) deallocate(pdat%dLml)
        if(allocated(pdat%d2Lml)) deallocate(pdat%d2Lml)
        if(allocated(pdat%Oml)) deallocate(pdat%Oml)
        if(allocated(pdat%Qml)) deallocate(pdat%Qml)
        if(allocated(pdat%Wml)) deallocate(pdat%Wml)
        if(allocated(pdat%A1m)) deallocate(pdat%A1m)
        if(allocated(pdat%A2m)) deallocate(pdat%A2m)
        if(allocated(pdat%A3m)) deallocate(pdat%A3m)
        if(allocated(pdat%A4m)) deallocate(pdat%A4m)
        if(allocated(pdat%A5m)) deallocate(pdat%A5m)
        if(allocated(pdat%A6m)) deallocate(pdat%A6m)
        if(allocated(pdat%A7m)) deallocate(pdat%A7m)
        if(allocated(pdat%A8m)) deallocate(pdat%A8m)
        if(allocated(pdat%A9m)) deallocate(pdat%A9m)
        if(allocated(pdat%A10m)) deallocate(pdat%A10m)
        if(allocated(pdat%B1m)) deallocate(pdat%B1m)
        if(allocated(pdat%B2m)) deallocate(pdat%B2m)
        if(allocated(pdat%B3m)) deallocate(pdat%B3m)
        if(allocated(pdat%B4m)) deallocate(pdat%B4m)
        if(allocated(pdat%B5m)) deallocate(pdat%B5m)
        if(allocated(pdat%B6m)) deallocate(pdat%B6m)
        if(allocated(pdat%B7m)) deallocate(pdat%B7m)
        if(allocated(pdat%B8m)) deallocate(pdat%B8m)
        if(allocated(pdat%B9m)) deallocate(pdat%B9m)
        if(allocated(pdat%B10m)) deallocate(pdat%B10m)

        pdat%pinesallocated = .false.

    end subroutine deallocatepines


    ! shpines - Pines nonsingular spherical harmonics evaluation
    subroutine shpines(Rbody, GM, pdat, degmax, dorder, cart, V, dV, d2V, d3V)

        implicit none

        ! Input/output
        real(ap),intent(in) :: RBody
        real(ap),intent(in) :: GM
        type(PinesData),intent(inout) :: pdat
        integer,intent(in) :: degmax
        integer,intent(in) :: dorder
        real(ap),intent(in) :: cart(3)
        real(ap),intent(out) :: V
        real(ap),intent(out) :: dV(3)
        real(ap),intent(out) :: d2V(3,3)
        real(ap),intent(out) :: d3V(3,3,3)

        ! Local
        real(ap) :: s, t, u ! Pines coordinates
        real(ap) :: r2
        real(ap) :: oneByr
        real(ap) :: oneByr2
        real(ap) :: oneByr3
        real(ap) :: RByr
        real(ap) :: a1, a2, a3, a4
        real(ap) :: a11, a13, a22, a24, a34, a12, a14, a23, a33, a44
        real(ap) :: a111, a113, a123, a143, a222, a224, a331, a333, a441, a443
        real(ap) :: a112, a114, a124, a221, a223, a243, a332, a334, a442, a444
        real(ap) :: Vx, Vy, Vz ! First derivatives
        real(ap) :: Vxx, Vxy, Vyy, Vxz, Vyz, Vzz  ! Second derivatives
        real(ap) :: Vxxx, Vxxy, Vxxz, Vyyx, Vyyy, Vyyz, Vzzx, Vzzy, Vzzz, Vxyz ! Third derivatives
        integer :: l, m ! Degree/order indices
        real(ap) :: knm, knmp1, knmp2 ! Normalization factors for Hemholtz derivatives
        real(ap) :: lam, coslam, sinlam
        real(ap) :: phi, cosphi
        integer :: idx

        ! DAC: Zero initialize outputs
        ! V = 0._ap
        ! dV = 0._ap
        ! d2V = 0._ap
        ! d3V = 0._ap

        pdat%rhol = 0.0_ap
        pdat%cosmlam = 0.0_ap
        pdat%sinmlam = 0.0_ap
        pdat%cosphipow = 0.0_ap
        pdat%Hml = 0.0_ap

        r2 = cart(1)*cart(1) + cart(2)*cart(2) + cart(3)*cart(3)
        oneByr2 = 1.0_ap/r2
        oneByr = sqrt(oneByr2)
        oneByr3 = oneByr*oneByr2
        RByr = oneByr*RBody
        s = cart(1)*oneByr
        t = cart(2)*oneByr
        u = cart(3)*oneByr

        ! Recursively evaluate rhol, Recursively evaluate cosmlam, sinmlam, cosphipow
        ! Note we can use l/m interchangably here
        pdat%rhol(1) = GM*oneByr
        lam = atan2(cart(2),cart(1))
        coslam = cos(lam)
        sinlam = sin(lam)
        pdat%cosmlam(1) = 1.0_ap
        pdat%sinmlam(1) = 0.0_ap
        phi = asin(cart(3)*oneByr)
        cosphi = cos(phi)
        pdat%cosphipow(1) = 1.0_ap
        do l=1,degmax
            pdat%rhol(l+1) = RByr*pdat%rhol(l)
            pdat%cosmlam(l+1) = coslam*pdat%cosmlam(l) - sinlam*pdat%sinmlam(l)
            pdat%sinmlam(l+1) = coslam*pdat%sinmlam(l) + sinlam*pdat%cosmlam(l)
            pdat%cosphipow(l+1) = cosphi*pdat%cosphipow(l)
        end do

        ! Evaluate Hemholtz polynomials
        ! Need to go beyond degmax for derivatives
        pdat%Hml(1,1) = 1.0_ap
        pdat%Hml(2,2) = sqrt(3.0_ap)*pdat%Hml(1,1) ! Sectorial
        pdat%Hml(1,2) = u*pdat%Hml(2,2) ! Sub-diagonal
        do l = 2,degmax+3
            pdat%Hml(l+1,l+1) = sqrt(0.5_ap*(2.0_ap*l + 1.0_ap)/l)*pdat%Hml(l,l) ! Sectorial
            pdat%Hml(l,l+1) = u*sqrt(2.0_ap*l)*pdat%Hml(l+1,l+1) ! Sub-diagonal
            do m = 0,l-2
                pdat%Hml(m+1,l+1) = u*pdat%gml_pc(m+1,l+1)*pdat%Hml(m+1,l) - pdat%hml_pc(m+1,l+1)*pdat%Hml(m+1,l-1)
            end do
        end do

        ! Table 13 lumped coefficients
        V = 0.0_ap
        idx = 0
        pdat%A1m = 0.0_ap; pdat%B1m = 0.0_ap; 
        do m = 0,degmax
            do l = m,degmax
                idx = idx + 1
                pdat%Cfactor(idx) = pdat%rhol(l+1)*pdat%Cvec(idx)
                pdat%Sfactor(idx) = pdat%rhol(l+1)*pdat%Svec(idx)

                pdat%A1m(m+1) = pdat%A1m(m+1) + pdat%Cfactor(idx)*pdat%Hml(m+1,l+1)
                pdat%B1m(m+1) = pdat%B1m(m+1) + pdat%Sfactor(idx)*pdat%Hml(m+1,l+1)
            end do
            V = V + (pdat%A1m(m+1)*pdat%cosmlam(m+1) + pdat%B1m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
        end do

        if (dorder > 0) then
            pdat%dHml = 0.0_ap
            pdat%A2m = 0.0_ap; pdat%B2m = 0.0_ap; 
            pdat%A3m = 0.0_ap; pdat%B3m = 0.0_ap; 

            ! Evaluate Hemholtz polynomial derivatives
            do l = 1,degmax
                knm = sqrt(0.5_ap*l*(l+1.0_ap))
                pdat%dHml(1,l+1) = knm*pdat%Hml(2,l+1)
                do m = 1,l
                    knm = sqrt((l-m)*(l+m+1.0_ap))
                    pdat%dHml(m+1,l+1) = knm*pdat%Hml(m+1+1,l+1)
                end do
            end do

            ! Table 14 Auxiliary functions
            do l=0,degmax
                do m=0,l
                    pdat%Lml(m+1,l+1) = (l+m+1.0_ap)*pdat%Hml(m+1,l+1) + u*pdat%dHml(m+1,l+1)
                end do
            end do

            ! Table 13 lumped coefficients
            idx = 0
            do m = 0,degmax
                do l = m,degmax
                    idx = idx + 1
                    pdat%A2m(m+1) = pdat%A2m(m+1) + pdat%Cfactor(idx)*pdat%dHml(m+1,l+1)
                    pdat%B2m(m+1) = pdat%B2m(m+1) + pdat%Sfactor(idx)*pdat%dHml(m+1,l+1)

                    pdat%A3m(m+1) = pdat%A3m(m+1) + pdat%Cfactor(idx)*pdat%Lml(m+1,l+1)
                    pdat%B3m(m+1) = pdat%B3m(m+1) + pdat%Sfactor(idx)*pdat%Lml(m+1,l+1)
                end do
            end do

            ! Table 10 Evaluate functions ai
            a1 = 0.0_ap; a2 = 0.0_ap; a3 = 0.0_ap; a4 = 0.0_ap
            m=0
            a3 = a3 + (pdat%A2m(1)*pdat%cosmlam(1) + pdat%B2m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a4 = a4 - (pdat%A3m(1)*pdat%cosmlam(1) + pdat%B3m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            do m = 1,degmax
                a1 = a1 + m*(pdat%A1m(m+1)*pdat%cosmlam(m) + pdat%B1m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a2 = a2 + m*(pdat%B1m(m+1)*pdat%cosmlam(m) - pdat%A1m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a3 = a3 + (pdat%A2m(m+1)*pdat%cosmlam(m+1) + pdat%B2m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
                a4 = a4 - (pdat%A3m(m+1)*pdat%cosmlam(m+1) + pdat%B3m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
            end do
            a1 = a1*oneByr; a2 = a2*oneByr; a3 = a3*oneByr; a4 = a4*oneByr

            Vx = a1 + s*a4
            Vy = a2 + t*a4
            Vz = a3 + u*a4

            dV(1) = Vx; dV(2) = Vy; dV(3) = Vz;

        end if

        if (dorder > 1) then
            pdat%d2Hml = 0.0_ap
            pdat%A4m = 0.0_ap; pdat%B4m = 0.0_ap; 
            pdat%A5m = 0.0_ap; pdat%B5m = 0.0_ap; 
            pdat%A6m = 0.0_ap; pdat%B6m = 0.0_ap; 

            ! Evaluate Hemholtz polynomial derivatives
            do l = 1,degmax
                knm = sqrt(0.5_ap*l*(l+1.0_ap))
                knmp1 = sqrt(max((l-1.0_ap)*(l+2.0_ap),0.0_ap))
                pdat%d2Hml(1,l+1) = knmp1*knm*pdat%Hml(3,l+1)
                do m = 1,l
                    knm = sqrt((l-m)*(l+m+1.0_ap))
                    knmp1 = sqrt(max((l-(m+1.0_ap))*(l+(m+1.0_ap)+1.0_ap),0.0_ap))
                    pdat%d2Hml(m+1,l+1) = knmp1*knm*pdat%Hml(m+1+2,l+1)
                end do
            end do

            ! Table 14 Auxiliary functions
            do l=0,degmax
                do m=0,l
                    pdat%dLml(m+1,l+1) = (l+m+2.0_ap)*pdat%dHml(m+1,l+1) + u*pdat%d2Hml(m+1,l+1)
                    pdat%Oml(m+1,l+1) = (l+m+1.0_ap)*(l+m+2.0_ap)*pdat%Hml(m+1,l+1) &
                        + 2.0_ap*u*(l+m+2.0_ap)*pdat%dHml(m+1,l+1) + u*u*pdat%d2Hml(m+1,l+1)
                end do
            end do

            ! Table 13 lumped coefficients
            idx = 0
            do m = 0,degmax
                do l = m,degmax
                    idx = idx + 1

                    pdat%A4m(m+1) = pdat%A4m(m+1) + pdat%Cfactor(idx)*pdat%d2Hml(m+1,l+1)
                    pdat%B4m(m+1) = pdat%B4m(m+1) + pdat%Sfactor(idx)*pdat%d2Hml(m+1,l+1)

                    pdat%A5m(m+1) = pdat%A5m(m+1) + pdat%Cfactor(idx)*pdat%dLml(m+1,l+1)
                    pdat%B5m(m+1) = pdat%B5m(m+1) + pdat%Sfactor(idx)*pdat%dLml(m+1,l+1)

                    pdat%A6m(m+1) = pdat%A6m(m+1) + pdat%Cfactor(idx)*pdat%Oml(m+1,l+1)
                    pdat%B6m(m+1) = pdat%B6m(m+1) + pdat%Sfactor(idx)*pdat%Oml(m+1,l+1)
                end do
            end do


            ! Table 11 Evaluate functions aij
            a11 = 0.0_ap; a12 = 0.0_ap; a13 = 0.0_ap; a14 = 0.0_ap
            a23 = 0.0_ap; a24 = 0.0_ap
            a33 = 0.0_ap; a34 = 0.0_ap
            a44 = 0.0_ap
            ! m=0
            a33 = a33 + (pdat%A4m(1)*pdat%cosmlam(1) + pdat%B4m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a34 = a34 - (pdat%A5m(1)*pdat%cosmlam(1) + pdat%B5m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a44 = a44 + (pdat%A6m(1)*pdat%cosmlam(1) + pdat%B6m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            ! m=1
            a13 = a13 + (pdat%A2m(2)*pdat%cosmlam(1) + pdat%B2m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a14 = a14 - (pdat%A3m(2)*pdat%cosmlam(1) + pdat%B3m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a23 = a23 + (pdat%B2m(2)*pdat%cosmlam(1) - pdat%A2m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a24 = a24 - (pdat%B3m(2)*pdat%cosmlam(1) - pdat%A3m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a33 = a33 + (pdat%A4m(2)*pdat%cosmlam(2) + pdat%B4m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a34 = a34 - (pdat%A5m(2)*pdat%cosmlam(2) + pdat%B5m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a44 = a44 + (pdat%A6m(2)*pdat%cosmlam(2) + pdat%B6m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            do m = 2,degmax
                a11 = a11 + m*(m-1)*(pdat%A1m(m+1)*pdat%cosmlam(m-1) + pdat%B1m(m+1)*pdat%sinmlam(m-1))*pdat%cosphipow(m-1)
                a12 = a12 + m*(m-1)*(pdat%B1m(m+1)*pdat%cosmlam(m-1) - pdat%A1m(m+1)*pdat%sinmlam(m-1))*pdat%cosphipow(m-1)
                a13 = a13 + m*(pdat%A2m(m+1)*pdat%cosmlam(m) + pdat%B2m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a14 = a14 - m*(pdat%A3m(m+1)*pdat%cosmlam(m) + pdat%B3m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a23 = a23 + m*(pdat%B2m(m+1)*pdat%cosmlam(m) - pdat%A2m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a24 = a24 - m*(pdat%B3m(m+1)*pdat%cosmlam(m) - pdat%A3m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a33 = a33 + (pdat%A4m(m+1)*pdat%cosmlam(m+1) + pdat%B4m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
                a34 = a34 - (pdat%A5m(m+1)*pdat%cosmlam(m+1) + pdat%B5m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
                a44 = a44 + (pdat%A6m(m+1)*pdat%cosmlam(m+1) + pdat%B6m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
            end do
            a11 = oneByr2*a11; a12 = oneByr2*a12; a13 = oneByr2*a13; a14 = oneByr2*a14
            a23 = oneByr2*a23; a24 = oneByr2*a24
            a33 = oneByr2*a33; a34 = oneByr2*a34
            a44 = oneByr2*a44
            a22 = -a11

            Vxx = a11 + 2.0_ap*s*a14 + a4*oneByr + s*s*a44 - s*s*a4*oneByr
            Vxy = a12 + s*t*a44 + s*a24 + t*a14 - s*t*a4*oneByr
            Vyy = a22 + 2.0_ap*t*a24 + a4*oneByr + t*t*a44 - t*t*a4*oneByr
            Vxz = a13 + s*u*a44 + s*a34 + u*a14 - s*u*a4*oneByr
            Vyz = a23 + t*u*a44 + t*a34 + u*a24 - t*u*a4*oneByr
            Vzz = a33 + 2.0_ap*u*a34 + a4*oneByr + u*u*a44 - u*u*a4*oneByr

            d2V(1,1) = Vxx; d2V(2,2) = Vyy; d2V(3,3) = Vzz
            d2V(1,2) = Vxy; d2V(2,1) = Vxy;
            d2V(1,3) = Vxz; d2V(3,1) = Vxz;
            d2V(2,3) = Vyz; d2V(3,2) = Vyz;

        end if

        if (dorder > 2) then
            pdat%d3Hml = 0.0_ap
            pdat%A7m = 0.0_ap; pdat%B7m = 0.0_ap; 
            pdat%A8m = 0.0_ap; pdat%B8m = 0.0_ap; 
            pdat%A9m = 0.0_ap; pdat%B9m = 0.0_ap; 
            pdat%A10m = 0.0_ap; pdat%B10m = 0.0_ap; 

            ! Evaluate Hemholtz polynomial derivatives
            do l = 1,degmax
                knm = sqrt(0.5_ap*l*(l+1.0_ap))
                knmp1 = sqrt(max((l-1.0_ap)*(l+2.0_ap),0.0_ap))
                knmp2 = sqrt(max((l-2.0_ap)*(l+3.0_ap),0.0_ap))
                pdat%d3Hml(1,l+1) = knmp2*knmp1*knm*pdat%Hml(4,l+1)
                do m = 1,l
                    knm = sqrt((l-m)*(l+m+1.0_ap))
                    knmp1 = sqrt(max((l-(m+1.0_ap))*(l+(m+1.0_ap)+1.0_ap),0.0_ap))
                    knmp2 = sqrt(max((l-(m+2.0_ap))*(l+(m+2.0_ap)+1.0_ap),0.0_ap))
                    pdat%d3Hml(m+1,l+1) = knmp2*knmp1*knm*pdat%Hml(m+1+3,l+1)
                end do
            end do

            ! Table 14 Auxiliary functions
            do l=0,degmax
                do m=0,l
                    pdat%d2Lml(m+1,l+1) = (l+m+3.0_ap)*pdat%d2Hml(m+1,l+1) + u*pdat%d3Hml(m+1,l+1)
                    pdat%Qml(m+1,l+1) = (l+m+2.0_ap)*(l+m+3.0_ap)*pdat%dHml(m+1,l+1) &
                        + 2.0_ap*u*(l+m+3.0_ap)*pdat%d2Hml(m+1,l+1) + u*u*pdat%d3Hml(m+1,l+1)
                    pdat%Wml(m+1,l+1) = (l+m+1.0_ap)*(l+m+2.0_ap)*(l+m+3.0_ap)*pdat%Hml(m+1,l+1) &
                        + 3.0_ap*u*(l+m+2.0_ap)*(l+m+3.0_ap)*pdat%dHml(m+1,l+1) & 
                        + 3.0_ap*(l+m+3.0_ap)*u*u*pdat%d2Hml(m+1,l+1) & 
                        + u*u*u*pdat%d3Hml(m+1,l+1)
                end do
            end do

            ! Table 13 lumped coefficients
            idx = 0
            do m = 0,degmax
                do l = m,degmax
                    idx = idx + 1

                    pdat%A7m(m+1) = pdat%A7m(m+1) + pdat%Cfactor(idx)*pdat%d3Hml(m+1,l+1)
                    pdat%B7m(m+1) = pdat%B7m(m+1) + pdat%Sfactor(idx)*pdat%d3Hml(m+1,l+1)

                    pdat%A8m(m+1) = pdat%A8m(m+1) + pdat%Cfactor(idx)*pdat%d2Lml(m+1,l+1)
                    pdat%B8m(m+1) = pdat%B8m(m+1) + pdat%Sfactor(idx)*pdat%d2Lml(m+1,l+1)

                    pdat%A9m(m+1) = pdat%A9m(m+1) + pdat%Cfactor(idx)*pdat%Qml(m+1,l+1)
                    pdat%B9m(m+1) = pdat%B9m(m+1) + pdat%Sfactor(idx)*pdat%Qml(m+1,l+1)

                    pdat%A10m(m+1) = pdat%A10m(m+1) + pdat%Cfactor(idx)*pdat%Wml(m+1,l+1)
                    pdat%B10m(m+1) = pdat%B10m(m+1) + pdat%Sfactor(idx)*pdat%Wml(m+1,l+1)

                end do
            end do

            ! Table 12 Evaluate functions aijk
            a111 = 0.0_ap; a113 = 0.0_ap; a114 = 0.0_ap
            a123 = 0.0_ap; a124 = 0.0_ap; a143 = 0.0_ap
            a222 = 0.0_ap; a243 = 0.0_ap
            a331 = 0.0_ap; a332 = 0.0_ap; a333 = 0.0_ap; a334 = 0.0_ap
            a441 = 0.0_ap; a442 = 0.0_ap; a443 = 0.0_ap; a444 = 0.0_ap
            ! m=0
            a333 = a333 + (pdat%A7m(1)*pdat%cosmlam(1) + pdat%B7m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a334 = a334 - (pdat%A8m(1)*pdat%cosmlam(1) + pdat%B8m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a443 = a443 + (pdat%A9m(1)*pdat%cosmlam(1) + pdat%B9m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a444 = a444 - (pdat%A10m(1)*pdat%cosmlam(1) + pdat%B10m(1)*pdat%sinmlam(1))*pdat%cosphipow(1)
            ! m=1
            a143 = a143 - (pdat%A5m(2)*pdat%cosmlam(1) + pdat%B5m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a243 = a243 - (pdat%B5m(2)*pdat%cosmlam(1) - pdat%A5m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a331 = a331 + (pdat%A4m(2)*pdat%cosmlam(1) + pdat%B4m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a332 = a332 + (pdat%B4m(2)*pdat%cosmlam(1) - pdat%A4m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a441 = a441 + (pdat%A6m(2)*pdat%cosmlam(1) + pdat%B6m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a442 = a442 + (pdat%B6m(2)*pdat%cosmlam(1) - pdat%A6m(2)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a444 = a444 - (pdat%A10m(2)*pdat%cosmlam(2) + pdat%B10m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a333 = a333 + (pdat%A7m(2)*pdat%cosmlam(2) + pdat%B7m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a334 = a334 - (pdat%A8m(2)*pdat%cosmlam(2) + pdat%B8m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a443 = a443 + (pdat%A9m(2)*pdat%cosmlam(2) + pdat%B9m(2)*pdat%sinmlam(2))*pdat%cosphipow(2)
            ! m=2
            a113 = a113 + 2.0_ap*(pdat%A2m(3)*pdat%cosmlam(1) + pdat%B2m(3)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a114 = a114 - 2.0_ap*(pdat%A3m(3)*pdat%cosmlam(1) + pdat%B3m(3)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a123 = a123 + 2.0_ap*(pdat%B2m(3)*pdat%cosmlam(1) - pdat%A2m(3)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a124 = a124 - 2.0_ap*(pdat%B3m(3)*pdat%cosmlam(1) - pdat%A3m(3)*pdat%sinmlam(1))*pdat%cosphipow(1)
            a143 = a143 - 2.0_ap*(pdat%A5m(3)*pdat%cosmlam(2) + pdat%B5m(3)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a243 = a243 - 2.0_ap*(pdat%B5m(3)*pdat%cosmlam(2) - pdat%A5m(3)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a331 = a331 + 2.0_ap*(pdat%A4m(3)*pdat%cosmlam(2) + pdat%B4m(3)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a332 = a332 + 2.0_ap*(pdat%B4m(3)*pdat%cosmlam(2) - pdat%A4m(3)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a441 = a441 + 2.0_ap*(pdat%A6m(3)*pdat%cosmlam(2) + pdat%B6m(3)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a442 = a442 + 2.0_ap*(pdat%B6m(3)*pdat%cosmlam(2) - pdat%A6m(3)*pdat%sinmlam(2))*pdat%cosphipow(2)
            a444 = a444 - (pdat%A10m(3)*pdat%cosmlam(3) + pdat%B10m(3)*pdat%sinmlam(3))*pdat%cosphipow(3)
            a333 = a333 + (pdat%A7m(3)*pdat%cosmlam(3) + pdat%B7m(3)*pdat%sinmlam(3))*pdat%cosphipow(3)
            a334 = a334 - (pdat%A8m(3)*pdat%cosmlam(3) + pdat%B8m(3)*pdat%sinmlam(3))*pdat%cosphipow(3)
            a443 = a443 + (pdat%A9m(3)*pdat%cosmlam(3) + pdat%B9m(3)*pdat%sinmlam(3))*pdat%cosphipow(3)
            do m = 3,degmax
                a111 = a111 + m*(m-1.0_ap)*(m-2.0_ap)*(pdat%A1m(m+1)*pdat%cosmlam(m-2) + pdat%B1m(m+1)*pdat%sinmlam(m-2))*pdat%cosphipow(m-2)
                a222 = a222 - m*(m-1.0_ap)*(m-2.0_ap)*(pdat%B1m(m+1)*pdat%cosmlam(m-2) - pdat%A1m(m+1)*pdat%sinmlam(m-2))*pdat%cosphipow(m-2)
                a113 = a113 + m*(m-1.0_ap)*(pdat%A2m(m+1)*pdat%cosmlam(m-1) + pdat%B2m(m+1)*pdat%sinmlam(m-1))*pdat%cosphipow(m-1)
                a114 = a114 - m*(m-1.0_ap)*(pdat%A3m(m+1)*pdat%cosmlam(m-1) + pdat%B3m(m+1)*pdat%sinmlam(m-1))*pdat%cosphipow(m-1)
                a123 = a123 + m*(m-1.0_ap)*(pdat%B2m(m+1)*pdat%cosmlam(m-1) - pdat%A2m(m+1)*pdat%sinmlam(m-1))*pdat%cosphipow(m-1)
                a124 = a124 - m*(m-1.0_ap)*(pdat%B3m(m+1)*pdat%cosmlam(m-1) - pdat%A3m(m+1)*pdat%sinmlam(m-1))*pdat%cosphipow(m-1)
                a143 = a143 - m*(pdat%A5m(m+1)*pdat%cosmlam(m) + pdat%B5m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a243 = a243 - m*(pdat%B5m(m+1)*pdat%cosmlam(m) - pdat%A5m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a331 = a331 + m*(pdat%A4m(m+1)*pdat%cosmlam(m) + pdat%B4m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a332 = a332 + m*(pdat%B4m(m+1)*pdat%cosmlam(m) - pdat%A4m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a441 = a441 + m*(pdat%A6m(m+1)*pdat%cosmlam(m) + pdat%B6m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a442 = a442 + m*(pdat%B6m(m+1)*pdat%cosmlam(m) - pdat%A6m(m+1)*pdat%sinmlam(m))*pdat%cosphipow(m)
                a444 = a444 - (pdat%A10m(m+1)*pdat%cosmlam(m+1) + pdat%B10m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
                a333 = a333 + (pdat%A7m(m+1)*pdat%cosmlam(m+1) + pdat%B7m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
                a334 = a334 - (pdat%A8m(m+1)*pdat%cosmlam(m+1) + pdat%B8m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
                a443 = a443 + (pdat%A9m(m+1)*pdat%cosmlam(m+1) + pdat%B9m(m+1)*pdat%sinmlam(m+1))*pdat%cosphipow(m+1)
            end do
            a111 = oneByr3*a111; a113 = oneByr3*a113; a114 = oneByr3*a114
            a123 = oneByr3*a123; a124 = oneByr3*a124; a143 = oneByr3*a143
            a222 = oneByr3*a222; a243 = oneByr3*a243
            a331 = oneByr3*a331; a332 = oneByr3*a332; a333 = oneByr3*a333; a334 = oneByr3*a334
            a441 = oneByr3*a441; a442 = oneByr3*a442; a443 = oneByr3*a443; a444 = oneByr3*a444
            a112 = -a222; a221 = -a111; a223 = -a113; a224 = -a114

            Vxxx = a111 + 3.0_ap*s*a114 + 3*s*s*a441 + s*s*s*a444 &
                + 3.0_ap*(1.0_ap - s*s)*(a14 + s*a44 - s*a4*oneByr)*oneByr
            Vxxy = a112 + 2.0_ap*s*a124 + s*s*a442 + t*a114 + 2.0_ap*t*s*a441 &
                + (1.0_ap - s*s)*a24*oneByr + t*(1.0_ap - 3.0_ap*s*s)*(a44 - a4*oneByr)*oneByr &
                + t*s*s*a444 - 2.0_ap*t*s*a14*oneByr
            Vxxz = a113 + 2.0_ap*s*a143 + s*s*a443 + u*a114 + 2.0_ap*u*s*a441 &
                + (1.0_ap - s*s)*a34*oneByr + u*(1.0_ap - 3.0_ap*s*s)*(a44 - a4*oneByr)*oneByr &
                + u*s*s*a444 - 2.0_ap*u*s*a14*oneByr
            Vyyx = a221 + 2.0_ap*t*a124 + t*t*a441 + s*a224 + 2.0_ap*t*s*a442 &
                + (1.0_ap - t*t)*a14*oneByr + s*(1.0_ap - 3.0_ap*t*t)*(a44 - a4*oneByr)*oneByr &
                + t*t*s*a444 - 2.0_ap*t*s*a24*oneByr
            Vyyy = a222 + 3.0_ap*t*a224 + 3.0_ap*t*t*a442 + t*t*t*a444 &
                + 3.0_ap*(1.0_ap - t*t)*(a24 + t*a44 - t*a4*oneByr)*oneByr
            Vyyz = a223 + 2.0_ap*t*a243 + t*t*a443 + u*a224 + 2.0_ap*t*u*a442 &
                + (1.0_ap - t*t)*a34*oneByr + u*(1.0_ap - 3.0_ap*t*t)*(a44 - a4*oneByr)*oneByr &
                + t*t*u*a444 - 2*t*u*a24*oneByr
            Vzzx = a331 + 2.0_ap*u*a143 + u*u*a441 + s*a334 + 2.0_ap*u*s*a443 &
                + (1.0_ap - u*u)*a14*oneByr + s*(1.0_ap - 3.0_ap*u*u)*(a44 - a4*oneByr)*oneByr &
                + u*u*s*a444 - 2.0_ap*u*s*a34*oneByr
            Vzzy = a332 + 2.0_ap*u*a243 + u*u*a442 + t*a334 + 2.0_ap*u*t*a443 &
                + (1.0_ap - u*u)*a24*oneByr + t*(1.0_ap - 3.0_ap*u*u)*(a44 - a4*oneByr)*oneByr &
                + u*u*t*a444 - 2.0_ap*u*t*a34*oneByr
            Vzzz = a333 + 3.0_ap*u*a334 + 3.0_ap*u*u*a443 + u*u*u*a444 &
                + 3.0_ap*(1.0_ap - u*u)*(a34 + u*a44 - u*a4*oneByr)*oneByr
            Vxyz = a123 + u*a124 + s*a243 + u*s*a442 + t*a143 + s*t*a443 + u*s*t*a444 &
                - s*t*a34*oneByr - 3.0_ap*u*s*t*a44*oneByr - u*s*a24*oneByr - u*t*a14*oneByr &
                + 3.0_ap*u*s*t*a4*oneByr2 + u*t*a441

            d3V(1,1,1) = Vxxx; d3V(2,2,2) = Vyyy; d3V(3,3,3) = Vzzz
            d3V(1,1,2) = Vxxy; d3V(2,1,1) = Vxxy; d3V(1,2,1) = Vxxy
            d3V(1,1,3) = Vxxz; d3V(3,1,1) = Vxxz; d3V(1,3,1) = Vxxz
            d3V(2,2,1) = Vyyx; d3V(1,2,2) = Vyyx; d3V(2,1,2) = Vyyx
            d3V(2,2,3) = Vyyz; d3V(3,2,2) = Vyyz; d3V(2,3,2) = Vyyz
            d3V(3,3,1) = Vzzx; d3V(1,3,3) = Vzzx; d3V(3,1,3) = Vzzx
            d3V(3,3,2) = Vzzy; d3V(2,3,3) = Vzzy; d3V(3,2,3) = Vzzy
            d3V(1,2,3) = Vxyz; d3V(3,1,2) = Vxyz; d3V(2,3,1) = Vxyz
            d3V(1,3,2) = Vxyz; d3V(2,1,3) = Vxyz; d3V(3,2,1) = Vxyz

        end if

    end subroutine shpines
end module tinysh
