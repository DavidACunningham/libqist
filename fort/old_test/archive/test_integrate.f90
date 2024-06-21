program test_integrate
    use genqist
    use, intrinsic :: iso_fortran_env, only: dp => real64, wp => real128
    implicit none
     ! TODO: SET UP VARIABLES
     real(dp), parameter :: pi_d=4._dp*atan(1._dp)
     real(dp), parameter :: mu=398600.435507_dp, &
                            LU=6378.137_dp, TU=sqrt(LU**3/mu)
     real(wp), parameter :: pi_q=4._wp*atan(1._wp)
     integer,  parameter :: points_per_orbit_integ=300, &
                          & num_orbits_integ=200, &
                          & tlen=points_per_orbit_integ*num_orbits_integ
     real(dp)            :: t_nondimensional(tlen), &
                            t_dimensional(tlen), &
                            UV(tlen), &
                            x0_pqw(6),x0_ijk(6), Yf_temp(7), dYf_temp(7,7), &
                            rotmat(6,6),&
                            d2Yf_temp(7,7,7), &
                            x(tlen,6), Yf(tlen,7), dYf(tlen,7,7), &
                            d2Yf(tlen,7,7,7), capz, &
                            elvec_dim(6), elvec_nd(6)
     integer             :: errFlag, i


     elvec_dim = [8000._dp, 0.05_dp, pi_d/3._dp, 0._dp, 0._dp, 0._dp]
     elvec_nd = elvec_dim
     elvec_nd(1) = elvec_dim(1)/LU
     x0_pqw = kep2cart(elvec_nd)
     rotmat = ijkfrompqw_rv(x0_pqw)
     x0_ijk = matmul(rotmat,x0_pqw)
     UV = [(i*2*pi_d/tlen, i=1,tlen)]


     ! TODO: COMPUTE TRUTH ORBIT USING PartialsKepDeltaE
     ! subroutine PartialsKepDeltaE(
     !                        x0,x,Yf,dYf,d2Yf,capZ,orderCase,errFlag)

     ! TODO: SAVE TRUTH REFERENCE ORBIT IN SPICE KERNEL

     ! TODO: INTEGRATE QIST MODEL CALLING SPICE KERNEL

     ! TODO: COMPUTE ERROR BETWEEN MODELS (SHOULD BE ZERO)

     contains
     function ijkfrompqw_rv(x) result(res)
        ! x is assumed to be a tuple of the form
        ! x = [r.T,v.T].T
        ! returns a matrix so the reverse transformation is just the transpose
        ! output gives a 6x6 matrix to operate on a [r,v] state vector.
        ! Intended use: v_ijk = matmul(ijkfrompqw,v_pqw) for some vector v
        ! The rotation matrix giving the coordinates of the new basis
        ! in terms of the old basis
        real(dp), intent(in) :: x(6)
        real(dp)             :: res(6,6)
        real(dp)             :: rvec(3), vvec(3), r, v, h(3), &
                              & what(3), evec(3), phat(3), qhat(3), rotmat(3,3)
        rvec = x(:3)
        vvec = x(4:)
        r = norm2(rvec)
        v = norm2(vvec)
        h = cross(rvec,vvec)
        what = h/norm2(h)
        evec = cross(vvec,h)/mu -rvec/r
        phat = evec/norm2(evec)
        qhat = cross(what,phat)
        rotmat(:,1) = phat
        rotmat(:,2) = qhat
        rotmat(:,3) = what
        res = 0._dp
        res(:3,:3) = rotmat
        res(4:,4:) = rotmat
     end function
     function kep2cart(x) result(res)
        ! x is assumed to be a tuple of the form
        ! x = (a, e, i, RAAN, argp, nu)
        ! outputs in the PQW frame by default
        real(dp), intent(in) :: x(6)
        real(dp)             :: res(6)
        real(dp)             :: p, r, vr, vt, rpqw(3), vpqw(3)
        p = x(1)*(1-x(2)**2)
        r = p/(1+x(2)*cos(x(6)))
        vr = sqrt(mu/p)*x(2)*sin(x(6))
        vt = sqrt(mu/p)*(1+x(2)*cos(x(6)))
        rpqw = matmul(rot3(-x(6)),[r,0._dp,0._dp])
        vpqw = matmul(rot3(-x(6)),[vr,vt,0._dp])
        res = [rpqw,vpqw]
     end function
     function rot1(ang) result(res)
         real(dp), intent(in) :: ang
         real(dp)             :: res(3,3)
         res(1,:) = [ 1._dp, 0._dp, 0._dp]
         res(2,:) = [0._dp, cos(ang), sin(ang)]
         res(3,:) = [0._dp, -sin(ang), cos(ang)]
     end function
     function rot2(ang)  result(res)
         real(dp), intent(in) :: ang
         real(dp)             :: res(3,3)
         res(1,:) = [cos(ang), 0._dp, -sin(ang)]
         res(2,:) = [0._dp,1._dp,0._dp]
         res(3,:) = [sin(ang), 0._dp, cos(ang)]

     end function
     function rot3(ang) result(res)
         real(dp), intent(in) :: ang
         real(dp)             :: res(3,3)
         res(1,:) =[cos(ang), sin(ang), 0._dp]
         res(2,:) =[-sin(ang),cos(ang), 0._dp]
         res(3,:) =[0._dp, 0._dp, 1._dp]
     end function
    function cross(a, b) result(res)
        real(dp), intent(in) :: a(3), b(3)
        real(dp)             :: res(3)
        res = [a(2)*b(3) - a(3)*b(2), &
               a(3)*b(1) - a(1)*b(3), &
               a(1)*b(2) - a(2)*b(1)]
    end function
end program
