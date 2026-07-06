!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module m_gf_noco
    use m_juDFT
    implicit none
    contains
    SUBROUTINE gf_noco_rotate(stars,qpw,vpw)
        use m_gf_types
        use m_gf_fft
        implicit none
        type(t_stars),intent(in)::stars
        complex,intent(INOUT)      :: qpw(:,:)
        complex,intent(INOUT)   :: vpw(:,:)

        integer :: jspin,n
        complex,allocatable :: v(:,:),q(:,:)

        allocate(v(27*stars%mx1*stars%mx2*stars%mx3,3))
        allocate(q(27*stars%mx1*stars%mx2*stars%mx3,3))

        !Put data on real-space
        DO jspin = 1, 3
            call gf_fft3d(q(:,jspin),qpw(:,jspin),stars,GF_FFT_TO_REAL_SPACE)
            if (jspin<3) &
                 call gf_fft3d(v(:,jspin),vpw(:,jspin),stars,GF_FFT_TO_REAL_SPACE)
        END DO

        !Loop over all real-space points
        DO n = 1, size(v,1)
            call priv_rotate(q(n,:),v(n,:))
        END DO

        !Put potential back in g-space
        DO jspin = 1, 3
            call gf_fft3d(v(:,jspin),vpw(:,jspin),stars,GF_FFT_TO_G_SPACE)
        END DO

        DEALLOCATE(q,v)

    END SUBROUTINE gf_noco_rotate

    SUBROUTINE priv_rotate(q,v)
        use m_constants
        implicit none
        complex,intent(in)    :: q(3)
        complex,intent(inout) :: v(3)


        real:: mx,my,mz,theta,phi
        real:: vup,vdown,veff,beff

        real,parameter:: eps= 1.0e-20
        real          :: pi

        pi=pimach()

         mx      =  2*real(q(3))
         my      = -2*aimag(q(3))
         mz      = real(q(1)-q(2))

         IF (abs(mz) .LE. eps) THEN
            theta = pi/2
         ELSEIF (mz .GE. 0.0) THEN
            theta = atan(sqrt(mx**2 + my**2)/mz)
         ELSE
            theta = atan(sqrt(mx**2 + my**2)/mz) + pi
         ENDIF

         IF (abs(mx) .LE. eps) THEN
            IF (abs(my) .LE. eps) THEN
               phi = 0.0
            ELSEIF (my .GE. 0.0) THEN
               phi = pi/2
            ELSE
               phi = -pi/2
            ENDIF
         ELSEIF (mx .GE. 0.0) THEN
            phi = atan(my/mx)
         ELSE
            IF (my .GE. 0.0) THEN
               phi = atan(my/mx) + pi
            ELSE
               phi = atan(my/mx) - pi
            ENDIF
         ENDIF

         vup   = real(v(1))
         vdown = real(v(2))
         veff  = (vup + vdown)/2.0
         beff  = (vup - vdown)/2.0

         v(1) = veff + beff*cos(theta)
        !  V_22
         v(2) = veff - beff*cos(theta)
        !  the real,imag part of V_21
         v(3) = cmplx(beff*sin(theta)*cos(phi),beff*sin(theta)*sin(phi))

    END SUBROUTINE priv_rotate
end module m_gf_noco
