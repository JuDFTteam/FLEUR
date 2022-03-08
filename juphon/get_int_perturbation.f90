!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------

MODULE m_get_int_perturbation
    USE m_juDFT
    USE m_fft3d
    USE m_types

    IMPLICIT NONE

CONTAINS

    SUBROUTINE get_int_local_perturbation(sym, stars, atoms, sphhar, &
                                        & input, den, den1, den1im, starsq)

        USE m_constants

        TYPE(t_input),  INTENT(IN)    :: input
        TYPE(t_sym),    INTENT(IN)    :: sym
        TYPE(t_stars),  INTENT(IN)    :: stars, starsq
        TYPE(t_sphhar), INTENT(IN)    :: sphhar
        TYPE(t_atoms),  INTENT(IN)    :: atoms
        TYPE(t_potden), INTENT(IN)    :: den
        TYPE(t_potden), INTENT(INOUT) :: den1, den1im

        INTEGER                       :: iden, jspin, ifft3
        INTEGER                       :: ityp, iri, imesh
        REAL                          :: rho_11, rho_22, rho_21r, rho_21i
        REAL                          :: mx, my, mz, m
        REAL                          :: rhotot, rho_up, rho_down, theta, phi
        COMPLEX                       :: m1, mx1, my1, mz1, n1, t1, p1, rho1_up, rho1_down

        REAL, ALLOCATABLE             :: ris(:,:), ris_real(:,:), ris_imag(:,:)
        REAL, ALLOCATABLE             :: fftwork(:)

        ifft3 = 27*stars%mx1*stars%mx2*stars%mx3

        !TODO: Make sure the indices for rho1 are 1,2,3,4 == n1,mx1,my1,mz1
        ALLOCATE (ris(ifft3,4),fftwork(ifft3))
        ALLOCATE (ris_real(ifft3,4),ris_imag(ifft3,4))

        DO iden = 1, 2
            CALL fft3d(ris(:,iden),      fftwork,          den%pw(:,iden),  stars,  +1)
            CALL fft3d(ris_real(:,iden), ris_imag(:,iden), den1%pw(:,iden), starsq, +1)
        END DO

        CALL fft3d(ris(:,3),      ris(:,4),      den%pw(:,3),  stars,  +1)
        CALL fft3d(ris_real(:,3), ris_imag(:,3), den1%pw(:,3), starsq, +1)
        CALL fft3d(ris_real(:,4), ris_imag(:,4), den1%pw(:,4), starsq, +1)

        DO imesh = 1, ifft3
            ! Get real space density matrix elements
            rho_11   = ris(imesh,1)
            rho_22   = ris(imesh,2)
            rho_21r  = ris(imesh,3)
            rho_21i  = ris(imesh,4)

            ! Calculate unperturbed magnetization density
            mx       =  2*rho_21r
            my       =  2*rho_21i !TODO: PLEASE get on the same page with this sign!
            mz       = rho_11 - rho_22
            m  = SQRT(mx**2 + my**2 + mz**2)

            ! Calculate perturbed total and magnetization density
            n1  = ris_real(imesh,1) + Imagunit * ris_imag(imesh,1)
            mx1 = ris_real(imesh,2) + Imagunit * ris_imag(imesh,2)
            my1 = ris_real(imesh,3) + Imagunit * ris_imag(imesh,3)
            mz1 = ris_real(imesh,4) + Imagunit * ris_imag(imesh,4)

            theta = den%theta_pw(imesh)
            phi = den%phi_pw(imesh)

            ! Calculate the perturbed absolute value of the magnetization density
            m1 = cmplx(0.0, 0.0)
            m1 = m1 + mx1 * SIN(theta) * COS(phi)
            m1 = m1 + my1 * SIN(theta) * SIN(phi)
            m1 = m1 + mz1 * COS(theta)

            ! Calculate the perturbed angles
            t1 = cmplx(0.0, 0.0)
            t1 = t1 + mx1 * COS(theta) * COS(phi) / m
            t1 = t1 + my1 * COS(theta) * SIN(phi) / m
            t1 = t1 - mz1 * SIN(theta)            / m

            p1 = cmplx(0.0, 0.0)
            p1 = p1 - mx1 * SIN(theta) * SIN(phi) / m
            p1 = p1 + my1 * SIN(theta) * COS(phi) / m


            rho1_up   = (n1 + m1)/2
            rho1_down = (n1 - m1)/2

            ris_real(imesh,1) =  REAL(rho1_up  )
            ris_real(imesh,2) =  REAL(rho1_down)
            ris_imag(imesh,1) = AIMAG(rho1_up  )
            ris_imag(imesh,2) = AIMAG(rho1_down)
            den1%theta_pw(imesh)   =  REAL(t1)
            den1%phi_pw(imesh)     =  REAL(p1)
            den1im%theta_pw(imesh) = AIMAG(t1)
            den1im%phi_pw(imesh)   = AIMAG(p1)
        END DO

        ! Fourier tranform the up- and down density perturbations back to reciprocal space:
        DO jspin = 1, input%jspins
            CALL fft3d(ris_real(:,jspin),ris_imag(:,jspin),den1%pw(:,jspin),starsq,-1)
        END DO
    END SUBROUTINE get_int_local_perturbation

    SUBROUTINE get_int_global_perturbation(stars,atoms,sym,input,den,den1,den1im,vTot,vTot1,starsq)
        TYPE(t_input), INTENT(IN)     :: input
        TYPE(t_sym),    INTENT(IN)    :: sym
        TYPE(t_stars),  INTENT(IN)    :: stars, starsq
        TYPE(t_atoms),  INTENT(IN)    :: atoms
        TYPE(t_potden), INTENT(IN)    :: den, den1, den1im, vTot
        TYPE(t_potden), INTENT(INOUT) :: vTot1

        INTEGER                       :: imeshpt, ipot, jspin
        INTEGER                       :: ifft3, i
        REAL                          :: theta, phi, v11, v22

        REAL, ALLOCATABLE             :: vis(:,:), vis_re(:,:), vis_im(:,:), vis2_re(:,:), vis2_im(:,:), v1re(:,:), v1im(:,:)
        REAL, ALLOCATABLE             :: fftwork(:)

        COMPLEX :: a11, a22, a21, a12, av11, av22, av21, av12
        COMPLEX :: t1, p1, v21, v12, v1up, v1down, v1, b1, v1mat11, v1mat22, v1mat21, v1mat12

        ifft3 = 27*stars%mx1*stars%mx2*stars%mx3
        IF (ifft3.NE.SIZE(den%theta_pw)) CALL judft_error("Wrong size of angles")

        ALLOCATE (vis(ifft3,4),vis_re(ifft3,4),vis_im(ifft3,4),fftwork(ifft3),vis2_re(ifft3,4),vis2_im(ifft3,4),v1re(ifft3,4),v1im(ifft3,4))

        DO jspin = 1, input%jspins
            CALL fft3d(vis(:, jspin), fftwork, vTot%pw(:, jspin), stars, +1)
            CALL fft3d(v1re(:, jspin), v1im(:, jspin), vTot1%pw(:, jspin), starsq, +1)
        END DO

        CALL fft3d(vis(:, 3), vis(:, 4), vTot%pw(:, 3), stars, +1)

        DO imeshpt = 1, ifft3
            v11   = vis(imeshpt, 1)
            v22   = vis(imeshpt, 2)
            v21   = vis(imeshpt, 3) + ImagUnit * vis(imeshpt, 4)
            v12   = vis(imeshpt, 3) - ImagUnit * vis(imeshpt, 4)

            theta = den%theta_pw(imeshpt)
            phi   = den%phi_pw(imeshpt)

            v1up   = v1re(imeshpt,1) + ImagUnit * v1im(imeshpt,1)
            v1down = v1re(imeshpt,2) + ImagUnit * v1im(imeshpt,2)

            t1 = den1%theta_pw(imeshpt) + ImagUnit * den1im%theta_pw(imeshpt)
            p1 = den1%phi_pw(imeshpt) + ImagUnit * den1im%phi_pw(imeshpt)

            v1 = (v1up + v1down) / 2.0
            b1 = (v1up - v1down) / 2.0

            a11 =      -ImagUnit      * p1 / 2.0
            a22 =       ImagUnit      * p1 / 2.0
            a21 =  EXP( ImagUnit*phi) * t1 / 2.0
            a12 = -EXP(-ImagUnit*phi) * t1 / 2.0

            v1mat11 = v1 + b1 * COS(theta)                    !11
            v1mat22 = v1 - b1 * COS(theta)                    !22
            v1mat21 =      b1 * SIN(theta)*EXP( Imagunit*phi) !21
            v1mat12 =      b1 * SIN(theta)*EXP(-Imagunit*phi) !12

            av11 = a11 * v11 + a12 * v21 !11
            av22 = a21 * v12 + a22 * v22 !22
            av21 = a21 * v11 + a22 * v21 !21
            av12 = a11 * v12 + a12 * v22 !12

            v1mat11 = v1mat11 + av11 + CONJG(av11)
            v1mat22 = v1mat22 + av22 + CONJG(av22)
            v1mat21 = v1mat21 + av21 + CONJG(av12)
            v1mat12 = v1mat12 + av12 + CONJG(av21)

            vis_re(imeshpt, 1) =  REAL(v1mat11)
            vis_re(imeshpt, 2) =  REAL(v1mat22)
            vis_re(imeshpt, 3) =  REAL(v1mat21)
            vis_re(imeshpt, 4) =  REAL(v1mat12)

            vis2_re(imeshpt, 1) =  REAL(v1mat11 * stars%ufft(imeshpt-1) + v11 * starsq%ufft(imeshpt-1))
            vis2_re(imeshpt, 2) =  REAL(v1mat22 * stars%ufft(imeshpt-1) + v22 * starsq%ufft(imeshpt-1))
            vis2_re(imeshpt, 3) =  REAL(v1mat21 * stars%ufft(imeshpt-1) + v21 * starsq%ufft(imeshpt-1))
            vis2_re(imeshpt, 4) =  REAL(v1mat12 * stars%ufft(imeshpt-1) + v12 * starsq%ufft(imeshpt-1))

            vis_im(imeshpt, 1) = AIMAG(v1mat11)
            vis_im(imeshpt, 2) = AIMAG(v1mat22)
            vis_im(imeshpt, 3) = AIMAG(v1mat21)
            vis_im(imeshpt, 4) = AIMAG(v1mat12)

            vis2_im(imeshpt, 1) = AIMAG(v1mat11 * stars%ufft(imeshpt-1) + v11 * starsq%ufft(imeshpt-1))
            vis2_im(imeshpt, 2) = AIMAG(v1mat22 * stars%ufft(imeshpt-1) + v22 * starsq%ufft(imeshpt-1))
            vis2_im(imeshpt, 3) = AIMAG(v1mat21 * stars%ufft(imeshpt-1) + v21 * starsq%ufft(imeshpt-1))
            vis2_im(imeshpt, 4) = AIMAG(v1mat12 * stars%ufft(imeshpt-1) + v12 * starsq%ufft(imeshpt-1))

        END DO

        DO ipot = 1, 4
            CALL fft3d(vis_re(:, ipot),  vis_im(:, ipot),  vTot1%pw(1, ipot),   starsq, -1)
            CALL fft3d(vis2_re(:, ipot), vis2_im(:, ipot), vTot1%pw_w(1, ipot), starsq, -1)
        END DO

    END SUBROUTINE get_int_global_perturbation
END MODULE m_get_int_perturbation
