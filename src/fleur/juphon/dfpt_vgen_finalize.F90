!--------------------------------------------------------------------------------
! Copyright (c) 2022 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_dfpt_vgen_finalize
   USE m_juDFT
   USE m_xcBfield
   USE m_plot
   USE m_constants
   USE m_lattHarmsSphHarmsConv

CONTAINS

   SUBROUTINE dfpt_vgen_finalize(fmpi,atoms,stars,sym,juphon,noco,nococonv,input,sphhar,vTot,vTot1,vTot1imag,denRot,den1Rot,den1imRot,starsq,killcont)
      !! Collinear case: put V1Theta+VTheta1 into V1%pw_w together
      !! Noco case: Correctly rotate back the potential into a 2x2 matrix (TODO)
        USE m_types
        USE m_constants
        USE m_get_int_perturbation
        USE m_get_mt_perturbation
        USE m_fft3d
        !USE m_rotate_mt_den_tofrom_local

        IMPLICIT NONE

        TYPE(t_mpi),      INTENT(IN)    :: fmpi
        TYPE(t_noco),     INTENT(IN)    :: noco
        TYPE(t_nococonv), INTENT(IN)    :: nococonv
        TYPE(t_sym),      INTENT(IN)    :: sym
        TYPE(t_juphon),   INTENT(IN)    :: juphon
        TYPE(t_stars),    INTENT(IN)    :: stars, starsq
        TYPE(t_atoms),    INTENT(IN)    :: atoms
        TYPE(t_input),    INTENT(IN)    :: input
        TYPE(t_sphhar),   INTENT(IN)    :: sphhar
        TYPE(t_potden),   INTENT(IN)    :: vTot
        TYPE(t_potden),   INTENT(INOUT) :: vTot1, vTot1imag, denRot, den1Rot, den1imRot
        INTEGER, INTENT(IN) :: killcont(2)

        INTEGER                         :: i, js, ifft3, iPhonon
        REAL,    ALLOCATABLE :: fftwork(:), vre(:), v1re(:), v1im(:)
        COMPLEX, ALLOCATABLE :: v1full(:)

        iPhonon = 0
        IF (juphon%l_phonon) iPhonon = 1

        ifft3 = 27*stars%mx1*stars%mx2*stars%mx3 !TODO: What if starsq/=stars in that regard?

        ALLOCATE(vre(ifft3),v1re(ifft3),v1im(ifft3),v1full(ifft3),fftwork(ifft3))

        vTot1%pw_w = CMPLX(0.0,0.0)

        IF (.NOT.noco%l_noco) THEN
            DO js = 1, input%jspins
               vre = 0.0; v1re = 0.0; v1im = 0.0
               CALL fft3d(v1re, v1im, vTot1%pw(:, js), starsq, +1)
               CALL fft3d(vre, fftwork, vTot%pw(:, js), stars, +1)
               v1full = killcont(1)*CMPLX(v1re,v1im) * stars%ufft + iPhonon * killcont(2) * vre * starsq%ufft1!-q
               v1re =  REAL(v1full)
               v1im = AIMAG(v1full)
               CALL fft3d(v1re, v1im, vTot1%pw_w(:, js), starsq, -1)
            END DO
        ELSE IF(noco%l_noco) THEN
            CALL get_int_global_perturbation(stars,atoms,sym,input,denRot,den1Rot,den1imRot,vTot,vTot1,starsq)
            IF (any(noco%l_unrestrictMT)) CALL get_mt_global_perturbation(atoms,sphhar,sym,denRot,den1Rot,den1imRot,noco,vTot,vTot1,vTot1imag)
        END IF

    END SUBROUTINE dfpt_vgen_finalize

END MODULE m_dfpt_vgen_finalize
