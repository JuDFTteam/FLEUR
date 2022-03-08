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

   SUBROUTINE dfpt_vgen_finalize(fmpi,atoms,stars,sym,noco,nococonv,input,sphhar,vTot,vTot1,vTot1imag,denRot,den1Rot,den1imRot,starsq)

        USE m_types
        USE m_constants
        USE m_get_int_perturbation
        !USE m_rotate_mt_den_tofrom_local

        IMPLICIT NONE

        TYPE(t_mpi),      INTENT(IN)    :: fmpi
        TYPE(t_noco),     INTENT(IN)    :: noco
        TYPE(t_nococonv), INTENT(INOUT) :: nococonv
        TYPE(t_sym),      INTENT(IN)    :: sym
        TYPE(t_stars),    INTENT(IN)    :: stars, starsq
        TYPE(t_atoms),    INTENT(IN)    :: atoms
        TYPE(t_input),    INTENT(IN)    :: input
        TYPE(t_sphhar),   INTENT(IN)    :: sphhar
        TYPE(t_potden),   INTENT(INOUT) :: vTot, vTot1, vTot1imag, denRot, den1Rot, den1imRot

        INTEGER                         :: i, js

        IF (.NOT.noco%l_noco) THEN
            ! TODO: Do we need this for dfpt??
            DO js=1,SIZE(vtot1%pw_w, 2)
                DO i=1,stars%ng3
                    vTot1%pw_w(i,js)=vTot1%pw_w(i,js) / starsq%nstr(i)
                END DO
            END DO
        ELSE IF(noco%l_noco) THEN
            CALL get_int_global_perturbation(stars,atoms,sym,input,denRot,den1Rot,den1imRot,vTot,vTot1,starsq)
            !IF (any(noco%l_unrestrictMT)) CALL rotate_mt_den_from_local(atoms,sphhar,sym,denRot,noco,vtot)
        END IF

    END SUBROUTINE dfpt_vgen_finalize

END MODULE m_dfpt_vgen_finalize
