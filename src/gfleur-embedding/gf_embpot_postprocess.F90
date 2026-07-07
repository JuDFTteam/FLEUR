!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions 
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
MODULE m_gf_embpot_postprocess
!Created by wortmann on Sep 10, 2012
    IMPLICIT NONE
    CONTAINS
    SUBROUTINE gf_embpot_spectral(jspin,nk,gfinp,layers,lapw,lapw_gf,cell,kpts)
    !This subroutine calculates the green function at the embedding plane
    !It uses the embedding potential of the layer=1,side=1 as input
    !Additionally, the embedding potential of layer=layers%num_layers,side=2 can be used
    !A surface can be simulated by calling the appropriate vacuum routines
    !The surface gr
    USE m_gf_types
    USE m_gf_embedding
    USE m_gf_math
    USE m_gf_vacuum
    USE m_gf_energies
    USE m_gf_io2dmat
    IMPLICIT NONE
    INTEGER,INTENT(IN)       :: jspin,nk
    TYPE(t_embinp),INTENT(IN) :: gfinp
    TYPE(t_layers),INTENT(IN):: layers
    TYPE(t_lapw),INTENT(IN)  :: lapw
      TYPE(t_lapw_gf),INTENT(IN) :: lapw_gf
    TYPE(t_cell),INTENT(IN)  :: cell
    TYPE(t_kpts),INTENT(IN)  :: kpts

    COMPLEX,DIMENSION(lapw_gf%nv2_tot,lapw_gf%nv2_tot)::sig1,sig2,imag_matrix
    INTEGER :: n,en

    sig2=0.0
    imag_matrix=0.0

    IF (abs(gfinp%imag_broad)>epsilon(1.0)) THEN
      DO n=1,size(sig2,1)
         imag_matrix(n,n)=cmplx(0.0,gfinp%imag_broad)
      ENDDO
    ENDIF

    DO en=1,gf_noen()
        CALL gf_GETEMB2(sig1,1,1,en,nk,jspin,lapw,lapw_gf)
        IF (gf_io2dmatFID(IO2D_EMB,layers%num_layers,2,IO2D_READ)>0) CALL gf_GETEMB2(sig2,layers%num_layers,2,en,nk,jspin,lapw,lapw_gf)
        IF (gfinp%l_simplevacuum) THEN
            CALL gf_simple_vac(en,nk,lapw,lapw_gf,kpts,cell,gfinp%vacuum_energy,sig2,gfinp%efield)
        ENDIF
        sig1=sig1+sig2+imag_matrix

        sig1=mat_inverse(sig1)
        CALL gf_write2dmat(IO2D_gmat,1,1,en,nk,jspin,lapw_gf,sig1)
    ENDDO
    END SUBROUTINE gf_embpot_spectral
END MODULE m_gf_embpot_postprocess
