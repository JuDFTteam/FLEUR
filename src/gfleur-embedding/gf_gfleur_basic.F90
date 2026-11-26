!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_gfleur_basic
      use m_juDFT
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE gf_gfleur_basic(jspin,layer,nk,embinp,layers,ld,gmpi,   &
     &     bk,lapw,lapw_gf,pot_aux)
!*************************************************
!     Set up of the Hamiltonian (always).
!     In the spectral mode: Diagonalization of the
!     Hamiltonian.
!     If the number of layers is larger than one
!     the transfer-matrix of the layer is computed
!     in addition.
!     Frank Freimuth, November 2007
!
!     Port to current FLEUR: works on the per-layer t_gflayer container
!     (ld) and calls the modern gf_hs1lapw (hsmt + gf_hs_int on dense
!     t_mat). LDA+U/SOC and the duplicated v25 datatypes are gone.
!*************************************************
      USE m_gf_types
      USE m_gf_hs1lapw
      USE m_gf_spectrum
      USE m_gf_get_spectrum
      USE m_gf_gfcn
      USE m_gf_tmat, ONLY  : gf_tmat_layers
      USE m_gf_energies, ONLY: gf_noen
      IMPLICIT NONE
      INTEGER,INTENT(IN)            :: jspin
      INTEGER,INTENT(IN)            :: layer
      INTEGER,INTENT(IN)            :: nk
      TYPE(t_embinp),INTENT(IN)     :: embinp
      TYPE(t_layers),INTENT(IN)     :: layers
      TYPE(t_gflayer),INTENT(INOUT) :: ld
      TYPE(t_gfmpi),INTENT(INOUT)   :: gmpi
                                                    ! Kpt list
      REAL, INTENT(IN)              :: Bk(:,:)
      TYPE(t_lapw),INTENT(INOUT)    :: lapw
      TYPE(t_lapw_gf),INTENT(INOUT) :: lapw_gf
      COMPLEX,INTENT(IN)            :: pot_aux(:,:)

      COMPLEX, ALLOCATABLE :: g(:,:)
      COMPLEX, ALLOCATABLE :: gij(:,:,:)
      INTEGER :: en
      LOGICAL :: l_writespectrum

      !<-- FLAPW Hamiltonian + overlap of the layer at this (k,spin)
      CALL timestart("gf_hs1lapw")
      CALL gf_hs1lapw(jspin,layer,nk,embinp,ld,gmpi,bk(1:3,nk),lapw,      &
     &                lapw_gf,ld%vTot%pw)
      CALL timestop("gf_hs1lapw")
      !>

      IF (embinp%l_spectral) THEN
         l_writespectrum = .FALSE.
         CALL timestart("gf_get_spectrum")
         CALL gf_spectrum_clean()
         CALL gf_get_spectrum(layer,jspin,embinp,ld%fi%cell,lapw,lapw_gf, &
     &        .FALSE.,embinp%l_fullgreen,embinp%l_nogno,                  &
     &        embinp%l_nohelpregion,l_writespectrum,nk,.FALSE.)
         CALL timestop("gf_get_spectrum")
      ENDIF

      !<-- calculate the transfer matrix of a stack of layers
      IF (layers%num_layers > 1) THEN
         ALLOCATE( gij(2*lapw_gf%nv2_tot,2*lapw_gf%nv2_tot,1) )
         DO en = 1,gf_noen()
            IF (embinp%napw(layer) > -1) THEN
               !<-- standard calculation of the retarded GF
               CALL timestart("gf_gfcn")
               CALL gf_GFCN(.TRUE.,layer,en,nk,jspin,ld%fi%cell,lapw,     &
     &              lapw_gf,embinp,ld%fi%noco%l_noco,ld%fi%sym%invs,      &
     &              .FALSE.,g,gij,                                        &
     &              real_energy=embinp%l_CBS.OR.embinp%curr.NE.0)
               CALL timestop("gf_gfcn")
            ELSE
               CALL juDFT_error("BUG:vacuum in gfleur_basic")
            ENDIF
            CALL timestart("gf_tmat")
            IF (pot_aux(1,jspin) /= pot_aux(2,jspin)) CALL                &
     &         juDFT_error("Cannot treat different aux_potentials")
            CALL gf_TMAT_layers(layer,layers,en,nk,jspin,gmpi,lapw,       &
     &           lapw_gf,embinp,pot_aux(1,jspin),gij(:,:,1))
            CALL timestop("gf_tmat")
         ENDDO
         DEALLOCATE( gij )
      ENDIF
      !>
      END SUBROUTINE gf_gfleur_basic
      END
