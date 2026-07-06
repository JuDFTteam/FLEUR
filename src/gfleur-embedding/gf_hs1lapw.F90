!--------------------------------------------------------------------------------
! Copyright (c) 2026 Peter Grünberg Institut, Forschungszentrum Jülich, Germany
! This file is part of FLEUR and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
      MODULE m_gf_hs1lapw
      use m_juDFT
      IMPLICIT NONE
      PRIVATE

      PUBLIC gf_hs1lapw
      CONTAINS
      !<-- S: gf_hs1lapw()

      SUBROUTINE gf_hs1lapw(jspin,layer,nk,embinp,ld,gmpi,bk,lapw,      &
     &                      lapw_gf,vpw)
!*****************************************************************
!     DESC:This subroutine generates the FLAPW-Hamiltonian and the
!     Overlap Matrix by calling the appropriate FLEUR-subroutines.
!
!     In the port to current FLEUR:
!     - the interstitial part comes from gf_hs_int (GF step function
!       on the full G-difference box, GF potential vpw)
!     - the MT part comes from the standard hsmt working on the
!       tlmplm/usdus data prepared by gf_tlmplminit (local_ham)
!     - the matrices are dense t_mat objects; the assembled H,S are
!       stored Hermitian-completed in gf_hsdata
!     Daniel Wortmann, Sat Mar  2 13:46:30 2002
!*****************************************************************

      USE m_gf_types
      USE m_gf_hsdata
      USE m_gf_hs_int
      USE m_hsmt
      USE m_gf_spectrum,ONLY: gf_spectrum_clean
      use m_gf_stepsanaly, only: gf_initstepsanaly, gf_stepf_nohelpregion
      IMPLICIT NONE
      !<--Arguments
      INTEGER,INTENT(IN)            :: jspin
      INTEGER,INTENT(IN)            :: layer
      INTEGER,INTENT(IN)            :: nk
      TYPE(t_embinp),INTENT(IN)     :: embinp
      TYPE(t_gflayer),INTENT(INOUT) :: ld
      TYPE(t_gfmpi),INTENT(INOUT)   :: gmpi
                                                    ! Kpt
      REAL, INTENT(IN)              :: Bk(3)
      TYPE(t_lapw),INTENT(INOUT)    :: lapw
      TYPE(t_lapw_gf),INTENT(INOUT) :: lapw_gf
      !the GF interstitial potential (warped, zero in aux volumes)
      COMPLEX,INTENT(IN)            :: vpw(:,:)

      !>
      !<-- Locals
      INTEGER                   :: nbas,err
      COMPLEX,ALLOCATABLE       :: ustep(:,:,:)
      CLASS(t_mat),ALLOCATABLE  :: hmat(:,:),smat(:,:)

      IF (ld%fi%noco%l_noco) CALL juDFT_error(                           &
     &     "noco not yet supported in the gfleur port",                  &
     &     calledby="gf_hs1lapw")

      nbas = lapw%nv(jspin) + ld%fi%atoms%nlotot

      !<-- the GF step function on the full G-difference box
      ALLOCATE(ustep(-ld%stars%mx1:ld%stars%mx1,                         &
     &               -ld%stars%mx2:ld%stars%mx2,                         &
     &        -2*embinp%napw(layer):2*embinp%napw(layer)),stat=err)
      if(err/=0) then
        write(*,*) "Step function:",layer
        write(*,*) ld%stars%mx1,ld%stars%mx2,embinp%napw(layer)
        call juDFT_error("Could not allocate ustep",calledby="gf_hs1lapw")
      endif
      IF(.NOT.embinp%l_nohelpregion)THEN
         CALL juDFT_error("l_nohelpregion required",calledby              &
     &        ="gf_hs1lapw")
      ELSE
         CALL gf_initstepsanaly(ld%stars,embinp%napw(layer))
         CALL gf_stepf_nohelpregion(layer,ld%stars%mx1                   &
     &        ,ld%stars%mx2,embinp%napw(layer),ustep)
      ENDIF
      !>

      IF ( embinp%l_solwil ) THEN
         CALL juDFT_error("Soler Wiliams Basis not supported anymore !")
      ENDIF

      CALL  gf_spectrum_clean()

      !<-- allocate the dense matrices (always complex)
      ALLOCATE(t_mat::smat(1,1),hmat(1,1))
      CALL smat(1,1)%init(.FALSE.,nbas,nbas)
      CALL hmat(1,1)%init(smat(1,1))
      !>

      !<-- interstitial contribution to H+S
      CALL timestart("gf_hsint")
      CALL gf_hs_int(lapw,ld%stars,embinp%napw(layer),ustep,jspin,       &
     &               gmpi%fmpi,bk,ld%fi%cell%bbmat,vpw,                  &
     &               hmat(1,1),smat(1,1))
      CALL timestop("gf_hsint")
      !>
      !<-- MT-contribution to H and S
      CALL timestart("hssphn")
      CALL hsmt(ld%fi%atoms,ld%fi%sym,ld%enpara,jspin,ld%fi%input,       &
     &          gmpi%fmpi,ld%fi%noco,ld%nococonv,ld%fi%cell,lapw,        &
     &          ld%usdus,ld%tlmplm,smat,hmat)
      CALL timestop("hssphn")
      !>

      !<-- store the assembled matrices (Hermitian completion inside)
      CALL gf_storeHS(jspin,nk,layer,hmat(1,1)%data_c,smat(1,1)%data_c)
      IF ((embinp%l_charge.or.embinp%l_dos.OR.embinp%l_writehs.OR.       &
     &     .NOT.embinp%l_spectral).AND.gmpi%fmpi%n_rank == 0)            &
     &     CALL gf_writeHS(embinp%l_savemem)
      !>

      END SUBROUTINE

      !>

      END
